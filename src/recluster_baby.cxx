#include "recluster_baby.hpp"

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "fastjet/ClusterSequence.hh"

#include "timer.hpp"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

using namespace std;
using namespace fastjet;

namespace{
  string input_file_path = "";
  string output_file_path = "";
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  if(input_file_path == "") ERROR("Must provide input file with -i");
  if(output_file_path == "") ERROR("Must provide output file with -o");

  auto pos = input_file_path.rfind("/");
  string file_name = input_file_path;
  if(pos != string::npos){
    file_name = input_file_path.substr(pos+1);
  }

  TFile input_file(input_file_path.c_str(), "read");
  if(!input_file.IsOpen()) ERROR("Could not open input file "+input_file_path);
  TFile output_file(output_file_path.c_str(), "recreate");
  if(!output_file.IsOpen()) ERROR("Could not open output file "+output_file_path);
  output_file.cd();

  TTree *in_treeglobal = static_cast<TTree*>(input_file.Get("treeglobal"));
  if(in_treeglobal == nullptr) ERROR("No treeglobal in input file");

  TTree *out_treeglobal = in_treeglobal->CloneTree(0);
  in_treeglobal->CopyAddresses(out_treeglobal);

  int num_entries = in_treeglobal->GetEntries();
  for(int event = 0; event < num_entries; ++event){
    in_treeglobal->GetEntry(event);
    out_treeglobal->Fill();
  }
  out_treeglobal->Write();
  
  TTree *in_tree = static_cast<TTree*>(input_file.Get("tree"));
  if(in_tree == nullptr) ERROR("No tree in input file");

  TTree *out_tree = in_tree->CloneTree(0);
  in_tree->CopyAddresses(out_tree);

  float mj14 = 0., mj14_original = 0., mj14_validation = 0.;
  vector<float> *jets_pt = nullptr, *jets_eta = nullptr,
    *jets_phi = nullptr, *jets_m = nullptr;
  vector<bool> *jets_islep = nullptr;
  TBranch *b_jets_pt = nullptr, *b_jets_eta = nullptr,
    *b_jets_phi = nullptr, *b_jets_m = nullptr, *b_jets_islep = nullptr;

  in_tree->SetBranchAddress("jets_pt", &jets_pt, &b_jets_pt);
  in_tree->SetBranchAddress("jets_eta", &jets_eta, &b_jets_eta);
  in_tree->SetBranchAddress("jets_phi", &jets_phi, &b_jets_phi);
  in_tree->SetBranchAddress("jets_m", &jets_m, &b_jets_m);
  in_tree->SetBranchAddress("jets_islep", &jets_islep, &b_jets_islep);
  
  in_tree->SetBranchAddress("mj14", &mj14_original);
  out_tree->Branch("mj14_original", &mj14_original);
  out_tree->SetBranchAddress("mj14", &mj14);
  out_tree->Branch("mj14_validation", &mj14_validation);

  num_entries = in_tree->GetEntries();
  TLorentzVector lv;
  PseudoJet pj;
  vector<PseudoJet> jets_with_lep, jets_no_lep;
  vector<PseudoJet> fjets_with_lep, fjets_no_lep;
  JetDefinition jd(antikt_algorithm, 1.4);
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int event = 0; event < num_entries; ++event){
    timer.Iterate();
    in_tree->GetEntry(event);

    jets_with_lep.clear();
    jets_no_lep.clear();
    
    for(size_t i = 0; i < jets_pt->size(); ++i){
      bool pass_jet = jets_pt->at(i) > 30. && fabs(jets_eta->at(i)) <= 2.4;
      bool pass_lep = jets_islep->at(i);
      if(!(pass_jet || pass_lep)) continue;
      lv.SetPtEtaPhiM(jets_pt->at(i),
                      jets_eta->at(i),
                      jets_phi->at(i),
                      jets_m->at(i));
      pj = PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E());

      jets_with_lep.push_back(pj);
      if(pass_jet && !pass_lep) jets_no_lep.push_back(pj);
    }

    ClusterSequence cs_with_lep(jets_with_lep, jd) , cs_no_lep(jets_no_lep, jd);
    fjets_with_lep = cs_with_lep.inclusive_jets();
    fjets_no_lep = cs_no_lep.inclusive_jets();

    mj14 = 0.;
    for(const auto &fjet: fjets_no_lep){
      mj14 += fjet.m();
    }
    mj14_validation = 0.;
    for(const auto &fjet: fjets_with_lep){
      mj14_validation += fjet.m();
    }

    out_tree->Fill();
  }
  out_tree->Write();

  output_file.Close();
  input_file.Close();
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"input_file", no_argument, 0, 'i'},
      {"output_file", no_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "i:o:", long_options, &option_index);
    
    if( opt == -1) break;
    
    string optname;
    switch(opt){
    case 'i':
      input_file_path = optarg;
      break;
    case 'o':
      output_file_path = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
      }else{
        ERROR("Bad option! Found option name "+string(optname));
      }
      break;
    default:
      ERROR("Bad option! getopt_long returned character code "+to_string(opt));
      break;
    }
  }
}

