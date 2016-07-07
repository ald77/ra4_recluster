#include "mismeasure_leptons.hpp"

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "fastjet/ClusterSequence.hh"

#include "timer.hpp"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

using namespace std;
using namespace fastjet;

namespace{
  string input_file_path = "";
  string output_file_path = "";
  set<float> p_vals = {0.001, 0.0001, 0.00001};
  set<float> a_vals = {1., 5., 10.};
  set<float> b_vals = {0., 100., 400.};

  TRandom3 prng(3556329863);

  int FindLep(int lep_type, float new_pt, float new_miniso,
              float pt, float eta, float phi,
              const vector<float> &pts,
              const vector<float> &etas,
              const vector<float> &phis,
              const vector<float> &ids){
    for(size_t i = 0; i < pts.size(); ++i){
      if(pt == pts.at(i)
         && eta == etas.at(i)
         && phi == phis.at(i)
         && abs(lep_type) == fabs(ids.at(i))){
        return i;
      }
    }
    if(new_pt>20.
       && ((lep_type==11&&new_miniso<0.1)
           || (lep_type==13&&new_miniso<0.2))){
      return pts.size();
    }
    return -1;
  }

  float Area(float pt){
    float radius;
    if(pt<50.){
      radius =  0.2;
    }else if(pt<200.){
      radius =  10./pt;
    }else{
      radius =  0.05;
    }
    return radius*radius;
  }

  int FindJet(float eta, float phi, bool islep_only,
              const vector<float> &etas,
              const vector<float> &phis,
              const vector<bool> &isleps){
    float min_delta_r = 9999.;
    int best_jet = -1;
    for(size_t i = 0; i < etas.size(); ++i){
      if(islep_only && !isleps.at(i)) continue;
      float delta_r = hypot(TVector2::Phi_mpi_pi(phi-phis.at(i)),
                            eta-etas.at(i));
      if(delta_r < min_delta_r){
        min_delta_r = delta_r;
        best_jet = i;
      }
    }
    return best_jet;
  }

  float MT(float pta, float phia, float ptb, float phib){
    return sqrt(2.*pta*ptb*(1.-cos(phib-phia)));
  }
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

  float mt = 0., ht = 0., met = 0., met_phi = 0., mj14 = 0.;
  vector<float> *jets_pt = nullptr, *jets_eta = nullptr,
    *jets_phi = nullptr, *jets_m = nullptr, *jets_csv = nullptr;
  vector<bool> * jets_islep = nullptr;
  int nels = 0, nmus = 0, nleps = 0;
  int njets = 0, nbm = 0;
  vector<float> *mus_pt = nullptr, *mus_eta = nullptr,
    *mus_phi = nullptr, *mus_miniso = nullptr;
  vector<float> *els_pt = nullptr, *els_eta = nullptr,
    *els_phi = nullptr, *els_miniso = nullptr;
  vector<float> *leps_pt = nullptr, *leps_eta = nullptr,
    *leps_phi = nullptr, *leps_id = nullptr;
  TBranch *b_jets_pt = nullptr, *b_jets_eta = nullptr, *b_jets_phi = nullptr,
    *b_jets_m = nullptr, *b_jets_csv = nullptr, *b_jets_islep = nullptr,
    *b_leps_pt = nullptr, *b_leps_eta = nullptr,
    *b_leps_phi = nullptr, *b_leps_id = nullptr,
    *b_mus_pt = nullptr, *b_mus_eta = nullptr,
    *b_mus_phi = nullptr, *b_mus_miniso = nullptr,
    *b_els_pt = nullptr, *b_els_eta = nullptr,
    *b_els_phi = nullptr, *b_els_miniso = nullptr;

  in_tree->SetBranchAddress("jets_pt", &jets_pt, &b_jets_pt);
  in_tree->SetBranchAddress("jets_eta", &jets_eta, &b_jets_eta);
  in_tree->SetBranchAddress("jets_phi", &jets_phi, &b_jets_phi);
  in_tree->SetBranchAddress("jets_m", &jets_m, &b_jets_m);
  in_tree->SetBranchAddress("jets_csv", &jets_csv, &b_jets_csv);
  in_tree->SetBranchAddress("jets_islep", &jets_islep, &b_jets_islep);
  in_tree->SetBranchAddress("leps_pt", &leps_pt, &b_leps_pt);
  in_tree->SetBranchAddress("leps_eta", &leps_eta, &b_leps_eta);
  in_tree->SetBranchAddress("leps_phi", &leps_phi, &b_leps_phi);
  in_tree->SetBranchAddress("leps_id", &leps_id, &b_leps_id);
  in_tree->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
  in_tree->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
  in_tree->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
  in_tree->SetBranchAddress("mus_miniso", &mus_miniso, &b_mus_miniso);
  in_tree->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
  in_tree->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
  in_tree->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
  in_tree->SetBranchAddress("els_miniso", &els_miniso, &b_els_miniso);
  in_tree->SetBranchAddress("njets", &njets);
  in_tree->SetBranchAddress("nleps", &nleps);
  in_tree->SetBranchAddress("nmus", &nmus);
  in_tree->SetBranchAddress("nels", &nels);
  in_tree->SetBranchAddress("nbm", &nbm);
  in_tree->SetBranchAddress("mt", &mt);
  in_tree->SetBranchAddress("ht", &ht);
  in_tree->SetBranchAddress("met", &met);
  in_tree->SetBranchAddress("met_phi", &met_phi);
  in_tree->SetBranchAddress("mj14", &mj14);

  TTree *out_tree = in_tree->CloneTree(0);
  in_tree->CopyAddresses(out_tree);

  vector<int> *mm_jet_index = nullptr,
    *mm_lep_index = nullptr,
    *mm_mu_index = nullptr,
    *mm_el_index = nullptr,
    *mm_njets = nullptr,
    *mm_nleps = nullptr,
    *mm_nmus = nullptr,
    *mm_nels = nullptr,
    *mm_nbm = nullptr;
  vector<float> *mm_jet_pt = nullptr,
    *mm_jet_eta = nullptr,
    *mm_jet_phi = nullptr,
    *mm_jet_m = nullptr,
    *mm_lep_pt = nullptr,
    *mm_mu_pt = nullptr,
    *mm_el_pt = nullptr,
    *mm_mu_miniso = nullptr,
    *mm_el_miniso = nullptr,
    *mm_mt = nullptr,
    *mm_ht = nullptr,
    *mm_met = nullptr,
    *mm_met_phi = nullptr,
    *mm_mj14_lep = nullptr,
    *mm_mj14_nolep = nullptr;
  vector<bool> *mm_jet_islep = nullptr;

  out_tree->Branch("mm_jet_index", &mm_jet_index);
  out_tree->Branch("mm_lep_index", &mm_lep_index);
  out_tree->Branch("mm_mu_index", &mm_mu_index);
  out_tree->Branch("mm_el_index", &mm_el_index);
  out_tree->Branch("mm_njets", &mm_njets);
  out_tree->Branch("mm_nleps", &mm_nleps);
  out_tree->Branch("mm_nmus", &mm_nmus);
  out_tree->Branch("mm_nels", &mm_nels);
  out_tree->Branch("mm_jet_pt", &mm_jet_pt);
  out_tree->Branch("mm_jet_eta", &mm_jet_eta);
  out_tree->Branch("mm_jet_phi", &mm_jet_phi);
  out_tree->Branch("mm_jet_m", &mm_jet_m);
  out_tree->Branch("mm_jet_islep", &mm_jet_islep);
  out_tree->Branch("mm_lep_pt", &mm_lep_pt);
  out_tree->Branch("mm_mu_pt", &mm_mu_pt);
  out_tree->Branch("mm_el_pt", &mm_el_pt);
  out_tree->Branch("mm_mu_miniso", &mm_mu_miniso);
  out_tree->Branch("mm_el_miniso", &mm_el_miniso);
  out_tree->Branch("mm_mt", &mm_mt);
  out_tree->Branch("mm_ht", &mm_ht);
  out_tree->Branch("mm_nbm", &mm_nbm);
  out_tree->Branch("mm_met", &mm_met);
  out_tree->Branch("mm_met_phi", &mm_met_phi);
  out_tree->Branch("mm_mj14_lep", &mm_mj14_lep);
  out_tree->Branch("mm_mj14_nolep", &mm_mj14_nolep);

  size_t num_configs = p_vals.size() * a_vals.size() * b_vals.size();
  mm_jet_index->resize(num_configs);
  mm_lep_index->resize(num_configs);
  mm_mu_index->resize(num_configs);
  mm_el_index->resize(num_configs);
  mm_njets->resize(num_configs);
  mm_nleps->resize(num_configs);
  mm_nmus->resize(num_configs);
  mm_nels->resize(num_configs);
  mm_jet_pt->resize(num_configs);
  mm_jet_eta->resize(num_configs);
  mm_jet_phi->resize(num_configs);
  mm_jet_m->resize(num_configs);
  mm_jet_islep->resize(num_configs);
  mm_lep_pt->resize(num_configs);
  mm_mu_pt->resize(num_configs);
  mm_el_pt->resize(num_configs);
  mm_mu_miniso->resize(num_configs);
  mm_el_miniso->resize(num_configs);
  mm_mt->resize(num_configs);
  mm_ht->resize(num_configs);
  mm_nbm->resize(num_configs);
  mm_met->resize(num_configs);
  mm_met_phi->resize(num_configs);
  mm_mj14_lep->resize(num_configs);
  mm_mj14_nolep->resize(num_configs);

  //Caching
  vector<int> vint;
  TLorentzVector lv;
  PseudoJet pj;
  vector<PseudoJet> jets_with_lep, jets_no_lep;
  vector<PseudoJet> fjets_with_lep, fjets_no_lep;
  JetDefinition jd(antikt_algorithm, 1.4);

  num_entries = in_tree->GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int event = 0; event < num_entries; ++event){
    timer.Iterate();
    in_tree->GetEntry(event);

    size_t iconfig = 0;
    for(auto p: p_vals){
      for(auto a: a_vals){
	for(auto b: b_vals){
	  ++iconfig;

	  //Pick a lepton to modify
	  int ntot = mus_pt->size() + els_pt->size();
	  vint.clear();
	  for(int i = 0; i < ntot; ++i){
	    if(prng.Uniform() < p) vint.push_back(i);
	  }
	  int mod_lep = -1;
	  if(vint.size()>=1){
	    mod_lep = vint.at(prng.Integer(vint.size()));
	  }

	  if(mod_lep < 0){
	    //Event unchanged
	    mm_jet_index->at(iconfig) = -1;
	    mm_lep_index->at(iconfig) = -1;
	    mm_mu_index->at(iconfig) = -1;
	    mm_el_index->at(iconfig) = -1;
	    mm_njets->at(iconfig) = njets;
	    mm_nleps->at(iconfig) = nleps;
	    mm_nmus->at(iconfig) = nmus;
	    mm_nels->at(iconfig) = nels;
	    mm_jet_pt->at(iconfig) = -1.;
	    mm_jet_eta->at(iconfig) = 0.;
	    mm_jet_phi->at(iconfig) = 0.;
	    mm_jet_m->at(iconfig) = -1.;
	    mm_jet_islep->at(iconfig) = false;
	    mm_lep_pt->at(iconfig) = -1.;
	    mm_mu_pt->at(iconfig) = -1.;
	    mm_el_pt->at(iconfig) = -1.;
	    mm_mu_miniso->at(iconfig) = -1.;
	    mm_el_miniso->at(iconfig) = -1.;
	    mm_mt->at(iconfig) = mt;
	    mm_ht->at(iconfig) = ht;
	    mm_nbm->at(iconfig) = nbm;
	    mm_met->at(iconfig) = met;
	    mm_met_phi->at(iconfig) = met_phi;
	    mm_mj14_lep->at(iconfig) = mj14;
	    mm_mj14_nolep->at(iconfig) = 0.;

	    jets_no_lep.clear();
    	    for(size_t i = 0; i < jets_pt->size(); ++i){
	      bool pass_jet = jets_pt->at(i) > 30. && fabs(jets_eta->at(i)) <= 2.4;
	      bool pass_lep = jets_islep->at(i);
	      if(!pass_jet || pass_lep) continue;
	      lv.SetPtEtaPhiM(jets_pt->at(i),
			      jets_eta->at(i),
			      jets_phi->at(i),
			      jets_m->at(i));
	      pj = PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E());
	      jets_no_lep.push_back(pj);
	    }

	    ClusterSequence cs_no_lep(jets_no_lep, jd);
	    fjets_no_lep = cs_no_lep.inclusive_jets();
	    for(const auto &fjet: fjets_no_lep){
	      mm_mj14_nolep->at(iconfig) += fjet.m();
	    }
	  }else{
            //Event changed

            //Pick a lepton
            int lep_type;
            float old_pt, eta, phi, old_miniso;
            if(mod_lep < static_cast<int>(mus_pt->size())){
              lep_type = 13;
              old_pt = mus_pt->at(mod_lep);
              eta = mus_eta->at(mod_lep);
              phi = mus_phi->at(mod_lep);
              old_miniso = mus_miniso->at(mod_lep);
              mm_mu_index->at(iconfig) = mod_lep;
              mm_el_index->at(iconfig) = -1;
            }else{
              mod_lep -= mus_pt->size();
              lep_type = 11;
              old_pt = els_pt->at(mod_lep);
              eta = els_eta->at(mod_lep);
              phi = els_phi->at(mod_lep);
              old_miniso = els_miniso->at(mod_lep);
              mm_mu_index->at(iconfig) = -1;
              mm_el_index->at(iconfig) = mod_lep;
            }

            //Modify the lepton
            float new_pt = a*old_pt+b;
            TLorentzVector new_lep, old_lep;
            new_lep.SetPtEtaPhiM(new_pt, eta, phi, 0.);
            old_lep.SetPtEtaPhiM(old_pt, eta, phi, 0.);
            float new_miniso = old_miniso * (old_pt/new_pt) * (Area(new_pt)/Area(old_pt));
            int lep_index = FindLep(lep_type, new_pt, new_miniso,
                                    old_pt, eta, phi,
                                    *leps_pt, *leps_eta, *leps_phi, *leps_id);
            mm_lep_index->at(iconfig) = lep_index;
            if(lep_index >= 0){
              mm_nleps->at(iconfig) = nleps+1;
              mm_lep_pt->at(iconfig) = new_pt;
            }else{
              mm_nleps->at(iconfig) = nleps;
              mm_lep_pt->at(iconfig) = -1;
            }
            if(lep_type==11){
              mm_el_pt->at(iconfig) = new_pt;
              mm_el_miniso->at(iconfig) = new_miniso;
              if(lep_index >= 0) mm_nels->at(iconfig) = nels+1;
              else mm_nels->at(iconfig) = nels;
              mm_mu_pt->at(iconfig) = -1.;
              mm_mu_miniso->at(iconfig) = -1.;
              mm_nmus->at(iconfig) = nmus;
            }else{
              mm_el_pt->at(iconfig) = -1.;
              mm_el_miniso->at(iconfig) = -1.;
              mm_nels->at(iconfig) = nels;
              mm_mu_pt->at(iconfig) = new_pt;
              mm_mu_miniso->at(iconfig) = new_miniso;
              if(lep_index >= 0) mm_nmus->at(iconfig) = nmus+1;
              else mm_nmus->at(iconfig) = nmus;
            }

            //Modify associated jet
            int jet_index = FindJet(eta, phi, lep_index>=0 && lep_index<static_cast<int>(leps_pt->size()),
                                    *jets_eta, *jets_phi, *jets_islep);
            mm_jet_index->at(iconfig) = jet_index;
            TLorentzVector new_jet, old_jet;
            old_jet.SetPtEtaPhiM(jets_pt->at(jet_index),
                                 jets_eta->at(jet_index),
                                 jets_phi->at(jet_index),
                                 jets_m->at(jet_index));
            new_jet = old_jet + (new_lep-old_lep);
            mm_jet_pt->at(iconfig) = new_jet.Pt();
            mm_jet_eta->at(iconfig) = new_jet.Eta();
            mm_jet_phi->at(iconfig) = new_jet.Phi();
            mm_jet_m->at(iconfig) = new_jet.M();
            mm_jet_islep->at(iconfig) = lep_index >= 0;
            if(lep_index >= 0 && !jets_islep->at(jet_index)){
              //Jet became a lepton
              mm_njets->at(iconfig) = nleps-1;
              if(jets_csv->at(jet_index) > 0.8) mm_nbm->at(iconfig) = nbm-1;
              else mm_nbm->at(iconfig) = nbm;
            }else{
              //Jet is not a lepton or already was a lepton
              mm_njets->at(iconfig) = nleps;
              mm_nbm->at(iconfig) = nbm;
            }

            //Recompute HT and MJ
            mm_ht->at(iconfig) = 0.;
            for(size_t i = 0; i < jets_pt->size(); ++i){
              if(static_cast<int>(i) != jet_index){
                lv.SetPtEtaPhiM(jets_pt->at(i), jets_eta->at(i),
                                jets_phi->at(i), jets_m->at(i));
              }else{
                lv = new_jet;
              }
              bool pass_jet = lv.Pt()>30. && fabs(lv.Eta())<=2.4;
              bool pass_lep = jets_islep->at(i);
              if(!(pass_jet || pass_lep)) continue;
              pj = PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E());

              jets_with_lep.push_back(pj);
              if(pass_jet && !pass_lep){
                jets_no_lep.push_back(pj);
                mm_ht->at(iconfig) += lv.Pt();
              }
            }
            ClusterSequence cs_with_lep(jets_with_lep, jd), cs_no_lep(jets_no_lep, jd);
            fjets_with_lep = cs_with_lep.inclusive_jets();
            fjets_no_lep = cs_no_lep.inclusive_jets();
            mm_mj14_lep->at(iconfig) = 0.;
            for(const auto &fjet: fjets_with_lep) mm_mj14_lep->at(iconfig) += fjet.m();
            mm_mj14_nolep->at(iconfig) = 0.;
            for(const auto &fjet: fjets_no_lep) mm_mj14_nolep->at(iconfig) += fjet.m();

            //Modify MET and mT
            float met_x = met*cos(met_phi) + (old_lep.Px()-new_lep.Px());
            float met_y = met*sin(met_phi) + (old_lep.Py()-new_lep.Py());
            mm_met->at(iconfig) = hypot(met_x, met_y);
            mm_met_phi->at(iconfig) = atan2(met_y, met_x);
            mm_mt->at(iconfig) = MT(new_lep.Pt(), new_lep.Phi(),
                                    mm_met->at(iconfig), mm_met_phi->at(iconfig));
	  }
	}
      }
    }
  }
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

