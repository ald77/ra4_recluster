#! /usr/bin/env python

import argparse
import os
import subprocess

def FullPath(path):
    return os.path.abspath(os.path.expanduser(path))

def Recluster(file_dir):
    full_path = FullPath(file_dir)
    if os.path.isdir(full_path):
        contents = os.listdir(full_path)
        for fd in contents:
            Recluster(os.path.join(full_path, fd))
    elif os.path.isfile(full_path):
        base_path = "/net/cms2/cms2r0/babymaker/babies"
        out_dir = os.path.join(base_path, "reclustered")
        rel_path = os.path.relpath(full_path, base_path)
        out_path = os.path.join(out_dir, rel_path)
        if not os.path.exists(os.path.dirname(out_path)):
            os.makedirs(os.path.dirname(out_path))
        if not os.path.exists(out_path):
            subprocess.call(["JobSubmit.csh","./run/wrapper.sh","./run/recluster_baby.exe","-i",full_path,"-o",out_path])
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reclusters jets for set of files",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", nargs="*", default=[], metavar="INPUT_DIRS", help="List of files and directories to recluster")
    args = parser.parse_args()

    for file_dir in args.input:
        Recluster(file_dir)
