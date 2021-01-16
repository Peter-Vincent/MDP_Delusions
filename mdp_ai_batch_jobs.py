import os
import sys
import stat
import numpy as np
import subprocess
## Set the overall task parameters
parent_dir = os.getcwd()
parent_folder_name = parent_dir + '/' + sys.argv[1]
os.mkdir(parent_folder_name)
###
num_trials = 256
init_path  = "switch_sequence.mat"
num_sims   = sys.argv[2]
a_conting_var  = np.linspace(0.9,0.9,1)
e_dirich_var      = np.linspace(np.log(600),np.log(600),1)
affect_param_var  = np.linspace(0,0,1)
## Open and write to the file we will submit
submission_num = 0
num_a_conting = len(a_conting_var)
num_e_dirich  = len(e_dirich_var)
num_affect_param = len(affect_param_var)
for a_conting in a_conting_var:
    for e_dirich_exp in e_dirich_var:
        e_dirich = np.exp(e_dirich_exp)
        for affect_param in affect_param_var:                
            f = open("submission.txt","a")
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name=batch_{}\n".format(submission_num))
            f.write("#SBATCH --ntasks=1\n")
            f.write("#SBATCH --mem=8gb\n")
            f.write("#SBATCH --time=05:00:00\n")
            f.write("#SBATCH --output=output_file_{}\n".format(submission_num))
            f.write("\n")
            f.write("module load MATLAB/2016b\n")
            f.write("matlab -nodisplay -nosplash -r MDP_Delusions_harness({},'{}',{},{},{},{},'{}')".format(
                num_trials,init_path,num_sims,a_conting,e_dirich,affect_param,parent_folder_name))

            st = os.stat("submission.txt")
            os.chmod("submission.txt",st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

# os.remove("submission.txt")