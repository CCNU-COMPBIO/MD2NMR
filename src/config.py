#this is a config file for NMR relaxation rate calculation
#Written by Tiejun Wei Jan 2023, Subsequently modified and updated by Houfang Zhang and Yunhui Peng in July 2025
#This script and the main software is subject to the MIT liscence

### DEFINE PARAMETERS ###

""" 
1. Number of split to truncate the MD trajectory. The averaged N-H vector and its standard deviation will be used to assist scipy.curve_fit to gain a better fitting for the correlation function.
2. Tau_max in nanosecond. This is the threshold value to truncate the correlation function. Normally select 1.8ns or 6.5ns based on result from Musselman et al. 2010
3. The magnetic field strength. This is corresponds to the experiment result. Field Strength = 18.8 for 800 MHz [Teslas], 19.9 for 850 MHz
4. residue start and end index in list format. [[start_resid, end_resid], [start_resid_2, end_resid_2]]. Calculations will be based on residues in this range (including the start and end residues), excluding PRO residues since they don't have N-H on the backbone.
note this is 0-indexed resid. usually in pdb files residues are 1-indexed.
note for nucleosome we have two copies of tails, input like: [[0, 37], [136,172]]
5. Working Directory, where the trajectories and topology files are loaded # change as need. default is './data/'
6. Output Directory # change as need. default is './data/'
7. Timestep per frame (ps), note this is the actual time interval between frames in the input trajectory not the simulation step.
8. Number of exponentials; note these are used in the multi-exponential fitting of the autocorrelation function.
9. Align, note whether to remove global tumbling of the nucleosome core particle (NCP) by aligning frames to the histone fold domains.
10. Tumbling time constant in nanoseconds. If the system is a nucleosome, the default is 163.4 ns based on Rabdano et al., 2021
11. fit_log. Boolean flag indicating whether to fit the autocorrelation function using logarithmically spaced sampling points. 
    If fit_log = True, the correlation function will be sampled at log-spaced intervals instead of linear intervals.
    n_fit_log_point specifies the total number of logarithmic sampling points.
    If n_fit_log_point = None, it will be automatically determined as:
        n_fit_log_point = int(np.clip(np.sqrt(n_total), 20, 200)),
    where n_total is the total number of correlation time points available.
"""

n_split = 10 
tau_max = 50
B0 = 14.1 
resid = [[0, 37]] #[[0, 37], [136,172]]
use_chain_ids = True #True or False
chain_ids = ["E"]
wd = './data/'  
od = './data/'
dt = 40
n_exp = 3

traj_segment = [2000,50000]  #[start_snap,end_snap] or None
traj_stride = 1

use_align = True
tumbling_time = None 

fit_log = True
n_fit_log_point = None

prefix_list = ['WT_rw_run1_2000ns_40ps']#for batch mode.
