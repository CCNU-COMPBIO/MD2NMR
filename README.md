# MD2NMR
Written by TW in Jan 2023, Queen's University Panchenko's Lab.
For any enquiries please contact anna.panchenko@queensu.ca. back up email: t.wei@queensu.ca 

Subsequently modified and updated by Houfang Zhang and Yunhui Peng in July 2025, Central China Normal University. For any enquiries please contact yunhuipeng@ccnu.edu.cn

## Description
MD2NMR is a tool for calculating NMR relaxation observables (R1/R2/NOE/T1/T2/Tau_c) directly from MD trajectories. Initially written for calculations regarding nucleosome simulations but can be extended for other proteins/complexes. This software is subject to MIT license, with some functions imported from other repo, as noted individually in script comment section.

## Dependencies
The required packages are listed in 'requirements.txt'.
To create a new environment using anaconda: (replace myenv as appropriate)
conda create --name myenv --file requirements.txt

or use pip:
pip install virtualenv #create virtual environment
virtualenv test_env
source test_env/bin/activate
pip install numpy==2.2.6 pandas==2.3.3 scikit-learn==1.7.2 scipy==1.15.3 MDAnalysis==2.9.0 matplotlib==3.4.3 

## Usage
For single file mode (basic usage):
python md2nmr.py -t {$topology_file$} -y {$trajectory_file_path$}

Before usage, please check the `config.py` file and make sure the parameters are suitable for the calculation. You can change the prefix list and working directory/output directory as needed. Double-check that the magnetic field strength corresponds to the experiment result: B0 in `config.py` file.

### Data Download Protocol
To download the testing MD trajectory, use the following command to download it under the `./data` directory:

python ./tests/download_data.py


### Testing Protocol
To test MD2NMR, run the following command after downloading the testing MD trajectory:
python ./src/md2nmr.py -t H3.pdb -y H3_1.xtc

The result should be in the output directory. Then, change the `$use_chain_ids$` in `config.py` to `"use_chain_ids = False"` and run the following command to test on the other trajectory:

python3 ./src/md2nmr.py -t WT_rw_run1.pdb -y WT_rw_run1_2000ns_40ps.xtc 


### Test on your own trajectory.
Follow the steps below:
1. Check the magnetic field parameters. In the config.py the mag field unit is Tesla.
2. Check if the trajectory file is long enough. Usually an accurate calculation will need at least 200ns. In the config.py, make sure your traj length is longer than the number of split of trajectory(n_split) times the cutoff length (tau_max). e.g. for a 1000ns trajectory, recommend setting is: n_split = 50; tau_max = 1.8. See the config file for detailed explaination on parameters.
3. Check your topology file, usually this is a pdb file. Note that for different MD packages (AMBER/GROMACS/NAMD) the standard output format may be different. If your PDB file is seperated in different chains, set the change the `$use_chain_ids$` in `config.py` to `"use_chain_ids = True"` otherwise set to False (default).
4. make sure the working directory `"wd"` and the output directory `"od"` is set and you will see a prompt by the software the trajectory is loaded. 

### Batch Mode Protocol
For batch mode, check the `config.py` file and modify the `prefix_list` variable in it. Note that the batch mode will generate results for all traj/topo under working directory with satisfied prefix. Then, in the command line, run:

python ./src/md2nmr.py --batch=True
