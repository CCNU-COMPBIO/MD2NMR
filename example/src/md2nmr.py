import numpy as np
import pandas as pd
import config
import argparse
import os
import shutil

from math import sqrt
from util import *

###### define global params ######

wd = config.wd
od = config.od

prefix_list = config.prefix_list

###### define constants ######

H_gyro = 2*np.pi*42.57748*1e6     ## Gyromagnetic Ratio: Hydrogen ([rad]/[s][T]) 
N_gyro = -2*np.pi*4.317267*1e6    ## Gyromagnetic Ratio: Nitrogen ([rad]/[s][T])
B0 = config.B0                         ## Field Strength = 18.8 for 800MHz[Teslas] 19.9 for 850MHz  14.1 catherine
vn = 4.317267 * B0 * 1e6          ##resonance frequency of 15N in Hz, note this change with magnetic field B0.

## Need 5 Frequencies: ## J[0], J[wH], J[wN], J[wH-wN], J[wH+wN]
Larmor1H = H_gyro*B0              ## Larmor Frequency: Hydrogen ([rad]/[s])
Larmor15N = N_gyro*B0             ## Larmor Frequency: Hydrogen ([rad]/[s])
omDiff = Larmor1H - Larmor15N    ## Diff in Larmor Frequencies of Spin IS
omSum  = Larmor1H + Larmor15N    ## Sum of Larmor Frequencies of Spin IS

mu_0 = 4*np.pi*1e-7    ; ## Permeability of Free Space: ([H]/[m]) X
hbar = 1.0545718e-34  ; ## Reduced Plank's constant: [J] * [s] = [kg] * [m^2] * [s^-1]

R_NH = 1.02e-10                     ## distance between N-H atoms in Angstroms
dSigmaN = -170e-6               ##  CSA of the S-spin atom N = -170 ppm

FDD = (1./4.)*np.power((mu_0*hbar*H_gyro*N_gyro)/(4*np.pi*np.power(R_NH,3)),2) 
FCSA = (1.0/3.0)*(Larmor15N**2)*(dSigmaN**2) #(2.0/15.0)*(Larmor15N**2)*(dSigmaN**2)        ## CSA factor 

###### wrap up function ######
def gen_nmr_relaxation(FTOPN, FMDN, prefix):
    NH_Res_Fname = od + prefix + '_NH_Res.csv'
    CtOutFname = od + prefix + '_NH_Ct.csv'
    dCtOutFname = od + prefix + '_NH_dCt.csv'

    CtInName = od + prefix + '_NH_Ct.csv'
    dCtInName = od + prefix + '_NH_dCt.csv'

    DF_name = od + prefix + '_Relaxtion.csv'
        
    ## Calculate the NHVecs; Can be adapted to loop over multiple trajectories using TRAJLIST_LOC
    NHVecs = []
    NHV,NHRes = calc_NHVecs(FMDN, FTOPN, config.resid) #normalized N-H vectors.
    NHVecs.append(NHV)
    #print(NHV.shape[0])
    NH_Res_df = pd.DataFrame(NHRes, columns=["NH_Res"])
    NH_Res_df.to_csv(NH_Res_Fname)
    
    dt = config.dt #we use 100 ps ## timestep of trajectories: (ps)
    print('Total frames of MD trajectory loaded: ', np.array(NHVecs).shape[1])
    n_split = config.n_split
    tau_split = np.array(NHVecs).shape[1]*dt/n_split #time (ps) of splited trajectories.

    ## Split the vecs based off the tau_split you want and the time step. 
    vecs_split = split_NHVecs(NHVecs, dt, tau_split) 
    #print('After spliting it', n_split, ' times, the N-H vector now has a shape:', vecs_split.shape)
    #print("shape of the N-H vector: [num_of_split, frames_in_each_trunk, num_of_calc_residues, xyz_coor]")
    
    ## Calculate the correlation functions and the standard deviation in the correlation function.
    ## Save the correlation functions in a dataframe and then to a csv file for later use.
    Ct, dCt = calc_Ct(vecs_split)
    CtDF = pd.DataFrame(Ct, index = np.arange(1, Ct.shape[0]+1)*dt/1000) #this is the calculated correlation function dataframe.
    dCtDF = pd.DataFrame(dCt, index = np.arange(1, dCt.shape[0]+1)*dt/1000)
    CtDF.to_csv(CtOutFname)
    dCtDF.to_csv(dCtOutFname)
    
    CtDF = pd.read_csv(CtInName, index_col=0)
    dCtDF = pd.read_csv(dCtInName, index_col=0)
    #print(CtDF.shape)
        
    parameters_list = [4] ## A, tau_A, B, tau_B, (1-A-B), tau_C
    thresh=1.0 ## default
    FitDF = fitCorrF(CtDF, dCtDF, config.tau_max, parameters_list, thresh)
    
    ## Calculate spectral density from the FitDF by calling the J_direct_transform function for each of the 5 frequencies.
    ## Loop over the rows of the FitDF dataframe from fitCorrF function and calcuate the spectral densities.
    ## Save the spectral densities to a dictionary and append to a list.
    Jarr = [] # array of dictionaries.

    for i,fit in FitDF.iterrows():
        c = fit[['C_a','C_b','C_g','C_d']].values
        t = fit[['tau_a','tau_b','tau_g','tau_d']].values
        Jdict = {'0':0, '1H':0,'15N':0,'Sum':0,'Diff':0} 
        J0 = J_direct_transform(0, c, t)
        JH = J_direct_transform(Larmor1H, c, t)
        JN = J_direct_transform(Larmor15N, c, t)
        JSum = J_direct_transform(omSum,  c, t)
        JDiff = J_direct_transform(omDiff,  c, t)
        Jdict['1H'] = JH ; Jdict['15N'] = JN; Jdict['0'] = J0; 
        Jdict['Sum'] = JSum; Jdict['Diff'] = JDiff;
        Jarr.append(Jdict)
        
    ## Calculate NMR relaxation parameters for each residue by calling calc_NMR_relax 
    ## Save the T1, T2 and NOE parameters to a dataframe
    NMRRelaxDF = pd.DataFrame(np.zeros((len(Jarr),7)),index=range(1,len(Jarr)+1), columns=['T1','T2','NOE','T1/T2','R1','R2','tau_c'])
    for index in range(1,len(Jarr)+1):
        r1, r2, noe = calc_NMR_Relax(Jarr[index-1], FDD, FCSA, H_gyro, N_gyro)
        NMRRelaxDF.loc[index,'T1'] = 1/r1; 
        NMRRelaxDF.loc[index,'T2'] = 1/r2; 
        NMRRelaxDF.loc[index,'NOE'] = noe;
        NMRRelaxDF.loc[index,'T1/T2'] = r2/r1;
        NMRRelaxDF.loc[index,'R1'] = r1; 
        NMRRelaxDF.loc[index,'R2'] = r2;
        #NMRRelaxDF.loc[index,'tau_c'] = 1/(4*np.pi*vn) * sqrt(6*(r2/r1)-7) * 1e9; #note the tau_c is in ns.
        try:
            # Attempt to calculate tau_c
            NMRRelaxDF.loc[index, 'tau_c'] = 1/(4*np.pi*vn) * np.sqrt(6*(r2/r1)-7) * 1e9
        except (ValueError, ZeroDivisionError) as e:
            # Handle mathematical errors (negative sqrt or division by zero)
            print(f"Error calculating tau_c: {e}. Reason: Expression sqrt(6*(r2/r1)-7) is invalid (r1={r1}, r2={r2})")
            NMRRelaxDF.loc[index, 'tau_c'] = 'NA'  # Set to 'NA'
        except Exception as e:
            # Catch other potential exceptions
            print(f"Unexpected error occurred while calculating tau_c: {e}")
            NMRRelaxDF.loc[index, 'tau_c'] = 'NA'
            
        
    NMRRelaxDF['Resname'] = FitDF['Resname'].values
    NMRRelaxDF['RESNUM'] = NMRRelaxDF['Resname'].str.extract('([0-9]+)',expand=False).astype('int')+1
    
    ## Merge the NMR relaxation dataframes with the FitDF dataframe
    FitRelaxDF = FitDF.merge(NMRRelaxDF, how='left', left_on='Resname',right_on='Resname').set_index(NMRRelaxDF.index)
    ## Save FitRelaxDF to a csv file
    FitRelaxDF.to_csv(DF_name)

def main(args):
    ### print out the config and args ###
    
    print("###############    Software starting now.    ###############")
    print("Input argument and config files currently using:")
    
    print(args)
    for key, value in config.__dict__.items():
        if not key.startswith("__"):
            print(f"{key}: {value}")
    
    print("###############    Calculation starting...    ###############")
    if args.batch == True:
        print('Calculation using batch style.')
        print('Following files will be used in calculation:')
        print([wd + file for file in config.prefix_list])
        
        for prefix in prefix_list:
            print('Generating NMR relaxation for: ', prefix)
            FTOPN = wd + prefix + '.pdb' #topology file
            FMDN = wd + prefix + '.xtc' #"trajectory file

            gen_nmr_relaxation(FTOPN, FMDN, prefix)
    elif args.batch == False:
        print('Calculation using single traj style.')
        
        FTOPN = wd + str(args.topo) #topology file
        FMDN = wd + str(args.traj) #"trajectory file
        
        prefix = args.traj[:-4] #define the prefix for output files.
            
        gen_nmr_relaxation(FTOPN, FMDN, prefix)
        
    
###### main ######

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = 'This is a software that calculates the NMR relaxation rate based on Molecular Dynamic simulations.')
    
    parser.add_argument('-t', 
                        '--topo', 
                        type = str, 
                        help = 'Name of topology files, note that the file type has to be supported by mdtraj. Example: -i xxx.pdb or -i xxx.prmtop')
    
    parser.add_argument('-y', 
                        '--traj', 
                        type = str, 
                        help = 'Name of MD trajectory files, note that the file type has to be supported by mdtraj. The MD trajectory has to be pre-processed from raw data, preferrably superimposed on the main body of simulation and autoimmaged given the PBC box. To do this normally the box dimension has to be rotated and only xtc file support tilted/rotated box info storage. Example: -i xxx.xtc')
    
    parser.add_argument('--batch', type = bool, 
                        help='Run in batch mode or not. If set to True, the software will process all files named "$prefix$.pdb" and "$prefix$.xtc" in the working directory. the $prefix$ can be set in config.py file.',
                        default=False)
    
    args = parser.parse_args() #pass the user input args
    
    main(args)
        
    print("###############    Calculation completed. Software exiting...    ###############")
