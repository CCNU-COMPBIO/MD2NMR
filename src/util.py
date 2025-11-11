import numpy as np
import pandas as pd
import MDAnalysis as mda
import config
import os
import tempfile
import warnings

from MDAnalysis.analysis import align,dssp
from scipy.optimize import curve_fit
from math import sqrt
from tqdm import tqdm



def detect_resid_offset(u):
    """Detect residue ID offset."""
    min_resid = min(r.resid for r in u.residues)
    return 0 if min_resid == 0 else 1


def gen_res_selection_MDA(u, resid, element, use_chain_ids=False, chain_ids=None):
    """Generate MDAnalysis atom selection string."""
    resid_offset = detect_resid_offset(u)
    res_range_list = []
    for start, end in resid:
        res_range_list.append(f"(resid {start + resid_offset}:{end + resid_offset})")
    res_range_sel = " or ".join(res_range_list)
    
    if use_chain_ids and chain_ids:
        chain_selection = " or ".join(f"segid {cid}" for cid in chain_ids)
        sel = f"({chain_selection}) and (name {element} and not resname PRO) and ({res_range_sel})"
    else:
        sel = f"name {element} and not resname PRO and ({res_range_sel})"
    return sel


def calc_NHVecs(traj_file, top_file, resid, start_snap=0, end_snap=None):
    """
    Uses MDAnalysis to load the trajectory and get the atomic indices and coordinates to calculate the correlation functions.
    For each, trajectory load the trajectory using MDAnalysis, get the atomic index for the the N-H atoms and calculate the vector between the two.
    Append the vector to the NHVecs list for all the trajectories. 
    NHVecs should return a list of shape: (# Trajectories, # Snapshots, # Residues w/N-H Vectors, 3)
    resid is a list of list. providing the residue start/end id.
    
    """

    # Determine frame range to process
    if config.traj_segment is None:
        start_frame = start_snap
        end_frame = end_snap
    else:
        start_frame = (
            config.traj_segment[0]
            if (len(config.traj_segment) > 0 and config.traj_segment[0] is not None)
            else start_snap
        )
        end_frame = (
            config.traj_segment[1]
            if (len(config.traj_segment) > 1 and config.traj_segment[1] is not None)
            else end_snap
        )

    stride = getattr(config, "traj_stride", 1)

    warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")
    warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.analysis.base")
    u = mda.Universe(top_file, traj_file)
    total_frames = len(u.trajectory)
    end_frame_actual = end_frame if end_frame is not None else total_frames
    end_frame_actual = min(end_frame_actual, total_frames)
    frame_range = slice(start_frame, end_frame_actual, stride)
    n_frames = len(range(start_frame, end_frame_actual, stride))

    print(
        f"Total frames in trajectory: {total_frames}"
    )
    print(
        f"Processing frames from {start_frame} to {end_frame_actual - 1} with stride={stride} ({n_frames} frames)."
    )

    if config.use_align:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
            temp_pdb_path = tmp_pdb.name
            u.atoms.write(temp_pdb_path)
        u.trajectory[start_frame]  # Load first frame in range as reference
        #ref_u = mda.Universe(top_file)
        ref_u = mda.Universe(temp_pdb_path)
        ref_u.atoms.positions = u.atoms.positions.copy()


        # Select atoms for alignment based on secondary structure
        ss = dssp.DSSP(ref_u).run().results.dssp[0]
        ss_resids = [res.resid for res, s in zip(ref_u.residues, ss) if s in ['H','G','I','E','B']]
        ref_atoms = ref_u.select_atoms("resid " + " ".join(map(str, ss_resids)) + " and name CA")
        mobile_atoms = u.select_atoms("resid " + " ".join(map(str, ss_resids)) + " and name CA")

        for ts in tqdm(u.trajectory[frame_range], total=n_frames, desc="Aligning by secondary structure"):
            R, t = align.rotation_matrix(mobile_atoms.positions, ref_atoms.positions)
            mobile_atoms.positions = (
                (mobile_atoms.positions - np.mean(mobile_atoms.positions, axis=0)) @ R.T
                + np.mean(ref_atoms.positions, axis=0)
            )
        print("DSSP-based alignment completed.")
        if os.path.exists(temp_pdb_path):
            os.remove(temp_pdb_path)
            print(f"Temporary file {temp_pdb_path} deleted.")
            
    N_sel = gen_res_selection_MDA(u, resid, "N", config.use_chain_ids, config.chain_ids)
    H_sel = gen_res_selection_MDA(u, resid, "H*", config.use_chain_ids, config.chain_ids)
    Nitrogen = u.select_atoms(N_sel)
    Hydrogen = u.select_atoms(H_sel)

    #print(f"N selection: {N_sel} --> {len(Nitrogen)} atoms")
    #print(f"H selection: {H_sel} --> {len(Hydrogen)} atoms")

    NH_pairs, NH_Res = [], []
    for N_atom in Nitrogen:
        H_candidates = [H for H in Hydrogen if H.resid == N_atom.resid and H.segid == N_atom.segid]
        if H_candidates:
            NH_pairs.append((N_atom, H_candidates[0]))
            NH_Res.append(f"{N_atom.resname}{N_atom.resid}")

    if not NH_pairs:
        raise ValueError("No N-H pairs found! Check residue selection and topology.")

    NHVecs_tmp = np.zeros((n_frames, len(NH_pairs), 3), dtype=np.float32)
    for i, ts in enumerate(
        tqdm(u.trajectory[frame_range], total=n_frames, desc="Computing N-H vectors")
    ):
        vecs = np.array([H.position - N.position for N, H in NH_pairs])
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        NHVecs_tmp[i] = vecs / norms

    return NHVecs_tmp, NH_Res

    
def split_NHVecs(nhvecs, dt, tau):

    nFiles = len(nhvecs) ## number of trajectories
    nFramesPerChunk = int(tau/dt) ###tau/timestep 
    used_frames = np.zeros(nFiles,dtype=int)
    remainingFrames = np.zeros(nFiles,dtype=int)
    for i in range(nFiles):
        nFrames = nhvecs[i].shape[0]
        used_frames[i] = int(nFrames/nFramesPerChunk)*nFramesPerChunk
        remainingFrames[i] = nFrames % nFramesPerChunk
    
    nFramesTot=int(used_frames.sum())
    out = np.zeros((nFramesTot,nhvecs[0].shape[1],nhvecs[0].shape[2]), dtype=nhvecs[0].dtype)
    start = 0
    for i in range(nFiles):
        end = int(start+used_frames[i])
        endv = int(used_frames[i])
        out[start:end,...] = nhvecs[i][0:endv,...]
        start = end
        
    sh = out.shape
    vecs = out.reshape((int(nFramesTot/nFramesPerChunk), nFramesPerChunk, sh[-2], sh[-1]))
    print('vec shape after split:', vecs.shape)
    
    return vecs

def calc_Ct(nhvecs):

    sh = nhvecs.shape
    nReplicates=sh[0] ; nDeltas=int(sh[1]/2) ; nResidues=sh[2]
    Ct  = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    dCt = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    
    for delta in range(1,1+nDeltas):
        nVals=sh[1]-delta
        # = = Create < vi.v'i > with dimensions (nRep, nFr, nRes, 3) -> (nRep, nFr, nRes) -> ( nRep, nRes ), then average across replicates with SEM.
        tmp = -0.5 + 1.5 * np.square( np.einsum( 'ijkl,ijkl->ijk', nhvecs[:,:-delta,...] , nhvecs[:,delta:,...] ) )
        tmp  = np.einsum( 'ijk->ik', tmp ) / nVals
        Ct[delta-1]  = np.mean( tmp, axis=0 )
        dCt[delta-1] = np.std( tmp, axis=0 ) / ( np.sqrt(nReplicates) - 1.0 ) #if nReplicates is 1, NH_dCt will not be defined (because variance cannot be calculated for a single sample).
    
    return Ct, dCt

def _bound_check(func, params):
    """
    
    Checks if the fit returns a sum of the amplitudes greater than 1.
    
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    if len(params) == 1:
        return False
    elif len(params) %2 == 0 :
        s = sum(params[0::2])
        return (s>1)
    else:
        s = params[0]+sum(params[1::2])
        return (s>1)

def calc_chi(y1, y2, dy=[]):
    """
    Calculates the chi^2 difference between the predicted model and the actual data. 
    
    LICENSE INFO:
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    #if dy != []:
    if dy is not None and len(dy) > 0:
        return np.sum( (y1-y2)**2.0/dy )/len(y1)
    else:
        return np.sum( (y1-y2)**2.0 )/len(y1)

def calc_r2(y1, y2):
    ss_res = np.sum((y1 - y2)**2)
    ss_tot = np.sum((y1 - np.mean(y1))**2)
    return 1 - (ss_res/ss_tot)
def func_exp_decay(t, *params):
    """

    params = [A1, tau1, A2, tau2, ..., tau_n]
    """
    n_exp = config.n_exp
    if n_exp < 2 or n_exp > 4:
        raise ValueError("Only 2, 3, or 4 exponential components are supported.")

    # First (n-1) amplitudes and their corresponding time constants
    A_list = list(params[0:2*(n_exp-1):2])
    tau_list = list(params[1:2*(n_exp-1):2])
    tau_last = params[-1]

    # The last amplitude is determined to ensure normalization
    A_last = 1.0 - sum(A_list)

    # Calculate exponential decay terms
    exp_terms = [A_list[i] * np.exp(-t / tau_list[i]) for i in range(n_exp-1)]
    exp_terms.append(A_last * np.exp(-t / tau_last))
    result = np.sum(exp_terms, axis=0)
    if getattr(config, "use_align", False):
        tau_ref = config.tumbling_time if config.tumbling_time is not None else 163.4
        result = result / np.exp(-t / tau_ref)
    return result
  
def exponential_fitting(x, y, dy=np.empty([]), tau_mem=50.0):
    """
    Fitting using exponential decay functions using scipy curve_fit.
    Return the calculated chi-squared error with calc_chi function.
    Return the calculated r-squared error with calc_r2 function.
    """
    n_exp = config.n_exp 
    if n_exp < 2 or n_exp > 4:
        raise ValueError("Only 2, 3, or 4 exponential components are supported.")

    # Determine tumbling time constant
    tumbling_time = config.tumbling_time if config.tumbling_time is not None else 163.4

    # Initial guesses: equal amplitudes and logarithmically spaced time constants
    b1_guess = [1.0 / n_exp] * n_exp
    t1_guess = np.logspace(np.log10(0.001), np.log10(round(tumbling_time, 0)), n_exp)

    # Construct parameter guess and bounds dynamically
    # params = [A1, tau1, A2, tau2, ..., tau_(n-1), tau_n]
    guess = []
    lower_bounds = []
    upper_bounds = []

    for i in range(n_exp - 1):
        # Amplitude A_i
        guess.append(b1_guess[i])
        lower_bounds.append(0.0)
        upper_bounds.append(1.0)
        # Time constant tau_i
        guess.append(t1_guess[i])
        lower_bounds.append(0.0)
        upper_bounds.append(tumbling_time if config.use_align else np.inf)

    # Final time constant tau_n
    guess.append(t1_guess[-1])
    lower_bounds.append(0.0)
    upper_bounds.append(tumbling_time if config.use_align else np.inf)

    guess = tuple(guess)
    bounds = (tuple(lower_bounds), tuple(upper_bounds))

    func = func_exp_decay

    # === Fitting ===
    if dy is not None and len(dy) > 0:
        popt, popv = curve_fit(func, x, y, p0=guess, sigma=dy, bounds=bounds, method='trf', loss='soft_l1') 
    else:
        popt, popv = curve_fit(func, x, y, p0=guess, bounds=bounds, loss='soft_l1')
    
    #calculate the approximated y using the fitted functions.
    #note this is for error calculations.
    ymodel=[ func(x[i], *popt) for i in range(len(x)) ] 
    
    bExceed=_bound_check(func, popt)
    if bExceed: #disable the sum(weight) = 1, use the global scaling, S^2 ???
        #print("!!! WARNING, curve fitting in do_LSstyle_fit returns a sum>1. !!!")
        #return 9999.99, popt, np.sqrt(np.diag(popv)), ymodel 
        return calc_chi(y, ymodel, dy), calc_r2(y, ymodel), popt, popv, ymodel
    else:
        return calc_chi(y, ymodel, dy), calc_r2(y, ymodel), popt, popv, ymodel



def fitCorrF(CorrDF, dCorrDF, tau_mem, pars_l, threshold=1.0):
    """
        Input Variables:
            CorrDF: Dataframe containing the correlation functions. Columns are the NH-bond vectors, rows are timesteps. 
            dCorrDF: Error in the correlation function at time t
            tau_mem: Cut-Off time to remove noise at the tail of the correlation function 
            pars_l : parameters list. 
            
        Main function to fit the correlation function. 
        Loops over all residues with N-H vectors and calculates the fit, appends the best fit from findbest_Expstyle_fits2.
        Passes the set of lists to fitstoDF to return a data frame of the best fits for each residue. 
        
        Takes the correlation function CorrDF and errors in the correlation function, maximum tau mem to cut correlation
    """
    NH_Res = CorrDF.columns
    chi_list=[]
    r2_list = []
    names_list=[]
    pars_list=[]
    errs_list=[]
    ymodel_list=[]
    covarMat_list = []
    
    #loop on all residues
    #loop on all residues
    for i in CorrDF.columns:
               
        tstop = np.where(CorrDF.index.values==tau_mem)[0][0]
        time_values = CorrDF.index.values[:tstop]
        n_total = len(time_values)
        #print(n_total)
        y_values = CorrDF[i].values[:tstop]
        dy_values = dCorrDF[i].values[:tstop]
        if np.all(np.isnan(dy_values)):
            dy_values = []
        if config.n_fit_log_point == None:
             n_lag_points = int(np.clip(np.sqrt(n_total), 20, 200))
        else: 
             n_lag_points = config.n_fit_log_point
        #print(n_lag_points)
        if getattr(config, "fit_log", False):  # if config.fit_log == True
            lag_index = np.unique(
                np.logspace(0, np.log10(n_total - 1), n_lag_points, endpoint=False).astype(int)
            )
        else:
            lag_index = np.arange(n_total)
        #print(lag_index)

        x = time_values[lag_index]
        y = y_values[lag_index]
        dy = dy_values[lag_index] if len(dy_values) > 0 else []

        chi, r2, pars, covarMat, ymodel = exponential_fitting(x, y, dy, tau_mem)
        
        A_tau_pars = pars
        
        names = ['C_a', 'tau_a', 'C_b', 'tau_b', 'C_g', 'tau_g', 'C_d', 'tau_d']
        errs = np.sqrt(np.diag(covarMat)) # this is the standard covariance error.
        
        #append the calculated params to list
        chi_list.append(chi)
        r2_list.append(r2)
        names_list.append(names)
        pars_list.append(pars)
        errs_list.append(errs)
        ymodel_list.append(ymodel)
        
    FitDF = fitstoDF(NH_Res, chi_list, r2_list, pars_list, errs_list, names_list)
    
    return FitDF

def fitstoDF(resnames, chi_list, r2_list, pars_list, errs_list, names_list):
    ## Set Up columns indices and names for the data frame
    """
    Function that takes the residue names, chi^2, r2, parameters, errors and names of the fits and returns a data frame
    of the parameters.
    """
    mparnames = ['C_a', 'tau_a', 'C_b', 'tau_b', 'C_g', 'tau_g', 'C_d', 'tau_d']
    mtau_names = np.array(mparnames)[1::2]
    mc_names = np.array(mparnames)[::2]
    colnames = np.array(['Resname','NumExp'])
    tau_errnames = np.array([[c,"{}_err".format(c)] for c in mtau_names]).flatten()
    mc_errnames = np.array([[c, "{}_err".format(c)] for c in mc_names]).flatten()
    colnames = np.hstack([colnames,mc_errnames])
    colnames = np.hstack([colnames,tau_errnames])
    colnames = np.hstack([colnames,np.array(['Chi_Fit'])])
    colnames = np.hstack([colnames,np.array(['R2_Fit'])])
    FitDF = pd.DataFrame(0.0, index=np.arange(len(pars_list)), columns=colnames, dtype=float)
    FitDF['Resname'] = resnames
    FitDF['Chi_Fit'] = chi_list
    FitDF['R2_Fit'] = r2_list
            
    
    for i in range(len(pars_list)):
        npar = len(pars_list[i])
        if (npar%2)==1:
            ccut = npar-2
            tau_f, terr = pars_list[i][1:ccut+1:2], errs_list[i][1:ccut+1:2]
            tau_f = np.hstack([tau_f, pars_list[i][-1]])
            terr = np.hstack([terr, errs_list[i][-1]])
            sort_tau = np.argsort(tau_f)
            coeff, cerr= pars_list[i][0:ccut:2], errs_list[i][0:ccut:2]
            Clast = 1; Clasterr = 0.0;
            for n,m in zip(coeff, cerr):
                Clast -= n
                Clasterr += m
            
            coeff = np.hstack([coeff, np.array(Clast)])
            cerr = np.hstack([cerr, np.array(Clasterr)])
    
            tne = np.array([[c,"{}_err".format(c)] for c in mparnames[1:npar+1:2]]).flatten()
            cne = np.array([[c, "{}_err".format(c)] for c in mparnames[0:npar:2]]).flatten()
                
        else:
            tau_f, terr = pars_list[i][1::2], errs_list[i][1::2] 
            coeff, cerr= pars_list[i][0::2], errs_list[i][0::2]
            sort_tau = np.argsort(tau_f)[::-1]
            tne = np.array([[c,"{}_err".format(c)] for c in names_list[i][1::2]]).flatten()
            cne = np.array([[c, "{}_err".format(c)] for c in names_list[i][0::2]]).flatten()
    
        NumExp=np.array(len(tau_f))
        tau_err = np.array([[t,e] for t,e in zip(tau_f[sort_tau],terr[sort_tau])]).flatten()
        c_err = np.array([[c,e] for c,e in zip(coeff[sort_tau], cerr[sort_tau])]).flatten()
        namesarr = np.hstack([np.array('NumExp'),cne,tne])
        valarr = np.hstack([NumExp,c_err,tau_err])
    
        FitDF.loc[i,namesarr] = valarr
        
    FitDF['AUC_a'] = FitDF.C_a*FitDF.tau_a; FitDF['AUC_b'] = FitDF.C_b*FitDF.tau_b; 
    FitDF['AUC_g'] = FitDF.C_g*FitDF.tau_g; FitDF['AUC_d'] = FitDF.C_d*FitDF.tau_d;
    FitDF['AUC_Total'] = FitDF[['AUC_a','AUC_b','AUC_g','AUC_d']].sum(axis=1)
    
    
    return FitDF            
            
            
def J_direct_transform(om, consts, taus):
    
    """
        Calculation of the spectral density from the parameters of the fit by direct fourier transform
        
        return the sum of estimated spectral density on point: x = om. 
        
        MIT License

    Copyright (c) 2019 Alan Hicks

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
        
    """
    
    ndecay=len(consts)-1 ; noms=1;###lnden(om)
    Jmat = np.zeros( (ndecay, noms ) ) #Jmat is in shape [2,1] 2 is the number of fitted decay func
    for i in range(ndecay):
        Jmat[i] = consts[i]*(taus[i]*1e-9)/(1 + np.power((taus[i]*1e-9)*(om),2.))
    
    return (0.89)*Jmat.sum(axis=0)

def calc_NMR_Relax(J, fdd, fcsa, gammaH, gammaN):
    """
        Function to calculate the R1, R2 and NOE from the spectral densities and the physical parameters for the 
        dipole-dipole and csa contributions, fdd and fcsa. 
    """
    R1 = fdd * (J['Diff'] + 3*J['15N'] + 6*J['Sum']) + fcsa * J['15N']
    
    R2 = (1./2. * fdd * (4*J['0'] + J['Diff'] + 3*J['15N'] + 6*J['1H']) 
          + (1./6.) * fcsa*(4*J['0'] + 3*J['15N']) )
    
    NOE = 1 + ((fdd*gammaH)/(gammaN*R1))*(6*J['Sum'] - J['Diff']) #changed constant NOE = 1 + ((fdd*gammaH)/(*gammaN*R1))*(6*J['Sum'] - J['Diff'])
    
    return R1, R2, NOE
