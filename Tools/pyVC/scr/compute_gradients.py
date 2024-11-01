# src/compute_gradients.py

import numpy as np
from external_programs import *
from utils import *
from calculations import *

def H_grad(eig_val_a, eig_vec_a, H_soc_dim, step, disp_dir, dof_dir, base_folder,soc_der):
    ovspinfree_positive = collect_ovspinfree(
            H_soc_dim, 
            disp_dir[0], 
            dof_dir, 
            base_folder)
    ovspinfree_negative = collect_ovspinfree(
            H_soc_dim, 
            disp_dir[1], 
            dof_dir, 
            base_folder)
    ovrotsoc_positive = collect_ov_rot_soc(
            H_soc_dim, 
            disp_dir[0], 
            dof_dir, 
            base_folder)
    ovrotsoc_negative = collect_ov_rot_soc(
            H_soc_dim, 
            disp_dir[1], 
            dof_dir, 
            base_folder)

    forces = get_forces(H_soc_dim, disp_dir, dof_dir, base_folder,step)

    nacs_spinfree = differentiate_array(ovspinfree_positive,ovspinfree_negative,step)
    k_so       = spinfree2soc_array(nacs_spinfree,eig_vec_a)
    k_u        = differentiate_array(ovrotsoc_positive,ovrotsoc_negative,step)
    print_H_grad(soc_der, k_so, k_u, forces, eig_val_a, H_soc_dim, step, disp_dir, dof_dir, base_folder)

def compute_equilibrium(input_dict, path_eq):

    calc_info = get_info_det_from_orca_output(
            input_dict['basename_a'],
            None,
            None,
            path_eq)

    # Get size of H_soc matrix and phase_a
    H_soc_dim=0
    for roots, mult in zip(calc_info['N_roots'],calc_info['multiplicities']):
        H_soc_dim+=roots*mult
    phase_a=[1.0]*H_soc_dim

    # Get H_soc equilibrium
    eig_val_a, eig_vec_a =  get_H_soc(
            H_soc_dim,
            input_dict['basename_a'],
            path_eq,
            input_dict['b_mod'],
            input_dict['b_dir'],
            phase_a,
            'rescale')

    print_en(eig_val_a,path_eq)

    return eig_val_a, eig_vec_a, H_soc_dim




def compute_dof(eig_val_a, eig_vec_a, input_dict, path_calc, path_eq):

    # Get atom orbitals overlap matrix and MOs
    NAO = get_AO_MOs_from_gbw(
            input_dict['basename_a'],
            input_dict['basename_b'], 
            input_dict['orca_path'], 
            path_calc,
            path_eq)

    # Extract determinant info
    calc_info = get_info_det_from_orca_output(
            input_dict['basename_a'],
            input_dict['basename_b'], 
            path_calc,
            path_eq)
#   print_info()... to do

    # Get spin-free overlap matrices
    for block in range(calc_info['N_blocks']):
        get_spinfree_ovl(
                input_dict['basename_a'], 
                input_dict['basename_b'], 
                block, 
                NAO, 
                calc_info, 
                input_dict['wfoverlap_path'], 
                input_dict['wf_memory'], 
                path_calc,
                path_eq)
    
    phase_b, ov_spinfree = extract_overlap_matrix(
            calc_info['N_blocks'],
            calc_info['N_roots'],
            calc_info['multiplicities'],
            path_calc)


    # Clean folder    
    if input_dict['clean'] == 'yes':
        clean(path_calc,calc_info['N_blocks'])
   
    # Get size of H_soc matrix
    H_soc_dim=0
    for roots, mult in zip(calc_info['N_roots'],calc_info['multiplicities']):
        H_soc_dim+=roots*mult

    # Get H_soc
    eig_val_b, eig_vec_b =  get_H_soc(
            H_soc_dim,     
            input_dict['basename_b'],
            path_calc,
            input_dict['b_mod'],
            input_dict['b_dir'],
            phase_b,
            None)  

    print_en(eig_val_b,path_calc)

    # Phase correction on the displaced geometry
    eig_vec_b_phase_corrected = phase_correction(
            eig_vec_a, 
            eig_vec_b, 
            eig_val_b, 
            ov_spinfree, 
            input_dict['ov_phase_correction'],
            input_dict['degeneracy_thresh'])

    # Compute U^T U
    get_SOC_overlap_matrix(eig_vec_a,eig_vec_b_phase_corrected,path_calc)

