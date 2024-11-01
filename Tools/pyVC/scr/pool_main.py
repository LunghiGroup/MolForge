import argparse
import os
import multiprocessing
import time
from utils import *
from compute_gradients import compute_equilibrium
from compute_gradients import compute_dof
from compute_gradients import H_grad



def run_compute_dof(args):
    disp, dof, eig_val_a, eig_vec_a, input_dict, base_folder, path_eq = args
    path = os.path.join(base_folder, disp, dof)
    compute_dof(eig_val_a, eig_vec_a, input_dict, path, path_eq)

def main(input_file):
    
    start_time = time.time()
    base_folder = os.getcwd()
    input_dict = read_input(base_folder, input_file)
    path_eq=os.path.join(base_folder,input_dict['eq_dir'])

    #H_soc at equilibrium
    eig_val_a, eig_vec_a , H_soc_dim=  compute_equilibrium(input_dict, path_eq)

    #Read dof_dir and disp_dir:
    with open(os.path.join(base_folder,input_dict['dof_dir']), 'r') as file:
      dof_dir = [line.strip() for line in file]
    disp_dir = input_dict['disp_dir'].split()

    #Set type of calculation

    if input_dict['calculation'] == 'overlap':
       nproc =  int(input_dict['nproc'])
       pool = multiprocessing.Pool(processes=nproc)

        
       args_list = [
            (disp, dof, eig_val_a, eig_vec_a, input_dict, base_folder, path_eq)
            for disp in disp_dir for dof in dof_dir
        ] 
       pool.map(run_compute_dof, args_list)

       pool.close()
       pool.join()


    elif input_dict['calculation'] == 'H_grad':          

      H_grad(eig_val_a, 
              eig_vec_a, 
              H_soc_dim, 
              input_dict['step'], 
              disp_dir, 
              dof_dir, 
              base_folder,
              input_dict['soc_der'])

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vibronic couplings")

    # Argument for input file
    parser.add_argument(
        "-i", "--input",
        help="The path to the input configuration file",
        required=True
    )
    
    args = parser.parse_args()

    # Pass the input file argument to main
    main(args.input)    
