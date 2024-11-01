# src/calculations.py

from external_programs import *
from utils import *
import os
import numpy as np



def print_H_grad(soc_der, k_so, k_u, forces, eig_val_a, H_soc_dim, step, disp_dir, dof_dir, base_folder):
    
    n_dof = len(dof_dir)
    H_grad = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=np.complex128)

    if soc_der == 'yes':
        H_grad = k_so + k_u
    else:
        H_grad = k_so

    for i in range(n_dof):
        for j in range(H_soc_dim):
            for k in range(H_soc_dim):
                H_grad[i,j,k] *= (eig_val_a[k] - eig_val_a[j])

    H_grad += forces

    out_real = os.path.join(base_folder,'grad_matrix_R.txt')
    out_imag = os.path.join(base_folder,'grad_matrix_C.txt')
    out_module = os.path.join(base_folder,'grad_matrix_module.txt')

    #Write real part
    with open(out_real,"w") as f:
        f.write(f'{n_dof}\n')
    for i in range(n_dof):
           with open(out_real,"a") as f:
            f.write(f'{i+1}\n')
            np.savetxt(f, H_grad[i].real, fmt="%.20f")
    #Write imaginary part
    with open(out_imag,"w") as f:
        f.write(f'{n_dof}\n')
    for i in range(n_dof):
           with open(out_imag,"a") as f:
            f.write(f'{i+1}\n')
            np.savetxt(f, H_grad[i].imag, fmt="%.20f")   
    #Write module
    with open(out_module,"w") as f:
        f.write(f'{n_dof}\n')
    for i in range(n_dof):
           with open(out_module,"a") as f:
            f.write(f'{i+1}\n')
            np.savetxt(f, np.abs(H_grad[i]), fmt="%.1f")    

def get_forces(H_soc_dim, disp_dir, dof_dir, base_folder,step):
    pos_en = collect_energies_array(H_soc_dim, disp_dir[0], dof_dir, base_folder)
    neg_en = collect_energies_array(H_soc_dim, disp_dir[1], dof_dir, base_folder)

    n_dof = len(dof_dir)

    forces = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=np.complex128)
    for i in range(n_dof):
        forces[i] = (pos_en[i] - neg_en[i]) / (2 * float(step))
    
    return forces



def differentiate_array(array_positive, array_negative, step):

    n_dof = (array_positive.shape)[0]
    H_soc_dim = (array_positive.shape)[1]

    diff_array = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=array_positive.dtype)
    for i in range(n_dof):
        diff_array[i] = (array_positive[i] - array_negative[i]) / (2 * float(step))

#    if isinstance(array_positive[0,0,0],np.float128):
#       diff_array = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=np.float128)
#       for i in range(n_dof):
#           diff_array[i] = (array_positive[i] - array_negative[i]) / (2 * float(disp))
#
#    elif isinstance(array_positive[0,0,0],np.complex128):
#       diff_array = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=np.complex128)
#       for i in range(n_dof):
#           diff_array[i] = (array_positive[i] - array_negative[i]) / (2 * float(disp))

    diff_array = zero_diag_of_array(diff_array)
    return diff_array

def spinfree2soc_array(array,soc_vec):

    n_dof = (array.shape)[0]
    H_soc_dim = (array.shape)[1]

    rot_array = np.zeros((n_dof, H_soc_dim, H_soc_dim), dtype=array.dtype)

    for i in range(n_dof):
        rot_array[i] = np.matmul(np.matmul(np.conjugate(soc_vec.T),array[i]),soc_vec)

    return rot_array



def extract_overlap_matrix(N_blocks,N_roots,multiplicities,path_calc):
    

    ovl_matrix = []
    for i in range(N_blocks):
        b = []
        string = 'block_'+str(i+1)
        filename = os.path.join(path_calc,string,'WFo.out')
        wf_out = open(filename).readlines()

        for line_n, line in enumerate(wf_out):
            if 'Orthonormalized overlap matrix <PsiA_i|PsiB_j>' in line:
                ov_line = line_n
                break

        for j in range(N_roots[i]):
            a=[]
            for k in range(N_roots[i]):
                a.append(float(wf_out[ov_line+2+j].split()[2+k])) 
            b.append(a)
        ovl_matrix.append(b)  
     

    ov = [] # spin-free overlap matrix n_roots*n_roots
    ov_large = [] # expanded overlap matrix over multiplicities (n_roots*mult)x(n_roots*mult)
       
    for i in range(N_blocks):
      for l in range(multiplicities[i]):  
        prev_entries = 0
        next_entries = 0
        for j in range(i):
            prev_entries += N_roots[j]*multiplicities[j]
        prev_entries += l*N_roots[i]
        for j in range(i+1,N_blocks):
            next_entries += N_roots[j]*multiplicities[j]
        next_entries += N_roots[i]*(multiplicities[i]-l-1)  
        for k in range(N_roots[i]):
            h = []
            for j in range(prev_entries):
                h.append(0.0)
            for j in range(N_roots[i]):
                h.append(ovl_matrix[i][k][j])
            for j in range(next_entries):
                h.append(0.0)
            ov_large.append(h)




    for i in range(N_blocks):
        prev_blocks = i
        for j in range(N_roots[i]):
         h = []
         for p in range(prev_blocks):
            for q in range(N_roots[p]):
             h.append(0.0)
         for k in range(N_roots[i]):
            h.append(float(ovl_matrix[i][j][k]))
         for p in range(i+1,N_blocks):
            for q in range(N_roots[p]):
              h.append(0.0) 
         ov.append(h)

    # Phase correction here 
    phase_b = []
    for i in range(len(ov)):
             if (float(ov[i][i])/abs(float(ov[i][i])))!= 1:
                 for j in range(len(ov)):
                     ov[j][i] = float(ov[j][i])*(-1) 

    for i in range(len(ov_large)):
             if (float(ov_large[i][i])/abs(float(ov_large[i][i])))!= 1:
                 phase_b.append(float(-1)) #phase: PHASE.dat
                 for j in range(len(ov_large)):
                     ov_large[j][i] = float(ov_large[j][i])*(-1)
             else:
                 phase_b.append(1.0)                 
                
                 
    string = ''
    for i in range(len(ov)):
        for j in range(len(ov)):
            string += '%s ' % (ov[i][j])
        string+='\n'    
    overlap_matrix_file=os.path.join(path_calc,'ovl_matrix.txt')
    writefile(overlap_matrix_file,string)


    string=''
    for i in range(len(ov_large)):
        for j in range(len(ov_large)):
            string+='%s ' % (ov_large[i][j])
        string+='\n'
    overlap_matrix_file=os.path.join(path_calc,'ovl_matrix_large.txt')
    writefile(overlap_matrix_file,string)


    return phase_b, ov_large


def phase_correction(eig_vec_a,eig_vec_b,eig_val_b,ov_spinfree,correction,degeneracy_thresh):

    if correction == 'yes':
        M = np.matmul(np.matmul(eig_vec_a.getH(),ov_spinfree),eig_vec_b)
    else:    
        M = np.matmul(eig_vec_a.getH(),eig_vec_b)


    # Critical step to determine degeneracies and treat complex phases
    for i in range(len(eig_val_b)):
        for j in range(len(eig_val_b)):
            if (eig_val_b[i]-eig_val_b[j]) < float(degeneracy_thresh) and (eig_val_b[i]-eig_val_b[j]) > -1*float(degeneracy_thresh):
                continue
            else:
                M[i,j] = 0 + 1j*0

    u,s,vh = np.linalg.svd(M,full_matrices=False)
    P = np.matmul(u,vh)
    eig_vec_b_phase_corrected = np.matmul(eig_vec_b,P.getH())

    return eig_vec_b_phase_corrected



def get_SOC_overlap_matrix(eig_vec_a,eig_vec_b,path_calc):

    ov_SOC_real_file=os.path.join(path_calc,'ov_matrix_SOCR.txt')
    ov_SOC_im_file=os.path.join(path_calc,'ov_matrix_SOCI.txt')

    SOC_ov  =  np.matmul(eig_vec_a.getH(),eig_vec_b)

#    string_r=''
#    string_i=''

#    for i in range(len(eig_vec_a)):
#        for j in range(len(eig_vec_a)):
#           string_r+=f'{SOC_ov[i,j].real} ' 
#           string_i+=f'{SOC_ov[i,j].imag} '
#        string_r+='\n'
#        string_i+='\n'

    for i in range(len(eig_vec_a)):
        if SOC_ov[i,i].real < 0.9:
            print(f'Warning: eigenstate {i} has overlap significantly smaller than 1. Adjust degeneracy_thresh.')

    np.savetxt(ov_SOC_real_file, SOC_ov.real, fmt='%.18e')
    np.savetxt(ov_SOC_im_file, SOC_ov.imag, fmt='%.18e')


#    writefile(ov_SOC_real_file,string_r)
#    writefile(ov_SOC_imm_file,string_i)



def get_H_soc(dim,basename,path_calc,b_mod,b_dir,phase,rescale):
    
    out=open(f'{os.path.join(path_calc,basename)}.out').readlines()


    # Get matrices from ORCA output
    H_SOC = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Lx = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Ly = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Lz = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Sx = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Sy = np.matrix(np.zeros((dim, dim), dtype = np.complex128))
    Sz = np.matrix(np.zeros((dim, dim), dtype = np.complex128))    
    
    N_rows = dim//6
    rest=dim-N_rows*6

    for iline,line in enumerate(out):
     if 'SOC MATRIX (A.U.)' in line:
         SOC=1
     if ('Real part:' in line) and SOC==1:
      for i in range(dim):
                for j in range(N_rows):
                    for k in range(6):
                          H_SOC[i,k+j*6] = float(out[iline+2+i+(dim+1)*j].replace("-", " -").split()[k+1])
                for k in range(rest):
                    H_SOC[i,k+(dim-rest)] = float(out[iline+2+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])
     if ('Image part:' in line) and SOC==1:
      for i in range(dim):
                for j in range(N_rows):
                    for k in range(6):
                          H_SOC[i,k+j*6] = H_SOC[i,k+j*6].real + 1j*float(out[iline+2+i+(dim+1)*j].replace("-", " -").split()[k+1])
                for k in range(rest):
                    H_SOC[i,k+(dim-rest)] = H_SOC[i,k+(dim-rest)].real + 1j* float(out[iline+2+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])
      SOC=0
     if 'LX MATRIX IN CI BASIS' in line:
      for i in range(dim):
               for j in range(N_rows):
                   for k in range(6):
                         Lx[i,k+j*6] = 0 + 1j*float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])
               for k in range(rest):
                   Lx[i,k+(dim-rest)] = 0 + 1j*float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])
     if 'LY MATRIX IN CI BASIS' in line:
      for i in range(dim):
              for j in range(N_rows):
                  for k in range(6):
                        Ly[i,k+j*6] = 0 + 1j*float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])
              for k in range(rest):
                  Ly[i,k+(dim-rest)] = 0 + 1j*float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])
     if 'LZ MATRIX IN CI BASIS' in line:
      for i in range(dim):
              for j in range(N_rows):
                  for k in range(6):
                        Lz[i,k+j*6] = 0 + 1j*float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])
              for k in range(rest):
                  Lz[i,k+(dim-rest)] = 0 + 1j*float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])
     if 'SX MATRIX IN CI BASIS' in line:
      for i in range(dim):
               for j in range(N_rows):
                   for k in range(6):
                         Sx[i,k+j*6] = float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])*0.5
               for k in range(rest):
                   Sx[i,k+(dim-rest)] = float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])*0.5
     if 'SY MATRIX IN CI BASIS' in line:
      for i in range(dim):
              for j in range(N_rows):
                  for k in range(6):
                        Sy[i,k+j*6] =  0 + 1j*float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])*0.5
              for k in range(rest):
                  Sy[i,k+(dim-rest)] = 0 + 1j*float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])*0.5
     if 'SZ MATRIX IN CI BASIS' in line:
      for i in range(dim):
              for j in range(N_rows):
                  for k in range(6):
                        Sz[i,k+j*6] = float(out[iline+5+i+(dim+1)*j].replace("-", " -").split()[k+1])*0.5
              for k in range(rest):
                  Sz[i,k+(dim-rest)] = float(out[iline+5+i+(N_rows)*(dim+1)].replace("-", " -").split()[k+1])*0.5



    #Phase correction
    for i in range(dim):
        for j in range(dim):
           H_SOC[i,j] = H_SOC[i,j]*phase[i]*phase[j]
           Lx[i,j] = Lx[i,j]*phase[i]*phase[j]
           Ly[i,j] = Ly[i,j]*phase[i]*phase[j]
           Lz[i,j] = Lz[i,j]*phase[i]*phase[j]



    # Build H_soc with Zeeman 
    conv_factor = 2.12719107865677e-6 
    gs = 2.002319304362

    B = float(b_mod)
    B_vec = [float(val) for val in b_dir.split()]

    H_SOC_B = H_SOC + conv_factor*B*((gs*Sx+Lx)*B_vec[0]+(gs*Sy+Ly)*B_vec[1]+(gs*Sz+Lz)*B_vec[2])

    # Diagonalise H_soc
    eig_val,eig_vec = np.linalg.eigh(H_SOC_B)
    if rescale == 'rescale':
      eig_val = (eig_val-eig_val[0])*219474.63
    else:
       eig_val *= 219474.63 

    return eig_val, eig_vec
