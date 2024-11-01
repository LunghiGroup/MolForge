# src/utils.py

import os
import sys
import configparser
from io import StringIO 
import numpy as np

def writefile(filename: str, content: str) -> None:
    """
    Write content to a file.

    Args:
        filename (str): The name of the file to write to.
        content (str): The content to write; must be a string.

    Raises:
        TypeError: If the content is not a string.
        IOError: If the file cannot be written to.
    """
    if not isinstance(content, str):
        raise TypeError(f"Content of type {type(content)} cannot be written to file!")

    try:
        with open(filename, 'w') as f:
            f.write(content)
    except IOError:
        print(f"Could not write to file {filename}!")
        sys.exit(14)


def init(A):
    mat=[]
    for i in range(A):
          mat.append([])
    return mat

def zero_diag_of_array(array_input):
    array_output = array_input
    dim = (array_output.shape)[0]
    if isinstance(array_output[0,0,0],np.float128):
     for i in range(dim):
        np.fill_diagonal(array_output[i], 0)
    elif isinstance(array_output[0,0,0],np.complex128):
     for i in range(dim):
        np.fill_diagonal(array_output[i], 0 + 1j*0)
        
    return array_output 


def reorder_dof(dof_dir):
    ordered_dof_dir = sorted(dof_dir, key=lambda x: (int(x.split('_')[0]), int(x.split('_')[1])))
    return ordered_dof_dir


def collect_energies_array(H_soc_dim, disp_dir, dof_dir, base_folder):

    en_array = np.zeros((len(dof_dir),H_soc_dim,H_soc_dim), dtype = np.complex128)
    for dof_number, dof in enumerate(reorder_dof(dof_dir)):
        path = os.path.join(base_folder,disp_dir,dof,'soc_en.txt')
        values = np.loadtxt(path, dtype=np.complex128)
        np.fill_diagonal(en_array[dof_number], values)

    return en_array

def collect_ovspinfree(H_soc_dim, disp_dir, dof_dir, base_folder):

    ov_spinfree = np.zeros((len(dof_dir),H_soc_dim,H_soc_dim), dtype = np.complex128)
    for dof_number, dof in enumerate(reorder_dof(dof_dir)):
        path = os.path.join(base_folder,disp_dir,dof,'ovl_matrix_large.txt')
        ov_spinfree[dof_number] = np.loadtxt(path, dtype=np.complex128)

    return ov_spinfree

def collect_ov_rot_soc(H_soc_dim, disp_dir, dof_dir, base_folder):

    ov_rot_soc = np.zeros((len(dof_dir),H_soc_dim,H_soc_dim), dtype = np.complex128)
    for dof_number, dof in enumerate(reorder_dof(dof_dir)):
        path_real = os.path.join(base_folder,disp_dir,dof,'ov_matrix_SOCR.txt')
        path_imag = os.path.join(base_folder,disp_dir,dof,'ov_matrix_SOCI.txt')
        real_part = np.loadtxt(path_real, dtype=np.complex128)
        imag_part = np.loadtxt(path_imag, dtype=np.complex128)
        ov_rot_soc[dof_number] = real_part + 1j*imag_part
    
    return ov_rot_soc

def get_info_det_from_orca_output(basenameA,basenameB,path_calc,path_eq):

    #Extract information from ORCA output. All info in file A are the same
    #as in file B except for det_block    

    outA=open(os.path.join(path_eq,f'{basenameA}.out')).readlines()

    if basenameB is not None:
      outB=open(os.path.join(path_calc,f'{basenameB}.out')).readlines()
    
    
    for line_number,line in enumerate(outA):
        if 'Determined orbital ranges:' in line:
            Aline=line_number
            break
        
    N_blocks=int(outA[Aline+8].split()[5]) # number of blocks (different multiplicities)    

    
    counter=0
    N_roots=[]
    multiplicities=[]
    tot_roots=0
    for i in range(N_blocks):
     N_roots.append(int(outA[Aline + 13 + counter].split()[2])) #number of roots per block
     multiplicities.append(int(outA[Aline + 10 + counter].split()[2])) #multiplicity of each block
     tot_roots += int(outA[Aline+ 13 +counter].split()[2])
     counter += 6 + N_roots[i]
     
    det_block_A = []
    for line_number,line in enumerate(outA):
        if 'Spin-Determinant CI Printing' in line:
            det_block_A.append(line_number+4) #starting output determinant for each block 

    if basenameB is not None: 
     det_block_B = []
     for line_number,line in enumerate(outB):
        if 'Spin-Determinant CI Printing' in line:
            det_block_B.append(line_number+4) #starting output determinant for each block    
    else:
     det_block_B = None   
 
    ndet=[]
    for i in range(N_blocks):
     fline=det_block_A[i]
     while True:
       fline+=1
       line=outA[fline]
       if line == "\n":
        break
     ndet.append(fline-det_block_A[i])  #number of determinants per block
 
    # Store information in dictionary

    calc_info = {
            'N_blocks' : N_blocks,
            'N_roots' : N_roots,
            'ndet' : ndet,
            'det_block_a' : det_block_A,
            'det_block_b' : det_block_B,
            'multiplicities' : multiplicities
            }


    return calc_info


def get_det_from_orca_output(iline,ndet,N_roots,basename,detfile,NAO,path_calc):

    out=open(f'{os.path.join(path_calc,basename)}.out').readlines()
    
    for line_number,line in enumerate(out):
        if 'Determined orbital ranges:' in line:
            Aline=line_number
            break

    #get information about active space
    N_internal=int(out[Aline+1].replace("("," ").replace("-"," ").split()[3])
 
    N_active=int(out[Aline+2].replace("("," ").replace("-"," ").split()[3])

    N_external=int(out[Aline+3].replace("("," ").replace("-"," ").split()[3])


    #fill ci_coeff with determinant information from ORCA output    
    fline=iline+ndet   
    ci_coeff=[]
    for i in range(iline,fline):
      a=[]
      occ=(out[i].split()[0]).replace("[","").replace("]","").replace("u","a").replace("d","b").replace("0","e").replace("2","d")
      if N_internal != 0:
        for j in range(0,N_internal):
         occ="d"+occ   
      if N_external != 0:
        for j in range(0,N_external):
         occ=occ+"e"   
      counter=0
      a.append(occ)
      for j in range(N_roots):
          a.append(out[i+counter].split()[1])
          counter+=ndet+3
      ci_coeff.append(a)    

    string ='%i %i %i\n' %(N_roots,NAO,ndet)
    for i in range(ndet):
      string+='%s ' % (ci_coeff[i][0])  
      for j in range(N_roots):
         string+='%s ' % (ci_coeff[i][j+1])
      string+='\n'
    writefile(detfile,string)


def get_overlap_matrix(N_blocks,N_roots,multiplicities):
    

    ovl_matrix=[]
    for i in range(N_blocks):
        b=[]
        stringg='block_'+str(i+1)
        filename=os.path.join(path_calc,stringg,'WFo.out')
        f=open(filename).readlines()
        ASline=-1
        while True:
         ASline+=1
         line=f[ASline]
         if 'Orthonormalized overlap matrix <PsiA_i|PsiB_j>' in line:
           break
        for j in range(N_roots[i]):
            a=[]
            for k in range(N_roots[i]):
                a.append(float(f[ASline+2+j].split()[2+k])) 
            b.append(a)
        ovl_matrix.append(b)  
     


    ov=[] # spin-free overlap matrix n_roots*n_roots
    ov_large=[] # expanded overlap matrix over multiplicities (n_roots*mult)x(n_roots*mult)
       
    for i in range(N_blocks):
      for l in range(multiplicities[i]):  
        prev_entries=0
        next_entries=0
        for j in range(i):
            prev_entries+=N_roots[j]*multiplicities[j]
        prev_entries+=l*N_roots[i]
        for j in range(i+1,N_blocks):
            next_entries+=N_roots[j]*multiplicities[j]
        next_entries+=N_roots[i]*(multiplicities[i]-l-1)  
        for k in range(N_roots[i]):
            h=[]
            for j in range(prev_entries):
                h.append(0.0)
            for j in range(N_roots[i]):
                h.append(ovl_matrix[i][k][j])
            for j in range(next_entries):
                h.append(0.0)
            ov_large.append(h)




    for i in range(N_blocks):
        prev_blocks=i
        for j in range(N_roots[i]):
         h=[]
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
    k=[] 
    for i in range(len(ov)):
             if (float(ov[i][i])/abs(float(ov[i][i])))!= 1:
                 for j in range(len(ov)):
                     ov[j][i]=float(ov[j][i])*(-1) 

    for i in range(len(ov_large)):
             if (float(ov_large[i][i])/abs(float(ov_large[i][i])))!= 1:
                 k.append(float(-1)) #phase: PHASE.dat
                 for j in range(len(ov_large)):
                     ov_large[j][i]=float(ov_large[j][i])*(-1)
             else:
                 k.append(1.0)

                 
                
                 
    string=''
    for i in range(len(ov)):
        for j in range(len(ov)):
            string+='%s ' % (ov[i][j])
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


    return k

######################################################## not used
def print_info(N_blocks, N_roots, N_det, multiplicities):

      tot_roots=0
      SOC_dim=0
      string='Calculation info \n \n'
      string+='Number of blocks = %s \n \n' % (N_blocks)
      for i in range(N_blocks):
          string+='Block_%s: \n' % (i+1)
          string+='Multiplicity = %s \n' % (multiplicities[i])
          string+='Number of roots = %s \n' % (N_roots[i])
          string+='Number of determinants = %s \n \n' % (N_det[i])
      for i in N_roots:
          tot_roots+=i
      string+='Size of the overlap matrix Block_1 x ... x Block_N = %i' % (tot_roots)          
      for i in range(N_blocks):
          SOC_dim+=int(N_roots[i])*int(multiplicities[i])
     
      info_file=os.path.join(path_calc,'info.txt')
      writefile(info_file,string) 
      return SOC_dim



def get_energies_SOC(filename,out_file):

    f=open(filename).readlines()

    en_vec=[]
    for i in f:
        en_vec.append(float(i))
 

    en=[]
    string=''
    for i in range(len(en_vec)):
         h=[]
         prev=i
         nex=len(en_vec)-1-i
         for j in range(prev):
           h.append(0.0)
           string+='0.0 '
         h.append(float(en_vec[i])) #! no conversion
         string+='%s ' % (float(en_vec[i])*219474.63)
         for j in range(nex):
           h.append(0.0)
           string+='0.0 '
         en.append(h)
         string += '\n'

    writefile(out_file,string)
    return en
################################################################

def read_matrix(name, path):
    
    matrix = np.matrix(np.loadtxt(os.path.join(path,name), dtype=np.float128))

    return matrix

def read_complex_matrix(name_real,name_imag, path_real,path_imag):

    matrix_real = np.matrix(np.loadtxt(os.path.join(path_real,name_real), dtype=np.float128))
  
    matrix_imag = np.matrix(np.loadtxt(os.path.join(path_imag,name_imag), dtype=np.float128))
   
    matrix = np.matrix(matrix_real, dtype = np.complex128) + 1j*np.matrix(matrix_imag, dtype = np.complex128)

    return matrix
    
def read_array(name,path):

    array = np.loadtxt(os.path.join(path,name), dtype=np.float128)

    return array 


def print_en(eig_val, path) -> None:

    np.savetxt(os.path.join(path,'soc_en.txt'), eig_val, fmt='%.18e') 



def clean(path_calc, N_blocks) -> None:
    
 if os.path.exists(os.path.join(path_calc,'mo_a')):
    os.remove(os.path.join(path_calc,'mo_a'))

 if os.path.exists(os.path.join(path_calc,'mo_b')):
    os.remove(os.path.join(path_calc,'mo_b'))

 if os.path.exists(os.path.join(path_calc,'AO_overl.mixed')):
    os.remove(os.path.join(path_calc,'AO_overl.mixed'))

 for block in range(N_blocks):
   if os.path.exists(os.path.join(path_calc,f'block_{block+1}','memlog')):
      os.remove(os.path.join(path_calc,f'block_{block+1}','memlog'))
   if os.path.exists(os.path.join(path_calc,f'block_{block+1}','det_a')):
      os.remove(os.path.join(path_calc,f'block_{block+1}','det_a'))
   if os.path.exists(os.path.join(path_calc,f'block_{block+1}','det_b')):
      os.remove(os.path.join(path_calc,f'block_{block+1}','det_b'))   

def read_input(path_calc, inp):

 # Initialize the parser
 config = configparser.ConfigParser()

 # Set default values for parameters 
 config['DEFAULT'] = {
    'calculation' : 'overlap', #overlap or H_grad 
    'basename_a' : 'geom1',
    'basename_b' : 'geom2',
    'wfoverlap_path' : '/home/lorenzo/Desktop/WFoverlap/wfoverlap-master/wfoverlap/source',
    'wf_memory' : '1000', # in MB
    'orca_path' : '/home/lorenzo/pkgs/orca',
    'molforge_path' : '/home/lorenzo/pkgs/MolForge_branch3_B/MolForge-Hdyn',
    'b_dir' : '0 0 1',
    'b_mod' : '0.0', # in Tesla
    'dof_dir' : 'list.txt',
    'disp_dir': 'positive negative',
    'step' : '0.01', # in Ang
    'eq_dir' : '0_0',
    'clean' : 'yes',
    'ov_phase_correction' : 'yes',
    'degeneracy_thresh' : '1e-5', # in cm-1
    'soc_der' : 'yes'
 }
    
 # Read input with fictitious Section key
 with open(os.path.join(path_calc,inp),'r') as file:
     content = file.read()
 content = '[DEFAULT_SECTION]\n' + content
 config.read_file(StringIO(content))

 # Prepare a dictionary to hold the input
 input_dict = {key: config.get('DEFAULT_SECTION', key) for key in config['DEFAULT_SECTION']}

 return input_dict







