# src/external_programs.py

import os
import subprocess as sp
import shutil
import re
from utils import *







def get_AO_MOs_from_gbw(basenameA, basenameB, ORCApath, path_calc, path_eq) -> int :

    nblock = 6 # matrix format in fragovl output    

    # run orca_fragovl
    proc = sp.run(
        [f'{ORCApath}/orca_fragovl', #command
         f'{path_eq}/{basenameA}.gbw', #arg1
         f'{path_calc}/{basenameB}.gbw'],  #arg2
        cwd=path_calc, 
        capture_output=True, 
        text=True)  
    
    out_fragovl = proc.stdout.split('\n')    

    # get size of matrices
    # all have Number of Atomic Orbitals (NAO) size
    for line in reversed(out_fragovl):
      s = line.split()
      if len(s) >= 1:
        NAO = int(line.split()[0])+1 # the +1 is because ORCA counts from 0      
        break

    # find start of AO, MO_A, and MO_B matrices     
    for line_number, line in enumerate(out_fragovl):
        if 'FRAGMENT-FRAGMENT OVERLAP MATRIX' in line:
            AO_line = line_number+3
        if 'FRAGMENT A MOs MATRIX' in line:
            MOA_line = line_number+3
        if 'FRAGMENT B MOs MATRIX' in line:
            MOB_line = line_number+3
            break
        
    #initialize matrices
    #AO overlap
    ao_ovl=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]

    #Molecular Orbitals fragment A coefficients    
    MO_A=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]

    #Molecular Orbitals fragment B coefficients 
    MO_B=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]    


    #Fill the matrices from the fragovl output. 
    #The column is the same for each matrix.
    #The row depends on the repective initial line.
    #The matrices are filled as the transpose form of the matrix in the
    # output. For MO is ok, but the AO are extracted in transpose format.
    #float_pattern is used to avoid formatting issues from ORCA
    float_pattern =  r'[-+]?\d*\.\d{12}|\d+'

    #todo: "def fix_malformed_entry(entry):" otherwise the use of findall for all the entries slows down a lot the calculation

    for x in range(NAO):
      for y in range(NAO):
        block = x//nblock
        col = x%nblock+1
        row_AO = block*(NAO+1)+y+AO_line
        row_MOA = block*(NAO+1)+y+MOA_line
        row_MOB = block*(NAO+1)+y+MOB_line
        ao_ovl[y][x] = float(re.findall(float_pattern, out_fragovl[row_AO])[col])
        MO_A[x][y]   = float(re.findall(float_pattern, out_fragovl[row_MOA])[col])
        MO_B[x][y]   = float(re.findall(float_pattern, out_fragovl[row_MOB])[col])
#        ao_ovl[y][x] = float( out_fragovl[row_AO].replace("-"," -").split()[col])
#        MO_A[x][y] = float( out_fragovl[row_MOA].replace("-"," -").split()[col])
#        MO_B[x][y] = float( out_fragovl[row_MOB].replace("-"," -").split()[col])
        
   
   


    # make mo_a, mo_b in MOLCAS format for WFOverlap
    #MO_A
    string_MOA = f'''2mocoef
header
 1
MO-coefficients from orca_fragovl
 1
 {NAO}   {NAO}
 a
mocoef
(*)
''' 
    x = 0
    for imo,mo in enumerate(MO_A):
        for i in mo:
            if x >= 3:
                string_MOA += '\n'
                x = 0
            string_MOA += f'{i: 6.12e} ' 
            x += 1
        if x > 0:
            string_MOA += '\n'
            x = 0

    string_MOA += 'orbocc\n(*)\n'
    x = 0
    for i in range(NAO):
        if x >= 3:
            string_MOA += '\n'
            x = 0
        string_MOA += f'{0.0: 6.12e} '
        x += 1
    
        
    #MO_B
    string_MOB = f'''2mocoef
header
 1
MO-coefficients from orca_fragovl
 1
 {NAO}   {NAO}
 a
mocoef
(*)
'''
    x=0
    for imo,mo in enumerate(MO_B):
        for i in mo:
            if x >= 3:
                string_MOB += '\n'
                x = 0
            string_MOB += f'{i: 6.12e} ' 
            x += 1
        if x > 0:
            string_MOB += '\n'
            x = 0

    string_MOB += 'orbocc\n(*)\n'
    x = 0
    for i in range(NAO):
        if x >= 3:
            string_MOB += '\n'
            x = 0
        string_MOB += f'{0.0: 6.12e} '
        x += 1

    # make AO mixed overlap file for WFOverlap
    string_AO = f'{NAO} {NAO}\n' 
    for i in range(0,NAO):
      for j in range(0,NAO):
          string_AO += f'{ao_ovl[i][j]: .15e} '          
      string_AO += '\n'

    writefile(os.path.join(path_calc,'AO_overl.mixed'),string_AO)
    writefile(os.path.join(path_calc,'mo_a'),string_MOA)
    writefile(os.path.join(path_calc,'mo_b'),string_MOB)

    return NAO
    
def get_SOC_coefficients_MolForge(l_max, SOC_dir, filename, N_roots, multiplicity,k):

    SOC_dim=0
    for i in range(len(N_roots)):
      SOC_dim+=int(N_roots[i])*int(multiplicity[i])
    os.system('if [ -d "%s" ]; then rm -Rf %s; fi' % (SOC_dir,SOC_dir))
    os.system('mkdir %s' % (SOC_dir))
    os.system('cp %s %s' % (filename, SOC_dir))
    os.system('cp ../B.dat %s' % (SOC_dir))
    os.chdir(SOC_dir)
    string=''
    file_k=os.path.join(SOC_dir,'PHASE.dat')
    for i in range(SOC_dim):
        string+='%s \n' % (k[i])
    writefile(file_k,string)    
    field=0.000000000000000000000000000000000000001
    os.system('%s/bin/get_SH_Orca.x %s %s %s %s >> SH.out' %  (MolForge_path, filename, JMULT, l_max, field))
    os.chdir(path_calc)
    real_file=os.path.join(path_calc,SOC_dir,'SOC_EIGVEC_R.dat')
    imm_file=os.path.join(path_calc,SOC_dir,'SOC_EIGVEC_C.dat')
    real=[]
    for i in range(SOC_dim):
          real.append([])
    imm=[]
    for i in range(SOC_dim):
          imm.append([])
    SOC_coeff=[]
    for i in range(SOC_dim):
          SOC_coeff.append([])
    f=open(real_file).readlines()
    g=open(imm_file).readlines()
    for j in range(SOC_dim):
        for i in range(SOC_dim):
                SOC_coeff[j].append(complex(float(f[j].split()[i]),float(g[j].split()[i])))      

    return SOC_coeff



def get_spinfree_ovl(basenameA, basenameB, block, NAO, calc_info, wfoverlap_path, wf_memory, path_calc, path_eq):

    #block directory
    block_dir=os.path.join(path_calc,'block_'+str(block+1))

    #remove if exits and create block_dir
    if os.path.exists(block_dir):
        shutil.rmtree(block_dir)
    os.makedirs(block_dir)

    #extract determinant information for fragment A and B in current block
    get_det_from_orca_output(
            calc_info['det_block_a'][block],
            calc_info['ndet'][block],
            calc_info['N_roots'][block],
            basenameA,
            os.path.join(block_dir,'det_a'),
            NAO,
            path_eq)

    get_det_from_orca_output(
            calc_info['det_block_b'][block],
            calc_info['ndet'][block],
            calc_info['N_roots'][block],
            basenameB,
            os.path.join(block_dir,'det_b'),
            NAO,
            path_calc)

    #Prepare WFoverlap input   
    with open(f'{block_dir}/WFo.in', 'w') as input:
      input.write(f'a_mo={path_calc}/mo_a \n'
                f'b_mo={path_calc}/mo_b \n'
                f'a_det={block_dir}/det_a \n'
                f'b_det={block_dir}/det_b \n'
                f'mix_aoovl={path_calc}/AO_overl.mixed')
    #Run WFoverlap  
    proc = sp.run(
         [f'{wfoverlap_path}/wfoverlap_ascii.x', #command
          '-m',f'{wf_memory}','-f',f'{block_dir}/WFo.in'], #arg
         stdout=open(f'{block_dir}/WFo.out','w'),
         stderr=open(f'{block_dir}/WFo.error','w'),
         cwd=block_dir 
         )



