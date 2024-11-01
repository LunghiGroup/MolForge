In this tutorial, we will compute the T1 relaxation time for the CrN(pyrdtc)2 compound using the theory developed in Sci. Adv. 8, eabn7880 (2022) and J. Chem. Theory Comput. 20, 323−332 (2024).
Steps:

1) Run Ab Initio Calculations
All ab initio input files, including displaced geometries, are provided. Perform calculations for both positive and negative displacements in the directories ./vibronic_couplings/positive/X_Y and ./vibronic_couplings/negative/X_Y, where X is the atom label and Y is the Cartesian displacement. The equilibrium geometry calculation is in ./vibronic_couplings/0_0.


2) Compute Overlaps Using PyVC
Use the Python code pyVC to compute all overlaps. Use the input file VC.inp and set "calculation : overlap". To run the code, use:

../pyVC/scr/main.py -i VC.inp

Ensure that paths to WFOVERLAP and ORCA are specified in "wfoverlap_path" and "orca_path" in the input file. Additionally, "b_dir" and "b_mod" specify the direction and magnitude of the external magnetic field; "step" defines the displacement used for numerical differentiation, and "dof.txt" lists the degrees of freedom to be computed. For further information on input keywords, refer to pyVC/scr/utils.py.

3) Extract Vibronic Coupling Matrix Elements
After computing the overlaps, extract the vibronic coupling matrix elements by setting "calculation : H_grad" and running again 

../pyVC/scr/main.py -i VC.inp

This will generate the gradient matrices, grad_matrix_R.txt (real part) and grad_matrix_C.txt (complex part).

4) Compute Dynamics
Move to the T1 directory to compute the dynamics. In this folder, ensure you have:

-    Harmonic force constants (FC2),
-    Energies of the first N excited states (eig.dat),
-    Vibronic couplings produced in the previous step (grad_matrix_R.txt, grad_matrix_C.txt).

These files are provided to reproduce the relaxation of CrN(pyrdtc)2 with a magnetic field along g_perp and eigenvalues from the NEVPT2(1,5) level of theory. Run the dynamics executable Hdyn.x as follows:

/Hdyn.x arg1 arg2 arg3 arg4 arg5

where:

    arg1: Temperature in Kelvin,
    arg2: Gaussian smearing in cm−1,
    arg3: Number of states used to compute the vibronic couplings,
    arg4: "0" for Orbach, "1" for Raman,
    arg5: Number of degrees of freedom (3 * Number of atoms).

For this case, to simulate  Raman relaxation mechanism at 100K with smearing of 30 cm-1: 

/Hdyn.x 100 30 10 1 102


5) Extract T1 Relaxation Time
The T1 relaxation time can be identified as the inverse of the first non-zero eigenvalue of the diagonalized Limbladian Matrix, printed in the output.
