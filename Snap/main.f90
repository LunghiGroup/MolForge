        program main
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use potential_class
        use SNAP_fit_class
        use coul_fit_class
        use VdW_class
        use SNAP_FF_class
        use Coul_FF_class
        use trajectory_class
        use minimizer_class
        use md_class

        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(SNAP_fit)                            :: SNAP
        type(coul_fit)                            :: COUL
        type(SNAP_FF),target                      :: FF_SNAP
        type(VdW_FF),target                       :: FF_VdW
        type(COUL_FF)                             :: FF_Coul
        type(MD)                                  :: Dinamica
        type(minimizer)                           :: Minimizzatore
        type(lammps_obj)                          :: frame
        

        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=120)                        :: geometry_file,energy_file,dipoles_file,shift_file,&
                                                     frame_file,forces_file,charges_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        real(kind=dbl)                            :: lambda_en,cutoff_en
        real(kind=dbl)                            :: lambda_dip,cutoff_dip
        integer                                   :: SNAP_twojmax_dip
        real(kind=dbl)                            :: SNAP_cutoff_dip

        logical                                   :: energy_flag,md_flag,VdW_flag,coulomb_flag,minim_flag,train,&
                                                     phonon_flag,train_SNAP,train_COUL
        logical                                   :: flag_forces,flag_energy
        double precision,dimension(:),allocatable :: energies,coul_energy,snap_energy
        double precision,dimension(:),allocatable :: forces,gradients
        double precision                          :: R_screen
        double precision,dimension(:),allocatable :: tot_charge
        character(len=5)                          :: set_type_en,set_type_dip
        character(len=10)                         :: keyword_din,keyword_min
        double precision,allocatable              :: force(:,:)
        double precision                          :: E_plus,E_minus,ave_atom
        integer                                   :: tot_atom
        integer                                   :: atom,direction,k
        double precision,dimension(:),allocatable :: coul_energy_plus,coul_energy_minus
        double precision                          :: temperature,timestep
        real(kind=dbl)                            :: T_in
        real(kind=dbl)                            :: delta
        logical                                   :: active_learning
        integer                                   :: tot_steps_md
        double precision,allocatable              :: vec(:),grad(:),grad_SNAP(:),grad_VdW(:)
        integer                                   :: idist
        integer,dimension(4)                      :: iseed
        real(kind=dbl)                            :: weight
        real(kind=dbl)                            :: error
        real(kind=dbl)                            :: val

        train_SNAP=.true.
        train_COUL=.false.
        
        !the Coulomb flag below is related to switching on the Coulomb potential for the SNAP fitting of energies and forces
        coulomb_flag=.true.
        VdW_flag=.true.

        md_flag=.false.
        minim_flag=.false.
        active_learning=.false.
        
        !lines to set the parameters for the SNAP fit of energies and forces
        SNAP%nconfig=19
        SNAP%cutoff=5.0d0
        SNAP%twojmax=11
        SNAP%weight=dsqrt(64.0d0*3.0)
        SNAP%lambda=0.1d0
        SNAP%set_type='TRAIN'
        SNAP%flag_energy=.true.
        SNAP%flag_forces=.true.

        !lines to set the parameters for the Coulomb FF for the SNAP fitting of energies and forces
        
        SNAP%FF_Coul%num_bisp=56 !THIS LINE IS REDUNDANT WITH THE ONE TWO LINES BELOW  
        SNAP%FF_Coul%tot_kinds=3 !THIS SHOULD NOT BE SET FROM HERE BUT OBTAINED IN THE SUBROUTINE THAT ADD OR SUBTRACT THE ENERGIES
        !AND FORCES
        SNAP_twojmax_dip=8
        SNAP_cutoff_dip=4.0d0
        SNAP%FF_Coul%r_screen=5.0d0*A_to_B

        !lines to set the parameters for the LAMMPS frame which will be the starting point for the dynamics or minimization
        frame%cutoff_en=4.0d0
        frame%twojmax_en=8
        frame%cutoff_dip=4.0d0
        frame%twojmax_dip=8

        !lines to tell the programs where are the files relevant to the simulation
        geometry_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_tr_AL_++"
        SNAP%energy_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/ener_tr_AL_++"
        SNAP%forces_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/grad_tr_AL_++"
        shift_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/shift_tr_AL_++"
        frame_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_start"
        COUL%dipoles_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/dipoles_tr_AL_++"
        COUL%charges_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/charges_tr_AL_++"

        !lines to set the parameters for the SNAP fit of dipoles  
        COUL%nconfig=19
        COUL%cutoff=4.0d0
        COUL%twojmax=8
        COUL%weight=dsqrt(1.0d0) ! at the moment this weight doesn't do anything
        COUL%lambda=10.0d0
        COUL%set_type='TRAIN'
        
        !lines to call the training of the dipoles
        if (train_COUL) then 
         call COUL%import_set(file_input=trim(geometry_file),len_file_inp=len_trim(geometry_file),type="DIPOLE")         
         call shift_geom(set_array=COUL%set,shift_file=shift_file)
         call COUL%import_labels
         
         call get_tot_kinds(COUL%set,tot_kinds)
         allocate(COUL%coeff_mask(tot_kinds))
         COUL%coeff_mask=.true.

         call COUL%build_matrix
         call COUL%build_target
         call COUL%LLS("DIPOLE")
         call COUL%predict_target
        end if
        
        if ((train_SNAP).or.(active_learning)) then
         call SNAP%import_set(file_input=trim(geometry_file),len_file_inp=len_trim(geometry_file),type="ENERGY")
         
         ! by setting the value of cutoff_dip and two_jmax_dip we do not affect the original parameter of the force field since the
         ! variables affecting it are set(i)%cutoff_en and set(i)%twojmax_en (for further reference see the import_set subroutine in
         ! the def_linear_model.f90 file

         do i=1,SNAP%nconfig
          SNAP%set(i)%cutoff_dip=SNAP_cutoff_dip
          SNAP%set(i)%twojmax_dip=SNAP_twojmax_dip
         end do
         call SNAP%import_labels
        if (VdW_flag) then 
          call SNAP%add_sub_VdW("sub")
        end if
        if (coulomb_flag) then

        call SNAP%FF_Coul%import_coeff
        call SNAP%FF_Coul%import_coeff_forces

          call SNAP%add_sub_Coul("sub")
        end if  
         call SNAP%build_matrix
         call SNAP%build_target
        
        if (train_SNAP) then
         call SNAP%LLS("ENERGY")
         call SNAP%predict_target
        else if (active_learning) then
         call SNAP%import_coeff
        end if

        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!MD block
        
        Dinamica%num_pot=2
        keyword_din="SNAP_VdW"
        iseed=[1469,2425,122,693]
        T_in=5.0d0
        Dinamica%T_bath=5.0d0
        Dinamica%step_size=1.0d0*41.49d0
        Dinamica%max_steps=1000
        Dinamica%ensemble='nvt'
        delta=1.5
        
        if (md_flag) then
        
        call import_lammps_obj(nconfig=1,file_input=trim(frame_file),len_file_inp=len_trim(frame_file),set_scalar=frame)       
        frame%mass=amu_to_emass*frame%mass
        !FF_SNAP and FF_VdW are optional arguments and as new potentials are introduced in the code, they can be added here
        call Dinamica%init_potentials(FF_SNAP,FF_VdW,frame,trim(frame_file))
        call Dinamica%import_linear_fits(active_learning,SNAP)
        call Dinamica%link_potentials(keyword_din,FF_SNAP,FF_VdW)
        call Dinamica%init(frame,iseed,T_in)
        call Dinamica%propagate(frame,active_learning,delta)
        
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!! MINIMIZATION BLOCK
        if (minim_flag) then
        
        Minimizzatore%num_pot=2
        delta=1.25
        keyword_min="SNAP_VdW"
        
        call import_lammps_obj(nconfig=1,file_input=trim(frame_file),len_file_inp=len_trim(frame_file),set_scalar=frame)
        call Minimizzatore%init_potentials(FF_SNAP,FF_VdW,frame,trim(frame_file))
        call Minimizzatore%import_linear_fits(active_learning,SNAP)
        call Minimizzatore%link_potentials(keyword_min,FF_SNAP,FF_VdW)
        call Minimizzatore%init(frame)
        call Minimizzatore%minimize(active_learning=active_learning,delta=delta)
        
        end if
        end program
