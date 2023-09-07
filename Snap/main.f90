        program main
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use potential_class
        use SNAP_fit_class
        use VdW_class
        use SNAP_FF_class
        use md_class
        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(SNAP_fit)                            :: SNAP
        type(SNAP_FF)                             :: FF_SNAP
        type(VdW_FF)                              :: FF_VdW
        type(MD)                                  :: Dinamica
        type(lammps_obj)                          :: frame
        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=120)                        :: geometry_file,energy_file,dipoles_file,shift_file,atom_string, &
                                                        frame_file,forces_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        real(kind=dbl)                            :: lambda_en,cutoff_en
        double precision                          :: lambda_dip,cutoff_dip
        logical                                   :: dipole_flag,energy_flag,md_flag,VdW_flag,coul_flag,minim_flag,train,&
                rampa_flag,phonon_flag,single_eval
        logical                                   :: flag_forces,flag_energy
        double precision,dimension(:),allocatable :: energies,coul_energy,snap_energy
        double precision,dimension(:),allocatable :: forces,gradients
        double precision                          :: R_screen
        double precision,dimension(:),allocatable :: tot_charge
        character(len=5)                          :: set_type_en,set_type_dip
        double precision,allocatable              :: force(:,:)
        double precision                          :: E_plus,E_minus,ave_atom
        integer                                   :: tot_atom
        integer                                   :: atom,direction,k
        double precision,dimension(:),allocatable :: coul_energy_plus,coul_energy_minus
        double precision                          :: temperature,timestep
        real(kind=dbl)                            :: T_in
        integer                                   :: tot_steps_md
        double precision,allocatable              :: vec(:),grad(:),grad_SNAP(:),grad_VdW(:)
        integer                                   :: idist
        integer,dimension(4)                      :: iseed
        real(kind=dbl)                            :: weight
        real(kind=dbl)                            :: error
        real(kind=dbl)                            :: val

        train=.true.
        VdW_flag=.true.
        md_flag=.true.

        SNAP%nconfig=19
        SNAP%cutoff=4.0d0
        SNAP%twojmax=8
        SNAP%weight=dsqrt(17.0d0*3.0)
        SNAP%lambda=0.1d0
        SNAP%set_type='TRAIN'
        
        frame%cutoff=4.0d0
        frame%twojmax=8
        
        geometry_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_tr_AL_++"
        SNAP%energy_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/ener_tr_AL_++"
        SNAP%forces_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/grad_tr_AL_++"
        frame_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_start"

        SNAP%flag_energy=.true.
        SNAP%flag_forces=.true.

        if (train) then
         call SNAP%import_set(file_input=trim(geometry_file),len_file_inp=len_trim(geometry_file))
         call SNAP%import_labels
        if (VdW_flag) then 
          call SNAP%add_sub_VdW("sub")
        end if  
         call SNAP%build_matrix
         call SNAP%build_target
         call SNAP%LLS
         !call SNAP%get_uncertainty(SNAP%set(9),.true.,error)
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!MD block
        
        Dinamica%num_pot=2
        iseed=[1469,2425,122,693]
        T_in=5.0d0
        Dinamica%T_bath=5.0d0
        Dinamica%step_size=1.0d0*41.49d0
        Dinamica%max_steps=4000
        Dinamica%ensemble='nvt'
        
        if (md_flag) then
        
        call import_lammps_obj(nconfig=1,file_input=trim(frame_file),len_file_inp=len_trim(frame_file),set_scalar=frame)       
        !transformazione unita' da kcal related to atomic units
        frame%mass=amu_to_emass*frame%mass
        !!!!!!!!!!!!!!!!!!!!!!
        
        call number_bispec(frame%twojmax,FF_SNAP%num_bisp)
        FF_SNAP%tot_kinds=frame%nkinds

        call FF_SNAP%import_coeff
        FF_SNAP%frame=frame
        call FF_SNAP%get_fgrad(vec,val,grad)
        write(*,*) grad
        if (VdW_flag) then
         call FF_VdW%import_coeff
         call FF_VdW%import_geo(frame_file)
        end if
        
        call Dinamica%link_potentials("SNAP_VdW",FF_SNAP,FF_VdW)
        call Dinamica%initialize_vel(frame,iseed,T_in)
        call Dinamica%propagate(frame)
        end if



        !!!!!!!!!!!!!!!!!!!!!!!!!!
        end program
