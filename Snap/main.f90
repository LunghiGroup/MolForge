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
        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(SNAP_fit)                            :: SNAP
        type(SNAP_FF),target                        :: FF_SNAP
        type(VdW_FF),target                              :: FF_VdW
        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=120)                        :: geometry_file,energy_file,dipoles_file,shift_file,atom_string, &
                                                        frame_file,forces_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        real(kind=dbl)                            :: lambda_en,cutoff_en
        double precision                          :: lambda_dip,cutoff_dip
        logical                                   :: dipole_flag,energy_flag,md_flag,VdW_flag,coul_flag,minim_flag, &
                train_ff,rampa_flag,phonon_flag,single_eval
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
        integer                                   :: tot_steps_md
        double precision                          :: val,step
        double precision,allocatable              :: vec(:),grad(:)
        double precision,allocatable              :: VdW_en(:),test_1(:,:),test_2(:,:)
        double precision,allocatable              :: grads(:,:)
        double precision                          :: edisp
        double precision                          :: test(3,3)
        integer                                   :: INFO
        integer                                   :: len_shift_file
        integer                                   :: idist
        integer,dimension(4)                      :: iseed
        integer                                   :: N
        double precision, allocatable             :: X_M(:), Y(:),X(:),mean(:),sigma(:)
        double precision                          :: y_tmp,y_tmp_2
        real(kind=dbl)                            :: weight
        real(kind=dbl)                            :: error

        train_ff=.true.
        VdW_flag=.true.

        SNAP%nconfig=19
        cutoff_en=4.5d0
        twojmax_en=12
        
        SNAP%weight=20.0d0
        SNAP%lambda=12.0d0
        SNAP%set_type='TRAIN'

        geometry_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_tr_AL_++"
        energy_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/ener_tr_AL_++"
        forces_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/grad_tr_AL_++"
        frame_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/geo_start"
        
        SNAP%twojmax= twojmax_en
        SNAP%flag_energy=.true.
        SNAP%flag_forces=.true.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !GENERATE FF POTENTIAL
        if (train_ff) then
   
         call import_lammps_obj_list(nconfig=SNAP%nconfig,file_input=trim(geometry_file),&
                len_file_inp=len_trim(geometry_file),set_array=SNAP%set)
        
        do i=1,SNAP%nconfig
         SNAP%set(i)%twojmax = twojmax_en
         SNAP%set(i)%cutoff  = cutoff_en
        end do 
!!!!!!!!questa parte andrebbe tutta in un import del training set
         if (SNAP%flag_energy) then

          allocate(SNAP%energies(size(SNAP%set)))

          open(100,file=trim(energy_file))

          do i=1,SNAP%nconfig
           read(100,*) SNAP%energies(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!questa andrebbe in una parte di aggiustamento delle energie in base a dispersione and
!stuff           
           if (VdW_flag) then

                   FF_VdW%frame=SNAP%set(i)
                   call FF_VdW%get_fval(vec,FF_VdW%energy)
                   SNAP%energies(i)=SNAP%energies(i)-FF_VdW%energy*Har_to_Kc
           end if
          
          end do

          close(100)

         end if
         
         if (SNAP%flag_forces) then
          call get_ave_atoms(SNAP%set,ave_atom,tot_atom)
          allocate(gradients(3*tot_atom))

          open(100,file=trim(forces_file))

           do j=1,3*tot_atom
            read(100,*) gradients(j)
           end do

          close(100)

          if (VdW_flag) then

            do i=1,SNAP%nconfig
             
             FF_VdW%frame=SNAP%set(i)
             call FF_VdW%get_fgrad(vec,FF_VdW%energy,FF_VdW%grad)


             gradients((i-1)*3*SNAP%set(i)%nats+1:i*3*SNAP%set(i)%nats)=&
             gradients((i-1)*3*SNAP%set(i)%nats+1:i*3*SNAP%set(i)%nats)-FF_VdW%grad*F_conv

           end do
          end if
        
         SNAP%forces=-gradients
        
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
         call SNAP%build_matrix
         call SNAP%build_target
         call SNAP%LLS_solve
         call SNAP%get_uncertainty(SNAP%set(9),.true.,error)
         write(*,*)error
         
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!MD block
        
        end program
