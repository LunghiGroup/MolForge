        program main
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use SNAP_fit_class
        use VdW_class
        use SNAP_FF_class
        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(SNAP_fit)                            :: SNAP
        type(SNAP_FF)                             :: FF_SNAP
        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=120)                        :: geometry_file,energy_file,dipoles_file,shift_file,atom_string, &
                                                        shift_geo_file,forces_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        real(kind=dbl)                            :: lambda_en,cutoff_en
        double precision                          :: lambda_dip,cutoff_dip
        logical                                   :: dipole_flag,energy_flag,md_flag,VdW_flag,coul_flag,minim_flag, &
                train_ff,rampa_flag,phonon_flag
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

        train_ff=.false.
        VdW_flag=.false.
        
        SNAP%nconfig=375
        cutoff_en=4.0d0
        twojmax_en=8
        
        SNAP%weight=1.0d0
        SNAP%lambda=0.0d0
        SNAP%set_type='TRAIN'

        geometry_file="/home/valeriobriganti/Desktop/Progetto_Phonons/Test_ravera/FitSnap/geo_test_ravera"
        energy_file="/home/valeriobriganti/Desktop/Progetto_Phonons/Test_ravera/FitSnap/ener_test_ravera"
        forces_file="/home/valeriobriganti/Desktop/MolForge_SNAP/Snap/test_files/grad_tr_AL_++"
       
        call import_lammps_obj_list(SNAP%set,SNAP%nconfig,trim(geometry_file),len_trim(geometry_file))     
       
        do i=1,SNAP%nconfig
         SNAP%set(i)%twojmax = twojmax_en
         SNAP%set(i)%cutoff  = cutoff_en
        end do


        SNAP%twojmax= twojmax_en
        
        !the coeffiecients that you set to be true are the ones that you put equal to zero
        
        call number_bispec(SNAP%twojmax,SNAP%num_bisp)
        call get_ave_atoms(SNAP%set,SNAP%ave_atom,SNAP%tot_atom)
        call get_tot_kinds(SNAP%set,SNAP%tot_kinds)
        
        allocate(SNAP%coeff_mask(SNAP%tot_kinds))
        
        SNAP%coeff_mask=.true.
        SNAP%coeff_mask(1)=.false.
        SNAP%flag_energy=.true.
        SNAP%flag_forces=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,SNAP%nconfig
       
         call SNAP%set(i)%initialize()
         call SNAP%set(i)%setup(SNAP%set(i)%nkinds)
         call SNAP%set(i)%get_desc()
         call SNAP%set(i)%get_der_desc()
         call SNAP%set(i)%finalize()

        end do
        !GENERATE FF POTENTIAL
        
        if (train_ff) then

         allocate(SNAP%energies(size(set)))

         open(100,file=trim(energy_file))

          do j=1,SNAP%nconfig
           read(100,*) SNAP%energies(j)
          end do

         close(100)
        

         if (SNAP%flag_forces) then

          allocate(gradients(3*SNAP%tot_atom))

          open(100,file=trim(forces_file))

           do j=1,3*SNAP%tot_atom
            read(100,*) gradients(j)
           end do

          close(100)

         end if
        
         if (VdW_flag) then

          do i=1,SNAP%nconfig
        
           call grimme_d3(SNAP%set(i))
         
           if (SNAP%flag_energy) then

           SNAP%energies(i)=SNAP%energies(i)-SNAP%set(i)%en_VdW
            
           end if
           
           if (SNAP%flag_forces) then

           allocate(SNAP%set(i)%grads_VdW(3,SNAP%set(i)%nats))

           do j=1,SNAP%set(i)%nats
            do k=1,3
             gradients(((i-1)*SNAP%set(i)%nats*3) +(j-1)*3+k)=&
             gradients(((i-1)*SNAP%set(i)%nats*3) +(j-1)*3+k)-SNAP%set(i)%grads_VdW(k,j)
            end do
           end do
           
           deallocate(SNAP%set(i)%grads_VdW)

           end if
        
          end do
         !SNAP%forces=-gradients       
         end if 
         call SNAP%fit

        end if
        
        FF_SNAP%frame=SNAP%set(1)
        FF_SNAP%num_bisp=SNAP%num_bisp
        FF_SNAP%tot_kinds=SNAP%tot_kinds
        call FF_SNAP%import

        call FF_SNAP%get_fval(vec,FF_SNAP%energy)
        call FF_SNAP%get_fgrad(vec,FF_SNAP%energy,FF_SNAP%grad)
        write(*,*) FF_SNAP%energy
        write(*,*) grad
        end program
