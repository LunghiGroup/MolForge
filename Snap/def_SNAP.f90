        module SNAP_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use lammps_class
        use max_class
        use kind_class
        use parameters_class
        use lapack_inverse
        implicit none

        contains

        subroutine generate_SNAP(set,num_bisp,energies,forces,lambda,ezero_constr,set_type,flag_energy,flag_forces,weight)
        implicit none

        type(lammps_obj),dimension(:),allocatable,intent(in)                  :: set
        real(kind=dbl), dimension(:,:), allocatable                           :: SNAP_matrix
        real(kind=dbl), dimension(:,:), allocatable                           :: C,F
        real(kind=dbl), dimension(:),allocatable                              :: b
        real(kind=dbl), dimension(:),allocatable                              :: ML_forces
        real(kind=dbl), dimension(:),allocatable,intent(inout)                :: energies,forces
        real(kind=dbl)                                                        :: ave_atom
        real(kind=dbl),intent(in)                                             :: weight
        character(len=5),intent(in)                                           :: set_type
        logical,dimension(:),allocatable,intent(in)                           :: ezero_constr
        real(kind=dbl),dimension(:),allocatable                               :: ML_energies
        integer                                                               :: tot_kinds,tot_atom,comp
        integer                                                               :: i,j,k,l,m
        integer,intent(in)                                                    :: num_bisp
        integer                                                               :: size_ref
        integer                                                               :: counter,count_kinds_ezero,ez_cons_rows
        integer                                                               :: tot_atoms
        integer                                                               :: start_cycle,end_cycle,help_counter
        integer                                                               :: start_snap_force
        integer                                                               :: ierror
        integer                                                               :: status
        character(len=80)                                                     :: err_msg
        character(len=1)                                                      :: TRANS
        integer                                                               :: MF,N,NRHS,LDA,LDB,LWORK,INFO
        real(kind=dbl),dimension(:),allocatable                               :: WORK
        real(kind=dbl),intent(in)                                             :: lambda
        real(kind=dbl), dimension(:,:), allocatable                           :: D,D_t,test,test_copy
        logical,intent(in)                                                    :: flag_forces,flag_energy
    
        call get_tot_kinds(set,tot_kinds)
        call get_ave_atoms(set,ave_atom,tot_atom)

        ez_cons_rows = count(ezero_constr)

        if ((flag_energy).and.(.not.flag_forces)) then

         size_ref=size(set) + ez_cons_rows

        end if

        if ((flag_forces).and.(.not.flag_energy)) then

         size_ref=3*tot_atom + ez_cons_rows

        end if

        if ((flag_forces).and.(flag_energy)) then

         size_ref=size(set)+3*tot_atom+ez_cons_rows

        end if

        if (lambda.ne.0.0) then

         size_ref=size_ref+(num_bisp-1)*tot_kinds

        end if
        
        allocate(SNAP_matrix(size_ref,num_bisp*tot_kinds))
        SNAP_matrix=0.0

        allocate(b(size_ref))
        b=0.0

        if ((flag_energy).and.(.not.flag_forces)) then

         b(1:size(set))=energies

        end if

        if ((flag_forces).and.(.not.flag_energy)) then

         b(1:3*tot_atom)=forces

        end if

        if ((flag_forces).and.(flag_energy)) then

        b(1:size(set))=weight*energies
        b(size(set)+1:size(set)+3*tot_atom)=forces

        end if

        if ((flag_energy).and.(.not.flag_forces)) then

         start_cycle=size(set)+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(set)

        end if

        if ((flag_forces).and.(.not.flag_energy)) then

         start_cycle=3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=3*tot_atom
        
        end if

        if ((flag_forces).and.(flag_energy)) then

         start_cycle=size(set)+3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(set)+3*tot_atom

        end if

        open(unit=1,file="vettore_target",action="write",iostat=ierror,iomsg=err_msg)
        
         do i=1,size_ref
          write(112,*) b(i)
         end do
        close(1)

        if (lambda.ne.0.0) then

         counter=0
         do j=start_cycle,end_cycle

          if (mod(j-help_counter,num_bisp-1)==1) then
           counter=counter+1
          end if

         SNAP_matrix(j,j-help_counter+counter) = sqrt(lambda)

         end do

        end if

        !tecnicamente questa matrice la stai inizializzando per righe che non e' il massimo in fortran, ma dovrebbe essere un minor
        !problem perche' non sara' lo step piu' time consuming

        if ((flag_energy).or.((flag_energy).and.(flag_forces))) then

         do j=1,size(set)
          do i=1,set(j)%nats
           do k=1,num_bisp
           
            SNAP_matrix(j,(set(j)%kind(i)-1)*num_bisp+k) = SNAP_matrix(j,(set(j)%kind(i)-1)*num_bisp+k) + set(j)%at_desc(i)%desc(k)
           
           end do
          end do
         end do

        allocate(C(size(set),num_bisp*tot_kinds))

        C=SNAP_matrix(1:size(set),:)
        SNAP_matrix(1:size(set),:)=SNAP_matrix(1:size(set),:)*weight
        
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FLAG E0 TO REGULARIZE THE FIRST COEFFICIENTS OF BISPECTRUM

        counter=0

        do j=1,tot_kinds
         if (ezero_constr(j)) then

        counter=counter+1
        SNAP_matrix(size_ref-ez_cons_rows+counter,(j-1)*num_bisp+1)=1

         end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CONSTRUCTING THE SNAP MATRIX OF FORCES

        if ((flag_forces).and.(.not.flag_energy)) then
         start_snap_force=0
        else if ((flag_forces).and.(flag_energy)) then
         start_snap_force=size(set)
        end if

        if (flag_forces) then

         do j=1,size(set)
          do i=1,set(j)%nkinds
           do l=1,set(j)%nats
            do comp=1,3
             do k=1,num_bisp-1

              SNAP_matrix(start_snap_force+(j-1)*set(j)%nats*3+(l-1)*3+comp,num_bisp*(i-1)+1+k)=&
                set(j)%der_at_desc(l)%desc((i-1)*((num_bisp-1)*3)+(num_bisp-1)*(comp-1)+k)
             
             end do
            end do
           end do
          end do
         end do

         allocate(F(3*tot_atom,num_bisp*tot_kinds))
         F=SNAP_matrix(start_snap_force+1:start_snap_force+3*tot_atom,:)
        
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SOLVING THE LINEAR LEAST SQUARES PROBLEM
        
        if (set_type=='TRAIN') then
        
         TRANS='N'
         M=size_ref
         LDA=size_ref
         N=tot_kinds*num_bisp
         NRHS=1
         LDB=max(M,N)
         LWORK=2*min(M,N)
         allocate(WORK(LWORK))

         call dgels(TRANS,M,N,NRHS,SNAP_matrix,LDA,b,LDB,WORK,LWORK,INFO)
         
         deallocate(WORK)
         deallocate(SNAP_matrix)

         if (info.ne.0) then
          write(*,*) 'Convergence issues: could not solve the linear least square problem (def_SNAP.F90/fit_energy)'
         end if

         open(11,file='snapcoeff_energy',action='write')
          do i=1,num_bisp*tot_kinds
          write(11,*) b(i)
          end do
         close(11)

         b=b(1:num_bisp*tot_kinds)

         !in case the set is for validation/test
         else

         deallocate(b)
         allocate(b(num_bisp*tot_kinds))

         open(11,file='snapcoeff_energy',action='read')
          do i=1,tot_kinds*num_bisp
           read(11,*) b(i)
          end do
         close(11)

        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CALCULATING THE PREDICTED ENERGIES

        if (flag_energy) then

         allocate(ML_energies(size(set)))
          ML_energies=matmul(C,b)
         deallocate(C)

        end if

        if (flag_forces) then

         allocate(ML_forces(3*tot_atom))
         ML_forces=matmul(F,b)         
         deallocate(F)
        
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRINTING TO FILES COEFFICIENTS OF FIT AND ENERGIES PREDICTED

        if (flag_energy) then

        open(11,file='energy_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML energy','             ','DFT energy','           ','Error'

        do i=1,size(set)
         write(11,*) ' ', i, ML_energies(i), energies(i),ML_energies(i)-energies(i)
        end do

         write(11,*)'# RMS=',(1/ave_atom)*(sqrt(sum((ML_energies-energies)**2)/size(set))), 'hartree/atom'
        close(11)

        end if

        if (flag_forces) then

        open(11,file='forces_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML force','             ','DFT force','           ','Error'

        do i=1,3*tot_atom
         write(11,*) ' ', i, ML_forces(i), forces(i),ML_forces(i)-forces(i)
        end do

        write(11,*)'# RMS=',(sqrt(sum((ML_forces(:)-forces(:))**2)&
        /(dble(3*tot_atom)))), 'a.u./bohr'

        close(11)

        end if
        deallocate(b)

        end subroutine generate_SNAP

        end module SNAP_class
