        module SNAP_fit_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use parameters_class
        use lapack_inverse
        implicit none

        type                                         :: SNAP_fit
        type(lammps_obj), allocatable                :: set(:)
        real(kind=dbl),allocatable                   :: matrix(:,:)
        real(kind=dbl),allocatable                   :: target(:)
        real(kind=dbl), allocatable                  :: energies(:)
        real(kind=dbl), allocatable                  :: forces(:)
        real(kind=dbl)                               :: weight
        real(kind=dbl)                               :: lambda
        real(kind=dbl), allocatable                  :: beta(:)
        integer                                      :: twojmax
        integer                                      :: num_bisp
        integer                                      :: nconfig
        logical, allocatable                         :: coeff_mask(:)
        logical                                      :: flag_energy
        logical                                      :: flag_forces
        character(len=5)                             :: set_type
        
        contains

        procedure                                    :: LLS_solve
        procedure                                    :: build_matrix

        end type SNAP_fit

        contains
        
        subroutine build_matrix(this)
        implicit none

        class(SNAP_fit)                                                       :: this
        integer                                                               :: counter,count_kinds_ezero,ez_cons_rows
        integer                                                               :: i,j,k,l,m,comp
        integer                                                               :: size_ref
        integer                                                               :: start_cycle,end_cycle,help_counter
        integer                                                               :: start_snap_force
        integer                                                               :: tot_atom
        integer                                                               :: tot_kinds
        real(kind=dbl)                                                        :: ave_atom
        do i=1,this%nconfig

         call this%set(i)%initialize()
         call this%set(i)%setup(this%set(i)%nkinds)
         if (this%flag_energy) then
          call this%set(i)%get_desc()
         end if
         if (this%flag_forces) then
          call this%set(i)%get_der_desc()
         end if
         call this%set(i)%finalize()

        end do

        call number_bispec(this%twojmax,this%num_bisp)
        call get_tot_kinds(this%set,tot_kinds)
        call get_ave_atoms(this%set,ave_atom,tot_atom)

        allocate(this%coeff_mask(tot_kinds))
        this%coeff_mask=.true.
        
        if (this%flag_energy) then
         this%coeff_mask(1)=.false.
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Calculation dimension of the matrix

        ez_cons_rows = count(this%coeff_mask)

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         size_ref=size(this%set) + ez_cons_rows

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         size_ref=3*tot_atom + ez_cons_rows

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         size_ref=size(this%set)+3*tot_atom+ez_cons_rows

        end if

        if (this%lambda.ne.0.0) then

         size_ref=size_ref+(this%num_bisp-1)*tot_kinds

        end if

        !!!!!!!!!!!!!!!!!!!!!!!! Allocation of the matrix and target vector

        allocate(this%matrix(size_ref,this%num_bisp*tot_kinds))
        this%matrix=0.0

        allocate(this%target(size_ref))
        this%target=0.0

        !!!!!!!!!!!!!!!!!!!!!!!! Initialization of matrix and target vector

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         this%target(1:size(this%set))=this%energies*this%weight

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         this%target(1:3*tot_atom)=this%forces

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then
         write(*,*)size(this%target,1)
         this%target(1:size(this%set))=this%weight*this%energies
         this%target(size(this%set)+1:size(this%set)+3*tot_atom)=this%forces

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! definizione variabili per fare i cicli

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         start_cycle=size(this%set)+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         start_cycle=3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=3*tot_atom

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         start_cycle=size(this%set)+3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)+3*tot_atom

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! qui scriviamo proprio la matrice di SNAP        

        if (this%lambda.ne.0.0) then

         counter=0
         do j=start_cycle,end_cycle

          if (mod(j-help_counter,this%num_bisp-1)==1) then
           counter=counter+1
          end if
        
        this%matrix(j,j-help_counter+counter) = sqrt(this%lambda)

         end do

        end if

        if ((this%flag_energy).or.((this%flag_energy).and.(this%flag_forces))) then

         do j=1,size(this%set)
          do i=1,this%set(j)%nats
           do k=1,this%num_bisp

            this%matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) = &
                    this%matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) + this%set(j)%at_desc(i)%desc(k)

           end do
          end do
         end do

        this%matrix(1:size(this%set),:)=this%matrix(1:size(this%set),:)*this%weight

        end if

        counter=0

        do j=1,tot_kinds
         if (this%coeff_mask(j)) then

        counter=counter+1
        this%matrix(size_ref-ez_cons_rows+counter,(j-1)*this%num_bisp+1)=1

         end if
        end do

        if ((this%flag_forces).and.(.not.this%flag_energy)) then
         start_snap_force=0
        else if ((this%flag_forces).and.(this%flag_energy)) then
         start_snap_force=size(this%set)
        end if

        if (this%flag_forces) then

        do j=1,size(this%set)
          do i=1,this%set(j)%nkinds
           do l=1,this%set(j)%nats
            do comp=1,3
             do k=1,this%num_bisp-1

              this%matrix(start_snap_force+(j-1)*this%set(j)%nats*3+(l-1)*3+comp,this%num_bisp*(i-1)+1+k)=&
                this%set(j)%der_at_desc(l)%desc((i-1)*((this%num_bisp-1)*3)+(this%num_bisp-1)*(comp-1)+k)

             end do
            end do
           end do
          end do
         end do

        end if
        
        do i=1,this%nconfig
         deallocate(this%set(i)%at_desc)
         deallocate(this%set(i)%der_at_desc)
        end do
        
        end subroutine build_matrix


        subroutine LLS_solve(this)
        implicit none
        
        class(SNAP_fit)                                                       :: this
        real(kind=dbl), dimension(:,:), allocatable                           :: C,F
        real(kind=dbl),dimension(:,:), allocatable                            :: temp_matrix
        real(kind=dbl), dimension(:), allocatable                             :: ML_forces
        real(kind=dbl),dimension(:),allocatable                               :: ML_energies
        character(len=1)                                                      :: TRANS
        integer                                                               :: i
        integer                                                               :: M
        integer                                                               :: MF,N,NRHS,LDA,LDB,LWORK,INFO
        real(kind=dbl)                                                        :: ave_atom
        integer                                                               :: tot_atom
        integer                                                               :: start_snap_force
        real(kind=dbl),dimension(:),allocatable                               :: WORK
        
        if (this%set_type == "TRAIN") then

         TRANS='N'
         M=size(this%matrix,1)
         LDA=size(this%matrix,1)
         N=size(this%matrix,2)
         NRHS=1
         LDB=max(M,N)
         LWORK=2*min(M,N)
         allocate(WORK(LWORK))
        
         allocate(temp_matrix(size(this%matrix,1),size(this%matrix,2)))
         temp_matrix=this%matrix
         call dgels(TRANS,M,N,NRHS,this%matrix,LDA,this%target,LDB,WORK,LWORK,INFO)
         
         deallocate(WORK)

         if (info.ne.0) then
          write(*,*) 'Convergence issues: could not solve the linear least square problem'
         end if

         open(11,file='snapcoeff_energy',action='write')
          do i=1,size(this%matrix,2)
          write(11,*) this%target(i)
          end do
         close(11)


         allocate(this%beta(size(this%matrix,2)))
         this%beta=this%target(1:size(this%matrix,2))

         else

         deallocate(this%target)
         allocate(this%beta(size(this%matrix,2)))

         open(11,file='snapcoeff_energy',action='read')
          do i=1,size(this%matrix,2)
           read(11,*) this%beta(i)
          end do
         close(11)

        end if

        if (this%flag_energy) then
         allocate(C(size(this%set),size(this%matrix,2)))
         
         C=temp_matrix(1:size(this%set),:)/this%weight
         allocate(ML_energies(size(this%set)))
          ML_energies=matmul(C,this%beta)
         deallocate(C)

        end if

        if (this%flag_forces) then

         call get_ave_atoms(this%set,ave_atom,tot_atom) 
         
         if (this%flag_energy) then
           start_snap_force=size(this%set)
         else 
           start_snap_force=0
         end if
        
        allocate(F(3*tot_atom,size(this%matrix,2)))
        F=temp_matrix(start_snap_force+1:start_snap_force+3*tot_atom,:)
        
        deallocate(this%matrix)

         allocate(ML_forces(3*tot_atom))
         ML_forces=matmul(F,this%beta)         
         deallocate(F)
        
        end if
        
        call get_ave_atoms(this%set,ave_atom,tot_atom)

        if (this%flag_energy) then

        open(11,file='energy_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML energy','             ','DFT energy','           ','Error'

        do i=1,size(this%set)
         write(11,*) ' ', i, ML_energies(i), this%energies(i),ML_energies(i)- this%energies(i)
        end do

         write(11,*)'# RMS=',(1/ave_atom)*(sqrt(sum((ML_energies-this%energies)**2)/&
                 size(this%set))), 'kcal/mol/atom'
        close(11)

        end if

        if (this%flag_forces) then

        open(11,file='forces_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML force','             ','DFT force','           ','Error'

        do i=1,3*tot_atom
         write(11,*) ' ', i, ML_forces(i), this%forces(i),ML_forces(i)-this%forces(i)
        end do

        write(11,*)'# RMS=',(sqrt(sum((ML_forces(:)-this%forces(:))**2)&
        /(dble(3*tot_atom)))), 'kcal/mol/angstrom'

        close(11)

        end if

        end subroutine LLS_solve

        end module SNAP_fit_class
