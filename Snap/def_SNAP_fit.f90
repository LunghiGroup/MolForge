        module SNAP_fit_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use atoms_class
        use lammps_class
        use parameters_class
        use lapack_inverse
        implicit none

        type                                         :: SNAP_fit
        type(lammps_obj), allocatable                :: set(:)
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
        real(kind=dbl)                               :: ave_atom
        integer                                      :: tot_atom
        integer                                      :: tot_kinds
        
        contains

        procedure                        :: fit

        end type SNAP_fit

        contains
        
        subroutine fit(this)
        implicit none
        
        class(SNAP_fit)                                                       :: this
        real(kind=dbl), dimension(:,:), allocatable                           :: SNAP_matrix
        real(kind=dbl), dimension(:,:), allocatable                           :: C,F
        real(kind=dbl), dimension(:), allocatable                             :: ref_values                      
        real(kind=dbl), dimension(:), allocatable                             :: ML_forces
        real(kind=dbl),dimension(:),allocatable                               :: ML_energies
        integer                                                               :: i,j,k,l,m,comp
        integer                                                               :: size_ref
        integer                                                               :: counter,count_kinds_ezero,ez_cons_rows
        integer                                                               :: start_cycle,end_cycle,help_counter
        integer                                                               :: start_snap_force
        integer                                                               :: ierror
        integer                                                               :: status
        character(len=80)                                                     :: err_msg
        character(len=1)                                                      :: TRANS
        integer                                                               :: MF,N,NRHS,LDA,LDB,LWORK,INFO
        real(kind=dbl),dimension(:),allocatable                               :: WORK
        
        ez_cons_rows = count(this%coeff_mask)

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         size_ref=size(this%set) + ez_cons_rows

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         size_ref=3*this%tot_atom + ez_cons_rows

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         size_ref=size(this%set)+3*this%tot_atom+ez_cons_rows

        end if

        if (this%lambda.ne.0.0) then

         size_ref=size_ref+(this%num_bisp-1)*this%tot_kinds

        end if
        
        allocate(SNAP_matrix(size_ref,this%num_bisp*this%tot_kinds))
        SNAP_matrix=0.0

        allocate(ref_values(size_ref))
        ref_values=0.0

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         ref_values(1:size(this%set))=this%energies

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         ref_values(1:3*this%tot_atom)=this%forces

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         ref_values(1:size(this%set))=this%weight*this%energies
         ref_values(size(this%set)+1:size(this%set)+3*this%tot_atom)=this%forces

        end if

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         start_cycle=size(this%set)+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         start_cycle=3*this%tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=3*this%tot_atom
        
        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         start_cycle=size(this%set)+3*this%tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)+3*this%tot_atom

        end if

        if (this%lambda.ne.0.0) then

         counter=0
         do j=start_cycle,end_cycle

          if (mod(j-help_counter,this%num_bisp-1)==1) then
           counter=counter+1
          end if

          SNAP_matrix(j,j-help_counter+counter) = sqrt(this%lambda)

         end do

        end if

        if ((this%flag_energy).or.((this%flag_energy).and.(this%flag_forces))) then

         do j=1,size(this%set)
          do i=1,this%set(j)%nats
           do k=1,this%num_bisp
           
            SNAP_matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) = &
                    SNAP_matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) + this%set(j)%at_desc(i)%desc(k)
           
           end do
          end do
         end do

        allocate(C(size(this%set),this%num_bisp*this%tot_kinds))

        C=SNAP_matrix(1:size(this%set),:)
        SNAP_matrix(1:size(this%set),:)=SNAP_matrix(1:size(this%set),:)*this%weight
        
        end if

        counter=0

        do j=1,this%tot_kinds
         if (this%coeff_mask(j)) then

        counter=counter+1
        SNAP_matrix(size_ref-ez_cons_rows+counter,(j-1)*this%num_bisp+1)=1

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

              SNAP_matrix(start_snap_force+(j-1)*this%set(j)%nats*3+(l-1)*3+comp,this%num_bisp*(i-1)+1+k)=&
                this%set(j)%der_at_desc(l)%desc((i-1)*((this%num_bisp-1)*3)+(this%num_bisp-1)*(comp-1)+k)
             
             end do
            end do
           end do
          end do
         end do

         allocate(F(3*this%tot_atom,this%num_bisp*this%tot_kinds))
         F=SNAP_matrix(start_snap_force+1:start_snap_force+3*this%tot_atom,:)
        
        end if

        if (this%set_type =='TRAIN') then
        
         TRANS='N'
         M=size_ref
         LDA=size_ref
         N=this%tot_kinds*this%num_bisp
         NRHS=1
         LDB=max(M,N)
         LWORK=2*min(M,N)
         allocate(WORK(LWORK))

         call dgels(TRANS,M,N,NRHS,SNAP_matrix,LDA,ref_values,LDB,WORK,LWORK,INFO)
         
         deallocate(WORK)
         deallocate(SNAP_matrix)

         if (info.ne.0) then
          write(*,*) 'Convergence issues: could not solve the linear least square problem (def_SNAP.F90/fit_energy)'
         end if

         open(11,file='snapcoeff_energy',action='write')
          do i=1,this%num_bisp*this%tot_kinds
          write(11,*) ref_values(i)
          end do
         close(11)


         allocate(this%beta(this%num_bisp*this%tot_kinds))
         this%beta=ref_values(1:this%num_bisp*this%tot_kinds)

         else

         deallocate(ref_values)
         allocate(this%beta(this%num_bisp*this%tot_kinds))

         open(11,file='snapcoeff_energy',action='read')
          do i=1,this%tot_kinds*this%num_bisp
           read(11,*) this%beta(i)
          end do
         close(11)

        end if

        if (this%flag_energy) then

         allocate(ML_energies(size(this%set)))
          ML_energies=matmul(C,this%beta)
         deallocate(C)

        end if

        if (this%flag_forces) then

         allocate(ML_forces(3*this%tot_atom))
         ML_forces=matmul(F,this%beta)         
         deallocate(F)
        
        end if
        

        if (this%flag_energy) then

        open(11,file='energy_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML energy','             ','DFT energy','           ','Error'

        do i=1,size(this%set)
         write(11,*) ' ', i, ML_energies(i), this%energies(i),ML_energies(i)- this%energies(i)
        end do

         write(11,*)'# RMS=',(1/this%ave_atom)*(sqrt(sum((ML_energies-this%energies)**2)/size(this%set))), 'kcal/mol/atom'
        close(11)

        end if

        if (this%flag_forces) then

        open(11,file='forces_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML force','             ','DFT force','           ','Error'

        do i=1,3*this%tot_atom
         write(11,*) ' ', i, ML_forces(i), this%forces(i),ML_forces(i)-this%forces(i)
        end do

        write(11,*)'# RMS=',(sqrt(sum((ML_forces(:)-this%forces(:))**2)&
        /(dble(3*this%tot_atom)))), 'kcal/mol/angstrom'

        close(11)

        end if

        end subroutine fit

        end module SNAP_fit_class
