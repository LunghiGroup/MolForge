        module coul_fit_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use parameters_class
        use lapack_inverse
        use VdW_class
        use linear_model_class
        implicit none

        type,extends(linear_model)                   :: coul_fit
        real(kind=dbl), allocatable                  :: dipoles(:)
        real(kind=dbl), allocatable                  :: charges(:)
        character(len=120)                           :: dipoles_file
        character(len=120)                           :: charges_file
        character(len=120)                           :: shift_file

        contains

        procedure                                    :: import_labels => import_dipoles_charges_shift
        procedure                                    :: import_coeff => import_coeff_SNAP_dip
        procedure                                    :: predict_target => predict_dipoles_charges
        procedure                                    :: build_matrix => build_matrix_SNAP_dipoles
        procedure                                    :: build_target => build_target_dipoles
        procedure                                    :: get_uncertainty => get_uncertainty_dipoles

        end type coul_fit

        contains

        subroutine import_dipoles_charges_shift(this)
        implicit none
        class(coul_fit)                      :: this
        integer                              :: tot_atom
        real(kind=dbl)                       :: ave_atom
        integer                              :: i,j

        allocate(this%dipoles(3*size(this%set)))
        allocate(this%charges(size(this%set)))

        open(1,file=trim(this%dipoles_file))

        do j=1,3*size(this%set)
         read(1,*) this%dipoles(j)
        end do

        close(1)
        
        open(1,file=trim(this%charges_file))

        do j=1,size(this%set)
         read(1,*) this%charges(j)
        end do

        close(1)

        end subroutine import_dipoles_charges_shift

        subroutine import_coeff_SNAP_dip(this)
        implicit none
        class(coul_fit)                                         :: this
        integer                                                 :: i
        integer                                                 :: counter
        allocate(this%beta(size(this%matrix,2)))

        open(11,file='snapcoeff_dipoles',action='read')
        ! setting this loop to go up to the number of columns of the matrix makes necessary to define it before calling this
        ! subroutine        
        do i=1,size(this%matrix,2)
         read(11,*) this%beta(i)
        end do
        close(11)
        
        end subroutine import_coeff_SNAP_dip

        subroutine build_matrix_SNAP_dipoles(this)
        implicit none
        
        class(coul_fit)                                         :: this
        integer                                                 :: tot_atom
        integer                                                 :: tot_kinds
        integer                                                 :: i,j,k,m
        integer                                                 :: ez_cons_rows
        integer                                                 :: size_ref
        integer                                                 :: counter
        real(kind=dbl)                                          :: ave_atom


        do i=1,this%nconfig
         
         call this%set(i)%initialize()
         call this%set(i)%setup(this%set(i)%nkinds)
         call this%set(i)%get_desc("DIPOLE")
         call this%set(i)%finalize()

        end do
        
        call number_bispec(this%twojmax,this%num_bisp)
        call get_tot_kinds(this%set,tot_kinds)
        call get_ave_atoms(this%set,ave_atom,tot_atom)
        ez_cons_rows=count(this%coeff_mask)

        if (this%lambda.ne.0.0) then
         size_ref=3*size(this%set)+(this%num_bisp-1)*tot_kinds+size(this%set)+ez_cons_rows
        else
         size_ref=3*size(this%set)+size(this%set)+ez_cons_rows
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!Allocation of the matrix

        allocate(this%matrix(size_ref,this%num_bisp*tot_kinds))
        this%matrix=0.0
        do j=1,size(this%set)
         do m=1,3
          do i=1,this%set(j)%nats
           do k=1,this%num_bisp

            this%matrix(3*(j-1)+m,(this%set(j)%kind(i)-1)*this%num_bisp+k) = this%matrix(3*(j-1)+ m&
            ,(this%set(j)%kind(i)-1)*this%num_bisp+k)+ this%set(j)%at_desc_dip(i)%desc(k)*this%set(j)%x(i,m)*A_to_B

           end do
          end do
         end do
        end do

        !converting the coordinates from angstrom to bohr.
        !this%matrix=A_to_B*this%matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (this%lambda.ne.0.0) then
         counter=0

         do j=3*size(this%set)+1,size_ref-size(this%set)-ez_cons_rows
         
          if (mod(j-3*size(this%set),this%num_bisp-1)==1) then
           counter=counter+1
          end if
          this%matrix(j,j-3*size(this%set)+counter) = sqrt(this%lambda)
         
         end do

        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!here we define the rows corresponding to conservation of charge
        counter=0

        do j=size_ref-size(this%set)-ez_cons_rows+1,size_ref-ez_cons_rows
         counter=counter+1
         do i=1,this%set(counter)%nats
          do k=1,this%num_bisp

           this%matrix(j,(this%set(counter)%kind(i)-1)*this%num_bisp+k) = this%matrix&
           (j,(this%set(counter)%kind(i)-1)*this%num_bisp+k) +this%set(counter)%at_desc_dip(i)%desc(k)

          end do
         end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!here we constrain to 0 a set of coefficients

        counter=0
        
        do j=1,tot_kinds
         if (this%coeff_mask(j)) then
          counter=counter+1
          this%matrix(size_ref - ez_cons_rows + counter,(j-1)*this%num_bisp + 1)= 1
         end if
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!deallocation of descriptors

        do i=1,this%nconfig
         deallocate(this%set(i)%at_desc_dip)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end subroutine build_matrix_SNAP_dipoles

        subroutine build_target_dipoles(this)
        implicit none
        class(coul_fit)         :: this
        integer                 :: tot_kinds
        integer                 :: j
        integer                 :: tot_atom
        real(kind=dbl)          :: ave_atom

        call get_tot_kinds(this%set,tot_kinds)
        
        allocate(this%target(size(this%matrix,1)))
        this%target=0.0
        
        this%target(1:3*size(this%set))=this%dipoles
        this%target(3*size(this%set)+(this%num_bisp-1)*tot_kinds+1:3*size(this%set)+(this%num_bisp-1)&
                *tot_kinds+size(this%set))=this%charges
        
        open(222,file='target',position="append")

        do j=1,size(this%matrix,1)
         write(222,*) this%target(j)
        end do

        close(222)

        end subroutine build_target_dipoles
        
        subroutine predict_dipoles_charges(this)
        implicit none
        class(coul_fit)                                                       :: this
        real(kind=dbl), dimension(:,:), allocatable                           :: C
        real(kind=dbl), dimension(:), allocatable                             :: ML_dipoles
        real(kind=dbl)                                                        :: ave_atom
        integer                                                               :: i 
        integer                                                               :: tot_atom

        allocate(C(3*size(this%set),size(this%matrix,2)))
        C=this%matrix(1:3*size(this%set),:)
        allocate(ML_dipoles(3*size(this%set)))
        ML_dipoles=matmul(C,this%beta)
        deallocate(C)

        !call get_ave_atoms(this%set,ave_atom,tot_atom)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        open(11,file='dipoles_rms_fpp.dat',action='write')

        write(11,*) '#','      ','Num config','       ','ML dipole','             ','DFT dipole','           ','Error'

        do i=1,3*size(this%set)
         write(11,*) ' ', i, ML_dipoles(i), this%dipoles(i),ML_dipoles(i)-this%dipoles(i)
        end do

        write(11,*)'# RMS=',sqrt(sum((ML_dipoles-this%dipoles)**2)/(3*size(this%set))), 'a.u. (dipoles)'
        close(11)

        deallocate(this%dipoles,ML_dipoles)
        
        end subroutine predict_dipoles_charges

        subroutine get_uncertainty_dipoles(this,frame,calc_sz_flag,error)
        implicit none
        class(coul_fit)                           :: this
        type(lammps_obj)                          :: frame
        logical,intent(in)                        :: calc_sz_flag
        real(kind=dbl),intent(out)                :: error
        end subroutine get_uncertainty_dipoles
        
        end module coul_fit_class
