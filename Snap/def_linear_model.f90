        module linear_model_class
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use parameters_class
        use lapack_inverse
        use VdW_class
        implicit none 
        
        type,abstract                                :: linear_model
        type(lammps_obj), allocatable                :: set(:)
        real(kind=dbl),allocatable                   :: matrix(:,:)
        real(kind=dbl),allocatable                   :: target(:)
        real(kind=dbl)                               :: weight
        real(kind=dbl)                               :: lambda
        real(kind=dbl), allocatable                  :: beta(:)
        real(kind=dbl)                               :: s_z
        real(kind=dbl)                               :: cutoff
        integer                                      :: twojmax
        integer                                      :: num_bisp
        integer                                      :: nconfig
        logical, allocatable                         :: coeff_mask(:)
        character(len=5)                             :: set_type
        contains
        procedure                                    :: import_set
        procedure(import_labels),deferred            :: import_labels
        procedure(import_coeff),deferred             :: import_coeff
        procedure                                    :: LLS
        procedure(predict_target),deferred           :: predict_target
        procedure(build_matrix),deferred             :: build_matrix
        procedure(build_target),deferred             :: build_target
        procedure(get_uncertainty),deferred          :: get_uncertainty
        
        end type linear_model
        
        abstract interface
        
        subroutine import_labels(this)
        import linear_model
        implicit none
        class(linear_model)     :: this
        end subroutine

        subroutine import_coeff(this)
        import linear_model
        implicit none
        class(linear_model)     :: this
        end subroutine

        subroutine predict_target(this)
        import linear_model
        implicit none
        class(linear_model)     :: this
        end subroutine

        subroutine build_matrix(this)
        import linear_model
        implicit none
        class(linear_model)     :: this
        end subroutine

        subroutine build_target(this)
        import linear_model
        implicit none
        class(linear_model)     :: this
        end subroutine

        subroutine get_uncertainty(this,frame,calc_sz_flag,error)
        use lammps_class
        use kind_class
        import linear_model
        implicit none
        class(linear_model)                                     :: this
        type(lammps_obj)                                        :: frame
        logical,intent(in)                                      :: calc_sz_flag
        real(kind=dbl),intent(out)                              :: error
        end subroutine

        end interface

        contains

        subroutine import_set(this,file_input,len_file_inp,type)
        implicit none
        class(linear_model)                       :: this
        integer                                   :: nconfig, i,j,nats,ntypes
        character(len=100),allocatable            :: tmp(:,:)
        character(len=100),dimension(10)          :: tmp_cell_nkinds
        character(len=100)                        :: filename
        character(len=*)                          :: file_input
        character(len=*)                          :: type
        integer,intent(in)                        :: len_file_inp
        character(len=150)                        :: file_inp

        call import_lammps_obj(nconfig=this%nconfig,file_input=file_input,len_file_inp=len_file_inp,set_array=this%set)
        
        if (type=="ENERGY") then
         do i=1,this%nconfig

          this%set(i)%twojmax_en = this%twojmax
          this%set(i)%cutoff_en  = this%cutoff

         end do
        else if (type=="DIPOLE") then
         do i=1,this%nconfig

          this%set(i)%twojmax_dip = this%twojmax
          this%set(i)%cutoff_dip  = this%cutoff

         end do
        end if

        end subroutine import_set
        
        subroutine LLS(this,fitting_quantity)
        implicit none
        class(linear_model)                                                   :: this
        real(kind=dbl),dimension(:,:), allocatable                            :: temp_matrix
        real(kind=dbl), dimension(:), allocatable                             :: temp_target
        character(len=1)                                                      :: TRANS
        character(len=*)                                                      :: fitting_quantity
        integer                                                               :: i
        integer                                                               :: M
        integer                                                               :: MF,N,NRHS,LDA,LDB,LWORK,INFO
        real(kind=dbl),dimension(:),allocatable                               :: WORK
        
        if (this%set_type=="TRAIN") then
         TRANS='N'
         M=size(this%matrix,1)
         LDA=size(this%matrix,1)
         N=size(this%matrix,2)
         NRHS=1
         LDB=max(M,N)
         LWORK=2*min(M,N)
         allocate(WORK(LWORK))
        
         allocate(temp_matrix(size(this%matrix,1),size(this%matrix,2)))
         allocate(temp_target(size(this%matrix,1)))
         temp_matrix=this%matrix
         temp_target=this%target
         call dgels(TRANS,M,N,NRHS,temp_matrix,LDA,temp_target,LDB,WORK,LWORK,INFO)
         
         deallocate(WORK)
         deallocate(temp_matrix)
         if (info.ne.0) then
          write(*,*) 'Convergence issues: could not solve the linear least square problem'
         end if
        
         if (fitting_quantity=='ENERGY') then
          open(11,file='snapcoeff_energy',action='write')
           do i=1,size(this%matrix,2)
            write(11,*) temp_target(i)
           end do
          close(11)
         else if (fitting_quantity=='DIPOLE') then
          open(11,file='snapcoeff_dipoles',action='write')
           do i=1,size(this%matrix,2)
            write(11,*) temp_target(i)
           end do
          close(11)
         end if

         allocate(this%beta(size(this%matrix,2)))
         this%beta=temp_target(1:size(this%matrix,2))

         else

         call this%import_coeff

        end if
        
        end subroutine

        end module linear_model_class
