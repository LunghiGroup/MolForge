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
        procedure(LLS),deferred                      :: LLS
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

        subroutine LLS(this)
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

        subroutine import_set(this,file_input,len_file_inp)
        implicit none
        class(linear_model)                       :: this
        integer                                   :: nconfig, i,j,nats,ntypes
        character(len=100),allocatable            :: tmp(:,:)
        character(len=100),dimension(10)          :: tmp_cell_nkinds
        character(len=100)                        :: filename
        character(len=*)                          :: file_input
        integer,intent(in)                        :: len_file_inp
        character(len=150)                        :: file_inp

        call import_lammps_obj(nconfig=this%nconfig,file_input=file_input,len_file_inp=len_file_inp,set_array=this%set)

        do i=1,this%nconfig

         this%set(i)%twojmax = this%twojmax
         this%set(i)%cutoff  = this%cutoff

        end do

        end subroutine import_set
        

        end module linear_model_class
