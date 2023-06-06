        module general_types_class
               
        type :: sarc
         double precision                   :: k(3)
         complex(8), allocatable            :: v(:)
         integer, allocatable               :: vi(:)
        end type sarc

        type :: sub_space
         integer, allocatable               :: kind(:)
         integer                            :: max_ex
         double precision                   :: dist_max
        end type sub_space

        type :: vector_dbl
         double precision, allocatable :: v(:)
        end type vector_dbl

        type :: vector_int
         integer, allocatable :: v(:)
        end type vector_int

        type :: vector_cmplx
         complex(8), allocatable :: v(:)
        end type vector_cmplx

        type :: mat_cmplx
         complex(8), allocatable :: mat(:,:)
         contains
         procedure      :: delete => delete_cmplx_mat
        end type mat_cmplx

        contains

        subroutine delete_cmplx_mat(this)
        implicit none
        class(mat_cmplx)  :: this

         if(allocated(this%mat)) deallocate(this%mat)

        return
        end subroutine delete_cmplx_mat


        end module  general_types_class
