        module target_functions_class
        implicit none
         
         type, abstract :: target_function 
          contains
          procedure(get_fval), deferred    :: get_fval
          procedure(get_fgrad), deferred   :: get_fgrad
         end type target_function

         interface
          subroutine get_fval(this,vec,val)
          import target_function
          implicit none
          class(target_function)           :: this
          double precision, allocatable    :: vec(:)
          double precision                 :: val
          end subroutine get_fval
         end interface

         interface
          subroutine get_fgrad(this,vec,val,grad)
          import target_function
          implicit none
          class(target_function)           :: this
          double precision                 :: val
          double precision, allocatable    :: grad(:),vec(:)
          end subroutine get_fgrad
         end interface

        contains

        end module target_functions_class
