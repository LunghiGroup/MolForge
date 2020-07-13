        module sna_class
        implicit none

         type sna
          double precision, allocatable  :: B(:)
          integer                        :: nparams
          double precision, allocatable  :: inp(:)
          contains
          procedure                      :: get_output
         end type sna

        contains

         subroutine get_output(this,inp,val)
         implicit none
         class(sna)                     :: this
         double precision, allocatable  :: inp(:)
         double precision               :: val
         integer                        :: i

          if(allocated(this%inp)) deallocate(this%inp)
          this%inp=inp
          
          val=0.0d0
          do i=1,this%nparams
           val=val+this%B(i)*inp(i)           
          enddo 

         return
         end subroutine get_output
        
        end module sna_class
