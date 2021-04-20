        module dist_class
        implicit none
   
         type dist1D
          double precision, allocatable       :: dist(:)
          double precision                    :: sigma
          double precision                    :: step
          integer                             :: nsteps
          integer                             :: nelem
          integer                             :: type_smear=1 ! 0=Lorentzian 1=Gaussian
          contains
          procedure                           :: alloc_dist
          procedure                           :: delete_dist
          procedure                           :: update_dist
          procedure                           :: norm_dist
          procedure                           :: get_val
         end type dist1D

        contains

          subroutine alloc_dist(this)
          implicit none
          class(dist1d)  :: this           
           allocate(this%dist(this%nsteps))
           this%dist=0.0d0
           this%nelem=0
          return
          end subroutine alloc_dist


          subroutine delete_dist(this)
          implicit none
          class(dist1d)  :: this            
           if(allocated(this%dist)) deallocate(this%dist)
          return
          end subroutine delete_dist


          subroutine update_dist(this,val,coeff)
          use units_parms
          implicit none
          class(dist1d)  :: this
          double precision :: val,coeff
          integer          :: k,l,v
          
           k=NINT(val/this%step)
           l=NINT(this%sigma/this%step)

           do v=-3*l,3*l
            if((v+k).gt.0 .and. (v+k).lt.this%nsteps)then
             this%dist(k+v)=this%dist(k+v)+coeff*deltaG(DBLE(v),DBLE(l))
            endif
           enddo
           
           this%nelem=this%nelem+1

          return
          end subroutine update_dist

          subroutine norm_dist(this)
          implicit none
          class(dist1d)  :: this
           this%dist=this%dist/this%nelem          
          return
          end subroutine norm_dist


          function get_val(this,ener) result(val)
          implicit none
          class(dist1D)        :: this
          double precision     :: val,ener
          integer              :: k

           val=0.0d0
           k=NINT(ener/this%step)
           if(k.le.0 .or. k.gt.this%nsteps) return
           val=this%dist(k)

          return
          end function get_val

        end module dist_class
