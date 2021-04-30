        module dist_class
        implicit none
   
         type dist1D
          double precision, allocatable       :: dist(:)
          double precision                    :: sigma
          double precision                    :: step
          double precision                    :: shift=0.00
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

         type dist2D
          double precision, allocatable       :: dist(:,:)
          double precision                    :: sigma1,sigma2
          double precision                    :: step1,step2
          double precision                    :: shift1=0.00
          double precision                    :: shift2=0.00
          integer                             :: nsteps1,nsteps2
          integer                             :: nelem
          integer                             :: type_smear=1 ! 0=Lorentzian 1=Gaussian
          contains
          procedure                           :: alloc_dist => alloc_dist_2
          procedure                           :: delete_dist => delete_dist_2
          procedure                           :: update_dist => update_dist_2
          procedure                           :: norm_dist => norm_dist_2
          procedure                           :: get_val => get_val_2
         end type dist2D

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

          subroutine alloc_dist_2(this)
          implicit none
          class(dist2d)  :: this           
           allocate(this%dist(this%nsteps1,this%nsteps2))
           this%dist=0.0d0
           this%nelem=0
          return
          end subroutine alloc_dist_2

          subroutine delete_dist_2(this)
          implicit none
          class(dist2d)  :: this            
           if(allocated(this%dist)) deallocate(this%dist)
          return
          end subroutine delete_dist_2


          subroutine update_dist(this,val,coeff)
          use units_parms
          implicit none
          class(dist1d)  :: this
          double precision :: val,coeff
          integer          :: k,l,v
          
           val=val-this%shift

           k=NINT(val/this%step)
           l=NINT(this%sigma/this%step)

           do v=-3*l,3*l
            if((v+k).gt.0 .and. (v+k).lt.this%nsteps)then
             this%dist(k+v)=this%dist(k+v)+coeff*delta(DBLE(v),DBLE(l))
            endif
           enddo
           
           this%nelem=this%nelem+1

          return
          end subroutine update_dist

          subroutine update_dist_2(this,val1,val2,coeff)
          use units_parms
          implicit none
          class(dist2d)    :: this
          double precision :: val1,val2,coeff
          integer          :: k1,l1,k2,l2,v1,v2
          
           val1=val1-this%shift1
           k1=NINT(val1/this%step1)
           l1=NINT(this%sigma1/this%step1)

           val2=val2-this%shift2
           k2=NINT(val2/this%step2)
           l2=NINT(this%sigma2/this%step2)

           do v1=-3*l1,3*l1
            if((v1+k1).gt.0 .and. (v1+k1).lt.this%nsteps1)then
             do v2=-3*l2,3*l2
              if((v2+k2).gt.0 .and. (v2+k2).lt.this%nsteps2)then
               this%dist(k1+v1,k2+v2)=this%dist(k1+v1,k2+v2)+coeff*delta(DBLE(v1),DBLE(l1))*delta(DBLE(v2),DBLE(l2))
              endif
             enddo
            endif
           enddo
           
           this%nelem=this%nelem+1

          return
          end subroutine update_dist_2


          subroutine norm_dist(this)
          implicit none
          class(dist1d)  :: this
           this%dist=this%dist/this%nelem          
          return
          end subroutine norm_dist

          subroutine norm_dist_2(this)
          implicit none
          class(dist2d)  :: this
           this%dist=this%dist/this%nelem          
          return
          end subroutine norm_dist_2

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

          function get_val_2(this,ener1,ener2) result(val)
          implicit none
          class(dist2D)        :: this
          double precision     :: val,ener1,ener2
          integer              :: k1,k2

           val=0.0d0
           k1=NINT(ener1/this%step1)
           k2=NINT(ener2/this%step2)
           if(k1.le.0 .or. k1.gt.this%nsteps1) return
           if(k2.le.0 .or. k2.gt.this%nsteps2) return
           val=this%dist(k1,k2)

          return
          end function get_val_2

        end module dist_class
