        module SNAP_FF_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use lammps_class
        use max_class
        use kind_class
        use parameters_class
        use lapack_inverse
        use potential_class
        implicit none

        type, extends(potential)                :: SNAP_FF
        real(kind=dbl), allocatable             :: beta(:)
        integer                                 :: tot_kinds
        integer                                 :: num_bisp

        contains
        
        procedure                        :: import_coeff => import_SNAP_coeff
        procedure                        :: get_fval  => get_SNAP_energy
        procedure                        :: get_fgrad => get_SNAP_force
        
        end type SNAP_FF
        
        contains

        subroutine import_SNAP_coeff(this)
        implicit none
        class(SNAP_FF)                                :: this
        integer                                       :: i

        allocate(this%beta(this%num_bisp*this%tot_kinds))
        open(222,file='snapcoeff_energy',action='read')
         do i=1,this%num_bisp*this%tot_kinds
          read(222,*) this%beta(i)
         end do
        close(222)

        end subroutine import_SNAP_coeff

        subroutine get_SNAP_energy(this,vec,val)
        implicit none
        class(SNAP_FF)                                  :: this
        real(kind=dbl)                                  :: val
        real(kind=dbl), allocatable                     :: vec(:)
        integer                                         :: i,j,tot_kinds

        val=0.0
        
        call this%frame%initialize()
        call this%frame%setup(this%frame%nkinds)
        call this%frame%get_desc()
        call this%frame%finalize()

        do i=1,this%frame%nats
         do j=1,this%num_bisp

          val=val+this%beta((this%frame%kind(i)-1)*this%num_bisp+j)*this%frame%at_desc(i)%desc(j)

         end do
         deallocate(this%frame%at_desc(i)%desc)
        end do

        deallocate(this%frame%at_desc)
        
        end subroutine get_SNAP_energy

        subroutine get_SNAP_force(this,vec,val,grad)
        implicit none
        class(SNAP_FF)                  :: this
        real(kind=dbl)                  :: val
        integer                         :: i,j,k,m
        real(kind=dbl), allocatable     :: vec(:),grad(:)

        call this%frame%initialize()
        call this%frame%setup(this%frame%nkinds)
        call this%frame%get_der_desc()
        call this%frame%finalize()

        allocate(grad(this%frame%nats*3))
        grad=0.0

        do i=1,this%frame%nats
         do j=1,this%frame%nkinds
          do m=1,3
           do k=1,this%num_bisp-1

            grad((i-1)*3+m)=grad((i-1)*3+m)+this%frame%der_at_desc(i)%desc((j-1)*3*(this%num_bisp-1)&
            +(m-1)*(this%num_bisp-1)+k)*this%beta((this%frame%kind(j)-1)*this%num_bisp+k+1)

           end do
          end do
         end do
         deallocate(this%frame%der_at_desc(i)%desc)
        end do       
         
        deallocate(this%frame%der_at_desc)
        end subroutine get_SNAP_force

        end module SNAP_FF_class
