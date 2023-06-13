        module md_class
        use potential_class
        use target_functions_class
        implicit none

        type, extends(atoms_group) :: md
         
         integer :: max_steps
         integer :: step_size
         integer :: iter
         type(pot) :: potential
         contains
         procedure :: propagate
         procedure :: velocity_verlet
         procedure :: bussi

        end type md

        contains

        subroutine propagate
        implicit none
        class(md), intent(inout) :: this
        integer                  :: iter
         
         do while (this%iter < this%max_steps)

          call this%velocity_verlet
          call this%bussi

         end do

        end subroutine propagate

        subroutine velocity_verlet
        implicit none
        class(md), intent(inout) :: this
        double precision         :: energy
        integer                  :: i,j
         
         if (this%iter==1) then 
          
          call this%potential%get_fgrad(posits,energy,forces)
         
          do i=1,this%object_lammps%nats
           do j=1,3

            this%acc(i,j)=forces((i-1)*3+j)/this%mass(this%kind(i))

           end do
          end do

         end if

         !half-step velocities

         this%v=this%v+0.5d0*acc*this%step_size
         
         !
         !full step positions
         !

         this%x=this%x+this%v*this%step_size
         
         !

         call this%potential%get_fgrad(posits,energy,forces)

         do i=1,this%object_lammps%nats
          do j=1,3

           this%acc(i,j)=forces((i-1)*3+j)/this%mass(this%kind(i))

          end do
         end do

         !full-step velocities

         this%v=this%v+0.5d0*this%acc*this%step_size

        end subroutine velocity_verlet

        subroutine bussi

        end subroutine bussi

        end module md_class
