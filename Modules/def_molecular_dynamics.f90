        module md_class
        use SNAP_FF_class
        use lammps_class
        use rescale
        use target_functions_class
        implicit none

        type                    :: MD
         
        integer                 :: max_steps
        integer                 :: step_size
        character(len=3)        :: ensemble
        type(potential)         :: pot
        character

        end type md
        
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
         
        allocate(vec(this%potential%frame%nats*3))
       
        do while (this%iter < this%max_steps)
         call this%velocity_verlet
         if ((this%ensemble)=='nvt') call this%bussi
        end do

        end subroutine propagate

        subroutine velocity_verlet
        implicit none
        class(md), intent(inout) :: this
        integer                  :: i,j
        

        do i=1,this%potential%frame%nats
         do j=1,3

          vec((i-1)*3+j)=this%potential%frame%x(i,j)

         end do
        end do
 
        
        if (this%iter==1) then 
          
         call this%potential%get_fgrad(vec,this%potential%energy,this%potential%grad)
         
          do i=1,this%potential%frame%nats
           do j=1,3

            acc(i,j)=-this%potential%grad((i-1)*3+j)/this%mass(this%kind(i))

           end do
          end do

         end if

         this%potential%frame%v=this%potential%frame%v+0.5d0*acc*this%step_size
         this%potential%frame%x=this%potential%frame%x+this%potential%frame%v*this%step_size
         
         call this%potential%get_fgrad(posits,energy,forces)

         do i=1,this%object_lammps%nats
          do j=1,3

           this%acc(i,j)=forces((i-1)*3+j)/this%mass(this%kind(i))

          end do
         end do

         this%v=this%v+0.5d0*this%acc*this%step_size

        end subroutine velocity_verlet

        end module md_class
