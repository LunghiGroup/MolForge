        module md_class
        use SNAP_FF_class
        use lammps_class
        use potential_class
        use SNAP_FF_class
        use VdW_class
        !use rescale
        use target_functions_class
        implicit none

        type                                    :: MD
         
        type(potential_list),allocatable        :: pot(:)
        integer                                 :: max_steps
        integer                                 :: step_size
        integer                                 :: num_pot
        integer                                 :: iter
        character(len=3)                        :: ensemble
        contains
        procedure                               :: initialize_vel        
        procedure                               :: link_potentials
        !procedure :: propagate
        !procedure :: velocity_verlet
        !procedure :: bussi
        
        end type MD

        contains
        
        subroutine initialize_vel(this,frame,iseed,T_in)
        implicit none
        class(md),intent(inout)                    :: this
        real(kind=dbl), intent(in)                 :: T_in
        type(lammps_obj)                           :: frame
        integer                                    :: idist=3
        integer                                    :: N=3
        integer                                    :: i
        integer,dimension(4)                       :: iseed
        real(kind=dbl),dimension(3)                :: X
        real(kind=dbl),dimension(3)                :: lin_mom_tot
        real(kind=dbl)                             :: E_kin
        real(kind=dbl)                             :: T

        allocate(frame%v(frame%nats,3))
        
        do i=1,frame%nats
         call dlarnv(idist,iseed,N,X)
         frame%v(i,:)=dsqrt(boltz*T_in/(frame%mass(frame%kind(i))))*X
        end do

        E_kin=0.0d0
        lin_mom_tot=0.0d0

        do i=1,frame%nats
         lin_mom_tot=lin_mom_tot + frame%v(i,:)*frame%mass(frame%kind(i))
        end do

        do i=1,frame%nats
         frame%v(i,:)=frame%v(i,:)-(lin_mom_tot/(dble(frame%nats)*frame%mass(frame%kind(i))))
        end do

        do i=1,frame%nats
         E_kin=E_kin+0.5d0*sum((frame%v(i,:)**2))*(frame%mass(frame%kind(i)))
        end do

        T=(2.0d0*E_kin)/((3.0d0*frame%nats-3.0d0)*boltz)
        frame%v(:,:)=dsqrt(T_in/T)*frame%v(:,:)

        end subroutine initialize_vel

        subroutine link_potentials(this,keyword,SNAP,VdW)
        implicit none
        class(md)                       :: this
        character(len=*)                :: keyword
        type(SNAP_FF),target            :: SNAP
        type(VdW_FF),optional,target    :: VdW
        integer                         :: i

        allocate(this%pot(this%num_pot))
        
        if (keyword=="SNAP") then
         this%pot(1)%item=>SNAP
         do i=2,this%num_pot
          this%pot(i)%item=>NULL()
         end do      
        else if ((keyword=="SNAP_VdW").and.(present(VdW))) then 
         this%pot(1)%item=>SNAP
         this%pot(2)%item=>VdW
         do i=3,this%num_pot
          this%pot(i)%item=>NULL()
         end do
        
        end if
        end subroutine link_potentials

        !subroutine propagate(this,frame)
        !implicit none
        !class(md), intent(inout) :: this
        !type(lammps_obj)         :: frame
         
        !do while (this%iter < this%max_steps)
        ! call this%velocity_verlet(frame)
         !if ((this%ensemble)=='nvt') call this%bussi
        !end do

        !end subroutine propagate

        !subroutine velocity_verlet(this,frame)
        !implicit none
        !class(md), intent(inout) :: this
        !integer                  :: i,j
        
        !if (this%iter==1) then 
         ! do j=1,this%num_pot
         !call this%potential(i)%get_fgrad(vec,this%potential%energy,this%potential%grad)
          !end do

         ! do i=1,this%potential%frame%nats
         !  do j=1,3

         !   acc(i,j)=-this%potential%grad((i-1)*3+j)/this%mass(this%kind(i))

         !  end do
         ! end do

        ! end if

        ! this%potential%frame%v=this%potential%frame%v+0.5d0*acc*this%step_size
        ! this%potential%frame%x=this%potential%frame%x+this%potential%frame%v*this%step_size
         
        ! call this%potential%get_fgrad(posits,energy,forces)

        ! do i=1,this%object_lammps%nats
        !  do j=1,3

         !  this%acc(i,j)=forces((i-1)*3+j)/this%mass(this%kind(i))

         ! end do
         !end do

         !this%v=this%v+0.5d0*this%acc*this%step_size

        !end subroutine velocity_verlet

        end module md_class
