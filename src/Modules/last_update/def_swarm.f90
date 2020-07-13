        module particles_swarm_class
        use target_functions_class
        implicit none

         type particle 
          integer                        :: nval
          double precision               :: ener
          double precision               :: omega=0.6d0
          double precision               :: phi_1=1.70d0
          double precision               :: phi_2=1.70d0   
          double precision               :: lrate=0.99d0
          double precision, allocatable  :: res(:)
          double precision, allocatable  :: val(:)
          double precision, allocatable  :: vel(:)
          double precision, allocatable  :: step(:)
          double precision, allocatable  :: max_val(:)
          double precision, allocatable  :: min_val(:)
          double precision, allocatable  :: best_loc_val(:)
          double precision, allocatable  :: best_glo_val(:)
          double precision, allocatable  :: best_loc_ener
          double precision, allocatable  :: best_glo_ener
          class(target_function), pointer :: target_f => null()
          contains 
          procedure :: init_par
          procedure :: init_val
          procedure :: init_vel
          procedure :: get_ener
          procedure :: propagate_par
          procedure :: update_par
         end type particle

         type particles_swarm
          integer                     :: iter
          integer                     :: npar=20
          integer                     :: best_par
          integer                     :: max_iter=1000
          type(particle), allocatable :: par(:)
          class(target_function), pointer :: target_f => null()
          logical                     :: print_val=.false.
          integer                     :: print_val_io=23
          character(len=100)        :: print_val_name='param_swarm.dat'
          contains
          procedure :: init_swarm
          procedure :: propagate_swarm
          procedure :: update_swarm
          procedure :: minimize
          procedure :: get_best_val
          procedure :: map_target_f
          procedure :: release_target_f
         end type particles_swarm
   
         contains

         subroutine map_target_f(this)
         implicit none
         class(particles_swarm)        :: this
         integer                       :: i
          do i=1,this%npar
           this%par(i)%target_f => this%target_f
          enddo
         return
         end subroutine map_target_f

         subroutine release_target_f(this)
         implicit none
         class(particles_swarm)        :: this
         integer                       :: i
          do i=1,this%npar
           this%par(i)%target_f => null()
          enddo
          this%target_f => null()
         return
         end subroutine release_target_f

         subroutine get_best_val(this,vec)
         implicit none
         class(particles_swarm)        :: this
         double precision, allocatable :: vec(:)

          vec=this%par(this%best_par)%val

         return
         end subroutine get_best_val

         subroutine minimize(this,max_iter)
         implicit none
         class(particles_swarm)    :: this
         integer                   :: i
         integer, optional         :: max_iter
          
          this%iter=0
          if(present(max_iter)) this%max_iter=max_iter

          if(this%print_val) &
                open(this%print_val_io,file=this%print_val_name)

          do while (this%iter.lt.this%max_iter)

           this%iter=this%iter+1
           call this%propagate_swarm()

           if(this%print_val) write(this%print_val_io,*) this%par(this%best_par)%val

           write(*,*) 'Swarm Iter: ',this%iter,this%best_par,this%par(1)%best_glo_ener

          enddo

          if(this%print_val) close(this%print_val_io)

!          write(*,*) '### Final Solution: ',this%par(this%best_par)%ener
!          do i=1,this%par(this%best_par)%nval
!           write(*,*) '### Optimized Vals: ',this%par(this%best_par)%val(i)
!          enddo

         return
         end subroutine minimize


         subroutine init_swarm(this,npar,nval,step,min_val,max_val)
         use random_numbers_class
         implicit none
         class(particles_swarm)                  :: this
         integer                                 :: nval,i,j
         integer, optional                       :: npar
         double precision, allocatable, optional :: step(:)
         double precision, allocatable           :: step1(:)
         double precision, allocatable, optional :: min_val(:),max_val(:)
         double precision, allocatable           :: min_val1(:),max_val1(:)

          call init_random_seed()

          if(present(npar)) this%npar=npar

          if(.not.present(step))then
           allocate(step1(nval))
           step1=0.6d0
          else
           step1 = step
          endif
          if(.not.present(min_val))then
           allocate(min_val1(nval))
           min_val1=-5.0d-1
          else
           min_val1 = min_val
          endif
          if(.not.present(max_val))then
           allocate(max_val1(nval))
           max_val1=5.0d-1
          else
           max_val1 = max_val
          endif

          allocate(this%par(npar))

          call this%map_target_f()    

          do i=1,npar
           this%par(i)%nval=nval
           call this%par(i)%init_par(nval,step1,min_val1,max_val1)
          enddo

          this%best_par=1
          do j=1,this%npar
           this%par(j)%best_glo_ener=this%par(1)%ener
           this%par(j)%best_glo_val=this%par(1)%val
          enddo
          do i=2,this%npar
           if(this%par(i)%ener.lt.this%par(i)%best_glo_ener)then
            this%best_par=i
            do j=1,this%npar
             this%par(j)%best_glo_ener=this%par(i)%ener
             this%par(j)%best_glo_val=this%par(i)%val
            enddo
           endif
          enddo

         return
         end subroutine init_swarm


         subroutine update_swarm(this)
         implicit none
         class(particles_swarm)    :: this
         integer         :: i,j

          do i=1,this%npar
           if(this%par(i)%ener.lt.this%par(i)%best_glo_ener)then
            this%best_par=i
            do j=1,this%npar
             this%par(j)%best_glo_ener=this%par(i)%ener
             this%par(j)%best_glo_val=this%par(i)%val
            enddo
           endif
          enddo
          
         return
         end subroutine update_swarm


         subroutine propagate_swarm(this)
         implicit none
         class(particles_swarm)    :: this
         integer         :: i

          do i=1,this%npar
           call this%par(i)%propagate_par()
          enddo
          call this%update_swarm()          
          
         return
         end subroutine propagate_swarm


         subroutine init_par(this,nval,step,min_val,max_val)
         use random_numbers_class
         implicit none
         class(particle)    :: this
         integer            :: nval,i
         double precision, allocatable :: val0(:),step(:)
         double precision, allocatable :: min_val(:),max_val(:)

          this%nval=nval
          allocate(this%val(nval))
          allocate(this%vel(nval))
          allocate(this%step(nval))
          allocate(this%res(nval))
          allocate(this%min_val(nval))
          allocate(this%max_val(nval))
          allocate(this%best_loc_val(nval))
          allocate(this%best_glo_val(nval))
          
          this%res=1.0d0
          this%step=step
          this%min_val=min_val
          this%max_val=max_val

          call this%init_val()
          call this%init_vel()
          call this%get_ener() 

          this%best_loc_val=this%val
          this%best_loc_ener=this%ener
          
         return
         end subroutine init_par


         subroutine init_val(this)
         use random_numbers_class
         implicit none
         class(particle)    :: this
         double precision   :: rand_num
         integer            :: i

          do i=1,this%nval
           call random_number(rand_num)              
           this%val=rand_num*(this%max_val-this%min_val)+this%min_val
          enddo

         return
         end subroutine init_val


         subroutine init_vel(this)
         use random_numbers_class
         implicit none
         class(particle)    :: this
         double precision   :: rand_num

          this%vel=0.0d0

         return
         end subroutine init_vel         


         subroutine propagate_par(this)
         use random_numbers_class
         implicit none
         class(particle)    :: this
         double precision   :: rand_num,phi1,phi2
         integer            :: i

          do i=1,this%nval
           
           call random_number(rand_num)              
           phi1=this%phi_1*rand_num
           call random_number(rand_num)              
           phi2=this%phi_2*rand_num

           this%vel(i)=this%omega*this%vel(i)                 &
                        +phi1*(this%best_glo_val(i)-this%val(i))  &
                        +phi2*(this%best_loc_val(i)-this%val(i))     


           if(this%vel(i).gt.this%step(i))then
             this%vel(i)=this%step(i)
           endif

           if(this%vel(i).lt.-this%step(i))then
             this%vel(i)=-this%step(i)
           endif

           this%val(i)=this%val(i)+this%vel(i)

!           if(this%val(i).ge.this%max_val(i))then
!             call random_number(rand_num)              
!             this%val(i)=this%max_val(i)-0.01*rand_num*this%step(i)
!           endif

!           if(this%val(i).le.this%min_val(i))then
!             call random_number(rand_num)              
!             this%val(i)=this%min_val(i)+0.01*rand_num*this%step(i)
!           endif

          enddo

          call this%update_par()

         return
         end subroutine propagate_par


         subroutine update_par(this)
         use random_numbers_class
         implicit none
         class(particle)    :: this

          call this%get_ener() 

          if(this%ener.lt.this%best_loc_ener)then           
           this%best_loc_val=this%val
           this%best_loc_ener=this%ener
          endif

         return
         end subroutine update_par

         subroutine get_ener(this)
         implicit none
         class(particle)      :: this

          call this%target_f%get_fval(this%val,this%ener)

         return
         end subroutine get_ener

        end module particles_swarm_class
