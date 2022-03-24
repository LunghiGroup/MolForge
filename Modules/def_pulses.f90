        module pulses_class
        use blacs_utils
        use lists_class
        implicit none

        type :: general_pulse
         logical                          :: sx
         logical                          :: dx
         double precision                 :: weight
         double precision                 :: n(3)
         integer                          :: spin
         double precision                 :: beta
         type(dist_cmplx_mat)             :: rot
         contains
         procedure              :: bcast => pulse_bcast         
        end type general_pulse

        type, extends(list) :: general_pulse_list
        contains
        procedure             :: add_node => add_pulse_node
        procedure             :: rd_node  => rd_pulse_node
        end type general_pulse_list

        contains

        subroutine pulse_bcast(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(general_pulse)    :: this
        integer                 :: nspins

         call mpi_bcast(this%sx,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(this%dx,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(this%weight,1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(this%n(1),1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(this%n(2),1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(this%n(3),1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(this%spin,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%beta,1,mpi_double_precision,0,mpi_comm_world,err)

        return
        end subroutine pulse_bcast

        subroutine set_pulse(pulse,type_pulse,spin)
        implicit none
        type(general_pulse),allocatable     :: pulse(:)
        integer                             :: spin
        integer                             :: i
        character(len=10)                   :: type_pulse

         if(allocated(pulse)) deallocate(pulse)

         select case (type_pulse)

          case ('PI')
           allocate(pulse(1))
           pulse(1)%dx=.true.  
           pulse(1)%sx=.true.  
           pulse(1)%weight=1.0d0
           pulse(1)%n(1)=0.0d0
           pulse(1)%n(2)=1.0d0
           pulse(1)%n(3)=0.0d0
           pulse(1)%spin=spin
           pulse(1)%beta=acos(-1.0d0)

          case ('PI2')
           allocate(pulse(1))
           pulse(1)%dx=.true.  
           pulse(1)%sx=.true.  
           pulse(1)%weight=1.0d0
           pulse(1)%spin=spin
           pulse(1)%n(1)=0.0d0
           pulse(1)%n(2)=1.0d0
           pulse(1)%n(3)=0.0d0
           pulse(1)%beta=acos(-1.0d0)/2.0d0

         end select

        return
        end subroutine set_pulse

        subroutine add_pulse_node(this_list,val)
        implicit none
        class(general_pulse_list)      :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: nspins

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (general_pulse)
            allocate(general_pulse::this_list%head%key)
            allocate(general_pulse::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (general_pulse)
            allocate(general_pulse::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (general_pulse)
           select type (arrow=>tmp_node%key)
            type is (general_pulse)
            arrow%dx=val%dx
            arrow%sx=val%sx
            arrow%weight=val%weight
            arrow%n=val%n
            arrow%spin=val%spin
            arrow%beta=val%beta
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_pulse_node

        subroutine rd_pulse_node(this,pulse)
        implicit none
        class(general_pulse_list)  :: this
        class(general_pulse)       :: pulse
        class(*),pointer           :: bho
        integer                    :: nspins

         select type (bho=>this%node%key)
          type is (general_pulse)
            pulse%dx=bho%dx
            pulse%sx=bho%sx
            pulse%weight=bho%weight
            pulse%spin=bho%spin
            pulse%beta=bho%beta
            pulse%n=bho%n
         end select

        return
        end subroutine rd_pulse_node



        end module pulses_class
