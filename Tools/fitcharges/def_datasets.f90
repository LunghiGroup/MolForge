        module fit_lammps_class
        implicit none

        TYPE data_set
         integer                           :: frames
         integer                           :: nats
         character (len=100)               :: inp_data
         character (len=100)               :: inp_ener
         character (len=100)               :: inp_forces
         character (len=100)               :: inp_traj
         double precision, allocatable     :: ener(:)
         character (len=5),allocatable     :: label(:)
         double precision, allocatable     :: x(:,:)
         double precision, allocatable     :: fx(:,:)
         double precision, allocatable     :: fy(:,:)
         double precision, allocatable     :: fz(:,:)
         double precision                  :: weight
        end TYPE data_set

        TYPE system
         type(data_file), allocatable   :: data(:)
         integer                        :: ndata
         character (len=100)            :: inp
         character (len=100)            :: inp_fix
         character (len=100)            :: inp_fit
         integer                        :: tot_frames
         integer                        :: npar2fit
         contains
         procedure                      :: read_sys
!         procedure                      :: delete
        END TYPE system

        TYPE kernel_kind
         integer                        :: nenvs=0
         double precision               :: sigma=0
         double precision, allocatable  :: B(:,:)
         contains
         procedure                      :: dealloc => dealloc_kernel_kind
        END TYPE kernel_kind

        TYPE kernel_global
         integer                        :: nkinds=0
         type(kernel_kind), allocatable :: K(:)
         contains
         procedure                      :: dealloc => dealloc_kernel_global
        END TYPE kernel_global

        
        contains


        subroutine dealloc_kernel_global(this)
        implicit none
        class(kernel_global)    :: this
        integer                 ::i
         if( allocated(this%K) ) then 
          do i=1,this%nkinds
           call this%K(i)%dealloc()
          enddo
          deallocate(this%K)
         endif
         this%nkinds=0
        return
        end subroutine dealloc_kernel_global

        subroutine dealloc_kernel_kind(this)
        implicit none
        class(kernel_kind)    :: this
         if( allocated(this%B) ) deallocate(this%B)
         this%nenvs=0
         this%sigma=0
        return
        end subroutine dealloc_kernel_kind

        subroutine read_sys(sys,data_file,fit_ener,fit_forces)
        implicit none
        class(system)           :: sys
        integer                 :: i,l,v,k,m,n
        character(len=100)      :: label,data_file
        logical                 :: fit_forces,fit_ener

         if(allocated(sys%data))then
          do i=1,sys%ndata 
           if(allocated(sys%data(i)%x)) deallocate(sys%data(i)%x)
           if(allocated(sys%data(i)%ener)) deallocate(sys%data(i)%ener)
           if(allocated(sys%data(i)%fx)) deallocate(sys%data(i)%fx)
           if(allocated(sys%data(i)%fy)) deallocate(sys%data(i)%fy)
           if(allocated(sys%data(i)%fz)) deallocate(sys%data(i)%fz)
          enddo 
          deallocate(sys%data)
         endif   

         open(8,file=trim(data_file))

         sys%tot_frames=0
         sys%npar2fit=0

         read(8,*) sys%ndata    

         allocate(sys%data(sys%ndata))

         do i=1,sys%ndata

          if(fit_forces .and. fit_ener)then
           read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj,sys%data(i)%inp_ener,&
                     sys%data(i)%inp_forces,sys%data(i)%weight
          endif

          if(fit_ener .and.  (.not. fit_forces) )then
           read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj, &
                     sys%data(i)%inp_ener,sys%data(i)%weight
          endif

          if((.not. fit_ener) .and. fit_forces )then
           read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj, &
                     sys%data(i)%inp_forces,sys%data(i)%weight
          endif

          open(12,file=sys%data(i)%inp_traj)

          if(fit_ener)then
           allocate(sys%data(i)%ener(sys%data(i)%frames))
           open(13,file=sys%data(i)%inp_ener)
          endif

          if(fit_forces)then
           open(14,file=sys%data(i)%inp_forces)
          endif

          do l=1,sys%data(i)%frames

           sys%tot_frames=sys%tot_frames+1

           read(12,*) sys%data(i)%nats
           read(12,*) 
           if(fit_forces)then
            read(14,*) 
            read(14,*)  
           endif

           if(.not.allocated(sys%data(i)%label))then
            allocate(sys%data(i)%label(sys%data(i)%nats))
           endif
           if(.not.allocated(sys%data(i)%x))then
            allocate(sys%data(i)%x(sys%data(i)%frames,3*sys%data(i)%nats))
           endif
           if(.not.allocated(sys%data(i)%fx) .and. fit_forces)then
            allocate(sys%data(i)%fx(sys%data(i)%frames,sys%data(i)%nats))
           endif
           if(.not.allocated(sys%data(i)%fy) .and. fit_forces)then
            allocate(sys%data(i)%fy(sys%data(i)%frames,sys%data(i)%nats))
           endif
           if(.not.allocated(sys%data(i)%fz) .and. fit_forces)then
            allocate(sys%data(i)%fz(sys%data(i)%frames,sys%data(i)%nats))
           endif

           v=1
           do k=1,sys%data(i)%nats             

            read(12,*) sys%data(i)%label(k),sys%data(i)%x(l,v),sys%data(i)%x(l,v+1),sys%data(i)%x(l,v+2)

            if(fit_forces)then
             read(14,*) sys%data(i)%fx(l,k),sys%data(i)%fy(l,k),sys%data(i)%fz(l,k)
             sys%npar2fit=sys%npar2fit+3
            endif
            v=v+3

           enddo

           if(fit_ener)then
            read(13,*) sys%data(i)%ener(l)
            sys%npar2fit=sys%npar2fit+1
           endif

          enddo

          close(12)
          if(fit_ener) close(13)
          if(fit_forces) close(14)

         enddo

         close(8)

        return
        end subroutine read_sys


        end module fit_lammps_class
