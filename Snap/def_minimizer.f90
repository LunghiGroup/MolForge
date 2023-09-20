        module minimizer_class
        use target_functions_class
        use potential_class
        use SNAP_fit_class
        use trajectory_class
        implicit none
        
        type,extends(trajectory)                :: minimizer
         real(kind=dbl),allocatable             :: vec(:)
         real(kind=dbl), allocatable            :: loc_lr(:)
         real(kind=dbl)                         :: max_grad=0.01e5
         real(kind=dbl)                         :: lr=0.0001d0
         integer                                :: max_iter=3000
         logical                                :: print_grad=.false.
         integer                                :: print_grad_io=22
         logical                                :: print_val=.false.
         integer                                :: print_val_io=23
         real(kind=dbl)                         :: eps=1.0e-7
         real(kind=dbl)                         :: beta1=0.9d0
         real(kind=dbl)                         :: beta2=0.999d0        
         contains
         procedure :: init => init_minimizer
         procedure :: minimize => minimize_adam
         procedure :: print_info
         end type minimizer

         contains
         
         subroutine init_minimizer(this,frame,loc_lr)
         implicit none
         class(minimizer)                       :: this
         type(lammps_obj)                       :: frame
         integer                                :: i,j
         real(kind=dbl),optional,allocatable    :: loc_lr(:)
        
         allocate(this%vec(frame%nats*3))
         do i=1,frame%nats
          do j=1,3

           this%vec((i-1)*3+j)=frame%x(i,j)

          end do
         end do
         
         if(present(loc_lr)) this%loc_lr=loc_lr

         end subroutine init_minimizer

         subroutine minimize_adam(this,max_iter,start_iter,active_learning,delta)
         implicit none
         class(minimizer)              :: this
         integer                       :: iter,iter0
         integer                       :: i,j,k
         integer, optional             :: max_iter,start_iter
         double precision              :: gradnorm
         double precision, allocatable :: gradres(:),gradres2(:)
         real(kind=dbl),allocatable    :: sum_grad(:),grad(:)
         real(kind=dbl)                :: sum_en
         real(kind=dbl)                :: val
         real(kind=dbl),optional       :: delta
         real(kind=dbl)                :: error
         logical                       :: calc_sz_flag
         logical                       :: active_learning
         
          iter0=0
          if(present(start_iter)) iter0=start_iter
          if(present(max_iter)) this%max_iter=max_iter
          if(this%print_grad) open(this%print_grad_io,file='grad.dat')
          if(this%print_val) open(this%print_val_io,file='param.dat')

          if(allocated(gradres)) deallocate(gradres)
          if(allocated(gradres2)) deallocate(gradres2)
          if(allocated(sum_grad)) deallocate(sum_grad)
          if(allocated(grad)) deallocate(grad)
          allocate(gradres(size(this%vec)))
          allocate(gradres2(size(this%vec)))
          allocate(sum_grad(size(this%vec)))

          gradres=0.0d0
          gradres2=0.0d0

          iter=1

          do while (iter.le.this%max_iter)
           sum_grad=0.0d0
           sum_en=0.0d0
           val=0.0

           do i=1,size(this%pot)
            call this%pot(i)%item%get_fgrad(this%vec,val,grad)
            call this%pot(i)%item%get_fval(this%vec,val)
            sum_grad=sum_grad+grad
            sum_en=sum_en+val
            deallocate(grad)
           end do

           if(this%print_grad) write(this%print_grad_io,*) sum_grad

           gradnorm=0.0d0
           do i=1,size(sum_grad)
            gradres(i)=this%beta1*gradres(i)+(1-this%beta1)*sum_grad(i)
            gradres2(i)=this%beta2*gradres2(i)+(1-this%beta2)*sum_grad(i)**2
            gradnorm=gradnorm+sum_grad(i)**2
           enddo
           

           if(this%print_val) write(this%print_val_io,*) this%vec
           if(this%print_grad) write(this%print_grad_io,*) sum_grad

           if(allocated(this%loc_lr))then
            this%vec=this%vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)*this%loc_lr
           else
            this%vec=this%vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)
           endif

           !here we need to update the frame coordinates inside the potentials
           do k=1,this%num_pot
            do i=1,this%pot(k)%item%frame%nats
             do j=1,3
             
              this%pot(k)%item%frame%x(i,j)=this%vec((i-1)*3+j)

             end do
            end do
           end do
          
           if (active_learning) then
           
            if (.not.present(delta)) then
             write(*,*) "Insert value for delta"
            end if

            if (iter==1) then
             calc_sz_flag=.true.
            else
             calc_sz_flag=.false.
            end if
           
            call this%linear_fit%get_uncertainty(this%pot(1)%item%frame,calc_sz_flag,error)
            if (error > delta*this%linear_fit%s_z) then
             write(*,*)error,delta*this%linear_fit%s_z, "Error above the threshold"
             call this%print_info(this%pot(1)%item%frame)
             stop
            end if
           end if
           
           call this%print_info(this%pot(1)%item%frame)

           iter=iter+1
 
           enddo
          
          
           if(this%print_grad) close(this%print_grad_io)
           if(this%print_val) close(this%print_val_io)
           if(allocated(gradres)) deallocate(gradres)
           if(allocated(gradres2)) deallocate(gradres2)
           if(allocated(sum_grad)) deallocate(sum_grad)
           if(allocated(grad)) deallocate(grad)
           end subroutine minimize_adam
           
           subroutine print_info(this,frame)
           implicit none
           class(minimizer)                :: this
           type(lammps_obj)                :: frame
           integer                         :: i

           open(111, file="traj_min_molforge.xyz", action="write",position='append')

           write(111,*)frame%nats
           write(111,*)'XXX'

           do i=1,frame%nats
            write(111,*) frame%label(frame%kind(i)),frame%x(i,:)
           end do

           close(111)

           end subroutine print_info

           end module minimizer_class
