        module minimizer_class
        use target_functions_class
        use potential_class
        implicit none

        type minimizer
         real(kind=dbl),allocatable             :: vec(:)
         double precision, allocatable          :: loc_lr(:)
         double precision                       :: max_grad=0.01e5
         double precision                       :: lr=0.0001d0
         integer                                :: max_iter=3000
         logical                                :: print_grad=.false.
         integer                                :: print_grad_io=22
         logical                                :: print_val=.false.
         integer                                :: print_val_io=23
         double precision                       :: eps=1.0e-7
         double precision                       :: beta1=0.9d0
         double precision                       :: beta2=0.999d0        
         type(potential_list),allocatable       :: pot(:)
         integer                                :: num_pot
         contains
         procedure :: init 
         procedure :: minimize => minimize_adam
         end type minimizer

         contains
         subroutine init(this,frame,loc_lr)
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
         end subroutine init

         subroutine minimize_adam(this,max_iter,start_iter)
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

           write(*,*) 'Grad Iter: ',iter+iter0,sqrt(gradnorm/size(sum_grad)),sum_en

           if(allocated(this%loc_lr))then
            this%vec=this%vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)*this%loc_lr
           else
            this%vec=this%vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)
           endif

           iter=iter+1
          
          do k=1,this%num_pot
           do i=1,this%pot(k)%item%frame%nats
            do j=1,3

             this%pot(k)%item%frame%x(i,j)=this%vec((i-1)*3+j)

            end do
           end do
          end do
          enddo

          if(this%print_grad) close(this%print_grad_io)
          if(this%print_val) close(this%print_val_io)
          if(allocated(gradres)) deallocate(gradres)
          if(allocated(gradres2)) deallocate(gradres2)
          if(allocated(sum_grad)) deallocate(sum_grad)
          if(allocated(grad)) deallocate(grad)
         end subroutine minimize_adam

         end module minimizer_class
