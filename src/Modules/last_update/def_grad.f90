        module gradmin_class
        use target_functions_class
        implicit none

         type gradient_descent
          integer                        :: nval
          double precision               :: ener          
          double precision, allocatable  :: val(:)
          double precision, allocatable  :: grad(:)
          double precision               :: max_grad=0.01e5
          double precision               :: lr=0.0001d0
          integer                        :: max_iter=3000
          logical                        :: print_grad=.false.
          integer                        :: print_grad_io=22
          logical                        :: print_val=.false.
          integer                        :: print_val_io=23
          class(target_function), pointer :: target_f => null()
          contains
          procedure :: init 
          procedure :: minimize => minimize_gd
          procedure :: get_ener 
          procedure :: get_grad
          procedure :: release_target_f
         end type gradient_descent

         type, extends(gradient_descent) :: adam
          double precision               :: eps=1.0e-7
          double precision               :: beta1=0.9d0
          double precision               :: beta2=0.999d0
          contains
          procedure :: minimize => minimize_adam
         end type adam

         contains

         subroutine release_target_f(this)
         implicit none
         class(gradient_descent)       :: this
         integer                       :: i
          this%target_f => null()
         return
         end subroutine release_target_f

         subroutine init(this,nval,vec)
         use random_numbers_class
         implicit none
         class(gradient_descent)                   :: this
         double precision              :: rand_num
         double precision, allocatable, optional :: vec(:)
         integer                       :: nval,i

          this%nval=nval
          allocate(this%val(this%nval))
          allocate(this%grad(this%nval))
          if(present(vec))then
           this%val=vec
          else
           do i=1,this%nval
            call random_number(rand_num) 
            this%val(i)=rand_num
           enddo
          endif
          this%grad=0.0d0

         return
         end subroutine init

         subroutine minimize_gd(this,max_iter)
         implicit none
         class(gradient_descent)       :: this
         integer                       :: iter,i
         integer, optional             :: max_iter
         double precision              :: gradnorm

          if(present(max_iter)) this%max_iter=max_iter
          if(this%print_grad) open(this%print_grad_io,file='grad.dat')
          if(this%print_val) open(this%print_val_io,file='param.dat')

          iter=1

          do while (iter.le.this%max_iter)          

           call this%target_f%get_fgrad(this%val,this%ener,this%grad) 

           gradnorm=0.0d0
           do i=1,size(this%grad)
!            if(this%grad(i).gt.this%max_grad) this%grad(i)=this%max_grad
!            if(this%grad(i).lt.this%max_grad) this%grad(i)=-this%max_grad
            gradnorm=gradnorm+this%grad(i)**2
           enddo

           if(this%print_grad) write(this%print_grad_io,*) this%grad
           if(this%print_val) write(this%print_val_io,*) this%val

           write(*,*) 'Grad Iter: ',iter,sqrt(gradnorm/size(this%grad)),this%ener

           this%val=this%val-this%lr*this%grad

           iter=iter+1
          enddo

          if(this%print_grad) close(this%print_grad_io)
          if(this%print_val) close(this%print_val_io)

         return
         end subroutine minimize_gd

         subroutine minimize_adam(this,max_iter)
         implicit none
         class(adam)                   :: this
         integer                       :: iter,i
         integer, optional             :: max_iter
         double precision              :: gradnorm
         double precision, allocatable :: gradres(:),gradres2(:)

          if(present(max_iter)) this%max_iter=max_iter
          if(this%print_grad) open(this%print_grad_io,file='grad.dat')
          if(this%print_val) open(this%print_val_io,file='param.dat')

          if(allocated(gradres)) deallocate(gradres)
          if(allocated(gradres2)) deallocate(gradres2)
          allocate(gradres(size(this%val)))
          allocate(gradres2(size(this%val)))
          gradres=0.0d0
          gradres2=0.0d0
          iter=1

          do while (iter.le.this%max_iter)          
           
           call this%target_f%get_fgrad(this%val,this%ener,this%grad)

           if(this%print_grad) write(this%print_grad_io,*) this%grad

           gradnorm=0.0d0
           do i=1,size(this%grad)
            gradres(i)=this%beta1*gradres(i)+(1-this%beta1)*this%grad(i)
            gradres2(i)=this%beta2*gradres2(i)+(1-this%beta2)*this%grad(i)**2
            gradnorm=gradnorm+this%grad(i)**2
           enddo

           if(this%print_val) write(this%print_val_io,*) this%val
           if(this%print_grad) write(this%print_grad_io,*) this%grad

           write(*,*) 'Grad Iter: ',iter,sqrt(gradnorm/size(this%grad)),this%ener

           this%val=this%val-this%lr*(gradres/(1-this%beta1**(iter+1)))& 
                /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)

           iter=iter+1
          enddo

          if(this%print_grad) close(this%print_grad_io)
          if(this%print_val) close(this%print_val_io)

         return
         end subroutine minimize_adam

         subroutine get_ener(this)
         implicit none
         class(gradient_descent)         :: this         
          call this%target_f%get_fval(this%val,this%ener)
         return
         end subroutine get_ener

         subroutine get_grad(this)
         implicit none
         class(gradient_descent)       :: this         
          call this%target_f%get_fgrad(this%val,this%ener,this%grad)
         return
         end subroutine get_grad

        end module gradmin_class
