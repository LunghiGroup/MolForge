        module nets_class
        implicit none

         type neuron
          integer                        :: nweights
          double precision, allocatable  :: weights(:)
          double precision               :: bias
          double precision               :: output
          double precision               :: grad_act
          character(len=10)              :: activ='relu'
          contains 
          procedure                      :: get_activity
         end type neuron

         type layer
          integer                        :: nneu
          type(neuron), allocatable      :: neu(:)
          contains
          procedure                      :: set_nodes
         end type layer

         type net
          integer                        :: nlayers     !! tot layers, including output one
          integer                        :: nneu
          integer                        :: nparams
          integer                        :: ninput
          integer                        :: noutput
          type(layer), allocatable       :: layers(:) 
          double precision, allocatable  :: inp(:)
          double precision, allocatable  :: gradinp(:,:)
          double precision, allocatable  :: mean_inp(:),sigma_inp(:)
          logical                        :: norm_inp=.false.
          contains
          procedure :: get_output
          procedure :: get_grad_nn
          procedure :: get_grad_weights
          procedure :: get_grad_bias
          procedure :: get_gradinp
          procedure :: backprop
          procedure :: set_topology
          procedure :: set_parameters
          procedure :: get_parameters
         end type net

         contains

         subroutine get_output(this,inps,outs)
         implicit none
         class(net)                    :: this
         double precision, allocatable :: outs(:),inps(:)
         integer                       :: i,j

          if(allocated(this%inp)) deallocate(this%inp)
          allocate(this%inp(this%ninput))

          if(this%norm_inp) inps=(inps-this%mean_inp)/this%sigma_inp
          this%inp=inps

          do i=1,this%nlayers

           if(i.gt.1)then
            if(allocated(inps)) deallocate(inps)
            allocate(inps(this%layers(i-1)%nneu))
            inps=outs
           endif        

           if(allocated(outs)) deallocate(outs)          
           allocate(outs(this%layers(i)%nneu))

           do j=1,this%layers(i)%nneu
            call this%layers(i)%neu(j)%get_activity(inps,outs(j))
           enddo

          enddo

         return
         end subroutine get_output

         subroutine get_gradinp(this)
         implicit none
         class(net)                    :: this
         double precision, allocatable :: Wmat(:)
         double precision, allocatable :: Wmat2(:,:),Wmat3(:)
         integer                       :: i,j,k,v,l,s,t
        
          if(allocated(this%gradinp)) deallocate(this%gradinp)
          allocate(this%gradinp(this%ninput,this%noutput))
          this%gradinp=0.0d0

          do l=this%noutput,1,-1

           if(this%nlayers.gt.1)then

            if(allocated(Wmat)) deallocate(Wmat)
            allocate(Wmat(this%layers(this%nlayers-1)%nneu))
            Wmat=0.0d0

            do j=this%layers(this%nlayers-1)%nneu,1,-1
             Wmat(j)=this%layers(this%nlayers)%neu(l)%weights(j)*this%layers(this%nlayers)%neu(l)%grad_act  
            enddo

           endif

           do i=this%nlayers-1,1,-1

            if(i.gt.1)then
             if(allocated(Wmat2)) deallocate(Wmat2)
             allocate(Wmat2(this%layers(i)%nneu,this%layers(i-1)%nneu))
             Wmat2=0.0d0
            endif

            do j=this%layers(i)%nneu,1,-1

             if(i.eq.1)then
              do k=this%ninput,1,-1
               this%gradinp(k,l)=Wmat(j)*this%layers(i)%neu(j)%weights(k)*this%layers(i)%neu(j)%grad_act
              enddo
             else
              do t=this%layers(i-1)%nneu,1,-1
               Wmat2(j,t)=this%layers(i)%neu(j)%weights(t)*this%layers(i)%neu(j)%grad_act
              enddo
             endif

            enddo ! j nodes

            if(i.gt.1)then
             if(allocated(Wmat3)) deallocate(Wmat3)
             allocate(Wmat3(this%layers(i-1)%nneu))
             Wmat3=0.0d0
             Wmat3=matmul(Wmat,Wmat2)
             deallocate(Wmat) 
             allocate(Wmat(this%layers(i-1)%nneu))
             Wmat=Wmat3            
            endif

           enddo ! i innerlayers                    

          enddo ! l outputs

          if(allocated(Wmat)) deallocate(Wmat)
          if(allocated(Wmat2)) deallocate(Wmat2)
          if(allocated(Wmat3)) deallocate(Wmat3)

         return
         end subroutine get_gradinp

         subroutine backprop(this,grad)
         implicit none
         class(net)                    :: this
         double precision, allocatable :: grad(:,:),Wmat(:)
         double precision, allocatable :: Wmat2(:,:),Wmat3(:)
         integer                       :: i,j,k,v,l,s,t
        
          if(allocated(grad)) deallocate(grad)
          allocate(grad(this%nparams,this%noutput))
          grad=0.0d0

          if(allocated(this%gradinp)) deallocate(this%gradinp)
          allocate(this%gradinp(this%ninput,this%noutput))
          this%gradinp=0.0d0

          do l=this%noutput,1,-1

           !!! Set the index v to the l-output node's params
           if(this%nlayers.eq.1)then
            v=this%nparams+(l-this%noutput)*(this%ninput+1)
           else
            v=this%nparams+(l-this%noutput)*(this%layers(this%nlayers-1)%nneu+1)
           endif
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !!! Compute the grad of output l wrt the l-bias
           grad(v,l)=this%layers(this%nlayers)%neu(l)%grad_act
           v=v-1
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !!! Compute the grad of output l wrt the l-weights
           if(this%nlayers.eq.1)then

            do j=this%ninput,1,-1
             grad(v,l)=this%layers(this%nlayers)%neu(l)%grad_act*this%inp(j)
             v=v-1   
            enddo

           else

            if(allocated(Wmat)) deallocate(Wmat)
            allocate(Wmat(this%layers(this%nlayers-1)%nneu))
            Wmat=0.0d0

            do j=this%layers(this%nlayers-1)%nneu,1,-1

             grad(v,l)=this%layers(this%nlayers-1)%neu(j)%output*this%layers(this%nlayers)%neu(l)%grad_act
             Wmat(j)=this%layers(this%nlayers)%neu(l)%weights(j)*this%layers(this%nlayers)%neu(l)%grad_act  
             v=v-1   

            enddo

           endif
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(this%nlayers.gt.1) v=this%nparams-this%noutput*(this%layers(this%nlayers-1)%nneu+1)  

           do i=this%nlayers-1,1,-1

            if(i.gt.1)then
             if(allocated(Wmat2)) deallocate(Wmat2)
             allocate(Wmat2(this%layers(i)%nneu,this%layers(i-1)%nneu))
             Wmat2=0.0d0
            endif

            do j=this%layers(i)%nneu,1,-1

             grad(v,l)=Wmat(j)*this%layers(i)%neu(j)%grad_act
             v=v-1

             if(i.eq.1)then

              do k=this%ninput,1,-1
               grad(v,l)=Wmat(j)*this%inp(k)*this%layers(i)%neu(j)%grad_act
               this%gradinp(k,l)=Wmat(j)*this%layers(i)%neu(j)%weights(k)*this%layers(i)%neu(j)%grad_act
               v=v-1
              enddo

             else

              do t=this%layers(i-1)%nneu,1,-1
               grad(v,l)=Wmat(j)*this%layers(i-1)%neu(t)%output*this%layers(i)%neu(j)%grad_act
               Wmat2(j,t)=this%layers(i)%neu(j)%weights(t)*this%layers(i)%neu(j)%grad_act
               v=v-1
              enddo

             endif

            enddo ! j nodes

            if(i.gt.1)then
             if(allocated(Wmat3)) deallocate(Wmat3)
             allocate(Wmat3(this%layers(i-1)%nneu))
             Wmat3=0.0d0
             Wmat3=matmul(Wmat,Wmat2)
             deallocate(Wmat) 
             allocate(Wmat(this%layers(i-1)%nneu))
             Wmat=Wmat3            
            endif

           enddo ! i innerlayers                    

          enddo ! l outputs

          if(allocated(Wmat)) deallocate(Wmat)
          if(allocated(Wmat2)) deallocate(Wmat2)
          if(allocated(Wmat3)) deallocate(Wmat3)

         return
         end subroutine backprop


         subroutine set_parameters(this,vec)
         implicit none
         class(net)                    :: this
         integer                       :: i,j,k,l
         double precision, allocatable :: vec(:)

          j=1
          do i=1,this%nlayers
           do k=1,this%layers(i)%nneu
            do l=1,this%layers(i)%neu(k)%nweights
             this%layers(i)%neu(k)%weights(l)=vec(j)
             j=j+1
            enddo
            this%layers(i)%neu(k)%bias=vec(j)
            j=j+1
           enddo
          enddo

         return
         end  subroutine set_parameters


         subroutine get_parameters(this,vec)
         implicit none
         class(net)                    :: this
         integer                       :: i,j,k,l
         double precision, allocatable :: vec(:)

          if (allocated(vec)) deallocate(vec)
          allocate(vec(this%nparams))

          j=1
          do i=1,this%nlayers
           do k=1,this%layers(i)%nneu
            do l=1,this%layers(i)%neu(k)%nweights
             vec(j)=this%layers(i)%neu(k)%weights(l)
             j=j+1
            enddo
            vec(j)=this%layers(i)%neu(k)%bias
            j=j+1
           enddo
          enddo

         return
         end  subroutine get_parameters
         
         subroutine set_topology(this,ninp,layers,nodes)
         implicit none
         class(net)                     :: this
         integer                        :: layers,i,j,ninp
         integer, allocatable           :: nodes(:)
         character(len=10), allocatable :: nodes_act(:)
         
          this%ninput=ninp
          this%noutput=nodes(layers)
          this%nlayers=layers

          allocate(this%layers(layers))

          this%nneu=sum(nodes)
          this%nparams=(this%ninput+1)*nodes(1)
          do i=2,layers
           this%nparams=this%nparams+(nodes(i-1)+1)*nodes(i)
          enddo

          call this%layers(1)%set_nodes(nodes(1),this%ninput) 

          do i=2,layers
           call this%layers(i)%set_nodes(nodes(i),nodes(i-1)) 
           ! add keyword in arguments to set activ type
          enddo

          do i=1,this%noutput
           this%layers(this%nlayers)%neu(i)%activ='none'
          enddo

         return
         end  subroutine set_topology

         subroutine set_nodes(this,nodes,ninp)
         implicit none
         class(layer)                    :: this
         integer                         :: nodes,ninp,i
         
          this%nneu=nodes
          allocate(this%neu(nodes))
          do i=1,nodes
           this%neu(i)%nweights=ninp
           allocate(this%neu(i)%weights(ninp))
           ! set activation funcion
          enddo

         return
         end subroutine set_nodes

         subroutine get_activity(this,v,val)
         implicit none
         class(neuron) :: this
         double precision, allocatable :: v(:)
         double precision              :: val,grad

          val=dot_product(this%weights,v)+this%bias

          if(trim(this%activ).eq.'none')then
           grad=1.0d0
          endif

          if(trim(this%activ).eq.'relu')then
           if(val.gt.0.0d0)then
            grad=1.0d0
           else
            grad=0.0d0
           endif
           val=max(val,0.0d0)
          endif

          if(trim(this%activ).eq.'sigma')then
           val=1/(1+exp(-1.0d0*val))
           grad=val*(1-val)
          endif

          this%output=val
          this%grad_act=grad

         return
         end subroutine get_activity

         subroutine get_grad_nn(this,grad)
         implicit none
         class(net)                    :: this
         double precision, allocatable :: grad(:,:)
         double precision              :: gradp
         integer                       :: i,j,k,v,l
                 
          if(allocated(grad)) deallocate(grad)
          allocate(grad(this%nparams,this%noutput))
          grad=0.0d0

          do l=1,this%noutput
           v=1
           do i=1,this%nlayers
            do j=1,this%layers(i)%nneu
             do k=1,this%layers(i)%neu(j)%nweights
              if(i.eq.this%nlayers .and. j.ne.l)then
               grad(v,l)=0.0d0
               v=v+1
              else
               call this%get_grad_weights(gradp,i,j,k,this%nlayers,l)
               grad(v,l)=gradp
               v=v+1
              endif
             enddo

             if(i.eq.this%nlayers .and. j.ne.l)then
              grad(v,l)=0.0d0
              v=v+1
             else
              call this%get_grad_bias(gradp,i,j,this%nlayers,l)
              grad(v,l)=gradp
              v=v+1
             endif

            enddo 
           enddo
          enddo

         return
         end subroutine get_grad_nn

         recursive subroutine get_grad_weights(this,gradp,layer,node,par,i,j)
         implicit none
         class(net)                    :: this
         double precision              :: val,gradp,new_gradp
         integer                       :: layer,node,i,j,par
         integer                       :: k,v

          gradp=0.0d0

          if(i.eq.layer)then
           val=this%layers(i)%neu(node)%grad_act 
           if(i.ne.1)then
            gradp=val*this%layers(i-1)%neu(par)%output
           else
            gradp=val*this%inp(par)
           endif
           return
          endif

          if(layer.eq.i-1)then
           val=this%layers(i)%neu(j)%grad_act
           call this%get_grad_weights(new_gradp,layer,node,par,i-1,node)
           gradp=gradp+new_gradp*val*this%layers(i)%neu(j)%weights(node)
          else
           val=this%layers(i)%neu(j)%grad_act
           do v=1,this%layers(i)%neu(j)%nweights
            call this%get_grad_weights(new_gradp,layer,node,par,i-1,v)
            gradp=gradp+new_gradp*val*this%layers(i)%neu(j)%weights(v)
           enddo
          endif

         return
         end subroutine get_grad_weights

         recursive subroutine get_grad_bias(this,gradp,layer,node,i,j)
         implicit none
         class(net)                    :: this
         double precision              :: val,gradp,new_gradp
         integer                       :: layer,node,i,j,par
         integer                       :: k,v

          gradp=0.0d0

          if(i.eq.layer)then
           gradp=this%layers(i)%neu(node)%grad_act
           return
          endif

          if(layer.eq.i-1)then
           val=this%layers(i)%neu(j)%grad_act
           call this%get_grad_bias(new_gradp,layer,node,i-1,node)
           gradp=gradp+new_gradp*val*this%layers(i)%neu(j)%weights(node)
          else
           val=this%layers(i)%neu(j)%grad_act
           do v=1,this%layers(i)%neu(j)%nweights
            call this%get_grad_bias(new_gradp,layer,node,i-1,v)
            gradp=gradp+new_gradp*val*this%layers(i)%neu(j)%weights(v)
           enddo
          endif

         return
         end subroutine get_grad_bias
        
         end module nets_class

