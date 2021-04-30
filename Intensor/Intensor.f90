        program NeuralNets
        use nets_class
        use random_numbers_class        
        use particles_swarm_class
        use external_functions_class
        implicit none
        type(particles_swarm)         :: swarm
        integer, allocatable          :: nodes(:)
        integer                       :: i,j,max_iter,npar,nnodes=1,ninp,iter,ii
        double precision              :: rand_num,rmse_tr,rmse_te,gradnorm,val,disp
        double precision, allocatable :: vec(:),step(:),min_val(:),max_val(:)
        double precision, allocatable :: inps(:),outs(:)
        double precision, allocatable :: grad(:),gradres(:),gradres2(:)
        double precision              :: beta1,beta2,alpha,eps,lr
        logical                       :: do_grad=.false.,do_swarm=.false.
        logical                       :: random_init=.true.,restart=.false.
        character(len=100)            :: word,restart_file
       
         do i=1,iargc()
          call getarg(i,word)

          select case (trim(word))

             case ('-adam')
                 do_grad=.true.

             case ('-swarm')
                 do_swarm=.true.

             case ('-deep')
                 call getarg(i+1,word)
                 read(word,*) nnodes
                 nnodes=nnodes+1
                 allocate(nodes(nnodes))
                 do j=1,nnodes-1
                  call getarg(i+1+j,word)
                  read(word,*) nodes(j)
                 enddo

             case ('-restart')
                 call getarg(i+1,word)
                 read(word,*) restart_file
                 random_init=.false.
                 restart=.true.

          end select

         enddo

!        Gen Train/Set

         call f%read_sets()
         f%batch_size=2000

!        Set Net Topology

         ninp=f%tr(1)%size_desc
         if(.not. allocated(nodes)) allocate(nodes(nnodes))
         nodes(nnodes)=f%ndim
         call f%mymodel%set_topology(ninp,nnodes,nodes)

         write(*,*) "N Parameters: ",f%mymodel%nparams,nodes

!        Set Initial Net Parameters
        
         allocate(vec(f%mymodel%nparams))
         if(restart)then
          open(13,file=trim(restart_file))
          read(13,*) ii
          if(f%mymodel%nparams.eq.ii)then
           write(*,*) 'Reading Parameters from: ',trim(restart_file)
           read(13,*) vec
           random_init=.false.
          else
           Write(*,*) 'Restart parameters conflict with topology'
           Write(*,*) 'Setting initial parameters to random'
           random_init=.true.
          endif
          close(13)
         endif
         
         if(random_init)then
          call init_random_seed()
          do i=1,f%mymodel%nparams
           call random_number(rand_num)
           vec(i)=rand_num
          enddo
         endif

         call f%mymodel%set_parameters(vec)

!        Set Swarm Parameters

         if(do_swarm)then

          npar=50

          allocate(step(f%mymodel%nparams))
          allocate(min_val(f%mymodel%nparams))
          allocate(max_val(f%mymodel%nparams))

          step=0.6d0
          min_val=-5d-1
          max_val=5d-1
          max_iter=1000

!         swarm%do_meta=.true.
          swarm%meta%ngauss=0
          swarm%meta%sigma=0.1
          swarm%meta%height=0.1
          allocate(swarm%meta%gauss(npar*max_iter))         

!        Net Training

          call swarm%init_swarm(npar,f%mymodel%nparams,step,min_val,max_val)
          call swarm%minimize(max_iter)
          vec=swarm%par(swarm%best_par)%val

         endif

!        Adam Minimization
         
         if(do_grad)then

          iter=1
          max_iter=3000
          allocate(gradres(size(vec)))
          allocate(gradres2(size(vec)))
          eps=1.0d-7
          beta1=0.9d0
          beta2=0.999d0
          lr=0.0005d0
          gradres=0.0d0
          gradres2=0.0d0

          do while (iter.le.max_iter)
           call f%get_fgrad(vec,grad,val)
           call f%get_fval_test(vec,rmse_te)

           gradres=beta1*gradres+(1-beta1)*grad
           gradres2=beta2*gradres2+(1-beta2)*grad**2
           
           gradnorm=0.0d0
           do i=1,size(grad)
            gradnorm=gradnorm+grad(i)**2
           enddo

           write(*,*) iter,val,sqrt(gradnorm/size(grad)),rmse_te

           vec=vec-lr*(gradres/(1-beta1**(iter+1)))/(sqrt(gradres2/(1-beta2**(iter+1)))+eps)

           iter=iter+1
          enddo

         endif

!        Output Parameters

         open(13,file='model_param.dat')
         write(13,*) size(vec)
         write(13,*) vec
         close(13)

!        Output Errors

         call f%get_fval(vec,rmse_tr)         
         call f%get_fval_test(vec,rmse_te)  
         call f%mymodel%set_parameters(vec)
         call f%out_results()

         close(12)
         write(*,*) 'RMSE Training and Test Sets: ',rmse_tr,rmse_te
 
        return
        end program NeuralNets
