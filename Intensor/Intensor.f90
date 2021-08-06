        program Intensor
        use atoms_class
        use particles_swarm_class
        use gradmin_class
        use mlmodel_trainer_class
        use sh_mlmodel_trainer_class
        use general_types_class 
        implicit none
        integer                         :: i,j
        integer                         :: t1,t2,rate
        integer                         :: t1_tot,t2_tot,rate_tot

        type sh_interpolator
         type(sh_mlmodel_trainer), pointer  :: model
         type(particles_swarm)           :: swarm
         type(adam)                      :: grad
        end type sh_interpolator

        type(particles_swarm)           :: swarm
        type(adam)                      :: grad
        type(mlmodel_trainer), pointer  :: model

        integer                         :: ninp,nlayers
        type(vector_int)                :: topo

        character(len=100)              :: input_file
        double precision, allocatable   :: vec(:),vecmax(:),vecmin(:),loc_lr(:) 
        logical, allocatable            :: fixval(:)
        double precision                :: val

        double precision, allocatable   :: L2val(:)
        integer, allocatable            :: L2id(:)

          call system_clock(t1_tot,rate_tot)

         ! Setup Training/Test Sets

          call system_clock(t1,rate)

          allocate(model)
          allocate(model%ML)                    

          input_file='input_datas'
          call model%read_sets(input_file)

          call system_clock(t2)
          write(*,*) 'Intensor set up the systems. Time:',real(t2-t1)/real(rate),'s'

         ! Setup ML model
          
          call system_clock(t1,rate)

          open(13,file='input_net')

          read(13,*) nlayers ! input layer + inner layers + output layer
          allocate(topo%v(nlayers-1))
          read(13,*) ninp,topo%v

          close(13)

          call model%set_nets(ninp,topo)
          deallocate(topo%v)

          call system_clock(t2)
          write(*,*) 'Intensor set up the model. Time:',real(t2-t1)/real(rate),'s'

          call system_clock(t1,rate)

!          allocate(L2id(model%ML%nparams-1))
!          allocate(L2val(model%ML%nparams-1))   
!          L2val=0.01d0
!          j=1
!          do i=1,model%ML%nparams-1
!           L2id(j)=i
!           j=j+1
!          enddo

!          call model%set_L2(L2id=L2id,L2val=L2val)

!          call model%set_L2(L2=0.00001d0)
          call model%std_sets()

          call system_clock(t2)
          write(*,*) 'Intensor set up L2 regularization and std sets. Time:',real(t2-t1)/real(rate),'s'
        
          allocate(vecmin(model%ML%nparams))
          allocate(vecmax(model%ML%nparams))
          vecmin=-0.5d0
          vecmax=0.5d0
          allocate(fixval(model%ML%nparams))
          fixval=.false.

          swarm%target_f => model

 !         call swarm%init_swarm(nval=model%ML%nparams,npar=25,min_val=vecmin,max_val=vecmax,fixval=fixval)
          call swarm%init_swarm(nval=model%ML%nparams,npar=25) 
!          swarm%print_val=.true.
          call swarm%minimize(max_iter=100)        
          call swarm%get_best_val(vec)
          call swarm%release_target_f()
          call model%get_fval(vec,val)

!          allocate(loc_lr(model%ML%nparams))
!          loc_lr=1.0d0
          grad%target_f => model
!         call grad%init(nval=model%ML%nparams,vec=vec,loc_lr=loc_lr)
          call grad%init(nval=model%ML%nparams,vec=vec)
!          grad%print_val=.true.
!          grad%print_grad=.true.
          grad%lr=1.0e-3        
          call grad%minimize(max_iter=5000,start_iter=2500)          
          call grad%release_target_f()

         ! Print Results

          call model%out_results()

          call system_clock(t2_tot)
          write(*,*) 'Intensor Simulation Concluded'
          write(*,*) 'Total Simulation Time:',real(t2_tot-t1_tot)/real(rate_tot),'s'

        return
        end program Intensor

