        program MolForge
        use atoms_class
        use particles_swarm_class
        use gradmin_class
        use ff_trainer_class
        use general_types_class 
        implicit none
        integer                       :: i,j,offset
        integer                       :: t1,t2,rate
        integer                       :: t1_tot,t2_tot,rate_tot

        type(particles_swarm)         :: swarm
!        type(gradient_descent)        :: grad
        type(adam)                    :: grad
        type(chi2), pointer           :: ff

        integer                       :: nnets
        integer, allocatable          :: types(:),ninp(:),nlayers(:)
        type(vector_int), allocatable :: topo(:) 

        character(len=100)            :: input_file
        double precision, allocatable :: vec(:),vecmax(:),vecmin(:),loc_lr(:) 
        logical, allocatable          :: fixval(:)
        double precision              :: val

        double precision, allocatable :: L2val(:)
        integer, allocatable          :: L2id(:)

          call system_clock(t1_tot,rate_tot)

         ! Setup Training/Test Sets

          call system_clock(t1,rate)

          allocate(ff)
          allocate(ff%FF)          
          ff%FF%do_local_ener=.true.
          ff%FF%do_coul_ener=.true.
          ff%FF%do_disp_ener=.false.

          input_file='input_datas'
          call ff%read_sets(input_file)

          call system_clock(t2)
          write(*,*) 'MolForge set up the systems. Time:',real(t2-t1)/real(rate),'s'

         ! Setup Force Field
          
          call system_clock(t1,rate)

          open(13,file='input_nets')
          read(13,*) nnets

          allocate(types(nnets))
          allocate(topo(nnets))
          allocate(ninp(nnets))
          allocate(nlayers(nnets))

          read(13,*) types
          read(13,*) nlayers

          do i=1,nnets 
           allocate(topo(i)%v(nlayers(i)-1))
           read(13,*) ninp(i),topo(i)%v
          enddo

          close(13)

          call ff%set_nets(nnets,types,ninp,topo)

          deallocate(ninp)
          deallocate(types)
          do i=1,size(topo)
           deallocate(topo(i)%v)
          enddo
          deallocate(topo)        

          call system_clock(t2)
          write(*,*) 'MolForge set up the force field. Time:',real(t2-t1)/real(rate),'s'

!          call system_clock(t1,rate)

!          allocate(L2id(ff%FF%nparams-1))
!          allocate(L2val(ff%FF%nparams-1))   
!          L2val=0.1d0
!          do i=1,size(L2id)
!           L2id(i)=i
!          enddo
!          call ff%set_L2(L2id=L2id,L2val=L2val)
!          call ff%ridge()
!          call ff%out_results()
!          stop

!          call system_clock(t2)
!          write(*,*) 'MolForge permormed Ridge regression. Time:',real(t2-t1)/real(rate),'s'

         ! Optimize Force Field

          call system_clock(t1,rate)

          allocate(L2id(ff%FF%nparams-1))
          allocate(L2val(ff%FF%nparams-1))   
          L2val=0.01d0
          j=1
          do i=1,ff%FF%local_nparams-1
           L2id(j)=i
           j=j+1
          enddo
          offset=ff%FF%local_nparams
          do i=1,ff%FF%coul_nparams
           L2id(j)=i+offset
           j=j+1
          enddo

          call ff%set_L2(L2id=L2id,L2val=L2val)

!          call ff%set_L2(L2=0.001d0)
!          call ff%std_sets()

          call system_clock(t2)
          write(*,*) 'MolForge set up L2 and std sets. Time:',real(t2-t1)/real(rate),'s'
         
          allocate(vecmin(ff%FF%nparams))
          allocate(vecmax(ff%FF%nparams))
          vecmin=-0.5d0
          vecmax=0.5d0
          vecmin(112)=-23915.318176339031
          vecmax(112)=-23915.318176339031
          vecmin(56)=0.0d0
          vecmax(56)=0.0d0
          vecmin(168)=0.0d0
          vecmax(168)=0.0d0
          vecmin(224)=0.0d0
          vecmax(224)=0.0d0
          allocate(fixval(ff%FF%nparams))
          fixval=.false.
          fixval(56)=.true.
          fixval(112)=.true.          
          fixval(168)=.true.          
          fixval(224)=.true.          

          swarm%target_f => ff
          call swarm%init_swarm(nval=ff%FF%nparams,npar=25,min_val=vecmin,max_val=vecmax,fixval=fixval)
!          call swarm%init_swarm(nval=ff%FF%nparams,npar=25) 
!          swarm%print_val=.true.
          call swarm%minimize(max_iter=1000)        
          call swarm%get_best_val(vec)
          call swarm%release_target_f()

          call ff%get_fval(vec,val)

          allocate(loc_lr(ff%FF%nparams))
          loc_lr=1.0d0
          loc_lr(56)=0.0d0
          loc_lr(112)=0.0d0
          loc_lr(168)=0.0d0
          loc_lr(224)=0.0d0
          grad%target_f => ff
          call grad%init(nval=ff%FF%nparams,vec=vec,loc_lr=loc_lr)
!!          call grad%init(nval=ff%FF%nparams)
!!          grad%print_val=.true.
!!          grad%print_grad=.true.
          grad%lr=1.0e-5
          call grad%minimize(max_iter=10000,start_iter=25000)          
          call grad%release_target_f()

         ! Print Results

          call ff%out_results()

          call system_clock(t2_tot)
          write(*,*) 'MolForge Simulation Concluded.'
          write(*,*) 'Total Simulation Time:',real(t2_tot-t1_tot)/real(rate_tot),'s'

        return
        end program MolForge

