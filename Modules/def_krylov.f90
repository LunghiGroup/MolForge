        module krylov_class
        use hilbert_dist_class
        use tensor_network_class
        use units_parms
        implicit none

        type                           :: krylov_parameters
         !---------------------------------------------------------------
         ! Krylov/Lanczos parameters.
         !---------------------------------------------------------------
         integer                       :: dimension=-1
         integer                       :: restart=-1
         integer                       :: max_restart=1000
         integer                       :: reorthogonalization=0
         logical                       :: full_basis=.true.
         logical                       :: adaptive_restart=.false.
         logical                       :: rebuild=.false.
         double precision              :: tol=1.0d-10
         double precision              :: breakdown_tol=1.0d-14
         double precision              :: orth_tol=1.0d-12
         !---------------------------------------------------------------
         ! Time-evolution parameters.
         !---------------------------------------------------------------
         double precision              :: tstart=0.0d0
         double precision              :: tfinal=0.0d0
         double precision              :: dt=0.01d0
         integer                       :: nsteps=0
         !---------------------------------------------------------------
         ! Tensor-network parameters.
         !---------------------------------------------------------------
         integer                       :: mps_bond_dimension_max=100
         integer                       :: mpo_bond_dimension_max=100
         double precision              :: mps_trunc_tol=1.0d-8
         double precision              :: mpo_trunc_tol=1.0d-8
         !---------------------------------------------------------------
         ! Output and observable parameters.
         !---------------------------------------------------------------
         integer                       :: output_stride=1
         logical                       :: store_states=.false.
         logical                       :: calc_energy=.false.
         logical                       :: calc_entropy=.false.
         logical                       :: report_truncation=.true.
         logical                       :: verbose=.false.
         logical                       :: save_checkpoint=.false.
         integer                       :: checkpoint_stride=1000
         contains
         procedure                     :: bcast_params
        end type krylov_parameters

        type                           :: krylov_space
         !---------------------------------------------------------------
         ! User-defined Krylov and tensor-network parameters.
         !---------------------------------------------------------------
         type(krylov_parameters)       :: parameters
         !---------------------------------------------------------------
         ! Hamiltonian represented as a Matrix Product Operator.
         !---------------------------------------------------------------
         type(mpo_operator)            :: H
         !---------------------------------------------------------------
         ! Spin operators represented as Matrix Product Operators.
         !---------------------------------------------------------------
         type(mpo_operator)            :: Sx
         type(mpo_operator)            :: Sy
         type(mpo_operator)            :: Sz
         !---------------------------------------------------------------
         ! Current propagated state.
         !---------------------------------------------------------------
         type(mps_state)               :: psi
         !---------------------------------------------------------------
         ! Krylov basis vectors
         !
         !     V(1), V(2), ..., V(m)
         !
         ! generated during the Lanczos iteration.
         !---------------------------------------------------------------
         type(mps_state), allocatable  :: V(:)
         !---------------------------------------------------------------
         ! Lanczos tridiagonal matrix coefficients.
         !---------------------------------------------------------------
         complex(8), allocatable       :: alpha(:)
         complex(8), allocatable       :: beta(:)
         !---------------------------------------------------------------
         ! Eigenvalues and eigenvectors of the Lanczos tridiagonal
         ! matrix.
         !
         ! These are retained during fixed-basis propagation so that
         ! the tridiagonal matrix does not need to be diagonalized
         ! at every timestep.
         !---------------------------------------------------------------
         double precision, allocatable :: eigenvalue(:)
         complex(8), allocatable       :: eigenvector(:,:)
         !---------------------------------------------------------------
         ! Number of Krylov vectors currently generated.
         !---------------------------------------------------------------
         integer                       :: current_dimension=0
         contains
         procedure                     :: initialize
         procedure                     :: initialize_state
         procedure                     :: initialize_H
         procedure                     :: initialize_basis
         procedure                     :: initialize_workspace
         procedure                     :: propagate
         procedure                     :: build_initial_vector
         procedure                     :: build_krylov_space
         procedure                     :: lanczos_step
         procedure                     :: orthogonalize
         procedure                     :: apply_H
         procedure                     :: diagonalize_tridiagonal
         procedure                     :: exponentiate_tridiagonal
         procedure                     :: reconstruct_state
         procedure                     :: construct_observables
         procedure                     :: compute_observables
         procedure                     :: finalize
        end type krylov_space

        contains

        subroutine bcast_params(this)
        use mpi
        use mpi_utils
        implicit none
        class(krylov_parameters)       :: this

        !-------------------------------------------------------------
        ! Broadcast all Krylov parameters from the master MPI process
        ! to all other processes.
        !-------------------------------------------------------------

        !-----------------------------------------------------
        ! Integer parameters
        !-----------------------------------------------------
        call MPI_BCAST(this%dimension,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%max_restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%nsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%mps_bond_dimension_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%mpo_bond_dimension_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%reorthogonalization,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%output_stride,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%checkpoint_stride,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)

        !-----------------------------------------------------
        ! Logical parameters
        !-----------------------------------------------------
        call MPI_BCAST(this%full_basis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%adaptive_restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%rebuild,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%store_states,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%calc_energy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%calc_entropy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%report_truncation,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%verbose,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%save_checkpoint,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err)

        !-----------------------------------------------------
        ! Double precision parameters
        !-----------------------------------------------------
        call MPI_BCAST(this%tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%breakdown_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%orth_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%tstart,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%tfinal,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)

        call MPI_BCAST(this%mps_trunc_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        call MPI_BCAST(this%mpo_trunc_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)

        end subroutine bcast_params

        subroutine initialize(this,spindy,initial_state_type,basis_state,product_state,sup_basis,sup_coeff)
        implicit none
        class(krylov_space)            :: this
        class(spins_hilbert)           :: spindy
        integer                        :: initial_state_type,basis_state(:),sup_basis(:,:)
        complex(8)                     :: product_state(:,:),sup_coeff(:)

         !---------------------------------------------------------------
         ! Initialise the complete Krylov propagator.
         !
         ! The Krylov parameters have already been transferred from
         ! control_variables into this%parameters by do_action.
         !
         ! The initialisation therefore consists only of constructing
         ! the physical and tensor-network objects required for the
         ! Lanczos propagation.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Validate Krylov and tensor-network parameters.
         !
         ! Fatal errors correspond to parameter combinations for which
         ! the calculation cannot be meaningfully performed.
         !
         ! Warnings identify unusual or potentially inefficient choices
         ! but do not prevent the calculation from proceeding.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Krylov dimension.
         !---------------------------------------------------------------
         if(this%parameters%dimension<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov dimension must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Krylov convergence tolerance.
         !---------------------------------------------------------------
         if(this%parameters%tol<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov tolerance must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Lanczos breakdown tolerance.
         !---------------------------------------------------------------
         if(this%parameters%breakdown_tol<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov breakdown tolerance must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Orthogonalisation tolerance.
         !---------------------------------------------------------------
         if(this%parameters%orth_tol<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov orthogonalisation tolerance must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Reorthogonalisation method.
         !
         !   0 = none
         !   1 = full
         !   2 = selective
         !---------------------------------------------------------------
         if(this%parameters%reorthogonalization<0.or.this%parameters%reorthogonalization>2)then
          if(mpi_id==0)then
           write(6,*) 'Err: Invalid Krylov reorthogonalisation option'
           write(6,*) '     Allowed values are 0 (NONE), 1 (FULL), 2 (SELECTIVE)'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Restart parameters.
         !---------------------------------------------------------------
         if(this%parameters%max_restart<0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Maximum number of Krylov restarts cannot be negative'
          endif
          stop
         endif

         if(this%parameters%restart>this%parameters%dimension)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov restart dimension cannot exceed Krylov dimension'
           write(6,*) '     Restart = ',this%parameters%restart
           write(6,*) '     KDIM    = ',this%parameters%dimension
          endif
          stop
         endif

         if(this%parameters%restart==this%parameters%dimension.and.this%parameters%restart>0)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Krylov restart dimension equals the full Krylov dimension'
           write(6,*) '         Restarting will therefore provide little reduction'
           write(6,*) '         in Krylov-space size'
          endif
         endif

         if(this%parameters%adaptive_restart.and.this%parameters%max_restart==0)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Adaptive Krylov restart is enabled but'
           write(6,*) '         maximum number of restarts is zero'
          endif
         endif

         !---------------------------------------------------------------
         ! Full-basis storage.
         !
         ! A full Krylov basis requires storing every MPS generated
         ! during the iteration and can therefore become expensive.
         !---------------------------------------------------------------
         if(this%parameters%full_basis.and.this%parameters%dimension>100)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Full Krylov basis requested with KDIM > 100'
           write(6,*) '         Memory requirements may become substantial'
          endif
         endif

         !---------------------------------------------------------------
         ! Tensor-network bond dimensions.
         !---------------------------------------------------------------
         if(this%parameters%mps_bond_dimension_max<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPS maximum bond dimension must be greater than zero'
          endif
          stop
         endif

         if(this%parameters%mpo_bond_dimension_max<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPO maximum bond dimension must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Tensor-network truncation tolerances.
         !---------------------------------------------------------------
         if(this%parameters%mps_trunc_tol<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPS truncation tolerance must be greater than zero'
          endif
          stop
         endif

         if(this%parameters%mpo_trunc_tol<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPO truncation tolerance must be greater than zero'
          endif
          stop
         endif

         if(this%parameters%mps_trunc_tol>1.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Warning: MPS truncation tolerance is greater than one'
           write(6,*) '         This may result in severe truncation'
          endif
         endif

         if(this%parameters%mpo_trunc_tol>1.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Warning: MPO truncation tolerance is greater than one'
           write(6,*) '         This may result in severe truncation'
          endif
         endif

         !---------------------------------------------------------------
         ! Time-evolution parameters.
         !---------------------------------------------------------------
         if(this%parameters%tfinal<this%parameters%tstart)then
          if(mpi_id==0)then
           write(6,*) 'Err: Final time is smaller than initial time'
           write(6,*) '     TSTART = ',this%parameters%tstart
           write(6,*) '     TFINAL = ',this%parameters%tfinal
          endif
          stop
         endif

         if(this%parameters%dt<=0.0d0.and.this%parameters%nsteps<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Time step must be greater than zero when NSTEP is not specified'
          endif
          stop
         endif

         if(this%parameters%nsteps<0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Number of time steps cannot be negative'
          endif
          stop
         endif

         if(this%parameters%nsteps>0.and.this%parameters%dt<=0.0d0)then
          if(mpi_id==0)then
           write(6,*) 'Err: DT must be greater than zero when NSTEP is specified'
          endif
          stop
         endif

         if(this%parameters%nsteps>0.and.this%parameters%tfinal>this%parameters%tstart)then
          if(mpi_id==0)then
           if(abs(this%parameters%tstart + &
                  this%parameters%nsteps*this%parameters%dt - &
                  this%parameters%tfinal) > &
              1.0d-10*max(1.0d0,abs(this%parameters%tfinal)))then

            write(6,*) 'Warning: NSTEP*DT does not exactly span the requested'
            write(6,*) '         time interval TSTART to TFINAL'
           endif
          endif
         endif

         !---------------------------------------------------------------
         ! Output parameters.
         !---------------------------------------------------------------
         if(this%parameters%output_stride<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Output stride must be greater than zero'
          endif
          stop
         endif

         if(this%parameters%checkpoint_stride<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Checkpoint stride must be greater than zero'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Informative warnings for potentially expensive calculations.
         !---------------------------------------------------------------
         if(this%parameters%calc_entropy.and.this%parameters%output_stride==1)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Entanglement entropy will be calculated'
           write(6,*) '         at every time step'
          endif
         endif

         if(this%parameters%store_states)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Storing all propagated MPS states may'
           write(6,*) '         require substantial memory'
          endif
         endif

         if(this%parameters%save_checkpoint.and.this%parameters%checkpoint_stride<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Checkpointing enabled with invalid checkpoint stride'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Numerical recommendation for Lanczos orthogonalisation.
         !---------------------------------------------------------------
         if(this%parameters%reorthogonalization==0.and.this%parameters%dimension>50)then
          if(mpi_id==0)then
           write(6,*) 'Warning: No Lanczos reorthogonalisation requested'
           write(6,*) '         with a relatively large Krylov dimension'
           write(6,*) '         Loss of orthogonality may become significant'
          endif
         endif

         if(this%parameters%reorthogonalization==2.and.this%parameters%dimension<10)then
          if(mpi_id==0)then
           write(6,*) 'Warning: Selective reorthogonalisation requested'
           write(6,*) '         for a very small Krylov dimension'
          endif
         endif

         !---------------------------------------------------------------
         ! Report the selected configuration when verbose output is
         ! requested.
         !---------------------------------------------------------------
         if(this%parameters%verbose.and.mpi_id==0)then
          write(6,*) ''
          write(6,*) '     Krylov configuration:'
          write(6,*) '       Krylov dimension       = ',this%parameters%dimension
          write(6,*) '       Krylov tolerance       = ',this%parameters%tol
          write(6,*) '       Breakdown tolerance    = ',this%parameters%breakdown_tol
          write(6,*) '       Orthogonalisation tol  = ',this%parameters%orth_tol
          write(6,*) '       Restart dimension      = ',this%parameters%restart
          write(6,*) '       Maximum restarts       = ',this%parameters%max_restart
          write(6,*) '       Reorthogonalisation    = ',this%parameters%reorthogonalization
          write(6,*) '       Full basis             = ',this%parameters%full_basis
          write(6,*) '       Adaptive restart       = ',this%parameters%adaptive_restart
          write(6,*) ''
          write(6,*) '     Tensor-network configuration:'
          write(6,*) '       Max MPS bond dimension     = ',this%parameters%mps_bond_dimension_max
          write(6,*) '       Max MPO bond dimension     = ',this%parameters%mpo_bond_dimension_max
          write(6,*) '       MPS truncation tol     = ',this%parameters%mps_trunc_tol
          write(6,*) '       MPO truncation tol     = ',this%parameters%mpo_trunc_tol
          write(6,*) ''
          write(6,*) '     Time evolution:'
          write(6,*) '       TSTART                 = ',this%parameters%tstart,'ps'
          write(6,*) '       TFINAL                 = ',this%parameters%tfinal,'ps'
          write(6,*) '       DT                     = ',this%parameters%dt,'ps'
          write(6,*) '       NSTEPS                 = ',this%parameters%nsteps
          write(6,*) ''
         endif

         !---------------------------------------------------------------
         ! Construct the initial MPS and Hamiltonian MPO on rank 0.
         !---------------------------------------------------------------
         if(mpi_id==0)then
          call this%initialize_state(spindy,initial_state_type,basis_state,product_state,sup_basis,sup_coeff)
          call this%initialize_H(spindy)
         endif

         !---------------------------------------------------------------
         ! Distribute the initial MPS and Hamiltonian MPO to all ranks.
         !---------------------------------------------------------------
         call this%psi%bcast_mps()
         call this%H%bcast_mpo()

         !---------------------------------------------------------------
         ! Initialise the Krylov basis and workspace independently
         ! on every MPI process.
         !---------------------------------------------------------------
         call this%initialize_basis()
         call this%initialize_workspace()

        return
        end subroutine initialize

        subroutine initialize_state(this,spindy,initial_state_type,basis_state,product_state,sup_basis,sup_coeff)
        implicit none
        class(krylov_space)            :: this
        class(spins_hilbert)           :: spindy
        integer                        :: initial_state_type,basis_state(:),sup_basis(:,:)
        complex(8)                     :: product_state(:,:),sup_coeff(:)

         !---------------------------------------------------------------
         ! Construct the initial Matrix Product State.
         !
         ! The tensor-network parameters are first copied into the MPS.
         ! The requested initial state is then constructed using one of
         ! the MPS initialisation routines.
         !
         ! The local spin basis for every spin kind is generated
         ! internally by the MPS initialisation routines.
         !---------------------------------------------------------------

         this%psi%bond_dimension_max = this%parameters%mps_bond_dimension_max
         this%psi%trunc_tol = this%parameters%mps_trunc_tol

         select case(initial_state_type)
          !--------------------------------------------------------------
          ! Computational basis state.
          !--------------------------------------------------------------
          case(0)
           call this%psi%init_basis(spindy,basis_state)

          !--------------------------------------------------------------
          ! Product state.
          !--------------------------------------------------------------
          case(1)
           call this%psi%init_product(spindy,product_state)

          !--------------------------------------------------------------
          ! General product of local superposition states.
          !--------------------------------------------------------------
          case(2)
           call this%psi%init_superposition(spindy,sup_basis,sup_coeff)

          case default
           write(6,*) 'Err: Unknown initial_state_type in initialize_state'
           stop

         end select

         !---------------------------------------------------------------
         ! The initial state should be normalised before entering the
         ! Lanczos iteration.
         !
         ! The MPS construction routines provide the state, while the
         ! Krylov initialisation establishes the normalisation convention.
         !---------------------------------------------------------------
         call this%psi%normalize()

        return
        end subroutine initialize_state

        subroutine initialize_H(this,spindy)
        implicit none
        class(krylov_space)            :: this
        class(spins_hilbert)           :: spindy

         !---------------------------------------------------------------
         ! Construct the Hamiltonian Matrix Product Operator.
         !
         ! The Hamiltonian MPO is assembled exactly from the spin
         ! Hamiltonian specified in the parent spins_hilbert object.
         !
         ! During construction all Hamiltonian terms are added exactly
         ! without truncation. After the complete Hamiltonian has been
         ! assembled, a single compression sweep is performed to reduce
         ! the MPO bond dimensions while maintaining the requested
         ! truncation accuracy.
         !
         ! The resulting compressed MPO is subsequently used throughout
         ! the Krylov propagation.
         !---------------------------------------------------------------

         this%H%bond_dimension_max = this%parameters%mpo_bond_dimension_max
         this%H%trunc_tol = this%parameters%mpo_trunc_tol

         call this%H%init_spin_H(spindy)

         !---------------------------------------------------------------
         ! Report the final Hamiltonian MPO bond dimension.
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          write(*,*)
          write(*,*) '     Hamiltonian MPO successfully constructed'
          write(*,*) '     MPO bond dimension : ',this%H%bond_dimension
          write(*,*) '     MPO truncation error : ',this%H%trunc_error
          write(*,*)
          flush(6)
         endif

        return
        end subroutine initialize_H

        subroutine initialize_basis(this)
        implicit none
        class(krylov_space)            :: this

         !---------------------------------------------------------------
         ! Allocate the Krylov basis and Lanczos coefficients.
         !
         ! Storage is allocated for the requested maximum Krylov
         ! dimension:
         !
         !     V(1), V(2), ..., V(m)
         !
         ! together with the diagonal and off-diagonal coefficients of
         ! the Lanczos tridiagonal matrix:
         !
         !     alpha(i)
         !     beta(i)
         !
         ! The MPS tensors contained in V are not constructed here.
         ! They are generated during the Lanczos iteration.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Remove any existing Krylov workspace.
         ! This allows initialize_basis to be called safely more than
         ! once if the Krylov object is reused.
         !---------------------------------------------------------------
         if(allocated(this%V)) then
          deallocate(this%V)
         endif

         if(allocated(this%alpha)) then
          deallocate(this%alpha)
         endif

         if(allocated(this%beta)) then
          deallocate(this%beta)
         endif

         !---------------------------------------------------------------
         ! Allocate the maximum Krylov basis.
         !
         ! Each element V(i) is an mps_state. The individual MPS tensors
         ! are allocated later when the corresponding Krylov vector is
         ! constructed.
         !---------------------------------------------------------------
         allocate(this%V(this%parameters%dimension))

         !---------------------------------------------------------------
         ! Allocate the Lanczos tridiagonal coefficients.
         !
         ! For an m-dimensional Krylov space:
         !
         !              alpha(1) beta(1)
         !                       alpha(2) beta(2)
         !                                ...
         !
         ! The final beta value is not required for an m-dimensional
         ! tridiagonal matrix, but allocating dimension elements keeps
         ! the indexing simple and is useful during the iteration.
         !---------------------------------------------------------------
         allocate(this%alpha(this%parameters%dimension))
         allocate(this%beta(this%parameters%dimension))

         !---------------------------------------------------------------
         ! Initialise all coefficients to zero.
         !---------------------------------------------------------------
         this%alpha = (0.0d0,0.0d0)
         this%beta = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! No Krylov vectors have been generated yet.
         !---------------------------------------------------------------
         this%current_dimension = 0

         !---------------------------------------------------------------
         ! Report the allocated Krylov workspace.
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          write(*,*)
          write(*,*) '     Krylov basis allocated'
          write(*,*) '     Maximum Krylov dimension : ',this%parameters%dimension
          flush(6)
         endif

        return
        end subroutine initialize_basis

        subroutine initialize_workspace(this)
        implicit none
        class(krylov_space)            :: this

         !---------------------------------------------------------------
         ! Initialise the Krylov workspace.
         !
         ! The temporary MPS objects required during the Lanczos
         ! iteration are created locally by the corresponding helper
         ! routines rather than being stored permanently in the Krylov
         ! object.
         !
         ! This routine therefore resets the persistent Krylov
         ! propagation state ready for a new evolution.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! No Krylov vectors have yet been generated.
         !---------------------------------------------------------------
         this%current_dimension = 0

         !---------------------------------------------------------------
         ! Reset the Lanczos coefficients.
         !
         ! This is useful if the Krylov object is reused for another
         ! propagation without reconstructing the entire object.
         !---------------------------------------------------------------
         if(allocated(this%alpha))then
          this%alpha = (0.0d0,0.0d0)
         endif

         if(allocated(this%beta))then
          this%beta = (0.0d0,0.0d0)
         endif

         !---------------------------------------------------------------
         ! Report successful workspace initialisation.
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          write(*,*) '     Krylov workspace initialised'
          flush(6)
         endif

        return
        end subroutine initialize_workspace

        subroutine propagate(this)
        implicit none
        class(krylov_space)            :: this
        integer                        :: istep
        integer                        :: nsteps
        double precision               :: dt,time
        complex(8), allocatable        :: y(:)
        integer                        :: t1,t2,rate

         !---------------------------------------------------------------
         ! Determine the number of propagation steps.
         !
         ! If nsteps has been explicitly specified, use it.
         ! Otherwise determine the number of steps from the requested
         ! time interval and timestep.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Performing Krylov Time Propagation'
          flush(6)
         endif

         nsteps = this%parameters%nsteps
         dt = this%parameters%dt

         if(nsteps<=0) then
          if(this%parameters%tfinal<=this%parameters%tstart) then
           if (mpi_id==0) then
            write(6,*) 'Err: Invalid time interval in Krylov propagation'
           endif
           stop
          endif

          nsteps = nint((this%parameters%tfinal-this%parameters%tstart)/dt)

          if(nsteps<=0) then
           if (mpi_id==0) then
            write(6,*) 'Err: Invalid number of Krylov propagation steps'
           endif
           stop
          endif
         endif

         !---------------------------------------------------------------
         ! Report propagation parameters.
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          write(*,*)
          write(*,*) '     Initial time        : ',this%parameters%tstart,'ps'
          write(*,*) '     Time step           : ',dt,'ps'
          write(*,*) '     Number of steps     : ',nsteps
          write(*,*) '     Maximum Krylov dim. : ',this%parameters%dimension
          write(*,*)
          flush(6)
         endif

         !---------------------------------------------------------------
         ! Calculate observables for the initial state.
         !
         ! This is timestep zero:
         !
         !     t = tstart
         !
         ! No propagation has yet taken place.
         !---------------------------------------------------------------
         call this%compute_observables(0,this%parameters%tstart)

         !---------------------------------------------------------------
         ! Time propagation loop.
         !
         ! At each timestep:
         !
         !     1. Construct the initial Krylov vector from psi if activated.
         !     2. Construct the Lanczos basis and tridiagonal
         !        representation of H.
         !     3. Exponentiate the projected Hamiltonian.
         !     4. Reconstruct the propagated MPS.
         !     5. Compress the propagated MPS according to the
         !        tensor-network truncation parameters.
         !     6. Calculate requested observables.
         !
         ! The physical state stored in this%psi is therefore the state
         ! at the beginning of each timestep and is replaced by the
         ! state at the end of that timestep.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Building the Krylov space before evolution.
         !---------------------------------------------------------------
         if(.not.(this%parameters%rebuild))then
          call this%build_initial_vector()
          call this%build_krylov_space()

          !-------------------------------------------------------------
          ! Diagonalise the projected Hamiltonian once.
          !-------------------------------------------------------------
          call this%diagonalize_tridiagonal()

          !-------------------------------------------------------------
          ! Initialise the projected state coefficients.
          !-------------------------------------------------------------
          if (allocated(y)) deallocate(y)
          allocate(y(this%current_dimension))

          y = (0.0d0,0.0d0)
          y(1) = (1.0d0,0.0d0)
         endif

         do istep=1,nsteps
          !-------------------------------------------------------------
          ! Rebuilding the Krylov subspace at each timestep if set.
          !-------------------------------------------------------------
          if(this%parameters%rebuild)then
           !-------------------------------------------------------------
           ! Construct V(1) from the current state.
           !-------------------------------------------------------------
           call this%build_initial_vector()

           !-------------------------------------------------------------
           ! Construct the Lanczos Krylov space.
           !
           ! The Lanczos routine determines the actual Krylov dimension
           ! reached during this timestep and stores the corresponding
           ! alpha and beta coefficients.
           !-------------------------------------------------------------
           call this%build_krylov_space()

           !-------------------------------------------------------------
           ! Diagonalize the tridiagonal matrix for the new Krylov space.
           !-------------------------------------------------------------
           call this%diagonalize_tridiagonal()
           if(allocated(y)) deallocate(y)
           allocate(y(this%current_dimension))
           y = (0.0d0,0.0d0)
           y(1) = (1.0d0,0.0d0)
          endif

          !-------------------------------------------------------------
          ! Propagate within the Krylov subspace.
          !-------------------------------------------------------------
          call this%exponentiate_tridiagonal(dt,y)

          !-------------------------------------------------------------
          ! Transform the propagated Krylov-space vector back into
          ! the full tensor-network Hilbert space.
          !-------------------------------------------------------------
          call this%reconstruct_state(y)

          !-------------------------------------------------------------
          ! Compress the resulting MPS.
          !
          ! The MPS class owns the compression operation, so the
          ! Krylov class should use that operation rather than
          ! duplicating the tensor-network compression machinery.
          !-------------------------------------------------------------
          call this%psi%compress()

          !-------------------------------------------------------------
          ! Re-normalize the MPS after compression.
          !-------------------------------------------------------------
          call this%psi%normalize()

          !--------------------------------------------------------------
          ! The state is now at timestep istep:
          !
          !     t = tstart + istep*dt
          !
          ! Calculate observables at this time.
          !
          ! The final timestep is always evaluated. For intermediate
          ! timesteps, output_stride controls the output frequency.
          !--------------------------------------------------------------
          time = this%parameters%tstart+dble(istep)*dt

          if(mod(istep,this%parameters%output_stride).eq.0.or.istep.eq.nsteps)then
           call this%compute_observables(istep,time)
          endif

          !-------------------------------------------------------------
          ! Optional progress reporting.
          !-------------------------------------------------------------
          if(this%parameters%verbose) then
           if(mpi_id.eq.0)then
            write(*,*) '     Krylov timestep : ',istep, &
                       ' / ',nsteps, &
                       '   dimension : ',this%current_dimension
            flush(6)
           endif
          endif
         enddo

         !---------------------------------------------------------------
         ! Final propagation report.
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine propagate

        subroutine build_initial_vector(this)
        implicit none
        class(krylov_space)            :: this

         !---------------------------------------------------------------
         ! Construct the first Lanczos basis vector.
         !
         ! The initial Krylov vector is the current physical state:
         !
         !     |v_1> = |psi> / ||psi||
         !
         ! Explicit normalisation is performed here even if the MPS
         ! was already normalised during initialisation. This protects
         ! the Lanczos procedure against accumulated numerical error
         ! when propagation is restarted.
         !---------------------------------------------------------------

         if(.not.allocated(this%V))then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov basis has not been allocated'
          endif
          stop
         endif

         if(this%parameters%dimension<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Invalid Krylov dimension'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Copy the current state into the first Krylov vector.
         !---------------------------------------------------------------
         this%V(1) = this%psi

         !---------------------------------------------------------------
         ! Normalise the first Krylov vector.
         !---------------------------------------------------------------
         call this%V(1)%normalize()

         !---------------------------------------------------------------
         ! Reset the Lanczos iteration state.
         !---------------------------------------------------------------
         this%current_dimension = 1

         this%alpha = (0.0d0,0.0d0)
         this%beta = (0.0d0,0.0d0)

        return
        end subroutine build_initial_vector

        subroutine build_krylov_space(this)
        implicit none
        class(krylov_space)            :: this
        integer                        :: m
        integer                        :: max_restart
        integer                        :: restart_count
        integer                        :: i
        logical                        :: converged
        logical                        :: breakdown
        double precision               :: beta_value
        double precision               :: residual

         !---------------------------------------------------------------
         ! Construct the Lanczos Krylov space.
         !
         ! The Krylov space is
         !
         !  K_m(H,psi) = span{|psi>, H|psi>, H^2|psi>, ..., H^(m-1)|psi>}.
         !
         ! The basis is constructed recursively using the Lanczos
         ! three-term recurrence.
         !
         ! The routine also handles:
         !
         !   * user-defined Krylov dimension
         !   * early convergence
         !   * Lanczos breakdown
         !   * optional restarting
         !   * optional adaptive restarting
         !   * optional full/selective reorthogonalisation
         !---------------------------------------------------------------

         m = this%parameters%dimension

         if(m<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Invalid Krylov dimension'
          endif
          stop
         endif

         max_restart = this%parameters%max_restart

         if(max_restart<0)then
          if(mpi_id==0)then
           write(6,*) 'Err: max_restart must be >= 0'
          endif
          stop
         endif

         restart_count = 0
         converged = .false.

         !---------------------------------------------------------------
         ! The basis has already been initialised with v_1.
         !---------------------------------------------------------------
         do
          !-------------------------------------------------------------
          ! Construct the remaining Lanczos vectors.
          !-------------------------------------------------------------
          do i=this%current_dimension,m
           if(i>1)then
            if(abs(this%beta(i-1))<=this%parameters%breakdown_tol)then
             breakdown = .true.
             exit
            endif
           endif

           call this%lanczos_step(i,breakdown)

           if(breakdown) exit

           !-----------------------------------------------------------
           ! The residual associated with the next Lanczos vector is
           ! proportional to beta(i).
           !
           ! A small beta therefore indicates that the Krylov space has
           ! effectively closed.
           !-----------------------------------------------------------

           residual = abs(this%beta(i))

           if(residual<=this%parameters%tol)then
            converged = .true.
            exit
           endif
          enddo

          !-------------------------------------------------------------
          ! Stop if the requested Krylov space has been constructed.
          !-------------------------------------------------------------
          if(converged) exit

          if(breakdown) exit

          if(this%current_dimension>=m) exit

          !-------------------------------------------------------------
          ! Restarting is only relevant if the current Krylov projection
          ! has not reached the desired convergence criterion.
          !-------------------------------------------------------------
          if(this%parameters%restart<=0) exit

          restart_count = restart_count+1

          if(restart_count>max_restart)then
           if(this%parameters%verbose)then
            if(mpi_id==0)then
             write(6,*) '     Maximum Krylov restarts reached'
            endif
           endif
           exit
          endif

          !-------------------------------------------------------------
          ! Restart from the current propagated state.
          !
          ! For the initial implementation this is deliberately simple:
          ! rebuild a fresh Krylov space from psi.
          !
          ! adaptive_restart can later modify the restart length based
          ! on convergence behaviour.
          !-------------------------------------------------------------
          if(this%parameters%adaptive_restart)then
           !------------------------------------------------------------
           ! Adaptive restarting will be implemented by modifying the
           ! effective Krylov dimension/restart length based on the
           ! observed convergence behaviour.
           !
           ! For now the current restart is retained.
           !------------------------------------------------------------
          endif
          call this%build_initial_vector()
         enddo

         if(this%parameters%verbose.and.mpi_id==0)then
          write(6,*) '     Krylov dimension constructed : ',this%current_dimension

          if(converged.and.mpi_id==0)then
           write(6,*) '     Krylov convergence achieved'
          endif

          if(breakdown.and.this%current_dimension==m.and.mpi_id==0)then
           write(6,*) '     Lanczos invariant subspace closure achieved'
          endif

          if(breakdown.and.this%current_dimension/=m.and.mpi_id==0)then
           write(6,*) '     Lanczos breakdown detected'
          endif
          flush(6)
         endif

        return
        end subroutine build_krylov_space

        subroutine lanczos_step(this,i,breakdown)
        implicit none
        class(krylov_space)            :: this
        integer                        :: i,j
        logical                        :: breakdown
        type(mps_state)                :: w
        complex(8)                     :: alpha_i
        complex(8)                     :: overlap
        double precision               :: beta_i

         !---------------------------------------------------------------
         ! Perform one Lanczos iteration.
         !
         ! Starting from |v_i>, calculate
         !
         !                 |w> = H |v_i>
         !
         ! and remove the components parallel to |v_i> and |v_(i-1)>.
         !
         ! The remaining norm defines beta_i and the next basis vector.
         !---------------------------------------------------------------
         breakdown = .false.

         !---------------------------------------------------------------
         ! Apply the Hamiltonian MPO.
         !---------------------------------------------------------------
         call this%apply_H(this%V(i),w)

         !---------------------------------------------------------------
         ! Calculate the diagonal Lanczos coefficient:
         !
         !              alpha_i = <v_i|H|v_i>.
         !---------------------------------------------------------------
         alpha_i = conjg(this%V(i)%inner_product(w))
         this%alpha(i) = alpha_i

         !---------------------------------------------------------------
         ! Remove alpha_i |v_i>.
         !---------------------------------------------------------------
         call w%lin_com((1.0d0,0.0d0),this%V(i),-alpha_i,w)

         !---------------------------------------------------------------
         ! Remove beta_(i-1) |v_(i-1)>.
         !---------------------------------------------------------------
         if(i>1)then
          call w%lin_com((1.0d0,0.0d0),this%V(i-1),-this%beta(i-1),w)
         endif

         !---------------------------------------------------------------
         ! Compress w to avoid large bond dimensions.
         !---------------------------------------------------------------
         call w%compress()

         !---------------------------------------------------------------
         ! Optional reorthogonalisation.
         !---------------------------------------------------------------
         if(this%parameters%reorthogonalization/=0)then
          call this%orthogonalize(w,i)
         endif

         !---------------------------------------------------------------
         ! Calculate beta_i.
         !---------------------------------------------------------------
         overlap = w%inner_product(w)
         beta_i = sqrt(max(0.0d0,real(overlap)))
         this%beta(i) = cmplx(beta_i,0.0d0,kind=8)

         !---------------------------------------------------------------
         ! A sufficiently small beta indicates exact or numerical
         ! termination of the Krylov sequence.
         !---------------------------------------------------------------
         if(beta_i<=this%parameters%tol)then
          this%current_dimension = i
          if(beta_i<=this%parameters%breakdown_tol)then
           breakdown = .true.
          endif
          return
         endif

         !---------------------------------------------------------------
         ! Construct the next Lanczos vector.
         !---------------------------------------------------------------
         if(i<this%parameters%dimension)then
          call w%rescale(cmplx(1.0d0/beta_i,0.0d0,8))
          this%V(i+1) = w
          this%current_dimension = i+1
         endif

        return
        end subroutine lanczos_step

        subroutine orthogonalize(this,w,i)
        implicit none
        class(krylov_space)            :: this
        type(mps_state)                :: w
        integer                        :: i
        integer                        :: j
        complex(8)                     :: overlap

         !---------------------------------------------------------------
         ! Reorthogonalise the candidate Lanczos vector against the
         ! previously generated Krylov basis.
         !
         ! reorthogonalization =
         !
         !     0 : no reorthogonalisation
         !     1 : full reorthogonalisation
         !     2 : selective reorthogonalisation
         !
         ! The selective implementation currently uses orth_tol as the
         ! threshold for deciding whether a basis vector has a significant
         ! overlap with the candidate.
         !---------------------------------------------------------------

         if(i<=1) return

         select case(this%parameters%reorthogonalization)
          !--------------------------------------------------------------
          ! No reorthogonalisation.
          !--------------------------------------------------------------
          case(0)
           return

          !--------------------------------------------------------------
          ! Full reorthogonalisation.
          !--------------------------------------------------------------
          case(1)
           do j=1,i-1
            overlap = conjg(this%V(j)%inner_product(w))
            if(abs(overlap)>this%parameters%orth_tol)then
             call w%lin_com((1.0d0,0.0d0),this%V(j),-overlap,w)
            endif
           enddo

          !--------------------------------------------------------------
          ! Selective reorthogonalisation.
          !
          ! Only basis vectors having an overlap above orth_tol are
          ! explicitly removed.
          !--------------------------------------------------------------
          case(2)
           do j=1,i-1
            overlap = conjg(this%V(j)%inner_product(w))

            if(abs(overlap) > this%parameters%orth_tol)then
             call w%lin_com((1.0d0,0.0d0),this%V(j),-overlap,w)
            endif
           enddo

          case default
           write(6,*) 'Err: Invalid Lanczos reorthogonalization option'
           stop

         end select

        return
        end subroutine orthogonalize

        subroutine apply_H(this,state,result)
        implicit none
        class(krylov_space)            :: this
        type(mps_state)                :: state
        type(mps_state)                :: result

         !---------------------------------------------------------------
         ! Apply the Hamiltonian MPO to an MPS.
         !
         !                 |result> = H|state>
         !
         ! The tensor-network contraction itself is handled by the MPS
         ! implementation.
         !---------------------------------------------------------------

         call state%apply_mpo(this%H,result)

        return
        end subroutine apply_H

        subroutine diagonalize_tridiagonal(this)
        implicit none
        class(krylov_space)            :: this
        integer                        :: m
        integer                        :: i
        integer                        :: info
        integer                        :: lwork
        complex(8), allocatable        :: T(:,:)
        complex(8), allocatable        :: work(:)
        double precision, allocatable  :: rwork(:)

         !---------------------------------------------------------------
         ! Diagonalise the Lanczos tridiagonal matrix.
         !
         ! The matrix is constructed from the Lanczos coefficients:
         !
         !                 T = Z Lambda Z^dagger
         !
         ! The eigenvalues and eigenvectors are stored in the Krylov
         ! object and reused during subsequent propagation steps.
         !---------------------------------------------------------------
         m = this%current_dimension

         if(m<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Cannot diagonalize empty Krylov space'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Remove any previous eigendecomposition.
         !---------------------------------------------------------------
         if(allocated(this%eigenvalue)) deallocate(this%eigenvalue)
         if(allocated(this%eigenvector)) deallocate(this%eigenvector)

         allocate(this%eigenvalue(m))
         allocate(this%eigenvector(m,m))

         allocate(T(m,m))
         allocate(rwork(max(1,3*m-2)))

         T = (0.0d0,0.0d0)

         do i=1,m
          T(i,i) = cmplx(real(this%alpha(i)),0.0d0,kind=8)

          if(i<m)then
           T(i,i+1) = this%beta(i)
           T(i+1,i) = conjg(this%beta(i))
          endif
         enddo

         this%eigenvector = T

         !---------------------------------------------------------------
         ! Workspace query for ZHEEV.
         !---------------------------------------------------------------
         lwork = -1
         allocate(work(1))

         call zheev('V','U',m,this%eigenvector,m,this%eigenvalue,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: ZHEEV workspace query failed, INFO=',info
          endif
          stop
         endif

         lwork = max(1,int(real(work(1))))
         deallocate(work)
         allocate(work(lwork))

         !---------------------------------------------------------------
         ! Diagonalise T.
         !
         ! On exit:
         !
         !     this%eigenvalue   = Lambda
         !     this%eigenvector  = Z
         !---------------------------------------------------------------
         call zheev('V','U',m,this%eigenvector,m,this%eigenvalue,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: ZHEEV failed, INFO=',info
          endif
          stop
         endif

         deallocate(T)
         deallocate(work)
         deallocate(rwork)

        return
        end subroutine diagonalize_tridiagonal

        subroutine exponentiate_tridiagonal(this,dt,y)
        implicit none
        class(krylov_space)            :: this
        double precision               :: dt
        complex(8), allocatable        :: y(:)
        integer                        :: m
        integer                        :: i
        complex(8), allocatable        :: coeff(:)
        complex(8), allocatable        :: ynew(:)
        complex(8)                     :: phase

         !---------------------------------------------------------------
         ! Exponentiate the projected Lanczos Hamiltonian.
         !
         ! The tridiagonal matrix is diagonalized as:
         !
         !                 T = Z Lambda Z^dagger
         !
         ! and therefore
         !
         !                 exp(-i T dt)
         !
         !              =  Z exp(-i Lambda dt) Z^dagger.
         !
         ! The projected coefficient vector is advanced according to
         !
         !                 y_new = Z exp(-i Lambda dt) Z^dagger y.
         !
         ! The incoming vector y is the current projected state.
         !---------------------------------------------------------------

         m = this%current_dimension

         if(m<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Cannot exponentiate empty Krylov space'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Check that the eigendecomposition is available.
         !---------------------------------------------------------------
         if(.not.allocated(this%eigenvalue).or..not.allocated(this%eigenvector))then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov tridiagonal matrix has not been diagonalized'
          endif
          stop
         endif

         if(size(this%eigenvalue)/=m.or.size(this%eigenvector,1)/=m.or.size(this%eigenvector,2)/=m)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov eigendecomposition dimension mismatch'
          endif
          stop
         endif

         if(.not.allocated(y))then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov coefficient vector is not allocated'
          endif
          stop
         endif

         if(size(y)/=m)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov coefficient dimension mismatch'
          endif
          stop
         endif

         allocate(coeff(m),ynew(m))

         !---------------------------------------------------------------
         ! Transform the current coefficients into the eigenbasis.
         !---------------------------------------------------------------
         coeff = matmul(conjg(transpose(this%eigenvector)),y)

         !---------------------------------------------------------------
         ! Apply the time-evolution phase factors.
         !---------------------------------------------------------------
         do i=1,m
          phase = exp(-(0.0d0,1.0d0)*this%eigenvalue(i)*dt*2.0d0*acos(-1.0d0)/hplank)
          coeff(i) = phase*coeff(i)
         enddo

         !---------------------------------------------------------------
         ! Transform back into the Lanczos basis.
         !---------------------------------------------------------------
         ynew = matmul(this%eigenvector,coeff)

         y = ynew

         deallocate(coeff)
         deallocate(ynew)

        return
        end subroutine exponentiate_tridiagonal

        subroutine reconstruct_state(this,y)
        implicit none
        class(krylov_space)            :: this
        complex(8)                     :: y(:)
        integer                        :: i
        type(mps_state)                :: result

         !---------------------------------------------------------------
         ! Reconstruct the propagated MPS from the Krylov basis:
         !
         !        |psi(t+dt)> = sum_j y_j |v_j>
         !
         ! The resulting MPS is subsequently compressed to control
         ! tensor-network bond growth.
         !---------------------------------------------------------------

         if(size(y)/=this%current_dimension)then
          if(mpi_id==0)then
           write(6,*) 'Err: Krylov coefficient dimension mismatch'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Start with the first Krylov vector.
         !---------------------------------------------------------------
         result = this%V(1)
         call result%rescale(y(1))

         !---------------------------------------------------------------
         ! Add the remaining Krylov vectors, compressing after
         ! each addition to prevent direct-sum bond growth.
         !---------------------------------------------------------------
         do i=2,this%current_dimension
          call result%lin_com((1.0d0,0.0d0),this%V(i),y(i),result)
          call result%compress()
         enddo

         !---------------------------------------------------------------
         ! Replace the current state.
         !---------------------------------------------------------------
         this%psi = result

        return
        end subroutine reconstruct_state

        subroutine construct_observables(this,nspins,spins)
        implicit none
        class(krylov_space)            :: this
        integer                        :: nspins,spins(:)
        integer                        :: i,j
        integer                        :: kind,d
        integer                        :: t1,t2,rate
        integer                        :: bond_dimension
        logical                        :: all_spins
        logical,allocatable            :: observable_site(:)

         !---------------------------------------------------------------
         ! Construct the Matrix Product Operators formed from the Sx, Sy
         ! and Sz operators of the input observable spins.
         !
         !        Sx_tot = sum_{i in observable spins} Sx_i
         !        Sy_tot = sum_{i in observable spins} Sy_i
         !        Sz_tot = sum_{i in observable spins} Sz_i
         !
         ! The MPOs are constructed by taking the spin operator matrices
         ! for each spin kind and adding them to the overall MPO at the
         ! chosen spin sites.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Spin Observable MPOs'
          flush(6)
         endif

         !---------------------------------------------------------------
         ! Check the number of observable spins.
         !---------------------------------------------------------------
         if(nspins.le.0)then
          if(mpi_id.eq.0)then
           write(6,*) 'Err: Number of observable spins must be greater than zero'
          endif
          stop
         endif

         if(size(spins).lt.nspins)then
          if(mpi_id.eq.0)then
           write(6,*) 'Err: Observable spin array is too small'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Determine whether the "all spins" condition has been
         ! requested.
         !---------------------------------------------------------------
         all_spins = .false.

         if(nspins.eq.1)then
          if(spins(1).eq.-1) all_spins = .true.
         endif

         !---------------------------------------------------------------
         ! Create a logical mask identifying the sites which contribute
         ! to the observable MPOs.
         !
         ! observable_site(i) = .true.
         !
         ! means that Sx(i), Sy(i) and Sz(i) are included.
         !---------------------------------------------------------------
         allocate(observable_site(this%psi%nsites))
         observable_site = .false.

         if(all_spins)then
          observable_site = .true.

         else
          do i=1,nspins
           if(spins(i).lt.1.or.spins(i).gt.this%psi%nsites)then
            if(mpi_id.eq.0)then
             write(6,*) 'Err: Invalid observable spin index : ', spins(i)
             write(6,*) '     Valid spin indices are 1 through ', this%psi%nsites
            endif
            stop
           endif
           observable_site(spins(i)) = .true.
          enddo
         endif

         !---------------------------------------------------------------
         ! Each observable MPO acts on the complete MPS, but represents
         ! only the sum of the corresponding spin operator over the
         ! selected observable sites.
         !---------------------------------------------------------------
         this%Sx%nsites = this%psi%nsites
         this%Sy%nsites = this%psi%nsites
         this%Sz%nsites = this%psi%nsites

         if (allocated(this%Sx%site)) deallocate(this%Sx%site)
         allocate(this%Sx%site(this%Sx%nsites))
         if (allocated(this%Sy%site)) deallocate(this%Sy%site)
         allocate(this%Sy%site(this%Sy%nsites))
         if (allocated(this%Sz%site)) deallocate(this%Sz%site)
         allocate(this%Sz%site(this%Sz%nsites))

         !---------------------------------------------------------------
         ! The observable MPOs use exactly the same local spin bases as
         ! the MPS.
         !---------------------------------------------------------------
         this%Sx%basis = this%psi%basis
         this%Sy%basis = this%psi%basis
         this%Sz%basis = this%psi%basis

         !---------------------------------------------------------------
         ! A sum of one-site operators has an exact MPO representation
         ! with bond dimension two:
         !
         !                 [ I   O_i ]
         !        W_i =    [ 0    I  ]
         !
         ! with appropriate boundary tensors.
         !
         ! For a one-site system the MPO has bond dimension one.
         !---------------------------------------------------------------
         if(this%psi%nsites.eq.1)then
          bond_dimension=1
         else
          bond_dimension=2
         endif

         this%Sx%bond_dimension = bond_dimension
         this%Sy%bond_dimension = bond_dimension
         this%Sz%bond_dimension = bond_dimension

         this%Sx%bond_dimension_max = bond_dimension
         this%Sy%bond_dimension_max = bond_dimension
         this%Sz%bond_dimension_max = bond_dimension

         this%Sx%trunc_error = 0.0d0
         this%Sy%trunc_error = 0.0d0
         this%Sz%trunc_error = 0.0d0

         this%Sx%is_hermitian = .true.
         this%Sy%is_hermitian = .true.
         this%Sz%is_hermitian = .true.

         !---------------------------------------------------------------
         ! Construct the local MPO tensors.
         !
         ! The MPO represents
         !
         !        O = sum_i O_i
         !
         ! where O_i is Sx, Sy or Sz if site i is selected, and zero
         ! otherwise.
         !
         ! For an interior site:
         !
         !              [ I   O_i ]
         !        W_i = [ 0    I  ].
         !
         ! The upper-right path inserts the local operator, while the
         ! diagonal identity path propagates through the chain.
         !---------------------------------------------------------------
         do i=1,this%psi%nsites
          kind = this%psi%site(i)%spin_kind
          d = this%psi%basis%kind(kind)%d

          !-------------------------------------------------------------
          ! Store common MPO tensor metadata.
          !-------------------------------------------------------------
          this%Sx%site(i)%site_index = i
          this%Sy%site(i)%site_index = i
          this%Sz%site(i)%site_index = i

          this%Sx%site(i)%spin_kind = kind
          this%Sy%site(i)%spin_kind = kind
          this%Sz%site(i)%spin_kind = kind

          this%Sx%site(i)%phys_in = d
          this%Sy%site(i)%phys_in = d
          this%Sz%site(i)%phys_in = d

          this%Sx%site(i)%phys_out = d
          this%Sy%site(i)%phys_out = d
          this%Sz%site(i)%phys_out = d

          !-------------------------------------------------------------
          ! One-site system.
          !-------------------------------------------------------------
          if(this%psi%nsites.eq.1)then
           allocate(this%Sx%site(i)%W(1,d,d,1))
           allocate(this%Sy%site(i)%W(1,d,d,1))
           allocate(this%Sz%site(i)%W(1,d,d,1))

           this%Sx%site(i)%W=(0.0d0,0.0d0)
           this%Sy%site(i)%W=(0.0d0,0.0d0)
           this%Sz%site(i)%W=(0.0d0,0.0d0)

           this%Sx%site(i)%left_dim=1
           this%Sx%site(i)%right_dim=1

           this%Sy%site(i)%left_dim=1
           this%Sy%site(i)%right_dim=1

           this%Sz%site(i)%left_dim=1
           this%Sz%site(i)%right_dim=1

           if(observable_site(i))then
            this%Sx%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sx)
            this%Sy%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sy)
            this%Sz%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sz)
           endif

          !-------------------------------------------------------------
          ! First site of a multi-site MPO.
          !
          !        W_1 = [ I  O_1 ]
          !-------------------------------------------------------------
          elseif(i.eq.1)then
           allocate(this%Sx%site(i)%W(1,d,d,2))
           allocate(this%Sy%site(i)%W(1,d,d,2))
           allocate(this%Sz%site(i)%W(1,d,d,2))

           this%Sx%site(i)%W=(0.0d0,0.0d0)
           this%Sy%site(i)%W=(0.0d0,0.0d0)
           this%Sz%site(i)%W=(0.0d0,0.0d0)

           this%Sx%site(i)%left_dim=1
           this%Sx%site(i)%right_dim=2

           this%Sy%site(i)%left_dim=1
           this%Sy%site(i)%right_dim=2

           this%Sz%site(i)%left_dim=1
           this%Sz%site(i)%right_dim=2

           !------------------------------------------------------------
           ! Identity propagation.
           !------------------------------------------------------------
           this%Sx%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sy%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sz%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I

           !------------------------------------------------------------
           ! Operator insertion.
           !
           ! A non-observable site has a zero insertion tensor.
           !------------------------------------------------------------
           if(observable_site(i))then
            this%Sx%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sx)
            this%Sy%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sy)
            this%Sz%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sz)
           endif

          !-------------------------------------------------------------
          ! Last site of a multi-site MPO.
          !
          !        W_N = [ O_N ]
          !              [ I   ]
          !-------------------------------------------------------------
          elseif(i.eq.this%psi%nsites)then
           allocate(this%Sx%site(i)%W(2,d,d,1))
           allocate(this%Sy%site(i)%W(2,d,d,1))
           allocate(this%Sz%site(i)%W(2,d,d,1))

           this%Sx%site(i)%W=(0.0d0,0.0d0)
           this%Sy%site(i)%W=(0.0d0,0.0d0)
           this%Sz%site(i)%W=(0.0d0,0.0d0)

           this%Sx%site(i)%left_dim=2
           this%Sx%site(i)%right_dim=1

           this%Sy%site(i)%left_dim=2
           this%Sy%site(i)%right_dim=1

           this%Sz%site(i)%left_dim=2
           this%Sz%site(i)%right_dim=1

           !------------------------------------------------------------
           ! Completed-operator path.
           !------------------------------------------------------------
           if(observable_site(i))then
            this%Sx%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sx)
            this%Sy%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sy)
            this%Sz%site(i)%W(1,:,:,1) = transpose(this%psi%basis%kind(kind)%Sz)
           endif

           !------------------------------------------------------------
           ! Identity path.
           !------------------------------------------------------------
           this%Sx%site(i)%W(2,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sy%site(i)%W(2,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sz%site(i)%W(2,:,:,1) = this%psi%basis%kind(kind)%I

          !-------------------------------------------------------------
          ! Interior site.
          !
          !             [ I   O_i ]
          !        W_i= [ 0    I  ]
          !-------------------------------------------------------------
          else
           allocate(this%Sx%site(i)%W(2,d,d,2))
           allocate(this%Sy%site(i)%W(2,d,d,2))
           allocate(this%Sz%site(i)%W(2,d,d,2))

           this%Sx%site(i)%W=(0.0d0,0.0d0)
           this%Sy%site(i)%W=(0.0d0,0.0d0)
           this%Sz%site(i)%W=(0.0d0,0.0d0)

           this%Sx%site(i)%left_dim=2
           this%Sx%site(i)%right_dim=2

           this%Sy%site(i)%left_dim=2
           this%Sy%site(i)%right_dim=2

           this%Sz%site(i)%left_dim=2
           this%Sz%site(i)%right_dim=2

           !------------------------------------------------------------
           ! Identity propagation on the two diagonal paths.
           !------------------------------------------------------------
           this%Sx%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sx%site(i)%W(2,:,:,2) = this%psi%basis%kind(kind)%I
           this%Sy%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sy%site(i)%W(2,:,:,2) = this%psi%basis%kind(kind)%I
           this%Sz%site(i)%W(1,:,:,1) = this%psi%basis%kind(kind)%I
           this%Sz%site(i)%W(2,:,:,2) = this%psi%basis%kind(kind)%I

           !------------------------------------------------------------
           ! Operator insertion path.
           !------------------------------------------------------------
           if(observable_site(i))then
            this%Sx%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sx)
            this%Sy%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sy)
            this%Sz%site(i)%W(1,:,:,2) = transpose(this%psi%basis%kind(kind)%Sz)
           endif
          endif
         enddo

         deallocate(observable_site)

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine construct_observables

        subroutine compute_observables(this,step_no,time)
        implicit none
        class(krylov_space)            :: this
        double precision               :: time
        integer                        :: step_no
        complex(8)                     :: mx
        complex(8)                     :: my
        complex(8)                     :: mz
        complex(8)                     :: norm
        integer                        :: unit_no
        integer                        :: ios

         !---------------------------------------------------------------
         ! Calculate expectation values of the spin observables.
         !
         ! The observable MPOs were constructed once before propagation:
         !
         !        Sx = sum_i Sx_i
         !        Sy = sum_i Sy_i
         !        Sz = sum_i Sz_i
         !
         ! where the sum runs over the requested observable spin sites.
         !
         ! For each observable:
         !
         !        Mx = <psi|S_x|psi>
         !        My = <psi|S_y|psi>
         !        Mz = <psi|S_z|psi>
         !
         ! The MPS norm is also calculated:
         !
         !        norm = <psi|psi>
         !
         ! The results are written to:
         !
         !        Mx_dynamics.dat
         !        My_dynamics.dat
         !        Mz_dynamics.dat
         !
         ! with columns:
         !
         !        step   time   expectation_value   norm
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Check that the observable MPOs have been constructed.
         !---------------------------------------------------------------
         if(.not.allocated(this%Sx%site))then
          if(mpi_id.eq.0)then
           write(6,*) 'Err: Sx observable MPO has not been constructed'
          endif
          stop
         endif

         if(.not.allocated(this%Sy%site))then
          if(mpi_id.eq.0)then
           write(6,*) 'Err: Sy observable MPO has not been constructed'
          endif
          stop
         endif

         if(.not.allocated(this%Sz%site))then
          if(mpi_id.eq.0)then
           write(6,*) 'Err: Sz observable MPO has not been constructed'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Calculate the norm.
         !---------------------------------------------------------------
         norm = this%psi%inner_product(this%psi)

         if(abs(norm)==0.0d0)then
          if(mpi_id.eq.0) then
           write(6,*) 'Err: Cannot calculate observables for zero MPS'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Calculate the physical expectation values.
         !---------------------------------------------------------------
         mx = this%psi%expectation_value(this%Sx)
         my = this%psi%expectation_value(this%Sy)
         mz = this%psi%expectation_value(this%Sz)

         !---------------------------------------------------------------
         ! Sx, Sy and Sz are Hermitian, so the expectation values should
         ! be real. Discard numerical roundoff in the imaginary parts.
         !---------------------------------------------------------------
         mx = cmplx(real(mx,kind=8),0.0d0,kind=8)
         my = cmplx(real(my,kind=8),0.0d0,kind=8)
         mz = cmplx(real(mz,kind=8),0.0d0,kind=8)

         !---------------------------------------------------------------
         ! Write from the MPI root process.
         !
         ! Output format:
         !
         !     timestep   time   M_alpha   norm
         !---------------------------------------------------------------
         if(mpi_id.eq.0)then
          inquire(file='Mx_dynamics.dat',number=unit_no)
          if(unit_no.eq.-1) open(11,file='Mx_dynamics.dat')

          inquire(file='My_dynamics.dat',number=unit_no)
          if(unit_no.eq.-1) open(12,file='My_dynamics.dat')

          inquire(file='Mz_dynamics.dat',number=unit_no)
          if(unit_no.eq.-1) open(13,file='Mz_dynamics.dat')

          write(11,*) step_no,time,real(mx,kind=8),real(norm,kind=8)
          flush(11)

          write(12,*) step_no,time,real(my,kind=8),real(norm,kind=8)
          flush(12)

          write(13,*) step_no,time,real(mz,kind=8),real(norm,kind=8)
          flush(13)
         endif

        return
        end subroutine compute_observables

        subroutine finalize(this)
        implicit none
        class(krylov_space)            :: this

         !---------------------------------------------------------------
         ! Release Krylov-specific workspace.
         !
         ! The object can subsequently be reused by calling initialize
         ! again.
         !---------------------------------------------------------------
         if(allocated(this%V))then
          deallocate(this%V)
         endif

         if(allocated(this%alpha))then
          deallocate(this%alpha)
         endif

         if(allocated(this%beta))then
          deallocate(this%beta)
         endif

         this%current_dimension = 0

         !---------------------------------------------------------------
         ! Reset the Hamiltonian and propagated state.
         !---------------------------------------------------------------
         this%H%nsites = 0
         this%psi%nsites = 0

         if(mpi_id.eq.0)then
          write(6,*) '     Krylov workspace finalised'
          flush(6)
         endif

        return
        end subroutine finalize

        end module krylov_class
