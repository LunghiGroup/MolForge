        module tensor_network_class
        use hilbert_dist_class
        implicit none

        !===============================================================
        ! Local spin basis for a single spin kind.
        !
        ! The spin operators are stored in the |m> eigenbasis of Sz.
        !===============================================================
        type                           :: spin_basis
         !---------------------------------------------------------------
         ! Twice the spin quantum number.
         !
         !     S2 = 2S
         !
         ! is stored to allow both integer and half-integer spins to be
         ! represented using an integer.
         !---------------------------------------------------------------
         integer                       :: S2=0
         !---------------------------------------------------------------
         ! Local Hilbert-space dimension.
         !
         !     d = 2S + 1
         !---------------------------------------------------------------
         integer                       :: d=0
         !---------------------------------------------------------------
         ! Magnetic quantum numbers
         !
         !     m = -S, ..., S.
         !---------------------------------------------------------------
         double precision, allocatable :: m(:)
         !---------------------------------------------------------------
         ! Local spin operators.
         !---------------------------------------------------------------
         complex(8), allocatable       :: I(:,:)
         complex(8), allocatable       :: Sx(:,:)
         complex(8), allocatable       :: Sy(:,:)
         complex(8), allocatable       :: Sz(:,:)
         complex(8), allocatable       :: Sp(:,:)
         complex(8), allocatable       :: Sm(:,:)
        end type

        !===============================================================
        ! Collection of local spin bases.
        !
        ! One basis is stored for each spin kind appearing in the spin
        ! system.
        !===============================================================
        type                           :: spin_basis_set
         integer                       :: nkinds = 0
         type(spin_basis), allocatable :: kind(:)
         contains
         procedure                     :: bcast_basis
        end type

        !===============================================================
        ! Matrix Product State tensor.
        !
        ! Each tensor has dimensions
        !
        !     A(left,physical,right)
        !
        ! where the virtual indices connect neighbouring tensors and the
        ! physical index corresponds to the local spin basis.
        !===============================================================
        type                           :: mps_tensor
         complex(8), allocatable       :: A(:,:,:)
         integer                       :: left_dim=0
         integer                       :: phys_dim=0
         integer                       :: right_dim=0
         integer                       :: site_index=0
         integer                       :: spin_kind=0
        end type

        !===============================================================
        ! Matrix Product State.
        !
        ! Represents a many-spin wavefunction in tensor-network form.
        !===============================================================
        type                           :: mps_state
         !---------------------------------------------------------------
         ! Local MPS tensors.
         !---------------------------------------------------------------
         type(mps_tensor), allocatable :: site(:)
         !---------------------------------------------------------------
         ! Number of physical spin sites.
         !---------------------------------------------------------------
         integer                       :: nsites=0
         !---------------------------------------------------------------
         ! Maximum permitted virtual bond dimension during compression.
         !---------------------------------------------------------------
         integer                       :: bond_dimension_max=500
         !---------------------------------------------------------------
         ! Largest bond dimension currently present in the MPS.
         !---------------------------------------------------------------
         integer                       :: bond_dimension=1
         !---------------------------------------------------------------
         ! Relative singular-value truncation tolerance.
         !---------------------------------------------------------------
         double precision              :: trunc_tol=1.0d-8
         !---------------------------------------------------------------
         ! Accumulated truncation error from the most recent compression.
         !---------------------------------------------------------------
         double precision              :: trunc_error=0.0d0
         !---------------------------------------------------------------
         ! Position of the orthogonality centre.
         !
         ! A value of zero indicates that the state is not currently in a
         ! canonical form.
         !---------------------------------------------------------------
         integer                       :: canonical_site=0
         !---------------------------------------------------------------
         ! Indicates whether the state has unit norm.
         !---------------------------------------------------------------
         logical                       :: is_normalized=.false.
         !---------------------------------------------------------------
         ! Local spin bases associated with the physical indices.
         !---------------------------------------------------------------
         type(spin_basis_set)          :: basis
         contains
         procedure                     :: bcast_mps
         procedure                     :: init_basis
         procedure                     :: init_product
         procedure                     :: init_superposition
         procedure                     :: rescale
         procedure                     :: inner_product
         procedure                     :: normalize
         procedure                     :: lin_com => lin_com_mps
         procedure                     :: compress => compress_mps
         procedure                     :: apply_mpo
         procedure                     :: expectation_value
        end type

        !===============================================================
        ! Matrix Product Operator tensor.
        !
        ! Each tensor has dimensions
        !
        !     W(left,physical_in,physical_out,right)
        !===============================================================
        type                           :: mpo_tensor
         complex(8), allocatable       :: W(:,:,:,:)
         integer                       :: left_dim=0
         integer                       :: right_dim=0
         integer                       :: phys_in=0
         integer                       :: phys_out=0
         integer                       :: site_index=0
         integer                       :: spin_kind=0
        end type

        !===============================================================
        ! Matrix Product Operator.
        !
        ! Represents an operator acting on an MPS, such as the spin
        ! Hamiltonian.
        !===============================================================
        type                           :: mpo_operator
         !---------------------------------------------------------------
         ! Local MPO tensors.
         !---------------------------------------------------------------
         type(mpo_tensor), allocatable :: site(:)
         !---------------------------------------------------------------
         ! Number of physical spin sites.
         !---------------------------------------------------------------
         integer                       :: nsites=0
         !---------------------------------------------------------------
         ! Maximum permitted virtual bond dimension during compression.
         !---------------------------------------------------------------
         integer                       :: bond_dimension_max=500
         !---------------------------------------------------------------
         ! Largest bond dimension currently present in the MPO.
         !---------------------------------------------------------------
         integer                       :: bond_dimension=1
         !---------------------------------------------------------------
         ! Relative singular-value truncation tolerance.
         !---------------------------------------------------------------
         double precision              :: trunc_tol=1.0d-8
         !---------------------------------------------------------------
         ! Accumulated truncation error from the most recent compression.
         !---------------------------------------------------------------
         double precision              :: trunc_error=0.0d0
         !---------------------------------------------------------------
         ! Indicates whether the operator is Hermitian.
         !---------------------------------------------------------------
         logical                       :: is_hermitian=.true.
         !---------------------------------------------------------------
         ! Local spin bases associated with the physical indices.
         !---------------------------------------------------------------
         type(spin_basis_set)          :: basis
         contains
         procedure                     :: bcast_mpo
         procedure                     :: lin_com => lin_com_mpo
         procedure                     :: compress => compress_mpo
         procedure                     :: compose_mpos
         procedure                     :: kronecker_product
         procedure                     :: init_spin_H
        end type

        contains

        subroutine bcast_basis(this)
        use mpi
        implicit none
        class(spin_basis_set)          :: this
        integer                        :: ierr
        integer                        :: k,d

         !---------------------------------------------------------------
         ! Broadcast the complete set of local spin bases.
         !
         ! Rank 0 owns the initially constructed basis. All other ranks
         ! allocate the required arrays and receive the basis matrices.
         !---------------------------------------------------------------

         call MPI_BCAST(this%nkinds,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

         if(mpi_id/=0)then
          if(allocated(this%kind)) deallocate(this%kind)
          allocate(this%kind(this%nkinds))
         endif

         do k=1,this%nkinds
          !-------------------------------------------------------------
          ! Broadcast the local spin dimensions.
          !-------------------------------------------------------------
          call MPI_BCAST(this%kind(k)%S2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          d = this%kind(k)%d

          !-------------------------------------------------------------
          ! Allocate local basis arrays on non-root processes.
          !-------------------------------------------------------------
          if(mpi_id/=0)then
           allocate(this%kind(k)%m(d))
           allocate(this%kind(k)%I(d,d))
           allocate(this%kind(k)%Sx(d,d))
           allocate(this%kind(k)%Sy(d,d))
           allocate(this%kind(k)%Sz(d,d))
           allocate(this%kind(k)%Sp(d,d))
           allocate(this%kind(k)%Sm(d,d))
          endif

          !-------------------------------------------------------------
          ! Broadcast the magnetic quantum numbers.
          !-------------------------------------------------------------
          call MPI_BCAST(this%kind(k)%m,d,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

          !-------------------------------------------------------------
          ! Broadcast the local spin operators.
          !-------------------------------------------------------------
          call MPI_BCAST(this%kind(k)%I,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%Sx,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%Sy,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%Sz,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%Sp,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%kind(k)%Sm,d*d,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         enddo

        return
        end subroutine bcast_basis

        subroutine bcast_mps(this)
        use mpi
        implicit none
        class(mps_state)               :: this
        integer                        :: ierr
        integer                        :: i
        integer                        :: dims(5)

         !---------------------------------------------------------------
         ! Broadcast a complete Matrix Product State.
         !
         ! Rank 0 owns the constructed MPS. All other ranks receive its
         ! metadata, local spin basis and tensor elements.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Broadcast general MPS metadata.
         !---------------------------------------------------------------
         call MPI_BCAST(this%nsites,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%bond_dimension_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%bond_dimension,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%canonical_site,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

         call MPI_BCAST(this%trunc_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%trunc_error,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

         call MPI_BCAST(this%is_normalized,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

         !---------------------------------------------------------------
         ! Allocate the site array on non-root processes.
         !---------------------------------------------------------------
         if(mpi_id/=0)then
          if(allocated(this%site)) deallocate(this%site)
          allocate(this%site(this%nsites))
         endif

         !---------------------------------------------------------------
         ! Broadcast the local spin basis.
         !---------------------------------------------------------------
         call this%basis%bcast_basis()

         !---------------------------------------------------------------
         ! Broadcast the individual MPS tensors.
         !
         ! dims contains:
         !   1 = left_dim
         !   2 = phys_dim
         !   3 = right_dim
         !   4 = site_index
         !   5 = spin_kind
         !---------------------------------------------------------------
         do i=1,this%nsites
          if(mpi_id==0)then
           dims(1) = this%site(i)%left_dim
           dims(2) = this%site(i)%phys_dim
           dims(3) = this%site(i)%right_dim
           dims(4) = this%site(i)%site_index
           dims(5) = this%site(i)%spin_kind
          endif

          call MPI_BCAST(dims,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          if(mpi_id/=0)then
           this%site(i)%left_dim = dims(1)
           this%site(i)%phys_dim = dims(2)
           this%site(i)%right_dim = dims(3)
           this%site(i)%site_index = dims(4)
           this%site(i)%spin_kind = dims(5)

           allocate(this%site(i)%A(dims(1),dims(2),dims(3)))
          endif

          !-------------------------------------------------------------
          ! Broadcast the complex MPS tensor elements.
          !-------------------------------------------------------------
          call MPI_BCAST(this%site(i)%A,dims(1)*dims(2)*dims(3),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         enddo

        return
        end subroutine bcast_mps

        subroutine bcast_mpo(this)
        use mpi
        implicit none
        class(mpo_operator)            :: this
        integer                        :: ierr
        integer                        :: i
        integer                        :: dims(6)

         !---------------------------------------------------------------
         ! Broadcast a complete Matrix Product Operator.
         !
         ! Rank 0 owns the constructed MPO. All other ranks receive its
         ! metadata, local spin basis and tensor elements.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Broadcast general MPO metadata.
         !---------------------------------------------------------------
         call MPI_BCAST(this%nsites,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%bond_dimension_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%bond_dimension,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

         call MPI_BCAST(this%trunc_tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(this%trunc_error,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

         call MPI_BCAST(this%is_hermitian,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

         !---------------------------------------------------------------
         ! Allocate the site array on non-root processes.
         !---------------------------------------------------------------
         if(mpi_id/=0)then
          if(allocated(this%site)) deallocate(this%site)
          allocate(this%site(this%nsites))
         endif

         !---------------------------------------------------------------
         ! Broadcast the local spin basis.
         !---------------------------------------------------------------
         call this%basis%bcast_basis()

         !---------------------------------------------------------------
         ! Broadcast the individual MPO tensors.
         !
         ! dims contains:
         !   1 = left_dim
         !   2 = phys_in
         !   3 = phys_out
         !   4 = right_dim
         !   5 = site_index
         !   6 = spin_kind
         !---------------------------------------------------------------
         do i=1,this%nsites
          if(mpi_id==0)then
           dims(1) = this%site(i)%left_dim
           dims(2) = this%site(i)%phys_in
           dims(3) = this%site(i)%phys_out
           dims(4) = this%site(i)%right_dim
           dims(5) = this%site(i)%site_index
           dims(6) = this%site(i)%spin_kind
          endif

          call MPI_BCAST(dims,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          if(mpi_id/=0)then
           this%site(i)%left_dim = dims(1)
           this%site(i)%phys_in = dims(2)
           this%site(i)%phys_out = dims(3)
           this%site(i)%right_dim = dims(4)
           this%site(i)%site_index = dims(5)
           this%site(i)%spin_kind = dims(6)

           allocate(this%site(i)%W(dims(1),dims(2),dims(3),dims(4)))
          endif

          !-------------------------------------------------------------
          ! Broadcast the complex MPO tensor elements.
          !-------------------------------------------------------------
          call MPI_BCAST(this%site(i)%W,dims(1)*dims(2)*dims(3)*dims(4),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         enddo

        return
        end subroutine bcast_mpo

        subroutine init_basis(this,spindy,basis_state)
        implicit none
        class(mps_state)               :: this
        class(spins_hilbert)           :: spindy
        integer                        :: basis_state(:),i,level,kind,d,S2,t1,t2,rate

         !---------------------------------------------------------------
         ! Construct a product-state Matrix Product State from a basis
         ! configuration.
         !
         ! Each physical spin site is initialised in one of its local
         ! Sz spin-projections.
         !
         ! The resulting MPS has bond dimension one and therefore
         ! represents an unentangled product state.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Basis State MPS'
          flush(6)
         endif

         !---------------------------------------------------------------
         ! Initialise MPS metadata.
         !---------------------------------------------------------------
         this%bond_dimension = 1
         this%canonical_site = 0
         this%trunc_error = 0.0d0
         this%is_normalized = .true.

         !---------------------------------------------------------------
         ! Allocate one tensor for each physical spin site.
         !---------------------------------------------------------------
         this%nsites = spindy%nspins

         if (allocated(this%site)) deallocate(this%site)
         allocate(this%site(this%nsites))

         if (allocated(this%basis%kind)) deallocate(this%basis%kind)

         this%basis%nkinds = spindy%nkinds
         allocate(this%basis%kind(spindy%nkinds))

         !---------------------------------------------------------------
         ! Construct the local spin basis for every spin kind.
         !
         ! Each basis contains the identity operator together with the
         ! Cartesian and ladder spin operators used throughout the tensor
         ! network library.
         !---------------------------------------------------------------
         do i=1,spindy%nkinds
          call init_spin_basis(this%basis%kind(i),nint(2.0d0*spindy%spin(i)))
         enddo

         !---------------------------------------------------------------
         ! Initialise each local MPS tensor.
         !
         ! The basis_state array stores twice the spin projection
         ! (2*m_S) in order to avoid half-integer arithmetic.
         ! This is converted into the corresponding local basis index.
         !---------------------------------------------------------------
         do i=1,this%nsites
          !!! Assumes kind labels are always ordered as 1,2,... !!!
          kind = spindy%kind(i)
          d = this%basis%kind(kind)%d
          S2 = this%basis%kind(kind)%S2

          !---------------------------------------------------------------
          ! Verify that the requested spin projection is valid for the
          ! local spin basis.
          !---------------------------------------------------------------
          if (abs(basis_state(i))>S2) then
           if (mpi_id==0) then
            write(6,*) 'Err: Invalid spin projection on site ',i
           endif
           stop
          endif

          if (mod(basis_state(i)+S2,2)/=0) then
           if(mpi_id==0) then
            write(6,*) 'Err: Invalid spin projection on site ',i
           endif
           stop
          endif

          this%site(i)%spin_kind = kind
          this%site(i)%site_index = i

          allocate(this%site(i)%A(1,d,1))
          this%site(i)%A = (0.0d0,0.0d0)

          this%site(i)%left_dim = 1
          this%site(i)%phys_dim = d
          this%site(i)%right_dim = 1

          level = (basis_state(i)+S2)/2+1

          this%site(i)%A(1,level,1) = (1.0d0,0.0d0)
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine init_basis

        subroutine init_product(this,spindy,product_state)
        implicit none
        class(mps_state)               :: this
        class(spins_hilbert)           :: spindy
        complex(8)                     :: product_state(:,:)
        integer                        :: i,kind,d,t1,t2,rate

         !---------------------------------------------------------------
         ! Construct a product-state Matrix Product State.
         !
         ! Each physical spin site is initialised in an arbitrary local
         ! spin state expressed in its Sz basis.
         !
         ! The resulting MPS has bond dimension one and therefore
         ! represents an unentangled product state.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Product State MPS'
          flush(6)
         endif

         !---------------------------------------------------------------
         ! Allocate one tensor for each physical spin site.
         !---------------------------------------------------------------
         this%nsites = spindy%nspins

         if (allocated(this%site)) deallocate(this%site)
         allocate(this%site(this%nsites))

         if (allocated(this%basis%kind)) deallocate(this%basis%kind)

         this%basis%nkinds = spindy%nkinds
         allocate(this%basis%kind(spindy%nkinds))

         !---------------------------------------------------------------
         ! Construct the local spin basis for every spin kind.
         !
         ! Each basis contains the identity operator together with the
         ! Cartesian and ladder spin operators used throughout the tensor
         ! network library.
         !---------------------------------------------------------------
         do i=1,spindy%nkinds
          call init_spin_basis(this%basis%kind(i),nint(2.0d0*spindy%spin(i)))
         enddo

         !---------------------------------------------------------------
         ! Verify that every supplied local product state is normalised.
         !
         ! One state vector is supplied for each spin kind and reused for
         ! every site belonging to that kind.
         !---------------------------------------------------------------
         do i=1,spindy%nkinds
          if (abs(sum(abs(product_state(:,i))**2)-1.0d0)>1.0d-6) then
           if (mpi_id==0) then
            write(6,*) 'Err: Product state for spin kind ',spindy%kind(i),' is not normalized'
           endif
           stop
          endif
         enddo

         !---------------------------------------------------------------
         ! Initialise each local MPS tensor.
         !
         ! The supplied local product-state amplitudes are copied directly
         ! into the physical dimension of each site tensor.
         !---------------------------------------------------------------
         do i=1,this%nsites
          kind = spindy%kind(i)
          d = this%basis%kind(kind)%d

          this%site(i)%spin_kind = kind
          this%site(i)%site_index = i

          allocate(this%site(i)%A(1,d,1))
          this%site(i)%A = (0.0d0,0.0d0)

          this%site(i)%left_dim = 1
          this%site(i)%phys_dim = d
          this%site(i)%right_dim = 1

          this%site(i)%A(1,:,1) = product_state(1:d,kind)
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine init_product

        subroutine init_superposition(this,spindy,sup_basis,sup_coeff)
        implicit none
        class(mps_state)               :: this
        class(spins_hilbert)           :: spindy
        integer                        :: sup_basis(:,:)
        complex(8)                     :: sup_coeff(:)
        integer                        :: i,j,kind,d,nstates,S2,level
        integer                        :: t1,t2,rate

         !---------------------------------------------------------------
         ! Construct a Matrix Product State representing a superposition
         ! of product states.
         !
         ! Each product state is specified by a basis configuration and
         ! a complex expansion coefficient.
         !
         ! The resulting MPS has an exact bond dimension equal to the
         ! number of product states in the superposition.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Superposition State MPS'
          flush(6)
         endif

         if(spindy%nspins<=0)then
          if(mpi_id==0) write(6,*) 'Err: Superposition requires at least one spin'
          stop
         endif

         this%nsites = spindy%nspins
         nstates = size(sup_coeff)

         if(nstates <= 0)then
          if(mpi_id == 0) write(6,*) 'Err: Superposition requires at least one basis state'
          stop
         endif

         !---------------------------------------------------------------
         ! Verify that the supplied basis configurations and expansion
         ! coefficients are dimensionally consistent.
         !---------------------------------------------------------------
         if(size(sup_basis,1)/=this%nsites) then
          if (mpi_id==0) then
           write(6,*) 'Err: Superposition basis has incorrect number of spins'
          endif
          stop
         endif

         if(size(sup_basis,2)/=nstates) then
          if (mpi_id==0) then
           write(6,*) 'Err: Superposition coefficient mismatch'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Verify that the superposition coefficients are normalised.
         !---------------------------------------------------------------
         if(abs(sum(abs(sup_coeff)**2)-1.0d0)>1.0d-4) then
          if (mpi_id==0)then
           write(6,*) 'Err: Superposition coefficients are not normalised'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Allocate one tensor for each physical spin site.
         !---------------------------------------------------------------
         if(allocated(this%site)) deallocate(this%site)
         allocate(this%site(this%nsites))

         if(allocated(this%basis%kind)) deallocate(this%basis%kind)

         this%basis%nkinds = spindy%nkinds
         allocate(this%basis%kind(this%basis%nkinds))

         !---------------------------------------------------------------
         ! Construct the local spin basis for every spin kind.
         !
         ! Each basis contains the identity operator together with the
         ! Cartesian and ladder spin operators used throughout the tensor
         ! network library.
         !---------------------------------------------------------------
         do i=1,spindy%nkinds
          call init_spin_basis(this%basis%kind(i),nint(2.0d0*spindy%spin(i)))
         enddo

         !---------------------------------------------------------------
         ! Construct the MPS tensors.
         !
         ! The virtual bond index labels the individual product states.
         ! Interior tensors are therefore block diagonal, while the final
         ! tensor stores the expansion coefficients.
         !---------------------------------------------------------------
         do i=1,this%nsites
          kind = spindy%kind(i)
          d = this%basis%kind(kind)%d
          S2 = this%basis%kind(kind)%S2

          if(kind<1.or.kind>spindy%nkinds)then
           if(mpi_id == 0) write(6,*) 'Err: Invalid spin kind at site ', i
           stop
          endif

          allocate(this%site(i)%A(merge(1,nstates,i==1),d,merge(1,nstates,i==this%nsites)))
          this%site(i)%A=(0.0d0,0.0d0)

          this%site(i)%left_dim = size(this%site(i)%A,1)
          this%site(i)%phys_dim = size(this%site(i)%A,2)
          this%site(i)%right_dim = size(this%site(i)%A,3)

          this%site(i)%site_index = i
          this%site(i)%spin_kind = kind

          do j=1,nstates
           if(abs(sup_basis(i,j))>S2) then
            if (mpi_id==0) then
             write(6,*) 'Err: Invalid spin projection at site ',i
            endif
            stop
           endif

           if(mod(sup_basis(i,j)+S2,2)/=0) then
            if (mpi_id==0) then
             write(6,*) 'Err: Invalid spin projection at site ',i
            endif
            stop
           endif

           level = (sup_basis(i,j)+S2)/2+1

           if(this%nsites==1) then
            this%site(i)%A(1,level,1) = this%site(i)%A(1,level,1)+sup_coeff(j)
           elseif(i==1) then
            this%site(i)%A(1,level,j)=1.0d0
           else if(i==this%nsites) then
            this%site(i)%A(j,level,1) = sup_coeff(j)
           else
            this%site(i)%A(j,level,j) = 1.0d0
           endif
          enddo
         enddo

         this%bond_dimension = 1
         do i=1,this%nsites
          this%bond_dimension = max(this%bond_dimension,this%site(i)%left_dim,this%site(i)%right_dim)
         enddo
         this%bond_dimension_max = max(this%bond_dimension_max,this%bond_dimension)
         this%canonical_site = 0
         this%trunc_error = 0.0d0
         this%is_normalized = .true.

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine init_superposition

        subroutine init_spin_basis(basis,S2)
        implicit none
        type(spin_basis)               :: basis
        integer                        :: S2,i
        double precision               :: S,m,coeff

         !---------------------------------------------------------------
         ! Construct the local spin basis for a single spin species.
         !
         ! The basis is ordered according to increasing Sz projection
         ! from -S to +S and stores the identity together with the
         ! Cartesian and ladder spin operators.
         !---------------------------------------------------------------

         basis%S2 = S2
         basis%d = S2+1

         S = 0.5d0*dble(S2)

         !---------------------------------------------------------------
         ! Allocate the basis vectors and local spin operators.
         !---------------------------------------------------------------
         if (allocated(basis%m)) deallocate(basis%m)
         if (allocated(basis%I)) deallocate(basis%I)
         if (allocated(basis%Sx)) deallocate(basis%Sx)
         if (allocated(basis%Sy)) deallocate(basis%Sy)
         if (allocated(basis%Sz)) deallocate(basis%Sz)
         if (allocated(basis%Sp)) deallocate(basis%Sp)
         if (allocated(basis%Sm)) deallocate(basis%Sm)

         allocate(basis%m(basis%d))
         allocate(basis%I(basis%d,basis%d))
         allocate(basis%Sx(basis%d,basis%d))
         allocate(basis%Sy(basis%d,basis%d))
         allocate(basis%Sz(basis%d,basis%d))
         allocate(basis%Sp(basis%d,basis%d))
         allocate(basis%Sm(basis%d,basis%d))

         basis%I = (0.0d0,0.0d0)
         basis%Sx = (0.0d0,0.0d0)
         basis%Sy = (0.0d0,0.0d0)
         basis%Sz = (0.0d0,0.0d0)
         basis%Sp = (0.0d0,0.0d0)
         basis%Sm = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Construct the basis states together with the identity and
         ! diagonal Sz operator.
         !---------------------------------------------------------------
         do i=1,basis%d
          m = -S+dble(i-1)
          basis%m(i) = m
          basis%I(i,i) = (1.0d0,0.0d0)
          basis%Sz(i,i) = cmplx(m,0.0d0,kind=8)
         enddo

         !---------------------------------------------------------------
         ! Construct the spin raising operator.
         !---------------------------------------------------------------
         do i=1,basis%d-1
          m = -S+dble(i-1)
          coeff = sqrt(S*(S+1.0d0)-m*(m+1.0d0))
          basis%Sp(i+1,i) = cmplx(coeff,0.0d0,kind=8)
         enddo

         !---------------------------------------------------------------
         ! Construct the spin lowering operator.
         !---------------------------------------------------------------
         do i=2,basis%d
          m = -S+dble(i-1)
          coeff = sqrt(S*(S+1.0d0)-m*(m-1.0d0))
          basis%Sm(i-1,i) = cmplx(coeff,0.0d0,kind=8)
         enddo

         !---------------------------------------------------------------
         ! Form the Cartesian spin operators from the ladder operators.
         !---------------------------------------------------------------
         basis%Sx = 0.5d0*(basis%Sp+basis%Sm)
         basis%Sy = -0.5d0*(0.0d0,1.0d0)*(basis%Sp-basis%Sm)

        return
        end subroutine init_spin_basis

        subroutine rescale(this,factor)
        implicit none
        class(mps_state)               :: this
        complex(8)                     :: factor

         !---------------------------------------------------------------
         ! Multiply a Matrix Product State by a complex scalar.
         !
         ! Since an overall scalar factor may be attached to any tensor in
         ! the network without changing the represented state, the factor
         ! is applied to the first site tensor only.
         !---------------------------------------------------------------

         this%site(1)%A = factor*this%site(1)%A

        return
        end subroutine rescale

        complex(8) function inner_product(this,other)
        implicit none
        class(mps_state)               :: this
        class(mps_state)               :: other
        complex(8), allocatable        :: env(:,:),new_env(:,:)
        integer                        :: i,l1,l2,r1,r2,s
        integer                        :: ld1,rd1,ld2,rd2

         !---------------------------------------------------------------
         ! Compute the inner product between two Matrix Product States.
         !
         ! The contraction proceeds sequentially from left to right by
         ! constructing the left environment tensor one site at a time.
         !
         ! Both MPS objects must describe the same physical Hilbert space.
         !
         ! The output inner product is
         !                <other|this>
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Check compatible number of sites
         !---------------------------------------------------------------
         if (this%nsites/=other%nsites) then
          if (mpi_id==0) then
           write(6,*) 'Err: MPS inner product size mismatch'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Sweep from left to right, updating the left environment by
         ! contracting one physical site at a time.
         !---------------------------------------------------------------
         allocate(env(1,1))
         env(1,1) = (1.0d0,0.0d0)

         do i=1,this%nsites
          ld1 = other%site(i)%left_dim
          rd1 = other%site(i)%right_dim

          ld2 = this%site(i)%left_dim
          rd2 = this%site(i)%right_dim

          if (other%site(i)%phys_dim/=this%site(i)%phys_dim) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPS physical dimensions do not match at site ',i
           endif
           stop
          endif

          !-------------------------------------------------------------
          ! Allocate the updated environment after contracting the
          ! current site tensors.
          !-------------------------------------------------------------
          allocate(new_env(rd1,rd2))
          new_env = (0.0d0,0.0d0)

          do r1=1,rd1
           do r2=1,rd2
            do l1=1,ld1
             do l2=1,ld2
              do s=1,this%site(i)%phys_dim
               new_env(r1,r2) = new_env(r1,r2)+conjg(other%site(i)%A(l1,s,r1))*env(l1,l2)*this%site(i)%A(l2,s,r2)
              enddo
             enddo
            enddo
           enddo
          enddo

          deallocate(env)
          call move_alloc(new_env,env)
         enddo

         !---------------------------------------------------------------
         ! After the final contraction the remaining scalar is the MPS
         ! inner product.
         !---------------------------------------------------------------
         inner_product = env(1,1)

         deallocate(env)

        return
        end function inner_product

        subroutine normalize(this)
        implicit none
        class(mps_state)               :: this
        complex(8)                     :: norm
        double precision               :: norm_real
        double precision               :: norm_imag

         !---------------------------------------------------------------
         ! Normalise a Matrix Product State.
         !
         ! The norm is evaluated by contracting the MPS with its Hermitian
         ! conjugate. The entire state is then rescaled so that
         !
         !     <ψ|ψ> = 1.
         !
         ! Since the overall scale of an MPS is arbitrary, the
         ! normalisation factor is applied only to the first site tensor.
         !---------------------------------------------------------------

         if(.not.allocated(this%site))then
          if(mpi_id==0)then
           write(6,*) 'Err: Cannot normalize an unallocated MPS'
          endif
          stop
         endif

         if(this%nsites<=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: Cannot normalize an empty MPS'
          endif
          stop
         endif

         norm = this%inner_product(this)

         norm_real = real(norm,kind=8)
         norm_imag = abs(aimag(norm))

         !---------------------------------------------------------------
         ! Check that the calculated norm is real and positive.
         !---------------------------------------------------------------
         if(norm_imag>1.0d-10*max(1.0d0,abs(norm_real)))then
          if(mpi_id==0)then
           write(6,*) 'Err: Significant imaginary component in MPS norm'
           write(6,*) 'Norm = ',norm
          endif
          stop
         endif

         if(norm_real<=tiny(1.0d0))then
          if(mpi_id==0)then
           write(6,*) 'Err: Cannot normalize zero MPS'
          endif
          stop
         endif

         call this%rescale(cmplx(1.0d0/sqrt(norm_real),0.0d0,kind=8))

         this%is_normalized = .true.

        return
        end subroutine normalize

        subroutine lin_com_mps(this,a,other,b,result)
        implicit none
        class(mps_state)               :: this
        type(mps_state)                :: other
        type(mps_state)                :: result
        type(mps_state)                :: tmp
        complex(8)                     :: a,b
        integer                        :: i,ld1,rd1,ld2,rd2,p

         !---------------------------------------------------------------
         ! Form a linear combination of two Matrix Product States.
         !
         ! The resulting state is
         !
         !     |ψ⟩ = a|ψ₁⟩ + b|ψ₂⟩.
         !
         ! The direct-sum construction is used on every virtual bond,
         ! producing an exact representation of the linear combination
         ! without truncation.
         !---------------------------------------------------------------

         if (this%nsites/=other%nsites) then
          if(mpi_id==0) then
           write(6,*) 'Err: MPS lin_com site mismatch'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Construct the direct-sum MPS representation in a temporary
         ! MPS.
         !
         ! The first and last tensors concatenate a single virtual bond,
         ! while interior tensors become block diagonal in the virtual
         ! indices.
         !---------------------------------------------------------------
         tmp%nsites = this%nsites
         tmp%basis = this%basis

         if(allocated(tmp%site)) deallocate(tmp%site)
         allocate(tmp%site(tmp%nsites))

         do i=1,this%nsites
          if (this%site(i)%spin_kind/=other%site(i)%spin_kind) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPS spin kinds do not match at site ',i
           endif
           stop
          endif

          if (this%site(i)%phys_dim/=other%site(i)%phys_dim) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPS physical dimensions do not match at site ',i
           endif
           stop
          endif

          ld1 = this%site(i)%left_dim
          rd1 = this%site(i)%right_dim

          ld2 = other%site(i)%left_dim
          rd2 = other%site(i)%right_dim

          tmp%site(i)%site_index = i
          tmp%site(i)%spin_kind = this%site(i)%spin_kind

          tmp%site(i)%phys_dim = this%site(i)%phys_dim

          !---------------------------------------------------------------
          ! One-site MPS.
          !
          ! There is no virtual bond. The linear combination is simply
          !
          !     |ψ⟩ = a|ψ₁⟩ + b|ψ₂⟩
          !
          ! represented by a single [1 x 1] virtual tensor.
          !---------------------------------------------------------------
          if(this%nsites==1)then
           allocate(tmp%site(1)%A(1,tmp%site(1)%phys_dim,1))

           tmp%site(1)%A = a*this%site(1)%A+b*other%site(1)%A

           tmp%site(1)%left_dim  = 1
           tmp%site(1)%right_dim = 1

          !-------------------------------------------------------------
          ! First site
          !-------------------------------------------------------------
          elseif(i==1)then
           allocate(tmp%site(i)%A(1,tmp%site(i)%phys_dim,rd1+rd2))
           tmp%site(i)%A=(0.0d0,0.0d0)

           do p=1,tmp%site(i)%phys_dim
            tmp%site(i)%A(1,p,1:rd1) = a*this%site(i)%A(1,p,:)
            tmp%site(i)%A(1,p,rd1+1:rd1+rd2) = b*other%site(i)%A(1,p,:)
           enddo

           tmp%site(i)%left_dim = 1
           tmp%site(i)%right_dim = rd1+rd2

          !-------------------------------------------------------------
          ! Last site
          !-------------------------------------------------------------
          elseif(i==this%nsites)then
           allocate(tmp%site(i)%A(ld1+ld2,tmp%site(i)%phys_dim,1))
           tmp%site(i)%A=(0.0d0,0.0d0)

           do p=1,tmp%site(i)%phys_dim
            tmp%site(i)%A(1:ld1,p,1) = this%site(i)%A(:,p,1)
            tmp%site(i)%A(ld1+1:ld1+ld2,p,1) = other%site(i)%A(:,p,1)
           enddo

           tmp%site(i)%left_dim = ld1+ld2
           tmp%site(i)%right_dim = 1

          !-------------------------------------------------------------
          ! Interior site
          !-------------------------------------------------------------
          else
           allocate(tmp%site(i)%A(ld1+ld2,tmp%site(i)%phys_dim,rd1+rd2))
           tmp%site(i)%A=(0.0d0,0.0d0)

           tmp%site(i)%A(1:ld1,:,1:rd1) = this%site(i)%A
           tmp%site(i)%A(ld1+1:ld1+ld2,:,rd1+1:rd1+rd2) = other%site(i)%A

           tmp%site(i)%left_dim = ld1+ld2
           tmp%site(i)%right_dim = rd1+rd2
          endif
         enddo

         !---------------------------------------------------------------
         ! Update the MPS metadata
         !---------------------------------------------------------------
         tmp%bond_dimension = 1
         do i=1,tmp%nsites
          tmp%bond_dimension = max(tmp%bond_dimension,tmp%site(i)%left_dim,tmp%site(i)%right_dim)
         end do
         tmp%bond_dimension_max = max(this%bond_dimension_max,other%bond_dimension_max)
         tmp%canonical_site = 0
         tmp%is_normalized = .false.

         !---------------------------------------------------------------
         ! Update the set result MPS
         !---------------------------------------------------------------
         result = tmp

        return
        end subroutine lin_com_mps

        subroutine lin_com_mpo(this,a,other,b,result)
        implicit none
        class(mpo_operator)            :: this
        type(mpo_operator)             :: other
        type(mpo_operator)             :: result
        type(mpo_operator)             :: tmp
        complex(8)                     :: a,b
        integer                        :: i,ld1,rd1,ld2,rd2,pin,pout

         !---------------------------------------------------------------
         ! Form a linear combination of two Matrix Product Operators.
         !
         ! The resulting operator is
         !
         !     O = a O_1 + b O_2.
         !
         ! The direct-sum construction is used on every virtual bond,
         ! producing an exact representation of the linear combination
         ! without truncation.
         !---------------------------------------------------------------

         if (this%nsites/=other%nsites) then
          if (mpi_id==0) then
           write(6,*) 'Err: MPO lin_com site mismatch'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Check that both MPOs use the same spin basis.
         !---------------------------------------------------------------
         if(this%basis%nkinds/=other%basis%nkinds)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPO basis mismatch in lin_com_mpo'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Construct the direct-sum MPO representation.
         !
         ! The first and last tensors concatenate a single virtual bond,
         ! while interior tensors become block diagonal in the virtual
         ! indices.
         !---------------------------------------------------------------
         tmp%nsites = this%nsites
         tmp%basis = this%basis

         if(allocated(tmp%site)) deallocate(tmp%site)
         allocate(tmp%site(tmp%nsites))

         do i=1,this%nsites
          ld1 = this%site(i)%left_dim
          rd1 = this%site(i)%right_dim

          ld2 = other%site(i)%left_dim
          rd2 = other%site(i)%right_dim

          if (this%site(i)%spin_kind/=other%site(i)%spin_kind) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO spin kinds do not match at site ',i
           endif
           stop
          endif

          if (this%site(i)%phys_in/=other%site(i)%phys_in) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO physical input dimensions do not match at site ',i
           endif
           stop
          endif

          if (this%site(i)%phys_out/=other%site(i)%phys_out) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO physical output dimensions do not match at site ',i
           endif
           stop
          endif

          tmp%site(i)%site_index = i
          tmp%site(i)%spin_kind = this%site(i)%spin_kind

          tmp%site(i)%phys_in = this%site(i)%phys_in
          tmp%site(i)%phys_out = this%site(i)%phys_out

          !---------------------------------------------------------------
          ! One-site MPO.
          !
          ! There is no virtual bond. The linear combination is simply
          !
          !     O = a O_1 + b O_2
          !
          ! represented by a single [1 x 1] virtual tensor.
          !---------------------------------------------------------------
          if(this%nsites==1)then
           allocate(tmp%site(i)%W(1,tmp%site(i)%phys_in,tmp%site(i)%phys_out,1))

           tmp%site(i)%W = a*this%site(i)%W+b*other%site(i)%W

           tmp%site(i)%left_dim  = 1
           tmp%site(i)%right_dim = 1

          !-------------------------------------------------------------
          ! First site
          !-------------------------------------------------------------
          elseif(i==1)then
           allocate(tmp%site(i)%W(1,tmp%site(i)%phys_in,tmp%site(i)%phys_out,rd1+rd2))
           tmp%site(i)%W=(0.0d0,0.0d0)

           do pin=1,tmp%site(i)%phys_in
            do pout=1,tmp%site(i)%phys_out
             tmp%site(i)%W(1,pin,pout,1:rd1) = a*this%site(i)%W(1,pin,pout,:)
             tmp%site(i)%W(1,pin,pout,rd1+1:rd1+rd2) = b*other%site(i)%W(1,pin,pout,:)
            enddo
           enddo

           tmp%site(i)%left_dim = 1
           tmp%site(i)%right_dim = rd1+rd2

          !-------------------------------------------------------------
          ! Last site
          !-------------------------------------------------------------
          elseif(i==this%nsites)then
           allocate(tmp%site(i)%W(ld1+ld2,tmp%site(i)%phys_in,tmp%site(i)%phys_out,1))
           tmp%site(i)%W=(0.0d0,0.0d0)

           do pin=1,tmp%site(i)%phys_in
            do pout=1,tmp%site(i)%phys_out
             tmp%site(i)%W(1:ld1,pin,pout,1) = this%site(i)%W(:,pin,pout,1)
             tmp%site(i)%W(ld1+1:ld1+ld2,pin,pout,1) = other%site(i)%W(:,pin,pout,1)
            enddo
           enddo

           tmp%site(i)%left_dim = ld1+ld2
           tmp%site(i)%right_dim = 1

          !-------------------------------------------------------------
          ! Interior site
          !-------------------------------------------------------------
          else
           allocate(tmp%site(i)%W(ld1+ld2,tmp%site(i)%phys_in,tmp%site(i)%phys_out,rd1+rd2))
           tmp%site(i)%W=(0.0d0,0.0d0)

           tmp%site(i)%W(1:ld1,:,:,1:rd1) = this%site(i)%W
           tmp%site(i)%W(ld1+1:ld1+ld2,:,:,rd1+1:rd1+rd2) = other%site(i)%W

           tmp%site(i)%left_dim = ld1+ld2
           tmp%site(i)%right_dim = rd1+rd2
          endif
         enddo

         !---------------------------------------------------------------
         ! Update the MPO metadata.
         !---------------------------------------------------------------
         tmp%bond_dimension = 1
         do i=1,tmp%nsites
          tmp%bond_dimension = max(tmp%bond_dimension,tmp%site(i)%left_dim,tmp%site(i)%right_dim)
         end do
         tmp%bond_dimension_max = max(this%bond_dimension_max,other%bond_dimension_max)
         if (this%is_hermitian.and.other%is_hermitian.and.aimag(a)==0.0d0.and.aimag(b)==0.0d0) then
          tmp%is_hermitian = .true.
         else
          tmp%is_hermitian = .false.
         endif
         tmp%trunc_tol = min(this%trunc_tol,other%trunc_tol)

         !---------------------------------------------------------------
         ! Update the set result MPO
         !---------------------------------------------------------------
         result = tmp

        return
        end subroutine lin_com_mpo

        subroutine compress_mps(this)
        implicit none
        class(mps_state)               :: this
        integer                        :: i,maxbond
        double precision               :: discarded

         !---------------------------------------------------------------
         ! Compress a Matrix Product State by truncating its virtual
         ! bond dimensions using singular value decompositions.
         !
         ! The truncation is performed sequentially from the first bond
         ! to the last. Singular values below the configured truncation
         ! threshold or beyond the maximum allowed bond dimension are
         ! discarded.
         !
         ! The total truncation error accumulated during the compression
         ! is stored in this%trunc_error.
         !---------------------------------------------------------------

         if (.not.allocated(this%site)) then
          if (mpi_id==0) then
           write(6,*) 'Err: Cannot compress unallocated MPS'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! First make the MPS right-canonical without truncation.
         !---------------------------------------------------------------
         call right_canonicalize_mps(this)

         !---------------------------------------------------------------
         ! Truncate each internal virtual bond from left to right.
         !
         ! Each SVD truncation updates the neighbouring site tensors and
         ! returns the contribution of the discarded singular values to
         ! the total truncation error.
         !---------------------------------------------------------------
         this%trunc_error = 0.0d0
         do i=1,this%nsites-1
          call svd_truncate_mps(this,i,discarded)
          this%trunc_error = this%trunc_error+discarded
         enddo

         !---------------------------------------------------------------
         ! Update the maximum virtual bond dimension after truncation.
         !
         ! Both left and right virtual dimensions are checked because
         ! truncation changes the bond dimensions of neighbouring tensors.
         !---------------------------------------------------------------
         maxbond=1
         do i=1,this%nsites
          maxbond = max(maxbond,this%site(i)%left_dim)
          maxbond = max(maxbond,this%site(i)%right_dim)
         enddo

         this%bond_dimension = maxbond
         this%canonical_site = this%nsites
         this%is_normalized = .false.

        return
        end subroutine compress_mps


        subroutine right_canonicalize_mps(this)
        implicit none
        class(mps_state)               :: this
        integer                        :: i
        integer                        :: ld,rd,phys
        integer                        :: m,n,k,info,lwork
        integer                        :: a,b,c,col
        complex(8), allocatable        :: theta(:,:)
        complex(8), allocatable        :: U(:,:)
        complex(8), allocatable        :: VT(:,:)
        complex(8), allocatable        :: work(:)
        complex(8), allocatable        :: newA(:,:,:)
        complex(8), allocatable        :: prevA(:,:,:)
        double precision, allocatable  :: S(:),rwork(:)

         !---------------------------------------------------------------
         ! Right-canonicalise sites from right to left.
         ! At each site:
         !   M(left,(phys,right)) = U*S*VT
         ! Store VT at this site and absorb U*S into the
         ! right bond of the previous site.
         !---------------------------------------------------------------

         do i=this%nsites,2,-1

          ld   = this%site(i)%left_dim
          rd   = this%site(i)%right_dim
          phys = this%site(i)%phys_dim

          if (this%site(i-1)%right_dim/=ld) then
           if (mpi_id==0) then
            write(6,*) 'Err: Invalid MPS bond in right canonicalisation'
            write(6,*) 'site ',i
           endif
           stop
          endif

          m = ld
          n = phys*rd
          k = min(m,n)

          allocate(theta(m,n))
          do a=1,ld
           do b=1,phys
            do c=1,rd
             col = (b-1)*rd+c
             theta(a,col) = this%site(i)%A(a,b,c)
            enddo
           enddo
          enddo

          allocate(U(m,k),VT(k,n),S(k),rwork(5*k))
          allocate(work(1))

          lwork=-1
          call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

          if (info/=0) then
           if (mpi_id==0) write(6,*) 'Err: MPS canonicalisation SVD workspace query failed',info
           stop
          endif

          lwork=max(1,int(real(work(1))))
          deallocate(work)
          allocate(work(lwork))

          call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

          if (info/=0) then
           if (mpi_id==0) write(6,*) 'Err: MPS right canonicalisation SVD failed',info
           stop
          endif

          !---------------------------------------------------------------
          ! The current tensor becomes VT, reshaped as (k,phys,rd).
          !---------------------------------------------------------------
          allocate(newA(k,phys,rd))
          do b=1,k
           do c=1,phys
            do a=1,rd
             col = (c-1)*rd+a
             newA(b,c,a) = VT(b,col)
            enddo
           enddo
          enddo

          !---------------------------------------------------------------
          ! Absorb U*S into the previous tensor.
          !---------------------------------------------------------------
          allocate(prevA(this%site(i-1)%left_dim,this%site(i-1)%phys_dim,k))
          prevA=(0.0d0,0.0d0)

          do a=1,this%site(i-1)%left_dim
           do b=1,this%site(i-1)%phys_dim
            do c=1,ld
             do col=1,k
              prevA(a,b,col) = prevA(a,b,col)+this%site(i-1)%A(a,b,c)*U(c,col)*S(col)
             enddo
            enddo
           enddo
          enddo

          deallocate(this%site(i)%A)
          call move_alloc(newA,this%site(i)%A)
          this%site(i)%left_dim = k

          deallocate(this%site(i-1)%A)
          call move_alloc(prevA,this%site(i-1)%A)
          this%site(i-1)%right_dim = k

          deallocate(theta,U,VT,S,rwork,work)
         enddo

        return
        end subroutine right_canonicalize_mps


        subroutine svd_truncate_mps(this,i,discarded)
        implicit none
        class(mps_state)               :: this
        integer                        :: i
        double precision               :: discarded
        integer                        :: ld1,rd1,ld2,rd2
        integer                        :: phys1,phys2
        integer                        :: m,n,k
        integer                        :: keep
        integer                        :: info
        integer                        :: lwork
        complex(8), allocatable        :: theta(:,:)
        complex(8), allocatable        :: U(:,:)
        complex(8), allocatable        :: VT(:,:)
        complex(8), allocatable        :: work(:)
        double precision, allocatable  :: S(:)
        double precision, allocatable  :: rwork(:)
        complex(8), allocatable        :: newA(:,:,:)
        complex(8), allocatable        :: newB(:,:,:)
        integer                        :: a,b,c,p,q,row,col
        double precision               :: norm

         !---------------------------------------------------------------
         ! Truncate the virtual bond between MPS sites i and i+1
         ! using an SVD.
         !
         ! The two neighbouring MPS tensors are first contracted
         ! over their shared virtual bond:
         !
         ! Theta(l,p,q,r) = sum_b A_i(l,p,b) A_{i+1}(b,q,r)
         !
         ! This is reshaped into a matrix
         !
         !    Theta_(l,p),(q,r)
         !
         ! and decomposed as
         !
         !    Theta = U S VT
         !
         ! The retained singular values are absorbed into the
         ! right tensor.
         !---------------------------------------------------------------

         discarded = 0.0d0

         !---------------------------------------------------------------
         ! Dimensions of the two tensors
         !---------------------------------------------------------------
         ld1   = this%site(i)%left_dim
         rd1   = this%site(i)%right_dim

         ld2   = this%site(i+1)%left_dim
         rd2   = this%site(i+1)%right_dim

         phys1 = this%site(i)%phys_dim
         phys2 = this%site(i+1)%phys_dim

         !---------------------------------------------------------------
         ! Check that the internal virtual bond matches.
         !---------------------------------------------------------------
         if (rd1/=ld2) then
          if (mpi_id==0) then
           write(6,*) 'Err: Invalid MPS bond in compression'
           write(6,*) 'site ',i
           write(6,*) 'right dimension = ',rd1
           write(6,*) 'next left dimension = ',ld2
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Matrix dimensions.
         !
         ! Rows    = left virtual index + physical index of site i
         ! Columns = physical index of site i+1 + right virtual index
         !---------------------------------------------------------------
         m = ld1*phys1
         n = phys2*rd2

         allocate(theta(m,n))
         theta = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Contract the two neighbouring MPS tensors.
         !---------------------------------------------------------------
         do a=1,ld1
          do p=1,phys1
           row = (a-1)*phys1+p
           do b=1,rd1
            do q=1,phys2
             col = (q-1)*rd2+1
             do c=1,rd2
              theta(row,col+c-1) = theta(row,col+c-1)+this%site(i)%A(a,p,b)*this%site(i+1)%A(b,q,c)
             enddo
            enddo
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Allocate reduced SVD arrays.
         !
         ! For ZGESVD with JOBU='S' and JOBVT='S':
         !
         ! U  = m x k
         ! S  = k
         ! VT = k x n
         !
         ! where k = min(m,n).
         !---------------------------------------------------------------
         k = min(m,n)

         allocate(U(m,k))
         allocate(VT(k,n))
         allocate(S(k))
         allocate(rwork(5*k))

         !---------------------------------------------------------------
         ! Workspace query for the reduced SVD.
         !---------------------------------------------------------------
         lwork = -1
         allocate(work(1))

         call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: ZGESVD workspace query failed'
           write(6,*) 'info = ',info
          endif
          stop
         endif

         lwork = max(1,int(real(work(1))))
         deallocate(work)
         allocate(work(lwork))

         !---------------------------------------------------------------
         ! Perform the reduced SVD.
         !---------------------------------------------------------------
         call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: SVD failed in MPS compression'
           write(6,*) 'info = ',info
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Determine the retained bond dimension.
         !
         ! Singular values are ordered from largest to smallest.
         !
         ! Retain the largest singular values such that the total
         ! discarded squared singular-value weight satisfies
         !
         !       sum(discarded S(a)^2) <= trunc_tol * norm
         !
         ! subject to the maximum permitted bond dimension.
         !
         ! The bond-dimension limit takes precedence if both conditions
         ! cannot be satisfied.
         !---------------------------------------------------------------
         keep = min(m,n)
         norm = sum(S*S)

         if(norm>0.0d0)then

          discarded = 0.0d0

          do while(keep>1)
           if(discarded+S(keep)*S(keep)<=this%trunc_tol*norm)then
            discarded = discarded+S(keep)*S(keep)
            keep = keep-1
           else
            exit
           endif
          enddo

         else
          keep = 1
         endif

         !---------------------------------------------------------------
         ! Enforce the maximum permitted bond dimension.
         !---------------------------------------------------------------
         keep = min(keep,this%bond_dimension_max)

         !---------------------------------------------------------------
         ! Calculate the actual discarded squared singular-value weight.
         !
         ! This includes any additional truncation imposed by the
         ! maximum bond dimension.
         !---------------------------------------------------------------
         discarded = 0.0d0

         do a=keep+1,min(m,n)
          discarded = discarded+S(a)*S(a)
         enddo

         !---------------------------------------------------------------
         ! Construct the new MPS tensors.
         !
         ! First tensor:
         !
         !    A'_i = U
         !
         ! Second tensor:
         !
         !    A'_{i+1} = S VT
         !---------------------------------------------------------------
         allocate(newA(ld1,phys1,keep))
         allocate(newB(keep,phys2,rd2))

         newA = (0.0d0,0.0d0)
         newB = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Reshape U into the first MPS tensor.
         !---------------------------------------------------------------
         do a=1,ld1
          do p=1,phys1
           row = (a-1)*phys1+p
           do b=1,keep
            newA(a,p,b) = U(row,b)
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Reshape S*VT into the second MPS tensor.
         !---------------------------------------------------------------
         do b=1,keep
          do q=1,phys2
           col = (q-1)*rd2+1
           do c=1,rd2
            newB(b,q,c) = S(b)*VT(b,col+c-1)
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Replace the original tensors.
         !---------------------------------------------------------------
         deallocate(this%site(i)%A)
         deallocate(this%site(i+1)%A)

         call move_alloc(newA,this%site(i)%A)
         call move_alloc(newB,this%site(i+1)%A)

         !---------------------------------------------------------------
         ! Update tensor dimensions.
         !---------------------------------------------------------------
         this%site(i)%right_dim = keep
         this%site(i+1)%left_dim = keep

         !---------------------------------------------------------------
         ! Release temporary arrays.
         !---------------------------------------------------------------
         deallocate(theta)
         deallocate(U)
         deallocate(VT)
         deallocate(S)
         deallocate(rwork)
         deallocate(work)

        return
        end subroutine svd_truncate_mps

        subroutine compress_mpo(this)
        implicit none
        class(mpo_operator)            :: this
        integer                        :: i,maxbond
        double precision               :: discarded

         !---------------------------------------------------------------
         ! Compress a Matrix Product Operator by truncating its virtual
         ! bond dimensions using singular value decompositions.
         !
         ! The truncation is performed sequentially from the first bond
         ! to the last. Singular values below the configured truncation
         ! threshold or beyond the maximum allowed bond dimension are
         ! discarded.
         !
         ! The total truncation error accumulated during the compression
         ! is stored in this%trunc_error.
         !---------------------------------------------------------------

         if (.not.allocated(this%site)) then
          if (mpi_id==0) then
           write(6,*) 'Err: Cannot compress unallocated MPO'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! First make the MPO right-canonical without truncating it.
         !---------------------------------------------------------------
         call right_canonicalize_mpo(this)

         !---------------------------------------------------------------
         ! Truncate each internal virtual bond from left to right.
         !
         ! Each SVD truncation updates the neighbouring site tensors and
         ! returns the contribution of the discarded singular values to
         ! the total truncation error.
         !---------------------------------------------------------------
         this%trunc_error = 0.0d0
         do i=1,this%nsites-1
          call svd_truncate_mpo(this,i,discarded)
          this%trunc_error = this%trunc_error+discarded
         enddo

         !---------------------------------------------------------------
         ! Update the maximum virtual bond dimension after truncation.
         !
         ! Both left and right virtual dimensions are checked because
         ! truncation changes the bond dimensions of neighbouring tensors.
         !---------------------------------------------------------------
         maxbond=1
         do i=1,this%nsites
          maxbond = max(maxbond,this%site(i)%left_dim)
          maxbond = max(maxbond,this%site(i)%right_dim)
         enddo

         this%bond_dimension = maxbond

        return
        end subroutine compress_mpo


        subroutine right_canonicalize_mpo(this)
        implicit none
        class(mpo_operator)            :: this
        integer                        :: i
        integer                        :: ld,rd,pin,pout
        integer                        :: m,n,k,info,lwork
        integer                        :: a,b,c,pi,po,col
        complex(8), allocatable        :: theta(:,:)
        complex(8), allocatable        :: U(:,:)
        complex(8), allocatable        :: VT(:,:)
        complex(8), allocatable        :: work(:)
        complex(8), allocatable        :: newW(:,:,:,:)
        complex(8), allocatable        :: prevW(:,:,:,:)
        double precision, allocatable  :: S(:),rwork(:)

         !--------------------------------------------------------------
         ! Right-canonicalise sites from right to left.
         ! At each site:
         !   M(left,(phys_in,phys_out,right)) = U*S*VT
         ! Store VT at this site and absorb U*S into the
         ! right bond of the previous site.
         !--------------------------------------------------------------

         do i=this%nsites,2,-1
          ld   = this%site(i)%left_dim
          rd   = this%site(i)%right_dim
          pin  = this%site(i)%phys_in
          pout = this%site(i)%phys_out

          if (this%site(i-1)%right_dim/=ld) then
           if (mpi_id==0) then
            write(6,*) 'Err: Invalid MPO bond in right canonicalisation'
            write(6,*) 'site ',i
           endif
           stop
          endif

          m = ld
          n = pin*pout*rd
          k = min(m,n)

          allocate(theta(m,n))
          do a=1,ld
           do pi=1,pin
            do po=1,pout
             do c=1,rd
              col = ((pi-1)*pout+(po-1))*rd+c
              theta(a,col) = this%site(i)%W(a,pi,po,c)
             enddo
            enddo
           enddo
          enddo

          allocate(U(m,k),VT(k,n),S(k),rwork(5*k))
          allocate(work(1))

          lwork=-1
          call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

          if (info/=0) then
           if (mpi_id==0) write(6,*) 'Err: MPO canonicalisation SVD workspace query failed',info
           stop
          endif

          lwork=max(1,int(real(work(1))))
          deallocate(work)
          allocate(work(lwork))

          call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

          if (info/=0) then
           if (mpi_id==0) write(6,*) 'Err: MPO right canonicalisation SVD failed',info
           stop
          endif

          !--------------------------------------------------------------
          ! The current tensor becomes VT, reshaped as (k,pin,pout,rd).
          !--------------------------------------------------------------
          allocate(newW(k,pin,pout,rd))
          do b=1,k
           do pi=1,pin
            do po=1,pout
             do c=1,rd
              col = ((pi-1)*pout+(po-1))*rd+c
              newW(b,pi,po,c) = VT(b,col)
             enddo
            enddo
           enddo
          enddo

          !--------------------------------------------------------------
          ! Absorb U*S into the previous tensor.
          !--------------------------------------------------------------
          allocate(prevW(this%site(i-1)%left_dim,this%site(i-1)%phys_in,this%site(i-1)%phys_out,k))
          prevW=(0.0d0,0.0d0)

          do a=1,this%site(i-1)%left_dim
           do pi=1,this%site(i-1)%phys_in
            do po=1,this%site(i-1)%phys_out
             do b=1,ld
              do c=1,k
               prevW(a,pi,po,c) = prevW(a,pi,po,c)+this%site(i-1)%W(a,pi,po,b)*U(b,c)*S(c)
              enddo
             enddo
            enddo
           enddo
          enddo

          deallocate(this%site(i)%W)
          call move_alloc(newW,this%site(i)%W)
          this%site(i)%left_dim = k

          deallocate(this%site(i-1)%W)
          call move_alloc(prevW,this%site(i-1)%W)
          this%site(i-1)%right_dim = k

          deallocate(theta,U,VT,S,rwork,work)
         enddo

        return
        end subroutine right_canonicalize_mpo

        subroutine svd_truncate_mpo(this,i,discarded)
        implicit none
        class(mpo_operator)            :: this
        integer                        :: i
        double precision               :: discarded
        integer                        :: ld1,rd1,ld2,rd2
        integer                        :: pin1,pout1,pin2,pout2
        integer                        :: m,n,k
        integer                        :: keep
        integer                        :: info
        integer                        :: lwork
        complex(8), allocatable        :: theta(:,:)
        complex(8), allocatable        :: U(:,:)
        complex(8), allocatable        :: VT(:,:)
        complex(8), allocatable        :: work(:)
        double precision, allocatable  :: S(:)
        double precision, allocatable  :: rwork(:)
        complex(8), allocatable        :: newW1(:,:,:,:)
        complex(8), allocatable        :: newW2(:,:,:,:)
        integer                        :: a,b,c
        integer                        :: pi,po,qi,qo
        integer                        :: row,col
        double precision               :: norm

         !---------------------------------------------------------------
         ! Truncate a single virtual bond of an MPO using an SVD.
         !
         ! The two neighbouring MPO tensors are contracted over
         ! their shared virtual bond:
         !
         !   Theta(l,pi,po,qi,qo,r) = sum_b W_i(l,pi,po,b) W_{i+1}(b,qi,qo,r)
         !
         ! This is reshaped into the matrix
         !
         !   row    = (l,pi,po)
         !   column = (qi,qo,r)
         !
         ! and decomposed as
         !
         !   Theta = U S VT
         !---------------------------------------------------------------
         discarded = 0.0d0

         !---------------------------------------------------------------
         ! Dimensions of the neighbouring tensors.
         !---------------------------------------------------------------
         ld1 = this%site(i)%left_dim
         rd1 = this%site(i)%right_dim

         ld2 = this%site(i+1)%left_dim
         rd2 = this%site(i+1)%right_dim

         pin1  = this%site(i)%phys_in
         pout1 = this%site(i)%phys_out

         pin2  = this%site(i+1)%phys_in
         pout2 = this%site(i+1)%phys_out

         !---------------------------------------------------------------
         ! Check the shared virtual bond.
         !---------------------------------------------------------------
         if (rd1/=ld2) then
          if (mpi_id==0) then
           write(6,*) 'Err: Invalid MPO bond in compression'
           write(6,*) 'site ',i
           write(6,*) 'right dimension = ',rd1
           write(6,*) 'next left dimension = ',ld2
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Dimensions of the SVD matrix.
         !
         ! Row space:
         !
         !    left virtual x input physical x output physical
         !
         ! Column space:
         !
         !    input physical x output physical x right virtual
         !---------------------------------------------------------------
         m = ld1*pin1*pout1
         n = pin2*pout2*rd2

         allocate(theta(m,n))
         theta = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Contract the two MPO tensors.
         !---------------------------------------------------------------
         do a=1,ld1
          do pi=1,pin1
           do po=1,pout1
            row = ((a-1)*pin1+(pi-1))*pout1+po
            do b=1,rd1
             do qi=1,pin2
              do qo=1,pout2
               col = ((qi-1)*pout2+(qo-1))*rd2+1
               do c=1,rd2
                theta(row,col+c-1) = theta(row,col+c-1)+this%site(i)%W(a,pi,po,b)*this%site(i+1)%W(b,qi,qo,c)
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Allocate reduced SVD arrays.
         !
         ! For ZGESVD with JOBU='S' and JOBVT='S':
         !
         ! U  = m x k
         ! S  = k
         ! VT = k x n
         !
         ! where k = min(m,n).
         !---------------------------------------------------------------
         k = min(m,n)

         allocate(U(m,k))
         allocate(VT(k,n))
         allocate(S(k))
         allocate(rwork(5*k))

         !---------------------------------------------------------------
         ! Workspace query for the reduced SVD.
         !---------------------------------------------------------------
         lwork = -1
         allocate(work(1))

         call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: ZGESVD workspace query failed'
           write(6,*) 'info = ',info
          endif
          stop
         endif

         lwork = max(1,int(real(work(1))))
         deallocate(work)
         allocate(work(lwork))

         !---------------------------------------------------------------
         ! Perform the reduced SVD.
         !---------------------------------------------------------------
         call zgesvd('S','S',m,n,theta,m,S,U,m,VT,k,work,lwork,rwork,info)

         if(info/=0)then
          if(mpi_id==0)then
           write(6,*) 'Err: SVD failed in MPO compression'
           write(6,*) 'info = ',info
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Determine the retained bond dimension.
         !
         ! Singular values are ordered from largest to smallest.
         !
         ! Retain the largest singular values such that the total
         ! discarded squared singular-value weight satisfies
         !
         !       sum(discarded S(a)^2) <= trunc_tol * norm
         !
         ! subject to the maximum permitted bond dimension.
         !
         ! The bond-dimension limit takes precedence if both conditions
         ! cannot be satisfied.
         !---------------------------------------------------------------
         keep = min(m,n)
         norm = sum(S*S)

         if(norm>0.0d0)then

          discarded = 0.0d0

          do while(keep>1)
           if(discarded+S(keep)*S(keep)<=this%trunc_tol*norm)then
            discarded = discarded+S(keep)*S(keep)
            keep = keep-1
           else
            exit
           endif
          enddo

         else
          keep = 1
         endif

         !---------------------------------------------------------------
         ! Enforce the maximum permitted bond dimension.
         !---------------------------------------------------------------
         keep = min(keep,this%bond_dimension_max)

         !---------------------------------------------------------------
         ! Calculate the actual discarded squared singular-value weight.
         !
         ! This includes any additional truncation imposed by the
         ! maximum bond dimension.
         !---------------------------------------------------------------
         discarded = 0.0d0

         do a=keep+1,min(m,n)
          discarded = discarded+S(a)*S(a)
         enddo

         !---------------------------------------------------------------
         ! Reconstruct the two MPO tensors.
         !
         ! First tensor:
         !
         !   W'_i = U
         !
         ! Second tensor:
         !
         !   W'_{i+1} = S VT
         !---------------------------------------------------------------
         allocate(newW1(ld1,pin1,pout1,keep))
         allocate(newW2(keep,pin2,pout2,rd2))

         newW1 = (0.0d0,0.0d0)
         newW2 = (0.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Reshape U into the first MPO tensor.
         !---------------------------------------------------------------
         do a=1,ld1
          do pi=1,pin1
           do po=1,pout1
            row = ((a-1)*pin1+(pi-1))*pout1+po
            do b=1,keep
             newW1(a,pi,po,b) = U(row,b)
            enddo
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Reshape S*VT into the second MPO tensor.
         !---------------------------------------------------------------
         do b=1,keep
          do qi=1,pin2
           do qo=1,pout2
            col = ((qi-1)*pout2+(qo-1))*rd2+1
            do c=1,rd2
             newW2(b,qi,qo,c) = S(b)*VT(b,col+c-1)
            enddo
           enddo
          enddo
         enddo

         !---------------------------------------------------------------
         ! Replace the original MPO tensors.
         !---------------------------------------------------------------
         deallocate(this%site(i)%W)
         deallocate(this%site(i+1)%W)

         call move_alloc(newW1,this%site(i)%W)
         call move_alloc(newW2,this%site(i+1)%W)

         !---------------------------------------------------------------
         ! Update the shared virtual bond dimension.
         !---------------------------------------------------------------
         this%site(i)%right_dim  = keep
         this%site(i+1)%left_dim = keep

         !---------------------------------------------------------------
         ! Release temporary arrays.
         !---------------------------------------------------------------
         deallocate(theta)
         deallocate(U)
         deallocate(VT)
         deallocate(S)
         deallocate(rwork)
         deallocate(work)

        return
        end subroutine svd_truncate_mpo

        subroutine apply_mpo(this,mpo,result)
        implicit none
        class(mps_state)        :: this
        class(mpo_operator)     :: mpo
        class(mps_state)        :: result
        integer                 :: i,a,b,c,d
        integer                 :: al,ar,wl,wr
        integer                 :: s_in,s_out,new_left,new_right,p_in,p_out

         !---------------------------------------------------------------
         ! Apply a Matrix Product Operator to a Matrix Product State.
         !
         ! The resulting MPS represents the state
         !
         !     |ψ'⟩ = W|ψ⟩.
         !
         ! The MPO and MPS virtual bonds are combined by a tensor-product
         ! construction, giving an exact representation of the result
         ! without truncation.
         !---------------------------------------------------------------

         if (this%nsites/=mpo%nsites) then
          if (mpi_id==0) then
           write(6,*) 'Err: MPS and MPO sizes do not match'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Allocate the output MPS.
         !
         ! The input state is not modified. Any existing tensors in the
         ! result are replaced by the newly constructed MPS.
         !---------------------------------------------------------------
         if(allocated(result%site)) deallocate(result%site)

         result%nsites = this%nsites

         allocate(result%site(result%nsites))

         result%basis = this%basis
         result%bond_dimension_max = this%bond_dimension_max

         !---------------------------------------------------------------
         ! Construct each output MPS tensor by contracting the physical
         ! input index of the MPO with the physical index of the MPS.
         !
         ! The output virtual indices combine the MPS and MPO virtual
         ! indices:
         !
         !     new_left  = mps_left  * mpo_left
         !     new_right = mps_right * mpo_right
         !
         ! The resulting physical index is the MPO output index.
         !---------------------------------------------------------------
         do i=1,this%nsites
          al=this%site(i)%left_dim
          ar=this%site(i)%right_dim

          wl=mpo%site(i)%left_dim
          wr=mpo%site(i)%right_dim

          p_in = mpo%site(i)%phys_in
          p_out = mpo%site(i)%phys_out

          if (p_in/=this%site(i)%phys_dim) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO input physical dimension mismatch on site ',i
           endif
           stop
          endif

          if (p_out/=this%site(i)%phys_dim) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO output physical dimension mismatch on site ',i
           endif
           stop
          endif

          new_left = al*wl
          new_right = ar*wr

          allocate(result%site(i)%A(new_left,p_out,new_right))
          result%site(i)%A = (0.0d0,0.0d0)

          !---------------------------------------------------------------
          ! Contract the MPO input index with the MPS physical index.
          !
          ! The combined virtual indices are ordered with the MPS index
          ! as the outer index and the MPO index as the inner index.
          !---------------------------------------------------------------
          do a=1,al
           do b=1,ar
            do c=1,wl
             do d=1,wr
              do s_in=1,p_in
               do s_out=1,p_out
                result%site(i)%A((a-1)*wl+c,s_out,(b-1)*wr+d) = result%site(i)%A((a-1)*wl+c,s_out,(b-1)*wr+d)+ &
                                                                mpo%site(i)%W(c,s_in,s_out,d)* &
                                                                this%site(i)%A(a,s_in,b)
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo

          !---------------------------------------------------------------
          ! Store the contracted tensor and update its dimensions.
          !---------------------------------------------------------------
          result%site(i)%left_dim = new_left
          result%site(i)%phys_dim = p_out
          result%site(i)%right_dim = new_right
          result%site(i)%site_index = i
          result%site(i)%spin_kind = this%site(i)%spin_kind
         enddo

         !---------------------------------------------------------------
         ! Update the MPS metadata.
         !
         ! Applying an MPO does not perform any truncation, so the
         ! resulting bond dimension is the exact bond dimension generated
         ! by the contraction.
         !---------------------------------------------------------------
         result%bond_dimension = 1
         do i=1,result%nsites
          result%bond_dimension = max(result%bond_dimension,result%site(i)%left_dim,result%site(i)%right_dim)
         end do
         !---------------------------------------------------------------
         ! The tensor contraction does not preserve a canonical form or
         ! normalization in general. No truncation has been performed.
         !---------------------------------------------------------------
         result%canonical_site = 0
         result%is_normalized = .false.
         result%trunc_error = 0.0d0

        return
        end subroutine apply_mpo

        subroutine compose_mpos(this,other,result)
        implicit none
        class(mpo_operator)     :: this,other,result
        integer                 :: i,a,b,c,d
        integer                 :: al,ar,bl,br
        integer                 :: new_left,new_right,s_in,s_mid,s_out
        complex(8), allocatable :: tmp(:,:,:,:)

         !---------------------------------------------------------------
         ! Compose two Matrix Product Operators by matrix multiplication.
         !
         ! The resulting operator is
         !
         !     W = W_this W_other,
         !
         ! so that other acts first and this acts second when the
         ! composed operator is applied to a state.
         !
         ! The composition is represented exactly by contracting the
         ! physical output index of other with the physical input index
         ! of this.
         !---------------------------------------------------------------

         if (this%nsites/=other%nsites) then
          if (mpi_id==0) then
           write(6,*) 'Err: MPO sizes do not match'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Allocate the output MPO.
         !
         ! The input MPOs are not modified. Any existing tensors in
         ! result are replaced by the newly constructed operator.
         !---------------------------------------------------------------
         if(allocated(result%site)) deallocate(result%site)

         result%nsites = this%nsites

         allocate(result%site(result%nsites))

         result%basis = this%basis
         result%bond_dimension_max = max(this%bond_dimension_max,other%bond_dimension_max)

         !---------------------------------------------------------------
         ! Construct each output MPO tensor by contracting the physical
         ! output index of other with the physical input index of this.
         !
         ! The output virtual indices combine the virtual indices of the
         ! two MPOs:
         !
         !     new_left  = this_left  * other_left
         !     new_right = this_right * other_right
         !---------------------------------------------------------------
         do i=1,this%nsites
          al = this%site(i)%left_dim
          ar = this%site(i)%right_dim

          bl = other%site(i)%left_dim
          br = other%site(i)%right_dim

          if (this%site(i)%phys_in/=other%site(i)%phys_out) then
           if (mpi_id==0) then
            write(6,*) 'Err: MPO composition physical dimension mismatch on site ',i
           endif
           stop
          endif

          new_left = al*bl
          new_right = ar*br

          allocate(tmp(new_left,other%site(i)%phys_in,this%site(i)%phys_out,new_right))
          tmp = (0.0d0,0.0d0)

          !---------------------------------------------------------------
          ! Contract the intermediate physical index.
          !
          ! The combined virtual indices are ordered with the virtual
          ! index of this as the outer index and the corresponding
          ! virtual index of other as the inner index.
          !---------------------------------------------------------------
          do a=1,al
           do b=1,ar
            do c=1,bl
             do d=1,br
              do s_in=1,other%site(i)%phys_in
               do s_mid=1,this%site(i)%phys_in
                do s_out=1,this%site(i)%phys_out
                 tmp((a-1)*bl+c,s_in,s_out,(b-1)*br+d) = tmp((a-1)*bl+c,s_in,s_out,(b-1)*br+d) + &
                                                         this%site(i)%W(a,s_mid,s_out,b)*other%site(i)%W(c,s_in,s_mid,d)
                enddo
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo

          !---------------------------------------------------------------
          ! Store the composed tensor and update its dimensions.
          !---------------------------------------------------------------
          allocate(result%site(i)%W(new_left,this%site(i)%phys_in,this%site(i)%phys_out,new_right))

          result%site(i)%W = tmp

          result%site(i)%left_dim  = new_left
          result%site(i)%right_dim = new_right
          result%site(i)%phys_in = other%site(i)%phys_in
          result%site(i)%phys_out = this%site(i)%phys_out

          result%site(i)%site_index = i
          result%site(i)%spin_kind = this%site(i)%spin_kind

          deallocate(tmp)
         enddo

         !---------------------------------------------------------------
         ! Update the MPO metadata.
         !
         ! Composition is performed without truncation, so the resulting
         ! bond dimension is the exact bond dimension generated by the
         ! tensor contraction.
         !---------------------------------------------------------------
         result%bond_dimension = 1
         do i=1,result%nsites
          result%bond_dimension = max(result%bond_dimension,result%site(i)%left_dim,result%site(i)%right_dim)
         enddo
         result%is_hermitian = .false.
         result%trunc_error = 0.0d0

        return
        end subroutine compose_mpos

        subroutine kronecker_product(this,other,result)
        implicit none
        class(mpo_operator)     :: this,other,result
        integer                 :: i,n
        integer                 :: ld,rd,pin,pout

         !---------------------------------------------------------------
         ! Construct the Kronecker product of two Matrix Product
         ! Operators.
         !
         ! The resulting MPO represents the tensor product
         !
         !     W = this ⊗ other.
         !
         ! The two operators remain on separate physical sites and their
         ! virtual bonds are not contracted.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Allocate the output MPO.
         !
         ! The input MPOs are not modified. Any existing tensors in
         ! result are replaced by copies of the two input MPOs.
         !---------------------------------------------------------------
         if(allocated(result%site)) deallocate(result%site)

         result%nsites = this%nsites+other%nsites
         allocate(result%site(result%nsites))

         !---------------------------------------------------------------
         ! Construct the combined local spin basis.
         !
         ! The Kronecker product combines the local bases of the two
         ! independent MPOs. The resulting spin kinds are assigned
         ! sequential labels starting from 1.
         !
         ! The local basis from the first MPO is copied first, followed
         ! by the local bases from the second MPO.
         !---------------------------------------------------------------
         result%basis%nkinds = this%basis%nkinds+other%basis%nkinds

         if(allocated(result%basis%kind)) deallocate(result%basis%kind)
         allocate(result%basis%kind(result%basis%nkinds))

         do i=1,this%basis%nkinds
          result%basis%kind(i) = this%basis%kind(i)
         enddo

         do i=1,other%basis%nkinds
          result%basis%kind(this%basis%nkinds+i) = other%basis%kind(i)
         enddo

         result%bond_dimension_max = max(this%bond_dimension_max,other%bond_dimension_max)
         result%bond_dimension = max(this%bond_dimension,other%bond_dimension)

         !---------------------------------------------------------------
         ! Copy the MPO metadata.
         !
         ! The two MPOs are disconnected, so the maximum and actual bond
         ! dimensions are determined independently by the two operators.
         !---------------------------------------------------------------
         result%bond_dimension_max = max(this%bond_dimension_max,other%bond_dimension_max)
         result%bond_dimension = max(this%bond_dimension,other%bond_dimension)

         !---------------------------------------------------------------
         ! Copy the tensors belonging to the 'this' MPO.
         !
         ! These sites retain their original virtual and physical
         ! dimensions and occupy the first part of the output chain.
         !---------------------------------------------------------------
         do i=1,this%nsites
          ld = this%site(i)%left_dim
          rd = this%site(i)%right_dim
          pin = this%site(i)%phys_in
          pout = this%site(i)%phys_out

          allocate(result%site(i)%W(ld,pin,pout,rd))

          result%site(i)%W = this%site(i)%W

          result%site(i)%left_dim = ld
          result%site(i)%right_dim = rd
          result%site(i)%phys_in = pin
          result%site(i)%phys_out = pout
          result%site(i)%site_index = i
          result%site(i)%spin_kind = this%site(i)%spin_kind
         enddo

         !---------------------------------------------------------------
         ! Copy the tensors belonging to the 'other' MPO.
         !
         ! These sites are appended after the sites of this MPO and their
         ! site indices are shifted accordingly.
         !---------------------------------------------------------------
         do i=1,other%nsites
          n = this%nsites+i

          ld = other%site(i)%left_dim
          rd = other%site(i)%right_dim
          pin = other%site(i)%phys_in
          pout = other%site(i)%phys_out

          allocate(result%site(n)%W(ld,pin,pout,rd))

          result%site(n)%W = other%site(i)%W

          result%site(n)%left_dim = ld
          result%site(n)%right_dim = rd
          result%site(n)%phys_in = pin
          result%site(n)%phys_out = pout
          result%site(n)%site_index = n
          result%site(n)%spin_kind = this%basis%nkinds+other%site(i)%spin_kind
         enddo

         !---------------------------------------------------------------
         ! No truncation is performed during the Kronecker product.
         ! Hermiticity is preserved when both input operators are
         ! Hermitian.
         !---------------------------------------------------------------
         result%trunc_error = 0.0d0
         result%is_hermitian = this%is_hermitian.and.other%is_hermitian

        return
        end subroutine kronecker_product

        function expectation_value(this,mpo) result(value)
        implicit none
        class(mps_state)               :: this
        class(mpo_operator)            :: mpo
        complex(8)                     :: value
        complex(8), allocatable        :: env(:,:,:)
        complex(8), allocatable        :: new_env(:,:,:)
        integer                        :: i
        integer                        :: lb,rb
        integer                        :: lw,rw
        integer                        :: lk,rk
        integer                        :: sin,sout
        integer                        :: ld,rd
        integer                        :: wld,wrd
        integer                        :: s

         !---------------------------------------------------------------
         ! Calculate the expectation value of an MPO operator.
         !
         !                 value = <psi|O|psi>
         !
         ! The MPS bra, MPO and MPS ket are contracted simultaneously.
         ! No intermediate MPS is constructed.
         !
         ! The environment indices are:
         !
         !     (bra virtual, MPO virtual, ket virtual)
         !---------------------------------------------------------------
         value = (0.0d0,0.0d0)

         if(this%nsites/=mpo%nsites)then
          if(mpi_id==0)then
           write(6,*) 'Err: MPS and MPO sizes do not match in expectation_value'
          endif
          stop
         endif

         if(.not.allocated(this%site).or..not.allocated(mpo%site))then
          if(mpi_id==0)then
           write(6,*) 'Err: Unallocated MPS or MPO in expectation_value'
          endif
          stop
         endif

         allocate(env(1,1,1))
         env = (1.0d0,0.0d0)

         !---------------------------------------------------------------
         ! Contract the MPS bra, MPO and MPS ket from left to right.
         !---------------------------------------------------------------
         do i=1,this%nsites
          ld = this%site(i)%left_dim
          rd = this%site(i)%right_dim

          wld = mpo%site(i)%left_dim
          wrd = mpo%site(i)%right_dim

          if(mpo%site(i)%phys_in/=this%site(i)%phys_dim.or. &
             mpo%site(i)%phys_out/=this%site(i)%phys_dim)then
           if(mpi_id==0)then
            write(6,*) 'Err: Physical dimension mismatch in expectation_value'
            write(6,*) 'Site ',i
           endif
           stop
          endif

          if(size(env,1)/=ld.or.size(env,2)/=wld.or.size(env,3)/=ld)then
           if(mpi_id==0)then
            write(6,*) 'Err: Environment dimension mismatch in expectation_value'
            write(6,*) 'Site ',i
           endif
           stop
          endif

          allocate(new_env(rd,wrd,rd))
          new_env = (0.0d0,0.0d0)

          do rb=1,rd
           do rw=1,wrd
            do rk=1,rd
             do lb=1,ld
              do lw=1,wld
               do lk=1,ld
                do sin=1,this%site(i)%phys_dim
                 do sout=1,this%site(i)%phys_dim
                  new_env(rb,rw,rk) = new_env(rb,rw,rk)+conjg(this%site(i)%A(lb,sout,rb))*&
                                      mpo%site(i)%W(lw,sin,sout,rw)*this%site(i)%A(lk,sin,rk)*&
                                      env(lb,lw,lk)
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo

          deallocate(env)
          call move_alloc(new_env,env)

         enddo

         !---------------------------------------------------------------
         ! The final environment contains the scalar expectation value.
         !---------------------------------------------------------------
         value = env(1,1,1)

         deallocate(env)

        return
        end function expectation_value

        subroutine init_spin_H(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i,kind,d,t1,t2,rate

         !---------------------------------------------------------------
         ! Construct a Matrix Product Operator representing the spin
         ! Hamiltonian of the supplied spin system.
         !
         ! The local spin basis is constructed for every spin kind and
         ! the Hamiltonian is assembled from its one-site and two-site
         ! physical interaction terms.
         !
         ! The resulting MPO is compressed after all Hamiltonian terms
         ! have been added.
         !---------------------------------------------------------------

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Spin Hamiltonian MPO'
          flush(6)
         endif

         !---------------------------------------------------------------
         ! Allocate one tensor for each physical spin site and initialise
         ! the MPO metadata.
         !
         ! The Hamiltonian is constructed before compression, so the
         ! initial bond dimension is set to one and subsequently enlarged
         ! as individual Hamiltonian terms are added.
         !---------------------------------------------------------------
         this%nsites = spindy%nspins

         !---------------------------------------------------------------
         ! Construct the local spin basis for every spin kind.
         !
         ! The resulting basis provides the physical operator matrices
         ! used by the Hamiltonian-term construction routines.
         !---------------------------------------------------------------
         if(allocated(this%basis%kind)) deallocate(this%basis%kind)

         this%basis%nkinds = spindy%nkinds

         allocate(this%basis%kind(this%basis%nkinds))

         do i=1,spindy%nkinds
          call init_spin_basis(this%basis%kind(i),nint(2.0d0*spindy%spin(i)))
         enddo

         !---------------------------------------------------------------
         ! Allocate one MPO tensor per physical site.
         !---------------------------------------------------------------
         if(allocated(this%site)) deallocate(this%site)
         allocate(this%site(this%nsites))

         !---------------------------------------------------------------
         ! Initialise the Hamiltonian accumulator to the ZERO operator.
         !
         ! Every site has a bond-dimension-1 tensor containing zero.
         !
         ! This represents
         !
         !        H = 0
         !
         ! before the individual Hamiltonian terms are added.
         !---------------------------------------------------------------
         this%bond_dimension = 1
         this%trunc_error = 0.0d0
         this%is_hermitian = .true.

         do i=1,this%nsites
          kind = spindy%kind(i)
          d = this%basis%kind(kind)%d

          this%site(i)%site_index = i
          this%site(i)%spin_kind = kind

          this%site(i)%left_dim = 1
          this%site(i)%right_dim = 1
          this%site(i)%phys_in = d
          this%site(i)%phys_out = d

          allocate(this%site(i)%W(1,d,d,1))
          this%site(i)%W = (0.0d0,0.0d0)
         enddo

         !---------------------------------------------------------------
         ! Add the physical Hamiltonian terms to the MPO.
         !
         ! Each construction routine adds its corresponding one-site or
         ! two-site operator contributions to the MPO representation.
         ! The individual terms are assembled before the final MPO
         ! compression.
         !---------------------------------------------------------------
         call add_zeeman_terms(this,spindy)

         call add_stevens_terms(this,spindy)

         call add_singleion_terms(this,spindy)

         call add_exchange_terms(this,spindy)

         call add_exchange_anisotropy_terms(this,spindy)

         call add_dipolar_terms(this,spindy)

         !---------------------------------------------------------------
         ! Compress the assembled MPO to the requested truncation
         ! tolerance and maximum bond dimension.
         !---------------------------------------------------------------
         call this%compress()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine init_spin_H

        subroutine add_zeeman_terms(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i,j,kind,g_index,d
        complex(8)              :: cx,cy,cz
        complex(8), allocatable :: op(:,:)

         !---------------------------------------------------------------
         ! Add the Zeeman interaction terms to the Hamiltonian MPO.
         !
         ! The Zeeman interaction is constructed for every spin whose
         ! spin kind has an associated g tensor. The resulting one-site
         ! operators are added directly to the Hamiltonian MPO.
         !---------------------------------------------------------------

         if(spindy%SH%nG.eq.0) return

         do i=1,spindy%nspins
          kind = spindy%kind(i)

          g_index = 0

         !---------------------------------------------------------------
         ! Construct the Zeeman operator for each physical spin site.
         !
         ! Each site is matched to the corresponding g tensor through its
         ! spin kind. Sites without an associated g tensor are skipped.
         !---------------------------------------------------------------
          do j=1,spindy%SH%nG
           if(spindy%SH%G(j)%kind.eq.kind)then
            g_index = j
            exit
           endif
          enddo

          if(g_index.eq.0) cycle

          !------------------------------------------------------------
          ! Form the local Zeeman operator
          !
          !     H = -μB B·g·S
          !
          ! expressed in the local spin basis.
          !------------------------------------------------------------
          cx = spindy%SH%G(g_index)%G(1,1)*spindy%Bfield(1) + &
               spindy%SH%G(g_index)%G(2,1)*spindy%Bfield(2) + &
               spindy%SH%G(g_index)%G(3,1)*spindy%Bfield(3)

          cy = spindy%SH%G(g_index)%G(1,2)*spindy%Bfield(1) + &
               spindy%SH%G(g_index)%G(2,2)*spindy%Bfield(2) + &
               spindy%SH%G(g_index)%G(3,2)*spindy%Bfield(3)

          cz = spindy%SH%G(g_index)%G(1,3)*spindy%Bfield(1) + &
               spindy%SH%G(g_index)%G(2,3)*spindy%Bfield(2) + &
               spindy%SH%G(g_index)%G(3,3)*spindy%Bfield(3)

          cx = -cx*spindy%bohr_mag(kind)
          cy = -cy*spindy%bohr_mag(kind)
          cz = -cz*spindy%bohr_mag(kind)

          d = this%basis%kind(kind)%d

          allocate(op(d,d))
          op = (0.0d0,0.0d0)

          op = cx*this%basis%kind(kind)%Sx + &
               cy*this%basis%kind(kind)%Sy + &
               cz*this%basis%kind(kind)%Sz

          !------------------------------------------------------------
          ! Add the local Zeeman operator to the Hamiltonian MPO.
          !------------------------------------------------------------
          call add_one_site_operator(this,i,op)

          deallocate(op)
         enddo

        return
        end subroutine add_zeeman_terms

        subroutine add_stevens_terms(this,spindy)
        use anglib
        use stevens_class
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i,j,k,q
        integer                 :: d
        integer                 :: T1,M1,T,M,T2,M2
        integer                 :: s1,s2
        integer                 :: kind
        double precision        :: spin
        double precision        :: q1,q2
        double precision        :: a
        complex(8)              :: mat_elem
        complex(8)              :: val_tmp
        complex(8), allocatable :: Omat(:,:)

         !---------------------------------------------------------------
         ! Add the Stevens crystal-field terms to the Hamiltonian MPO.
         !
         ! The crystal-field Hamiltonian is constructed from the supplied
         ! Stevens parameters for each spin kind. The corresponding
         ! Stevens operator matrices are evaluated in the local spin basis
         ! before being added as one-site operators to the Hamiltonian MPO.
         !---------------------------------------------------------------

         if(spindy%SH%nO.eq.0) return

         !---------------------------------------------------------------
         ! Construct the crystal-field operator for each physical spin.
         !
         ! All Stevens terms defined for the corresponding spin kind are
         ! accumulated into a single local operator matrix.
         !---------------------------------------------------------------
         do i=1,spindy%nspins
          kind = spindy%kind(i)
          spin = spindy%spin(kind)
          d = this%basis%kind(kind)%d

          T1 = nint(2.0d0*spin)
          T2 = T1

          allocate(Omat(d,d))
          Omat = (0.0d0,0.0d0)

          !------------------------------------------------------------
          ! Construct the local Stevens operator matrix.
          !
          ! The contribution from every Stevens parameter defined for
          ! this spin kind is accumulated in the local spin basis.
          !------------------------------------------------------------
          do j=1,spindy%SH%nO
           if(spindy%SH%O(j)%kind.ne.kind) cycle

           k = spindy%SH%O(j)%k
           T = 2*k

           do q=1,2*k+1
            q2 = dble(spindy%SH%O(j)%q(q))

            do s1=1,d
             q1 = this%basis%kind(kind)%m(s1)

             do s2=1,d
              M2 = nint(2.0d0*q2)
              M = nint(2.0d0*q2)
              M1 = nint(2.0d0*q1)

              !---------------------------------------------------------
              ! Evaluate the Stevens operator matrix element using the
              ! equivalent-operator formalism.
              !
              ! The required Clebsch-Gordan coefficients together with the
              ! spin-tensor normalisation factors are used to evaluate the
              ! matrix representation of O_k^q.
              !---------------------------------------------------------
              mat_elem = (0.0d0,0.0d0)

              if(spindy%SH%O(j)%q(q).eq.0)then
               ! Real tesseral component (q = 0)
               a = cleb(T2,M2,T,M,T1,M1)

               mat_elem = cmplx(a,0.0d0,8)

               call spintens(k,spin,a)

               mat_elem = mat_elem*a/sqrt(2.0d0*spin+1.0d0)

              elseif(spindy%SH%O(j)%q(q).gt.0)then
               ! Real tesseral component (q > 0)
               if(abs(nint(q2/2.0d0)-q2/2.0d0).lt.1.0d-5)then
                a = cleb(T2,M2,T,M,T1,M1)
                M = -M
                a = a+cleb(T2,M2,T,M,T1,M1)
                mat_elem = cmplx(a,0.0d0,8)

                call spintens(k,spin,a)

                mat_elem = mat_elem*a/sqrt(4.0d0*spin+2.0d0)
               else
                a = -cleb(T2,M2,T,M,T1,M1)
                M = -M
                a = a+cleb(T2,M2,T,M,T1,M1)
                mat_elem = cmplx(a,0.0d0,8)

                call spintens(k,spin,a)

                mat_elem = mat_elem*a/sqrt(4.0d0*spin+2.0d0)
               endif
              else
               ! Imaginary tesseral component (q < 0)
               if(abs(nint(q2/2.0d0)-q2/2.0d0).lt.1.0d-5)then
                a = cleb(T2,M2,T,M,T1,M1)
                M = -M
                a = a-cleb(T2,M2,T,M,T1,M1)
                mat_elem = cmplx(0.0d0,a,8)
                call spintens(k,spin,a)

                mat_elem = mat_elem*a/sqrt(4.0d0*spin+2.0d0)
               else
                a = cleb(T2,M2,T,M,T1,M1)
                M = -M
                a = a+cleb(T2,M2,T,M,T1,M1)
                mat_elem = cmplx(0.0d0,a,8)

                call spintens(k,spin,a)

                mat_elem = mat_elem*a/sqrt(4.0d0*spin+2.0d0)
               endif
              endif

              ! Apply the corresponding Stevens crystal-field parameter.
              val_tmp = mat_elem*spindy%SH%O(j)%B(q)

              Omat(s1,s2) = Omat(s1,s2)+val_tmp
             enddo
            enddo
           enddo
          enddo

          !---------------------------------------------------------------
          ! Add the crystal-field operator to the Hamiltonian MPO if at
          ! least one Stevens contribution was generated.
          !---------------------------------------------------------------
          if(any(abs(Omat)>0.0d0))then
           call add_one_site_operator(this,i,Omat)
          endif

          deallocate(Omat)
         enddo

        return
        end subroutine add_stevens_terms

        subroutine add_singleion_terms(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i,j
        integer                 :: kind
        integer                 :: d
        complex(8)              :: Dmat(3,3)
        complex(8), allocatable :: Sx(:,:)
        complex(8), allocatable :: Sy(:,:)
        complex(8), allocatable :: Sz(:,:)
        complex(8), allocatable :: Omat(:,:)

         !--------------------------------------------------------------
         ! Add the single-ion anisotropy terms to the Hamiltonian MPO.
         !
         !     H_DSI = S·D·S
         !           = Σ_{αβ} D_{αβ} S_α S_β,
         !
         ! where α,β ∈ {x,y,z}.
         !
         ! The D tensor is associated with a SPIN KIND, not a physical
         ! site. Therefore every physical site having the appropriate
         ! spin kind receives the corresponding DSI operator.
         !
         ! No secular approximation is made here. All nine Cartesian
         ! components of D are included.
         !
         ! The local spin operators are taken directly from
         !
         !     this%basis%kind(kind)%Sx
         !     this%basis%kind(kind)%Sy
         !     this%basis%kind(kind)%Sz
         !---------------------------------------------------------------

         if(spindy%SH%nDSI.eq.0) return

         do i=1,spindy%nspins
          kind = spindy%kind(i)
          d = this%basis%kind(kind)%d

          !-------------------------------------------------------------
          ! Allocate local spin operators.
          ! Make local copies of the spin operators used to construct
          ! the single-ion anisotropy operator.
          !-------------------------------------------------------------
          allocate(Sx(d,d))
          allocate(Sy(d,d))
          allocate(Sz(d,d))

          Sx = this%basis%kind(kind)%Sx
          Sy = this%basis%kind(kind)%Sy
          Sz = this%basis%kind(kind)%Sz

          allocate(Omat(d,d))
          Omat = (0.0d0,0.0d0)

          !-------------------------------------------------------------
          ! Add all DSI tensors belonging to this spin kind.
          ! Multiple DSI definitions for the same kind are therefore
          ! accumulated into one local operator before being added to
          ! the MPO.
          !-------------------------------------------------------------
          do j=1,spindy%SH%nDSI
           if(spindy%SH%DSI(j)%kind.ne.kind) cycle
           Dmat = spindy%SH%DSI(j)%D

           !-----------------------------------------------------------
           ! H_DSI = sum_{alpha,beta} D(alpha,beta) S_alpha S_beta
           !
           ! alpha,beta:
           !     1 -> x
           !     2 -> y
           !     3 -> z
           !-----------------------------------------------------------
           Omat = Omat &
                + Dmat(1,1)*matmul(Sx,Sx) &
                + Dmat(1,2)*matmul(Sx,Sy) &
                + Dmat(1,3)*matmul(Sx,Sz) &
                + Dmat(2,1)*matmul(Sy,Sx) &
                + Dmat(2,2)*matmul(Sy,Sy) &
                + Dmat(2,3)*matmul(Sy,Sz) &
                + Dmat(3,1)*matmul(Sz,Sx) &
                + Dmat(3,2)*matmul(Sz,Sy) &
                + Dmat(3,3)*matmul(Sz,Sz)
          enddo

          !-------------------------------------------------------------
          ! Insert the accumulated local anisotropy operator into the
          ! Hamiltonian MPO.
          !-------------------------------------------------------------
          call add_one_site_operator(this,i,Omat)

          deallocate(Sx)
          deallocate(Sy)
          deallocate(Sz)
          deallocate(Omat)
         enddo

        return
        end subroutine add_singleion_terms

        subroutine add_exchange_terms(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i,j,k
        integer                 :: kind1,kind2
        integer                 :: d1,d2
        complex(8)              :: Jval
        complex(8), allocatable :: op1(:,:),op2(:,:)

         !---------------------------------------------------------------
         ! Add the isotropic exchange interactions to the Hamiltonian MPO.
         !
         !     H_J = J S_i . S_j
         !         = J Sz_i Sz_j
         !         + J/2 S+_i S-_j
         !         + J/2 S-_i S+_j
         !
         ! Each exchange parameter is associated with a pair of spin
         ! kinds. The corresponding two-site operator is constructed in
         ! the direct-product basis and added to every matching pair of
         ! physical spins.
         !---------------------------------------------------------------

         if(spindy%SH%nJ.eq.0) return

         do k=1,spindy%SH%nJ
          kind1 = spindy%SH%J(k)%kind(1)
          kind2 = spindy%SH%J(k)%kind(2)
          Jval = spindy%SH%J(k)%J

          !-------------------------------------------------------------
          ! Loop over all unique pairs of physical spins.
          ! Restricting to i < j avoids double counting.
          !-------------------------------------------------------------
          do i=1,spindy%nspins-1
           do j=i+1,spindy%nspins

            !-----------------------------------------------------------
            ! Check whether the current pair of physical spins matches
            ! the spin kinds associated with this exchange interaction.
            ! Either ordering of different spin kinds is accepted.
            !-----------------------------------------------------------
            if(.not.((spindy%kind(i).eq.kind1.and.spindy%kind(j).eq.kind2) .or. &
                     (spindy%kind(i).eq.kind2.and.spindy%kind(j).eq.kind1))) then
             cycle
            endif

            d1 = this%basis%kind(spindy%kind(i))%d
            d2 = this%basis%kind(spindy%kind(j))%d

            !-----------------------------------------------------------
            ! Allocate the operator matrices acting on each spin.
            !-----------------------------------------------------------
            allocate(op1(d1,d1))
            allocate(op2(d2,d2))

            !-----------------------------------------------------------
            ! Add the longitudinal exchange contribution
            !
            !     J Sz_i Sz_j.
            !-----------------------------------------------------------
            op1 = Jval*this%basis%kind(spindy%kind(i))%Sz
            op2 = this%basis%kind(spindy%kind(j))%Sz

            call add_two_site_operator(this,i,j,op1,op2)

            !-----------------------------------------------------------
            ! Add the spin-flip contribution
            !
            !     (J/2) S+_i S-_j.
            !-----------------------------------------------------------
            op1 = 0.5d0*Jval*this%basis%kind(spindy%kind(i))%Sp
            op2 = this%basis%kind(spindy%kind(j))%Sm

            call add_two_site_operator(this,i,j,op1,op2)

            !-----------------------------------------------------------
            ! Add the reverse spin-flip contribution
            !
            !     (J/2) S-_i S+_j.
            !-----------------------------------------------------------
            op1 = 0.5d0*Jval*this%basis%kind(spindy%kind(i))%Sm
            op2 = this%basis%kind(spindy%kind(j))%Sp

            call add_two_site_operator(this,i,j,op1,op2)

            deallocate(op1)
            deallocate(op2)
           enddo
          enddo
         enddo

        return
        end subroutine add_exchange_terms

        subroutine add_exchange_anisotropy_terms(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i
        integer                 :: s1,s2
        integer                 :: kind1,kind2
        integer                 :: a,c,p,q
        integer                 :: d1,d2
        logical                 :: add_term
        complex(8)              :: coeff,D(3,3)
        complex(8), allocatable :: op1(:,:),op2(:,:)

         !---------------------------------------------------------------
         ! Add the exchange-anisotropy interactions to the Hamiltonian
         ! MPO.
         !
         ! For each exchange-anisotropy tensor
         !
         !     H_D2S = sum_{a,b=x,y,z} D(a,b) S_a(s1) S_b(s2),
         !
         ! the interaction is decomposed into its nine separable
         ! Cartesian products
         !
         !     D(a,b) S_a(s1) S_b(s2),
         !
         ! each of which is added individually to the Hamiltonian MPO.
         !
         ! Each D2S tensor is associated with a pair of spin kinds and
         ! is applied to every matching pair of physical spins whose
         ! separation satisfies
         !
         !     dist <= 0.5.
         !
         ! No secular approximation is made; all Cartesian spin
         ! components are retained.
         !---------------------------------------------------------------

         if(spindy%SH%nD2S.eq.0) return

         do s1=1,spindy%nspins-1
          kind1 = spindy%kind(s1)

          do s2=s1+1,spindy%nspins
           kind2 = spindy%kind(s2)

           !-----------------------------------------------------------
           ! Determine the corresponding physical spin position used
           ! by the MolForge distance array.
           !
           ! In the replicated-spin representation, s2 can refer to
           ! a spin in another copy.
           !-----------------------------------------------------------
           p = MOD((s2-1),spindy%nspins_pr)+1
           q = INT((s2-1)/spindy%nspins_pr)+1

           !-----------------------------------------------------------
           ! Check every D2S tensor defined in the input.
           !-----------------------------------------------------------
           do i=1,spindy%SH%nD2S
            add_term = .false.

            !---------------------------------------------------------
            ! The spin kinds must match the kinds specified by the
            ! D2S tensor.
            !---------------------------------------------------------
            if(kind1.eq.spindy%SH%D2S(i)%kind(1).and.kind2.eq.spindy%SH%D2S(i)%kind(2))then
             add_term = .true.

            elseif(kind1.eq.spindy%SH%D2S(i)%kind(2).and.kind2.eq.spindy%SH%D2S(i)%kind(1).and. &
                   spindy%SH%D2S(i)%kind(1).ne.spindy%SH%D2S(i)%kind(2))then
             add_term = .true.
            endif

            if(.not.add_term) cycle

            !---------------------------------------------------------
            ! The D2S tensor is only applicable when the two physical
            ! spins are sufficiently close.
            !---------------------------------------------------------
            if(spindy%dist(s1,1,p,q).gt.0.5d0) cycle

            !---------------------------------------------------------
            ! Local dimensions.
            !---------------------------------------------------------
            d1 = this%basis%kind(kind1)%d
            d2 = this%basis%kind(kind2)%d

            allocate(op1(d1,d1))
            allocate(op2(d2,d2))

            D = spindy%SH%D2S(i)%D

            !-----------------------------------------------------------
            ! Add each Cartesian product
            !
            !     D(a,b) S_a(s1) S_b(s2)
            !
            ! separately to the Hamiltonian MPO.
            !-----------------------------------------------------------
            do a=1,3
             do c=1,3
              coeff = D(a,c)

              if(abs(coeff).eq.0.0d0) cycle

              select case(a)
              case(1)
               op1 = coeff*this%basis%kind(kind1)%Sx
              case(2)
               op1 = coeff*this%basis%kind(kind1)%Sy
              case(3)
               op1 = coeff*this%basis%kind(kind1)%Sz
              end select

              select case(c)
              case(1)
               op2 = this%basis%kind(kind2)%Sx
              case(2)
               op2 = this%basis%kind(kind2)%Sy
              case(3)
               op2 = this%basis%kind(kind2)%Sz
              end select

              call add_two_site_operator(this,s1,s2,op1,op2)
             enddo
            enddo

            deallocate(op1)
            deallocate(op2)
           enddo
          enddo
         enddo

        return
        end subroutine add_exchange_anisotropy_terms

        subroutine add_dipolar_terms(this,spindy)
        implicit none
        class(mpo_operator)     :: this
        class(spins_hilbert)    :: spindy
        integer                 :: i
        integer                 :: s1,s2
        integer                 :: kind1,kind2
        integer                 :: d1,d2
        integer                 :: a,b
        complex(8)              :: coeff,D(3,3)
        complex(8), allocatable :: op1(:,:),op2(:,:)

         !---------------------------------------------------------------
         ! Add the spin-spin dipolar interactions to the Hamiltonian MPO.
         !
         ! The dipolar interaction tensors are generated by set_dipolar(),
         ! which computes the pair-specific interaction tensors from the
         ! supplied g tensors, Bohr magnetons and spin-spin separation.
         !
         ! Each entry SH%Ddip therefore corresponds to one physical spin
         ! pair together with its associated Cartesian dipolar tensor D.
         !
         ! The Hamiltonian contribution is
         !
         !     H_dip = sum_{a,b=x,y,z} D(a,b) S_a(s1) S_b(s2)
         !
         ! The interaction is decomposed into its nine separable
         ! Cartesian products
         !
         !     D(a,b) S_a(s1) S_b(s2),
         !
         ! each of which is added individually to the Hamiltonian MPO.
         !
         ! No secular approximation is made; all Cartesian spin
         ! components are retained.
         !---------------------------------------------------------------

         if(spindy%SH%nDdip.eq.0) return

         !---------------------------------------------------------------
         ! Loop over every dipolar interaction stored in the
         ! spin Hamiltonian.
         !---------------------------------------------------------------
         do i=1,spindy%SH%nDdip
          s1 = spindy%SH%Ddip(i)%kind(1)
          s2 = spindy%SH%Ddip(i)%kind(2)

          kind1 = spindy%kind(s1)
          kind2 = spindy%kind(s2)

          d1 = this%basis%kind(kind1)%d
          d2 = this%basis%kind(kind2)%d

          allocate(op1(d1,d1))
          allocate(op2(d2,d2))

          D = spindy%SH%Ddip(i)%D

          !-------------------------------------------------------------
          ! Add each Cartesian product
          !
          !     D(a,b) S_a(s1) S_b(s2)
          !
          ! separately to the Hamiltonian MPO.
          !-------------------------------------------------------------
          do a=1,3
           do b=1,3
            coeff = D(a,b)

            if(abs(coeff).eq.0.0d0) cycle

            !----------------------------------------------------------
            ! Select the Cartesian spin operator acting on the first
            ! spin.
            !----------------------------------------------------------
            select case(a)

            case(1)
             op1 = coeff*this%basis%kind(kind1)%Sx

            case(2)
             op1 = coeff*this%basis%kind(kind1)%Sy

            case(3)
             op1 = coeff*this%basis%kind(kind1)%Sz

            end select

            !----------------------------------------------------------
            ! Select the Cartesian spin operator acting on the second
            ! spin.
            !----------------------------------------------------------
            select case(b)

            case(1)
             op2 = this%basis%kind(kind2)%Sx

            case(2)
             op2 = this%basis%kind(kind2)%Sy

            case(3)
             op2 = this%basis%kind(kind2)%Sz

            end select

            call add_two_site_operator(this,s1,s2,op1,op2)
           enddo
          enddo

          deallocate(op1)
          deallocate(op2)
         enddo

        return
        end subroutine add_dipolar_terms

        subroutine add_one_site_operator(this,site,op)
        implicit none
        class(mpo_operator)     :: this
        integer                 :: site
        complex(8)              :: op(:,:)
        type(mpo_operator)      :: tmp,result
        integer                 :: i,d

         !---------------------------------------------------------------
         ! Add a one-site operator to the Hamiltonian MPO.
         !
         ! A temporary bond-dimension-1 MPO representing the supplied
         ! local operator is constructed and added to the Hamiltonian.
         ! Compression is intentionally deferred until the full
         ! Hamiltonian has been assembled.
         !---------------------------------------------------------------

         !---------------------------------------------------------------
         ! Check that the requested site is valid.
         !---------------------------------------------------------------
         if (site.lt.1.or.site.gt.this%nsites) then
          if (mpi_id.eq.0) then
           write(6,*) 'Err: Invalid site index in add_one_site_operator.'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Check that the supplied operator is square.
         !---------------------------------------------------------------
         if (size(op,1).ne.size(op,2)) then
          if (mpi_id.eq.0) then
           write(6,*) 'Err: Operator supplied to add_one_site_operator is not square.'
          endif
          stop
         endif

         d = this%basis%kind(this%site(site)%spin_kind)%d

         if(size(op,1).ne.d) then
          if (mpi_id.eq.0) then
           write(6,*) 'Err: Operator dimension does not match local Hilbert space.'
          endif
          stop
         endif

         !---------------------------------------------------------------
         ! Initialise a bond-dimension-1 MPO.
         !---------------------------------------------------------------
         tmp%nsites = this%nsites

         allocate(tmp%site(tmp%nsites))

         tmp%basis = this%basis
         tmp%bond_dimension_max = this%bond_dimension_max
         tmp%trunc_tol = this%trunc_tol
         tmp%bond_dimension = 1
         tmp%trunc_error = 0.0d0
         tmp%is_hermitian = .true.

         !---------------------------------------------------------------
         ! Construct the temporary MPO.
         !
         ! Every site is assigned the identity operator except the
         ! requested site, which receives the supplied operator.
         !---------------------------------------------------------------
         do i=1,tmp%nsites
          d = this%basis%kind(this%site(i)%spin_kind)%d

          allocate(tmp%site(i)%W(1,d,d,1))

          tmp%site(i)%left_dim   = 1
          tmp%site(i)%right_dim  = 1
          tmp%site(i)%phys_in    = d
          tmp%site(i)%phys_out   = d
          tmp%site(i)%site_index = i
          tmp%site(i)%spin_kind  = this%site(i)%spin_kind

          if(i.eq.site)then
           tmp%site(i)%W(1,:,:,1) = transpose(op)
          else
           tmp%site(i)%W(1,:,:,1) = this%basis%kind(this%site(i)%spin_kind)%I
          endif
         enddo

         !---------------------------------------------------------------
         ! Add the temporary MPO to the Hamiltonian.
         !---------------------------------------------------------------
         call this%lin_com((1.0d0,0.0d0),tmp,(1.0d0,0.0d0),result)

         !---------------------------------------------------------------
         ! Replace the current Hamiltonian with the summed MPO.
         !---------------------------------------------------------------
         if(allocated(this%site)) deallocate(this%site)

         this%nsites = result%nsites
         this%basis = result%basis
         this%bond_dimension = result%bond_dimension
         this%bond_dimension_max = result%bond_dimension_max
         this%trunc_error = result%trunc_error
         this%is_hermitian = result%is_hermitian

         call move_alloc(result%site,this%site)

        return
        end subroutine add_one_site_operator

        subroutine add_two_site_operator(this,s1,s2,op1,op2)
        implicit none
        class(mpo_operator)                 :: this
        integer                             :: s1,s2
        complex(8), intent(in)              :: op1(:,:),op2(:,:)
        type(mpo_operator)                  :: tmp,result
        integer                             :: i,d

         !---------------------------------------------------------------
         ! Add a separable two-site operator to the Hamiltonian MPO.
         !
         ! A temporary product MPO is constructed with
         !
         !     site s1 : op1
         !     site s2 : op2
         !     all other sites : identity
         !
         ! corresponding to the operator
         !
         !     I ... I ⊗ op1 ⊗ I ... I ⊗ op2 ⊗ I ... I.
         !
         ! The temporary MPO is then added to the Hamiltonian MPO using
         ! MPO linear combination. No compression is performed here;
         ! the complete Hamiltonian is compressed once after all
         ! Hamiltonian terms have been added.
         !---------------------------------------------------------------

         if (s1.eq.s2) then
          if (mpi_id.eq.0) then
           write(6,*) 'Err: add_two_site_operator requires two different sites'
          endif
          stop
         endif

         if(s1.lt.1.or.s1.gt.this%nsites.or.s2.lt.1.or.s2.gt.this%nsites) then
          if (mpi_id.eq.0) then
           write(6,*) 'Err: invalid site index in add_two_site_operator'
          endif
          stop
         endif

         tmp%nsites = this%nsites

         allocate(tmp%site(tmp%nsites))

         tmp%basis = this%basis
         tmp%bond_dimension_max = this%bond_dimension_max
         tmp%trunc_tol = this%trunc_tol
         tmp%bond_dimension = 1
         tmp%trunc_error = 0.0d0
         tmp%is_hermitian = .true.

         !---------------------------------------------------------------
         ! Construct the temporary product MPO.
         !
         ! The requested operators are placed on sites s1 and s2, while
         ! identity operators are placed on every other site.
         !---------------------------------------------------------------
         do i=1,tmp%nsites
          d = this%basis%kind(this%site(i)%spin_kind)%d

          allocate(tmp%site(i)%W(1,d,d,1))

          if(i.eq.s1)then
           tmp%site(i)%W(1,:,:,1) = transpose(op1)
          elseif(i.eq.s2)then
           tmp%site(i)%W(1,:,:,1) = transpose(op2)
          else
           tmp%site(i)%W(1,:,:,1) = this%basis%kind(this%site(i)%spin_kind)%I
          endif

          tmp%site(i)%left_dim   = 1
          tmp%site(i)%right_dim  = 1
          tmp%site(i)%phys_in    = d
          tmp%site(i)%phys_out   = d
          tmp%site(i)%site_index = i
          tmp%site(i)%spin_kind  = this%site(i)%spin_kind
         enddo

         !---------------------------------------------------------------
         ! Add the temporary MPO to the Hamiltonian MPO.
         !---------------------------------------------------------------
         call this%lin_com((1.0d0,0.0d0),tmp,(1.0d0,0.0d0),result)

         call move_alloc(result%site,this%site)

         this%basis = result%basis
         this%nsites = result%nsites
         this%bond_dimension = result%bond_dimension
         this%trunc_error = result%trunc_error
         this%is_hermitian = result%is_hermitian

         if(allocated(tmp%site)) deallocate(tmp%site)

        return
        end subroutine add_two_site_operator

        end module tensor_network_class
