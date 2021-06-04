        module spins_dist_rs_class         
        use spinham_class
        use spin_phonon_class
        use lattice_class
        implicit none

        type  :: SpinBath
         double precision     :: temp
         integer,allocatable  :: group(:)
        end type SpinBath
             
        type, extends(bravais_lattice) :: spins_group
         integer                            :: nspins=0
         integer                            :: nspins_pr=0
         integer                            :: nkinds=0
         integer, allocatable               :: kind(:)         
         double precision, allocatable      :: spin(:)
         double precision, allocatable      :: bohr_mag(:)
         double precision                   :: Bfield(3)=0.0d0
         double precision, allocatable      :: x(:,:)
         double precision, allocatable      :: dist(:,:,:,:)
         integer, allocatable               :: tr_map(:,:)
         double precision, allocatable      :: alpha0(:)
         double precision, allocatable      :: beta0(:)
         double precision, allocatable      :: gamma0(:)
         type(SpinBath), allocatable        :: spin_bath(:)
         type(SpinHamiltonian)              :: SH 
         type(SpinPhononHamiltonian)        :: SPH  
         type(SpinPhononHamiltonian)        :: SPH2
         contains
         procedure   ::  spin_bcast
         procedure   ::  dist_ij
         procedure   ::  cart2frac
         procedure   ::  frac2cart
!         procedure   ::  projSH
        end type spins_group

        contains


        subroutine spin_bcast(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_group)   :: this
        integer              :: i


         call this%lattice_bcast()
         call mpi_bcast(this%nspins,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nspins_pr,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nkinds,1,mpi_integer,0,mpi_comm_world,err)

         if(.not. allocated(this%kind)) allocate(this%kind(this%nspins))
         call mpi_bcast(this%kind,this%nspins,mpi_integer,0,mpi_comm_world,err)

         if(.not. allocated(this%bohr_mag)) allocate(this%bohr_mag(this%nkinds))
         call mpi_bcast(this%bohr_mag,this%nkinds,mpi_double_precision,0,mpi_comm_world,err)

         if(.not. allocated(this%spin)) allocate(this%spin(this%nkinds))
         call mpi_bcast(this%spin,this%nkinds,mpi_double_precision,0,mpi_comm_world,err)

         if(.not. allocated(this%x)) allocate(this%x(this%nspins,3))
         do i=1,this%nspins
          call mpi_bcast(this%x(i,:),3,mpi_double_precision,0,mpi_comm_world,err)
         enddo

         if(this%ntot.gt.1)then
          if(.not. allocated(this%tr_map)) allocate(this%tr_map(this%nspins,3))
          do i=1,this%nspins
           call mpi_bcast(this%tr_map(i,1),1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%tr_map(i,2),1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%tr_map(i,3),1,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif

         if(.not. allocated(this%alpha0)) allocate(this%alpha0(this%nspins))
         call mpi_bcast(this%alpha0,this%nspins,mpi_double_precision,0,mpi_comm_world,err)

         if(.not. allocated(this%beta0)) allocate(this%beta0(this%nspins))
         call mpi_bcast(this%beta0,this%nspins,mpi_double_precision,0,mpi_comm_world,err)

         if(.not. allocated(this%gamma0)) allocate(this%gamma0(this%nspins))
         call mpi_bcast(this%gamma0,this%nspins,mpi_double_precision,0,mpi_comm_world,err)

         this%alpha0=0.0d0
         this%beta0=0.0d0
         this%gamma0=0.0d0

         call mpi_bcast(this%Bfield,3,mpi_double_precision,0,mpi_comm_world,err)
         
        return
        end subroutine spin_bcast


        subroutine dist_ij(this)
        implicit none
        class(spins_group)      :: this
        integer                 :: i,j,celli,cellj,v1,v2,v
        double precision        :: c(3),a(3),b(3)

         call this%cart2frac()
         allocate(this%dist(this%nspins_pr,this%ntot,this%nspins_pr,this%ntot))

         do i=1,this%nspins_pr
          do j=1,this%nspins_pr
           do celli=1,this%ntot
            do cellj=1,this%ntot
             v1=(celli-1)*this%nspins_pr+i
             v2=(cellj-1)*this%nspins_pr+j
             do v=1,3
              a(v)=this%x(i,v)+this%rcell(celli,v)
              b(v)=this%x(j,v)+this%rcell(cellj,v)
             enddo
             c(1)=a(1)-b(1)
             c(1)=c(1)-nint(c(1)/dble(this%nx))*this%nx
             c(2)=a(2)-b(2)
             c(2)=c(2)-nint(c(2)/dble(this%ny))*this%ny
             c(3)=a(3)-b(3)
             c(3)=c(3)-nint(c(3)/dble(this%nz))*this%nz
             this%dist(i,celli,j,cellj)=this%dist_dir(c)
            enddo
           enddo
          enddo
         enddo

         call this%frac2cart()


        return
        end subroutine dist_ij

        subroutine frac2cart(this)
        implicit none
        class(spins_group) :: this
        integer            :: i
         do i=1,this%nspins
          this%x(i,:)=matmul(this%J,this%x(i,:))
         enddo
        return
        end subroutine frac2cart

        subroutine cart2frac(this)
        implicit none
        class(spins_group) :: this
        integer            :: i
         do i=1,this%nspins
          this%x(i,:)=matmul(this%Jinv,this%x(i,:))
         enddo
        return
        end subroutine cart2frac


        end module spins_dist_rs_class

