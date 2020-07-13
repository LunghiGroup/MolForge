        module phonons_class
        implicit none

         logical                        :: read_fc3=.false.
         logical                        :: do_brillouin_path=.false.
         logical                        :: do_brillouin_mesh=.true.
         logical                        :: boundary_scattering=.false.
         double precision, allocatable  :: temp(:)
         integer                        :: ntemps
         double precision               :: smear
         double precision               :: Length

        type kpoint
         double precision :: k(3)
         double precision :: weight
         double precision, allocatable       :: freq(:)
         double precision, allocatable       :: vel(:,:)
         double precision, allocatable       :: width(:,:)
         double precision, allocatable       :: krta(:,:,:,:)
         double precision, allocatable       :: kscf(:,:,:,:)
         complex(8),allocatable              :: hess(:,:)
         double precision, allocatable       :: proj(:,:)
         contains
         procedure                           :: diagD
         procedure                           :: projHess
        end type kpoint

        type dist1D
         double precision, allocatable       :: dist(:)
         double precision                    :: sigma
         double precision                    :: step
         integer                             :: nsteps
        end type dist1D
      
        type  brillouin 
         type(kpoint), allocatable      :: list(:) 
         type(kpoint), allocatable      :: path(:) 
         integer, allocatable           :: kinv(:)
         integer                        :: nx
         integer                        :: ny
         integer                        :: nz
         integer                        :: ntot
         integer                        :: nloc
         integer                        :: k_start
         integer, allocatable           :: mpi_conn(:)
         double precision, allocatable  :: krta(:,:,:)
         double precision, allocatable  :: kscf(:,:,:)
         type(dist1D)                   :: dos1p
         type(dist1D)                   :: dos2p
         type(dist1D), allocatable      :: pdos1p(:)
         contains 
         procedure        :: generate_mesh
         procedure        :: generate_path
         procedure        :: brillouin_bcast
         procedure        :: calc_disps
         procedure        :: calc_bands
         procedure        :: get_cart_disp
         procedure        :: calc_vel
         procedure        :: calc_linewidth
         procedure        :: calc_linewidth_sp
         procedure        :: calc_krta         
         procedure        :: calc_kscf
         procedure        :: remap_hess
         procedure        :: get_V3
         procedure        :: mpi_dist_kpoints
         procedure        :: calc_dos1p
         procedure        :: calc_dos2p
         procedure        :: smooth_bands
        end type brillouin

        type tetrahedron
         integer     :: indx(4)          
        end type tetrahedron

        contains

        subroutine gen_vars_bcast
        use mpi
        use mpi_utils
        implicit none 
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

         call MPI_BCAST(ntemps,1,mpi_integer,0,mpi_comm_world,err)
         if(.not.allocated(temp)) allocate(temp(ntemps))
         call mpi_bcast(temp,ntemps,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(smear,1,mpi_double_precision,0,mpi_comm_world,err)
         call MPI_BCAST(read_fc3,1,mpi_logical,0,mpi_comm_world,err)
         call MPI_BCAST(boundary_scattering,1,mpi_logical,0,mpi_comm_world,err)
         call MPI_BCAST(do_brillouin_path,1,mpi_logical,0,mpi_comm_world,err)
         call MPI_BCAST(do_brillouin_mesh,1,mpi_logical,0,mpi_comm_world,err)

        return
        end subroutine gen_vars_bcast

        subroutine mpi_dist_kpoints(this)
        use mpi
        use mpi_utils
        implicit none
        class(brillouin)       :: this
        integer                :: i,k,j,rest

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)       

        ! create groups of k point to be assigned to each mpi process
         
         this%nloc=this%ntot/mpi_nproc
         rest=this%ntot-this%nloc*mpi_nproc

         if(mpi_id.lt.rest)then
          this%nloc=this%nloc+1
         endif


         ! create an array to map k to mpi_id

         allocate(this%mpi_conn(this%ntot))

         do i=0,mpi_nproc-1
          if(i.lt.rest)then
           k=1
          else
           k=0
          endif
          j=i*int(this%ntot/mpi_nproc)+min(rest,i)+1
          k=int(this%ntot/mpi_nproc)+j+k-1
          this%mpi_conn(j:k)=i
         enddo

         ! index where to start for each process

         this%k_start=mpi_id*int(this%ntot/mpi_nproc)+min(rest,mpi_id)+1

        return
        end subroutine mpi_dist_kpoints

        subroutine calc_dos1p(this)
        use units_parms
        implicit none
        class(brillouin)       :: this
        integer                :: i,j,k,l,v
        double precision       :: norm

         allocate(this%dos1p%dist(this%dos1p%nsteps))

         this%dos1p%dist=0.0d0

         do j=1,this%ntot
          do i=1,size(this%list(1)%freq)

           k=NINT(this%list(j)%freq(i)/this%dos1p%step)
           l=NINT(this%dos1p%sigma/this%dos1p%step)

           do v=-3*l,3*l
            if((v+k).gt.0 .and. (v+k).lt.this%dos1p%nsteps)then
             this%dos1p%dist(k+v)=this%dos1p%dist(k+v)+delta(DBLE(v),DBLE(l))
            endif
           enddo

          enddo 
         enddo

         open(100,file='dos1p.dat')

         norm=0.0d0
         do i=1,this%dos1p%nsteps
          norm=norm+this%dos1p%dist(i)*this%dos1p%step
         enddo

         this%dos1p%dist=this%dos1p%dist/norm

         norm=0.0d0
         do i=1,this%dos1p%nsteps
          norm=norm+this%dos1p%step
          write(100,*) norm,this%dos1p%dist(i)
         enddo

         close(100)

        return
        end subroutine calc_dos1p

        subroutine calc_dos2p(this)
        use mpi
        use mpi_utils
        use units_parms
        implicit none
        class(brillouin)       :: this
        integer                :: i,j,k,l,v
        integer                :: i1,j1,i2,j2,k1,k2
        double precision       :: norm

         allocate(this%dos2p%dist(this%dos2p%nsteps))
         this%dos2p%dist=0.0d0

         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         j1=this%k_start

         do j=1,this%ntot
         do i1=1,size(this%list(1)%freq)

          if(this%list(j1)%freq(i1) .gt. 500.0d0) cycle
          if(this%list(j1)%freq(i1) .lt.0.0d0) cycle
          if(j1.eq.1 .and. i1.le.3) cycle

          do j2=1,this%ntot
          do i2=1,size(this%list(1)%freq)

           if(this%list(j2)%freq(i2) .gt. 500.0d0) cycle
           if(this%list(j2)%freq(i2) .lt.0.0d0) cycle
           if(j2.eq.1 .and. i2.le.3) cycle

           k1=(j1-1)*size(this%list(1)%freq)+i1
           k2=(j2-1)*size(this%list(1)%freq)+i2
           if(k1.lt.k2)cycle

           k=NINT(abs(this%list(j1)%freq(i1)-this%list(j2)%freq(i2))/this%dos2p%step)
           l=NINT(this%dos2p%sigma/this%dos2p%step)

           do v=-3*l,3*l
            if((v+k).gt.0 .and. (v+k).lt.this%dos2p%nsteps)then
            this%dos2p%dist(k+v)=this%dos2p%dist(k+v)+delta(DBLE(v),DBLE(l))*&
                ((bose(20.0d0,this%list(j1)%freq(i1))*&
                 (bose(20.0d0,this%list(j2)%freq(i2))+1.0d0))+&
                 (bose(20.0d0,this%list(j2)%freq(i2))*&
                 (bose(20.0d0,this%list(j1)%freq(i1))+1.0d0)))
            endif
           enddo

          enddo
          enddo
         enddo
         j1=j1+1
         enddo

         do i=1,this%dos2p%nsteps        
          norm=0.0d0
          call mpi_allreduce(this%dos2p%dist(i),norm,1,&
                mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,err)    
          this%dos2p%dist(i)=norm
         enddo

         if(mpi_id.eq.0)then

          open(100,file='dos2p.dat')

          norm=0.0d0
          do i=1,this%dos2p%nsteps
           norm=norm+this%dos2p%dist(i)*this%dos2p%step
          enddo
          this%dos2p%dist=this%dos2p%dist/norm

          norm=0.0d0
          do i=1,this%dos2p%nsteps
           norm=norm+this%dos2p%step
           write(100,*) norm,this%dos2p%dist(i)
          enddo

          close(100)

         endif

        return
        end subroutine calc_dos2p

        subroutine calc_disps(this,sys)
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,k,j,rest,s
        double precision       :: dist

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

         do i=1,size(this%path)
          call this%path(i)%diagD(sys) 
         enddo
        
         if(mpi_id.eq.0)then
          open(12,file='disp.dat')
          open(13,file='disp_modes.dat')
          dist=0.0d0
          do i=1,size(this%path)
           if(i.gt.1)then
            dist=dist+sqrt(sys%dist_rec(this%path(i)%k(:)-this%path(i-1)%k(:)))
           endif
           write(12,*) this%path(i)%k(:),dist,(this%path(i)%freq(s),s=1,size(this%path(i)%freq))
           do s=1,sys%nats*3
            write(13,*) 'K:',this%path(i)%k(:),'Mode:',s
            do j=1,sys%nats*3
             write(13,*) this%path(i)%k(:),dist,dble(this%path(i)%hess(j,s)),aimag(this%path(i)%hess(j,s))
            enddo
           enddo
          enddo
          close(12)
          close(13)
         endif

        return
        end subroutine calc_disps

        subroutine smooth_bands(this)
        use lapack_inverse
        implicit none
        class(brillouin)              :: this
        integer                       :: i,j,l,npt(3),inf,dimA,dimB,lwork
        double precision, allocatable :: A(:,:),B(:),work(:)

        
         npt=0
         do i=1,this%ntot
          do j=1,3        
           if(this%list(i)%freq(j).gt.2.0d0 .and. this%list(i)%freq(j).lt.7.0d0 )then
            npt(j)=npt(j)+1
           endif
          enddo
         enddo

         do j=1,3        

          l=1
          allocate(A(npt(j)+3,3))
          allocate(B(npt(j)+3))
          A=0.0d0
          B=0.0d0

          do i=1,this%ntot
           if(this%list(i)%freq(j).gt.2.0d0 .and. this%list(i)%freq(j).lt.7.0d0 )then
            A(l,1)=abs(this%list(i)%k(1))
            A(l,2)=abs(this%list(i)%k(2))
            A(l,3)=abs(this%list(i)%k(3))
            B(l)=this%list(i)%freq(j)
            write(*,*) l,i,this%list(i)%freq(j)
            l=l+1
           endif
          enddo

          do i=1,3 ! L2 norm
           A(npt(j)+i,i)=0.00001d0
           B(npt(j)+i)=0.0d0
          enddo

          dimA=size(A,2)
          dimB=size(B,1)

          lwork=dimB+64*dimB+1000
          allocate(work(lwork))

          call dgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf) 

          if(inf.ne.0)then
           write(*,*) 'dgels failed'
           stop
          endif
          deallocate(work)

          do i=1,this%ntot ! smooth
           if(this%list(i)%freq(j).lt.5.0d0 )then
            this%list(i)%freq(j)=0.0d0
            do l=1,3
             this%list(i)%freq(j)=this%list(i)%freq(j)+B(l)*abs(this%list(i)%k(l))
            enddo
           endif
          enddo

          deallocate(A)
          deallocate(B)

         enddo

        return
        end subroutine smooth_bands

        subroutine get_cart_disp(this,sys,ki,kj,hess)                
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        class(brillouin)             :: this
        type(atoms_group)            :: sys
        integer                      :: ki,kj,l,i,j
        double precision             :: mass1,coeff
        double complex, allocatable  :: hess(:)

         if(allocated(hess)) deallocate(hess)
         allocate(hess(size(this%list(ki)%hess,1)))
         hess=(0.0d0,0.0d0)

         l=1
         do i=1,sys%nats
          mass1=sys%mass(sys%kind(i))*1822.89d0
          coeff=0.5291772d0/sqrt(mass1*this%ntot*this%list(ki)%freq(kj)/219474.6313702d0)
          do j=1,3
           hess(l)=coeff*this%list(ki)%hess(l,kj)
           l=l+1
          enddo
         enddo       

        return
        end subroutine get_cart_disp

        subroutine calc_bands(this,sys)
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,k,j,rest

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

         j=this%k_start

         do i=1,this%nloc
          call this%list(j)%diagD(sys) 
          j=j+1
         enddo

         ! broadcast calculated bands 

         do i=1,this%ntot
          if(.not.allocated(this%list(i)%freq))then
           allocate(this%list(i)%freq(sys%nats*3))
!           allocate(this%list(i)%hess(sys%nats*3,sys%nats*3))
          endif
          call mpi_bcast(this%list(i)%freq,size(this%list(i)%freq),mpi_double_precision,this%mpi_conn(i),mpi_comm_world,err)
!          do j=1,size(this%list(i)%freq)           
!           call mpi_bcast(this%list(i)%hess(:,j),size(this%list(i)%freq),mpi_double_complex,this%mpi_conn(i),mpi_comm_world,err)
!          enddo
         enddo

         call this%smooth_bands()        
         this%list(1)%freq(1:3)=0.0d0           

        return
        end subroutine calc_bands

        
        subroutine remap_hess(this,kp,kpp,i,j,overlap)
        use units_parms
        implicit none
        class(brillouin)       :: this
        integer                :: i,j,v,kp,kpp,l,dims
        double precision       :: overlap
        double precision, allocatable :: prod(:)
        complex(8), allocatable :: prodc(:)

         dims=size(this%list(kp)%hess,1)
         allocate(prod(dims))
         allocate(prodc(dims))

         do v=1,dims          
          prodc(v)=dot_product(this%list(kp)%hess(:,i),this%list(kpp)%hess(:,v))
          prod(v)=conjg(prodc(v))*prodc(v)
         enddo

         j=maxloc(abs(prod),dim=1)
         overlap=abs(prod(j))

        return
        end subroutine remap_hess

        subroutine calc_vel(this,sys)
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,j,v,kp,kpp
        double precision       :: dist(3),overlap

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

        ! initialize group velocity

         do i=1,this%ntot
          allocate(this%list(i)%vel(sys%nats*3,3))
          this%list(i)%vel=0.0d0
         enddo

         if(mpi_id.eq.0) open(1212,file='vel_check')

        ! for each local k-point run over all the others and compute the
        ! velocity

         kp=this%k_start

         do v=1,this%nloc
          do kpp=1,this%ntot
                    
          dist(1)=this%list(kp)%k(1)-this%list(kpp)%k(1)
          dist(1)=dist(1)-nint(dist(1))
          dist(2)=this%list(kp)%k(2)-this%list(kpp)%k(2)
          dist(2)=dist(2)-nint(dist(2))
          dist(3)=this%list(kp)%k(3)-this%list(kpp)%k(3)
          dist(3)=dist(3)-nint(dist(3))

          if(abs(dist(2)).lt.1.0e-5  .and. &
             abs(dist(3)).lt.1.0e-5 ) then

           if(abs(dist(1)-1.0d0/this%nx).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j) write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,1)=this%list(kp)%vel(i,1)-this%list(kpp)%freq(j)/(2*1.0d0/this%nx)
            enddo
           endif
           if(abs(dist(1)+1.0d0/this%nx).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j) write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,1)=this%list(kp)%vel(i,1)+this%list(kpp)%freq(j)/(2*1.0d0/this%nx)
            enddo
           endif

          endif

          if(abs(dist(1)).lt.1.0e-5  .and. &
             abs(dist(3)).lt.1.0e-5 ) then

           if(abs(dist(2)-1.0d0/this%ny).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j)write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,2)=this%list(kp)%vel(i,2)-this%list(kpp)%freq(j)/(2*1.0d0/this%ny)
            enddo
           endif
           if(abs(dist(2)+1.0d0/this%ny).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j)write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,2)=this%list(kp)%vel(i,2)+this%list(kpp)%freq(j)/(2*1.0d0/this%ny)
            enddo
           endif

          endif

          if(abs(dist(1)).lt.1.0e-5  .and. &
             abs(dist(2)).lt.1.0e-5 ) then

           if(abs(dist(3)-1.0d0/this%nz).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j)write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,3)=this%list(kp)%vel(i,3)-this%list(kpp)%freq(j)/(2*1.0d0/this%nz)
            enddo
           endif
           if(abs(dist(3)+1.0d0/this%nz).lt.1.0e-6 )then
            do i=1,sys%nats*3
              call this%remap_hess(kp,kpp,i,j,overlap)
              if(i.ne.j)write(1212,*) kp,kpp,i,j,'+',overlap
              this%list(kp)%vel(i,3)=this%list(kp)%vel(i,3)+this%list(kpp)%freq(j)/(2*1.0d0/this%nz)
            enddo
           endif

          endif

          enddo
          kp=kp+1
         enddo        

         kp=this%k_start
         do v=1,this%nloc
          do i=1,sys%nats*3
           this%list(kp)%vel(i,:)=matmul(sys%J,this%list(kp)%vel(i,:))
          enddo
          kp=kp+1
         enddo


         ! broadcast calculated bands 

         do i=1,this%ntot
          do v=1,3
           call mpi_bcast(this%list(i)%vel(:,v),size(this%list(i)%vel,1),mpi_double_precision,this%mpi_conn(i),mpi_comm_world,err)
          enddo
         enddo
        
         if(mpi_id.eq.0)then 
           write(*,*) '   Group Velocities:'
          do i=1,this%ntot
            write(*,*) '      ',i,this%list(i)%k(:)
           do j=1,sys%nats*3
            write(*,*) '          ',j,this%list(i)%vel(j,:)
           enddo
          enddo
         endif


        return
        end subroutine calc_vel

        subroutine brillouin_bcast(this)
        use mpi
        use mpi_utils
        implicit none
        class(brillouin)     :: this

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)       

         call MPI_BCAST(this%nx,1,mpi_integer,0,mpi_comm_world,err) 
         call MPI_BCAST(this%ny,1,mpi_integer,0,mpi_comm_world,err) 
         call MPI_BCAST(this%nz,1,mpi_integer,0,mpi_comm_world,err) 
         call MPI_BCAST(this%ntot,1,mpi_integer,0,mpi_comm_world,err) 

         if(.not.allocated(this%list))then
          call this%generate_mesh(this%nx,this%ny,this%nz)
         endif
                    
        return
        end subroutine brillouin_bcast

        subroutine generate_path(this,npts,kf,ki)
        implicit none
        class(brillouin) :: this
        double precision :: step(3)
        integer          :: v,s,i
        double precision, allocatable :: ki(:,:),kf(:,:)
        integer, allocatable          :: npts(:)

         allocate(this%path(sum(npts)))
         v=1
         do s=1,size(npts)
          step(1)=(kf(1,s)-ki(1,s))/dble(npts(s))
          step(2)=(kf(2,s)-ki(2,s))/dble(npts(s))
          step(3)=(kf(3,s)-ki(3,s))/dble(npts(s))
          do i=1,npts(s)
           this%path(v)%k(1)=ki(1,s)+(i-1)*step(1)
           this%path(v)%k(2)=ki(2,s)+(i-1)*step(2)
           this%path(v)%k(3)=ki(3,s)+(i-1)*step(3)
           v=v+1
          enddo
         enddo
         
        return
        end subroutine generate_path

        subroutine generate_mesh(this,nx,ny,nz)
        implicit none
        class(brillouin) :: this
        double precision :: step(3)
        integer          :: i,j,k,s,v,nx,ny,nz

         this%nx=nx
         this%ny=ny
         this%nz=nz
         this%ntot=this%nx*this%ny*this%nz
         allocate(this%list(this%ntot))
         step(1)=1.0d0/this%nx
         step(2)=1.0d0/this%ny
         step(3)=1.0d0/this%nz

         v=1
         this%list(:)%weight=1.0d0/this%ntot

         do k=0,this%nz-1
          do s=0,this%ny-1
           do i=0,this%nx-1
            this%list(v)%k(1)=i*step(1)
            this%list(v)%k(2)=s*step(2)
            this%list(v)%k(3)=k*step(3)
            if(this%list(v)%k(1).gt.0.5d0)this%list(v)%k(1)=this%list(v)%k(1)-1.0d0
            if(this%list(v)%k(2).gt.0.5d0)this%list(v)%k(2)=this%list(v)%k(2)-1.0d0
            if(this%list(v)%k(3).gt.0.5d0)this%list(v)%k(3)=this%list(v)%k(3)-1.0d0
            v=v+1  
           enddo
          enddo
         enddo    


         allocate(this%kinv(this%ntot))
         this%kinv(:)=0

         do i=1,this%ntot
          do j=i,this%ntot

           if(abs(this%list(i)%k(1)+this%list(j)%k(1)).lt.1.0d-5 .and. &
              abs(this%list(i)%k(2)+this%list(j)%k(2)).lt.1.0d-5 .and. &
              abs(this%list(i)%k(3)+this%list(j)%k(3)).lt.1.0d-5 ) then

              this%kinv(i)=j
              this%kinv(j)=i

           endif

          enddo
          if(this%kinv(i).eq.0) this%kinv(i)=i
         enddo

         call this%mpi_dist_kpoints()

        return
        end subroutine generate_mesh

        subroutine diagD(ph,sys)
        use lattice_class
        use atoms_class
        implicit none
        class(kpoint)                 :: ph
        type(atoms_group)             :: sys
        integer                       :: v1,v2,i1,i2,s1,s2,l
        double precision              :: A,pi
        double precision, allocatable :: mass(:)

         pi=acos(-1.0d0)

         if(.not.allocated(ph%freq)) allocate(ph%freq(sys%nats*3))
         if(.not.allocated(ph%hess)) allocate(ph%hess(sys%nats*3,sys%nats*3))

         ph%freq=0.0d0
         ph%hess=(0.0d0,0.0d0)

         do l=1,sys%nx*sys%ny*sys%nz
        
          A=DOT_PRODUCT(sys%rcell(l,:),ph%k(:))

          do i1=1,sys%nats
          do i2=1,sys%nats
           do s1=1,3
           do s2=1,3

            v1=(i1-1)*3+s1
            v2=(i2-1)*3+s2

            ph%hess(v1,v2)=ph%hess(v1,v2)+sys%fcs2(l,i1,s1,i2,s2)*exp(CMPLX(0.0d0,2*pi*A,8))

           enddo
           enddo
          enddo
          enddo
         enddo

         allocate(mass(sys%nats))
         do i1=1,sys%nats
          mass(i1)=sys%mass(sys%kind(i1))
         enddo

         call diaghess(ph%hess,sys%nats,mass,ph%freq)

         deallocate(mass)

        return
        end subroutine diagD

        subroutine projHess(this,sys)
        use units_parms
        use atoms_class
        use proj_disp_class
        implicit none
        class(kpoint)                 :: this
        type(atoms_group)             :: sys
        type(molecule)                :: mol
        integer                       :: i1,s1,i,j,s,loc_nats
        double precision, allocatable :: mass(:),geo(:,:)
        double precision              :: mass1,coeff,norms(5)

         if(.not. allocated(this%proj))then
          allocate(this%proj(size(this%hess,2),3+sys%nkinds))
         endif

         this%proj=0.0d0
         loc_nats=sys%nats

         allocate(mass(loc_nats))
         do i1=1,loc_nats
          mass(i1)=sys%mass(sys%kind(i1))
         enddo

         allocate(geo(loc_nats,3))
         geo=sys%x(1:loc_nats,1:3)

         call mol%def_mol(geo,mass)
         geo=0.0d0

        ! def distorted geos 

         do s=1,size(this%hess,2)

          do j=1,loc_nats*3
           mass1=sys%mass(sys%kind((2+j)/3))
           coeff=bohr2ang/dsqrt(mass1*1822.89/219474.6313702)
           i1=(2+j)/3
           s1=mod(j-1,3)+1
           geo(i1,s1)=coeff*dble(this%hess(j,s))/sqrt(this%freq(s))
          enddo

          geo=geo+sys%x(1:loc_nats,:)

          call mol%def_mol_dist(geo)
          call mol%proj_disp()
          call mol%get_norms(norms)

          this%proj(s,1:3)=norms(2:4)/norms(1)
          if (norms(5).gt.1.0e-8) &
               write(*,*) 'mode ',s,' decomposition failed'

         enddo

         deallocate(mass)
         deallocate(geo)

         ! proj by kind

         do s=1,size(this%hess,2)
          do j=1,sys%nats
           i1=(2+j)/3
           i=sys%kind(i1)+3
           this%proj(s,i)=this%proj(s,i)+dble(this%hess(j,s)*conjg(this%hess(j,s)))
          enddo
         enddo

        return
        end subroutine projHess

        subroutine diaghess(hess,nat,mass,ener)
        use lapack_diag_simm
        implicit none
        integer v,s,val,resto,j,i,N,ialloc,t,nat,k
        double precision :: a,norm,mass(nat),mass2(3*nat),ener(3*nat)
        double complex   :: hess(3*nat,3*nat)
        
         N=nat*3

         i=1
         do j=1,nat
          mass2(i)=mass(j)
          mass2(i+1)=mass(j)
          mass2(i+2)=mass(j)
          i=i+3
         enddo
        
         do i=1,N
          do v=1,N
           hess(i,v)=hess(i,v)/dsqrt(mass2(i)*mass2(v))
          enddo           
         enddo

         call new_diag(N,hess,ener)

         do i=1,N               

          if(ener(i).gt.0.0d0)then
           ener(i)=sqrt(ener(i)/1822.89)    
!!!! 1822.89 turn the amu units of mass in a.u. (electron mass). Ener(i) alla fine è a.u.
!           do t=1,N  
!!!! the eigenvectors are transformed in cartesian displacements (angstrom) associated with the unit-less normal mode
!            hess(t,i)=hess(t,i)*0.5291772/dsqrt(mass2(t)*1822.89*ener(i))   
!!!! 1822.89 turn mass from amu to a.u. 
!!!! 0.5291 turn bohr in Ang
!           enddo

          ener(i)=ener(i)*219474.6313702

!!!! 1/(c[a.u.]*bohr2cm*2pi) gives ener as frequency in cm-1
!!!! o in modo equivalente
!!!! 1/( 3E8*1E2*0.036752*4.13566E-15) = 1/(c(m/s)*(m2cm)*(eV2Eh)*h(eVs)) turn ener(i) (t^{-1} [a.u.]) in cm^{-1}

          else
                
!           ener(i)=0.0d0
!           hess(:,i)=0.0d0

          endif
         enddo
       
        return
        end subroutine diaghess

        function get_V3(this,sys,q1,q2,q3,j1,j2,j3) result(V3)
        use lattice_class 
        use atoms_class
        implicit none
        class(brillouin)      :: this
        type(atoms_group)     :: sys
        complex(8)            :: V3
        double precision      :: A,B,pi,mass1,mass2,mass3,coeff
        integer               :: i1,i2,i3,j1,j2,j3,s2,s3,q1,q2,q3,i

         V3=(0.0d0,0.0d0)
         pi=acos(-1.0d0)

         do i=1,sys%fcs3%nfcs

          i1=sys%fcs3%nat(i,1)
          i2=sys%fcs3%nat(i,2)
          i3=sys%fcs3%nat(i,3)
          s2=sys%fcs3%cell(i,1)
          s3=sys%fcs3%cell(i,2)

          mass1=sys%mass(sys%kind((2+i1)/3))
          mass2=sys%mass(sys%kind((2+i2)/3))
          mass3=sys%mass(sys%kind((2+i3)/3))

          A=DOT_PRODUCT(sys%rcell(s2,:),this%list(q2)%k(:))
          B=DOT_PRODUCT(sys%rcell(s3,:),this%list(q3)%k(:))

          coeff=0.5291772/dsqrt(mass1*1822.89*this%list(q1)%freq(j1)/219474.6313702)
          coeff=coeff*0.5291772/dsqrt(mass2*1822.89*this%list(q2)%freq(j2)/219474.6313702)
          coeff=coeff*0.5291772/dsqrt(mass3*1822.89*this%list(q3)%freq(j3)/219474.6313702)

          V3=V3+coeff*sys%fcs3%val(i) &
               *this%list(q1)%hess(i1,j1)*this%list(q2)%hess(i2,j2)*this%list(q3)%hess(i3,j3)  &
               *exp(CMPLX(0.0d0,2*pi*(A+B),8))

         enddo

         V3=V3*219474.6313702/(0.5291772)**3

        return
        end function get_V3

        subroutine calc_linewidth_scf(this,sys,Dvet,Fvet) 
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: kp,kp1,kp2,kp3,i,i1,i2,l,j,v,s,jj
        double precision       :: ener0,V3sq,bose_fact,q2(3),width_sum
        double precision       :: Lamb1,Lamb2
        double precision, allocatable :: Dvet(:,:,:,:),Fvet(:,:,:,:)
        complex(8)             :: V3
        logical                :: plus,minus
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)       


         Dvet=0.0d0

        ! Calculation of V3 coefficients          
          
          kp=this%k_start
          do v=1,this%nloc
           do kp1=1,this%ntot
            do kp2=1,this%ntot

             plus=.false.
             minus=.false.
                
             q2(1)=this%list(kp1)%k(1)+this%list(kp2)%k(1)+this%list(kp)%k(1)
             q2(2)=this%list(kp1)%k(2)+this%list(kp2)%k(2)+this%list(kp)%k(2)
             q2(3)=this%list(kp1)%k(3)+this%list(kp2)%k(3)+this%list(kp)%k(3)

             if (abs(q2(1)-nint(q2(1))).lt.1.0D-6 .and.  &
                 abs(q2(2)-nint(q2(2))).lt.1.0D-6 .and.  &
                 abs(q2(3)-nint(q2(3))).lt.1.0D-6 ) then
                plus=.true.
                minus=.true.
             else
                cycle
             endif

!             q2(1)=-this%list(kp1)%k(1)-this%list(kp2)%k(1)+this%list(kp)%k(1)
!             q2(2)=-this%list(kp1)%k(2)-this%list(kp2)%k(2)+this%list(kp)%k(2)
!             q2(3)=-this%list(kp1)%k(3)-this%list(kp2)%k(3)+this%list(kp)%k(3)

!             if (abs(q2(1)-nint(q2(1))).lt.1.0D-6 .and.  &
!                 abs(q2(2)-nint(q2(2))).lt.1.0D-6 .and.  &
!                 abs(q2(3)-nint(q2(3))).lt.1.0D-6 ) plus=.true.

!             q2(1)=this%list(kp1)%k(1)-this%list(kp2)%k(1)+this%list(kp)%k(1)
!             q2(2)=this%list(kp1)%k(2)-this%list(kp2)%k(2)+this%list(kp)%k(2)
!             q2(3)=this%list(kp1)%k(3)-this%list(kp2)%k(3)+this%list(kp)%k(3)

!             if (abs(q2(1)-nint(q2(1))).lt.1.0D-6 .and.  &
!                 abs(q2(2)-nint(q2(2))).lt.1.0D-6 .and.  &
!                 abs(q2(3)-nint(q2(3))).lt.1.0D-6 ) minus=.true.
               
!              if( plus .or. minus) then

              do  i=1,3*sys%nats
              if (this%list(kp)%freq(i).lt.1.0d-6) cycle
              if (kp.eq.1 .and. i.le.3) cycle
               do  i1=1,3*sys%nats
               if (this%list(kp1)%freq(i1).lt.1.0d-6) cycle
               if (kp1.eq.1 .and. i1.le.3) cycle
                do  i2=1,3*sys%nats
                if (this%list(kp2)%freq(i2).lt.1.0d-6) cycle
                if (kp2.eq.1 .and. i2.le.3) cycle


                 Lamb1=this%list(kp1)%freq(i1)/this%list(kp)%freq(i)
                 Lamb2=this%list(kp2)%freq(i2)/this%list(kp)%freq(i)

                 ener0=-this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                   V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                   V3sq=CONJG(V3)*V3
                   do l=1,ntemps
                    bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1
                    do s=1,3
                     Dvet(kp,i,s,l)=Dvet(kp,i,s,l)+V3sq*bose_fact*delta(ener0,smear)* & 
                        (Lamb2*Fvet(kp2,i2,s,l)+Lamb1*Fvet(kp1,i1,s,l))
                     enddo
                   enddo
                 endif

                 ener0=this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1
                   do s=1,3
                    Dvet(kp,i,s,l)=Dvet(kp,i,s,l)+V3sq*bose_fact*delta(ener0,smear)* & 
                        (Lamb2*Fvet(kp2,i2,s,l)+Lamb1*Fvet(kp1,i1,s,l))
                   enddo
                  enddo
                 endif

                 ener0=this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then 
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=bose(temp(l),this%list(kp2)%freq(i2))-bose(temp(l),this%list(kp1)%freq(i1))
                   do s=1,3
                    Dvet(kp,i,s,l)=Dvet(kp,i,s,l)+V3sq*bose_fact*delta(ener0,smear)* & 
                        (-Lamb2*Fvet(kp2,i2,s,l)+Lamb1*Fvet(kp1,i1,s,l))
                   enddo
                  enddo
                 endif

                 ener0=-this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then 
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=-bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))
                   do s=1,3
                    Dvet(kp,i,s,l)=Dvet(kp,i,s,l)+V3sq*bose_fact*delta(ener0,smear)* & 
                        (Lamb2*Fvet(kp2,i2,s,l)-Lamb1*Fvet(kp1,i1,s,l))
                   enddo
                  enddo
                 endif

                enddo   ! i2 bands                
               enddo   ! i1 bands
              enddo   ! i bands

!             endif ! q2 condition

            enddo ! kp2
           enddo ! kp1
           kp=kp+1
          enddo ! kp


        !  reduction function to cumulate results from all nodes
        !

         do i=1,this%ntot
          do j=1,ntemps
           do jj=1,size(this%list(i)%width,1)
            do s=1,3
             width_sum=0.0d0
             call mpi_allreduce(Dvet(i,jj,s,j),width_sum,1, &
                  mpi_double_precision,MPI_SUM,mpi_comm_world,err)    
             Dvet(i,jj,s,j)=width_sum
            enddo
           enddo
          enddo
         enddo

         Dvet=Dvet*2*pi*pi/(16*this%ntot*hplank)
        !  final units are energy (cm-1)


        return
        end subroutine calc_linewidth_scf

        subroutine calc_kscf(this,sys)
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,k,j,s,t,l,kp
        integer                :: iter,max_iter
        double precision       :: bose_fact,mixing
        double precision, allocatable :: kold(:,:,:),Dvet(:,:,:,:),Fvet(:,:,:,:)

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)


         allocate(this%kscf(3,3,ntemps))
         this%kscf=0.0d0
         do i=1,this%ntot
          allocate(this%list(i)%kscf(sys%nats*3,3,3,ntemps))
          this%list(i)%kscf=0.0d0
         enddo

         allocate(kold(3,3,ntemps))
         allocate(Dvet(this%ntot,sys%nats*3,3,ntemps))
         allocate(Fvet(this%ntot,sys%nats*3,3,ntemps))
         kold=this%krta
         Fvet=0.0d0
         Dvet=0.0d0

         max_iter=1000
         mixing=1.0d0

         ! definire condizione iniziale equivalente a RTA

         do kp=1,this%ntot
          do i=1,3*sys%nats
           do l=1,ntemps
            do s=1,3
             Fvet(kp,i,s,l)=(hplank/(2*pi*this%list(kp)%width(i,l)))*this%list(kp)%vel(i,s)
            enddo
           enddo
          enddo
         enddo


         ! calculate kscf by iters on Fvet

         do iter=1,max_iter

          call calc_linewidth_scf(this,sys,Dvet,Fvet) 

          kp=this%k_start
          do k=1,this%nloc
           this%list(k)%kscf=0.0d0
           do i=1,3*sys%nats
            if (kp.eq.1 .and. i.le.3) cycle
             do l=1,ntemps
             if ( abs(this%list(kp)%width(i,l)) .lt. 1.0E-20) cycle   !!!
              do s=1,3
               Fvet(kp,i,s,l)=(1.0d0-mixing)*Fvet(kp,i,s,l)+   &
                mixing*((hplank/(2*pi*this%list(kp)%width(i,l)))*(this%list(kp)%vel(i,s)+Dvet(kp,i,s,l)))
              enddo
              bose_fact=bose(temp(l),this%list(kp)%freq(i))
              do s=1,3
               do t=1,3
                this%list(kp)%kscf(i,s,t,l)=this%list(kp)%kscf(i,s,t,l)+ &
                Fvet(kp,i,t,l)*this%list(kp)%vel(i,s)*                   &
                bose_fact*(bose_fact+1.0d0)*this%list(kp)%freq(i)**2
               enddo ! t
              enddo ! s
             enddo ! l
            enddo ! i
            do l=1,ntemps
             this%list(kp)%kscf(:,:,:,l)=0.1986d0*this%list(kp)%kscf(:,:,:,l)/    &
                                      (sys%vol*this%ntot*kboltz*(hplank**2)*(temp(l)**2))
            enddo
           kp=kp+1
          enddo ! k=kp


        !  broadcast single k contributions, compute on master node the
        !  total and broadcast back

          do i=1,this%ntot
           do k=1,3*sys%nats
            do s=1,3          
             do t=1,3          
              call mpi_bcast(this%list(i)%kscf(k,s,t,:),ntemps,mpi_double_precision,this%mpi_conn(i),mpi_comm_world,err)
             enddo
            enddo
           enddo
          enddo
 
         if(mpi_id.eq.0)then    

           this%kscf=0.0d0

           do i=1,this%ntot
            do k=1,3*sys%nats
             do s=1,3         
              do t=1,3    
               do l=1,ntemps
                this%kscf(s,t,l)=this%kscf(s,t,l)+this%list(i)%kscf(k,s,t,l)
               enddo
              enddo
             enddo
            enddo
           enddo

          endif

          do s=1,3
           do t=1,3
            call mpi_bcast(this%kscf(s,t,:),ntemps,mpi_double_precision,0,mpi_comm_world,err)
           enddo
          enddo
          
!!!!! condizione di convergenza e exit do

          if( maxval(abs(kold-this%kscf)).lt.1.0D-4) exit

          if(mpi_id.eq.0)then 

           write(*,*)   '   Iter=',iter,' Conv_thr=',maxval(abs(kold-this%kscf))
           write(*,*)   '   Total Full-BTE Thermal Conductivity:'

           do k=1,ntemps
            write(*,*) '            ',temp(k),this%kscf(:,:,k)
           enddo

          endif
          
          kold=this%kscf

         enddo


        !  print results

         if(mpi_id.eq.0)then 

          write(*,*)   '   Final Full-BTE Thermal Conductivity:'
          do k=1,ntemps
           write(*,*) '            ',temp(k),this%kscf(:,:,k)
          enddo

          open(11,file='thermal_cond_Full_BTE.dat')
          do l=1,ntemps
           write(11,*) temp(l),this%kscf(:,:,l)
          enddo
          close(11)

         endif


        !  deallocate stuff

         deallocate(kold)
         deallocate(Fvet)
         deallocate(Dvet)

        return
        end subroutine calc_kscf

        subroutine calc_krta(this,sys)
        use mpi
        use mpi_utils
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,k,j,s,t,l,kp
        double precision       :: bose_fact

         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

         allocate(this%krta(3,3,ntemps))
         this%krta=0.0d0
         do i=1,this%ntot
          allocate(this%list(i)%krta(sys%nats*3,3,3,ntemps))
          this%list(i)%krta=0.0d0
         enddo

         kp=this%k_start
         do k=1,this%nloc
          do i=1,3*sys%nats
           if (kp.eq.1 .and. i.le.3) cycle
            do l=1,ntemps
            if ( abs(this%list(kp)%width(i,l)) .lt. 1.0E-6) cycle
             bose_fact=bose(temp(l),this%list(kp)%freq(i))
             do s=1,3
              do t=1,3
               this%list(kp)%krta(i,s,t,l)=this%list(kp)%krta(i,s,t,l)+ &
              (hplank/(2*pi*this%list(kp)%width(i,l)))*                 &
               this%list(kp)%vel(i,s)*this%list(kp)%vel(i,t)*           &
               bose_fact*(bose_fact+1.0d0)*this%list(kp)%freq(i)**2
              enddo
             enddo
            enddo
           enddo
           do l=1,ntemps
            this%list(kp)%krta(:,:,:,l)=0.1986d0*this%list(kp)%krta(:,:,:,l)/    &
                                     (sys%vol*this%ntot*kboltz*(hplank**2)*(temp(l)**2))
          enddo
          kp=kp+1
         enddo


        !  broadcast single k contributions, compute on master node the
        !  total and broadcast back

         do i=1,this%ntot
          do k=1,3*sys%nats
           do s=1,3          
            do t=1,3          
             call mpi_bcast(this%list(i)%krta(k,s,t,:),ntemps,mpi_double_precision,this%mpi_conn(i),mpi_comm_world,err)
            enddo
           enddo
          enddo
         enddo

         if(mpi_id.eq.0)then
          do i=1,this%ntot
           do k=1,3*sys%nats
            do s=1,3         
             do t=1,3    
              do l=1,ntemps
               this%krta(s,t,l)=this%krta(s,t,l)+this%list(i)%krta(k,s,t,l)
              enddo
             enddo
            enddo
           enddo
          enddo
         endif

         do s=1,3
          do t=1,3
           call mpi_bcast(this%krta(s,t,:),ntemps,mpi_double_precision,0,mpi_comm_world,err)
          enddo
         enddo
                     

        !  print results

         if(mpi_id.eq.0)then 

!          write(*,*)   '   Partial RTA Thermal Conductivities:'
!          do i=1,this%ntot
!            write(*,*) '      ',i,this%list(i)%k(:)
!           do k=1,ntemps
!            write(*,*) '             ',temp(k),(this%list(i)%krta(j,:,:,k),j=1,3*sys%nats)
!           enddo
!          enddo
!          write(*,*) 
!          write(*,*) 

          write(*,*)   '   Total RTA Thermal Conductivity:'
          do k=1,ntemps
           write(*,*) '            ',temp(k),this%krta(:,:,k)
          enddo

          open(11,file='thermal_cond_RTA.dat')
          do l=1,ntemps
           write(11,*) temp(l),this%krta(:,:,l)
          enddo
          close(11)

         endif

        return
        end subroutine calc_krta

        subroutine calc_linewidth_sp(this,sys,kp,i)                
        use mpi
        use mpi_utils
        use lists_class
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        type(list)             :: list1,list2
        integer                :: nplet_loc,nplet_start,jj,nscatter
        integer                :: nplets,nplets_loc,nk,nstart,nloc
        integer                :: kp,kp1,kp2,kp3,i,i1,i2,l,j,v
        double precision       :: ener0,V3sq,bose_fact,q2(3),width_sum
        double precision       :: vel_tmp,coeff
        complex(8)             :: V3
        integer, allocatable   :: proc_grid(:),klist(:,:)
        integer, pointer       :: bho
           
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

        ! Calculation of V3 coefficients
 
        ! allocate width accordint to sys and ntemps

          if(.not. allocated(this%list(kp)%width))then
           allocate(this%list(kp)%width(sys%nats*3,ntemps))
          endif
          this%list(kp)%width(i,:)=0.0d0
          if(.not. allocated(this%list(kp)%hess))then
           call this%list(kp)%diagD(sys)
          endif

           call list1%init() 
           call list2%init()

           do kp1=1,this%ntot
            do kp2=kp1,this%ntot
                    
             q2(1)=this%list(kp1)%k(1)+this%list(kp2)%k(1)+this%list(kp)%k(1)
             q2(2)=this%list(kp1)%k(2)+this%list(kp2)%k(2)+this%list(kp)%k(2)
             q2(3)=this%list(kp1)%k(3)+this%list(kp2)%k(3)+this%list(kp)%k(3)

             if (abs(q2(1)-nint(q2(1))).lt.1.0D-6 .and.  &
                 abs(q2(2)-nint(q2(2))).lt.1.0D-6 .and.  &
                 abs(q2(3)-nint(q2(3))).lt.1.0D-6 ) then 
             else
              cycle
             endif

             call list1%add_node(kp1)
             call list2%add_node(kp2)
             
            enddo
           enddo      
           
           if(allocated(klist)) deallocate(klist)
           allocate(klist(list1%nelem,2))

           call list1%reboot()
           call list2%reboot()

           do v=1,list1%nelem           
            call list1%rd_val(klist(v,1))
            call list2%rd_val(klist(v,2))
            call list1%skip()
            call list2%skip()
           enddo

           call list1%delete()
           call list2%delete()

           if(allocated(proc_grid)) deallocate(proc_grid)
           call mpi_dist_nprocess(size(klist,1),nloc,nstart,proc_grid,mpi_comm_world)
           
           nk=nstart
           do v=1,nloc
             
             kp1=klist(nk,1)
             kp2=klist(nk,2)
                   
             if(kp1.eq.kp2) coeff=1.0d0
             if(kp1.ne.kp2) coeff=2.0d0

             if(kp1.ne.kp) call this%list(kp1)%diagD(sys) 
             if(kp2.ne.kp .and. kp2.ne.kp1) call this%list(kp2)%diagD(sys) 
                          
              do  i1=1,3*sys%nats
              if (this%list(kp1)%freq(i1).lt.1.0d-6) cycle
              if (kp1.eq.1 .and. i1.le.3) cycle
               do  i2=1,3*sys%nats
               if (this%list(kp2)%freq(i2).lt.1.0d-6) cycle
               if (kp2.eq.1 .and. i2.le.3) cycle

                ener0=-this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                if(abs(ener0).lt.smear*2)then
!                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
!                  V3sq=CONJG(V3)*V3
                 V3sq=1.0d0
                  do l=1,ntemps
                   bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1.0d0
                   this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)*coeff
                  enddo
                endif

                ener0=this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                if(abs(ener0).lt.smear*2)then
!                 V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
!                 V3sq=CONJG(V3)*V3
                 V3sq=1.0d0
                 do l=1,ntemps
                  bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1.0d0
                  this%list(kp)%width(i,l)=this%list(kp)%width(i,l)-V3sq*bose_fact*delta(ener0,smear)*coeff
                 enddo
                endif

                ener0=this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                if(abs(ener0).lt.smear*2)then
                 if(this%list(kp2)%freq(i2).gt.5*temp(size(temp)))cycle
                 if(this%list(kp1)%freq(i1).gt.5*temp(size(temp)))cycle
!                 V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
!                 V3sq=CONJG(V3)*V3
                 V3sq=1.0d0
                 do l=1,ntemps
                  bose_fact=bose(temp(l),this%list(kp2)%freq(i2))-bose(temp(l),this%list(kp1)%freq(i1))
                  this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)*coeff
                 enddo
                endif

                ener0=-this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                if(abs(ener0).lt.smear*2)then
                 if(this%list(kp2)%freq(i2).gt.5*temp(size(temp)))cycle
                 if(this%list(kp1)%freq(i1).gt.5*temp(size(temp)))cycle
!                 V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
!                 V3sq=CONJG(V3)*V3
                 V3sq=1.0d0
                 do l=1,ntemps
                  bose_fact=-bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))
                  this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)*coeff
                 enddo
                endif

               enddo   ! i2 bands         
              enddo   ! i1 bands
             
             if(allocated(this%list(kp1)%hess)) deallocate(this%list(kp1)%hess)
             if(allocated(this%list(kp2)%hess)) deallocate(this%list(kp2)%hess)

           nk=nk+1
          enddo ! kp1

          if(allocated(this%list(kp)%hess)) deallocate(this%list(kp)%hess)

          do j=1,ntemps
           width_sum=0.0d0
           call mpi_allreduce(this%list(kp)%width(i,j),width_sum,1,&
                mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,err)    
           this%list(kp)%width(i,j)=width_sum
          enddo

          this%list(kp)%width(i,:)=pi/(16*this%ntot)*this%list(kp)%width(i,:)
        !  final units are energy (cm-1)

        return
        end subroutine calc_linewidth_sp


        subroutine calc_linewidth(this,sys)                
        use mpi
        use mpi_utils
        use lists_class
        use lattice_class
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: nplet_loc,nplet_start,jj
        integer                :: nplets,nplets_loc
        integer                :: kp,kp1,kp2,kp3,i,i1,i2,l,j,v
        double precision       :: ener0,V3sq,bose_fact,q2(3),width_sum
        double precision       :: vel_tmp
        complex(8)             :: V3
        integer, allocatable   :: proc_grid(:)
        integer, pointer       :: bho
        logical                :: plus,minus
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)       

         if(.not.allocated(sys%fcs3))then
          if(mpi_id.eq.0) write(*,*) 'FC3 not allocated'
          return
         endif

        ! allocate width accordint to sys and ntemps

         do i=1,this%ntot
          allocate(this%list(i)%width(sys%nats*3,ntemps))
          this%list(i)%width=0.0d0
         enddo

        ! Calculation of V3 coefficients          
          
          kp=this%k_start
          do v=1,this%nloc
           do kp1=1,this%ntot
            do kp2=1,this%ntot

             q2(1)=this%list(kp1)%k(1)+this%list(kp2)%k(1)+this%list(kp)%k(1)
             q2(2)=this%list(kp1)%k(2)+this%list(kp2)%k(2)+this%list(kp)%k(2)
             q2(3)=this%list(kp1)%k(3)+this%list(kp2)%k(3)+this%list(kp)%k(3)

             if (abs(q2(1)-nint(q2(1))).lt.1.0D-6 .and.  &
                 abs(q2(2)-nint(q2(2))).lt.1.0D-6 .and.  &
                 abs(q2(3)-nint(q2(3))).lt.1.0D-6 ) then 
             else
              cycle
             endif
                                          
              do  i=1,3*sys%nats
              if (this%list(kp)%freq(i).lt.1.0d-6) cycle
              if (kp.eq.1 .and. i.le.3) cycle
               do  i1=1,3*sys%nats
               if (this%list(kp1)%freq(i1).lt.1.0d-6) cycle
               if (kp1.eq.1 .and. i1.le.3) cycle
                do  i2=1,3*sys%nats
                if (this%list(kp2)%freq(i2).lt.1.0d-6) cycle
                if (kp2.eq.1 .and. i2.le.3) cycle

                 ener0=-this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                   V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                   V3sq=CONJG(V3)*V3
                   do l=1,ntemps
                    bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1.0d0
                    this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)
                   enddo
                 endif

                 ener0=this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))+1.0d0
                   this%list(kp)%width(i,l)=this%list(kp)%width(i,l)-V3sq*bose_fact*delta(ener0,smear)
                  enddo
                 endif

                 ener0=this%list(kp2)%freq(i2)-this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=bose(temp(l),this%list(kp2)%freq(i2))-bose(temp(l),this%list(kp1)%freq(i1))
                   this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)
                  enddo
                 endif

                 ener0=-this%list(kp2)%freq(i2)+this%list(kp1)%freq(i1)+this%list(kp)%freq(i)
                 if(abs(ener0).lt.smear*3)then
                  V3=this%get_V3(sys,kp,kp1,kp2,i,i1,i2)
                  V3sq=CONJG(V3)*V3
                  do l=1,ntemps
                   bose_fact=-bose(temp(l),this%list(kp2)%freq(i2))+bose(temp(l),this%list(kp1)%freq(i1))
                   this%list(kp)%width(i,l)=this%list(kp)%width(i,l)+V3sq*bose_fact*delta(ener0,smear)
                  enddo
                 endif

                enddo   ! i2 bands                
               enddo   ! i1 bands
              enddo   ! i bands

            enddo ! kp2
           enddo ! kp1
           kp=kp+1
          enddo ! kp


        !  reduction function to cumulate results from all nodes
        !

         do i=1,this%ntot
          do j=1,ntemps
           do jj=1,size(this%list(i)%width,1)
            width_sum=0.0d0
            call mpi_allreduce(this%list(i)%width(jj,j),width_sum,1,&
                 mpi_double_precision,MPI_SUM,mpi_comm_world,err)    
            this%list(i)%width(jj,j)=width_sum
           enddo
          enddo
         enddo

         do j=1,this%ntot
          do jj=1,(sys%nats*3)
           this%list(j)%width(jj,:)=pi/(16*this%ntot)*this%list(j)%width(jj,:)
        !  final units are energy (cm-1)
          enddo
         enddo

         if(boundary_scattering)then
          do j=1,this%ntot
           do jj=1,(sys%nats*3)
            if (this%list(j)%freq(jj).lt.1.0d-6) cycle
            vel_tmp=sqrt(this%list(j)%vel(jj,1)**2+this%list(j)%vel(jj,2)**2)
            do i=1,ntemps
             this%list(j)%width(jj,i)=this%list(j)%width(jj,i)+(2.0d0*pi*vel_tmp/Length)
            enddo
!         !  final units are energy (cm-1)
           enddo
          enddo
         endif

        ! print results

         if(mpi_id.eq.0)then 
           write(*,*) '   Phonons Line-width:'
          do i=1,this%ntot
            write(*,*) '      ',i,this%list(i)%k(:)
           do j=1,size(this%list(i)%width,2)
            write(*,*) '          ',temp(j),this%list(i)%width(:,j)
           enddo
          enddo
         endif

         if(mpi_id.eq.0)then
          open(11,file='tau.dat')
          write(11,*) 'Temperature ','Bands '
          do kp=1,this%ntot
            write(11,*) '#### k:',this%list(kp)%k(:)
           do l=1,ntemps
            write(11,*) temp(l),hplank/2.0d0/pi/this%list(kp)%width(:,l)
           enddo
          enddo
          close(11)
         endif
        
        ! deallocate stuff

        return
        end subroutine calc_linewidth

        end module phonons_class
