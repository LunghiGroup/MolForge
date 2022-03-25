        module phonons_class
        implicit none

         logical                        :: read_fc3=.false.
         logical                        :: do_brillouin_path=.false.
         logical                        :: do_brillouin_mesh=.false.
         logical                        :: boundary_scattering=.false.
         double precision, allocatable  :: temp(:)
         integer                        :: ntemps
         integer                        :: type_smear=1 ! 1=Gaussian 0=Lorentzian
         double precision               :: smear
         double precision               :: Length


        type kpoint
         double precision :: k(3)
         double precision :: weight
         double precision, allocatable       :: freq(:)
         double precision, allocatable       :: vel(:,:)
         double precision, allocatable       :: width(:,:)
         double precision, allocatable       :: ir(:,:)
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
         logical                        :: effective_lt=.false.
         integer                        :: nx
         integer                        :: ny
         integer                        :: nz
         integer                        :: ntot
         integer                        :: nloc
         integer                        :: k_start
         integer, allocatable           :: mpi_conn(:)
         type(dist1D)                   :: dos1p
         type(dist1D)                   :: dos2p
         type(dist1D), allocatable      :: pdos1p(:)
         contains 
         procedure        :: generate_mesh
         procedure        :: generate_path
         procedure        :: brillouin_bcast
         procedure        :: calc_disps
         procedure        :: calc_IR
         procedure        :: calc_bands
         procedure        :: get_cart_disp
         procedure        :: calc_vel
         procedure        :: remap_hess
         procedure        :: mpi_dist_kpoints
         procedure        :: calc_dos1p
         procedure        :: calc_dos2p
         procedure        :: smooth_bands
         procedure        :: print_cart_disp
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
         call mpi_bcast(type_smear,1,mpi_integer,0,mpi_comm_world,err)
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

           if(this%list(j)%freq(i).le.0.0d0) cycle
           if(j.eq.1 .and. i.le.3) cycle

           k=NINT(this%list(j)%freq(i)/this%dos1p%step)
           l=NINT(this%dos1p%sigma/this%dos1p%step)

           do v=-3*l,3*l
            if((v+k).gt.0 .and. (v+k).lt.this%dos1p%nsteps)then
             this%dos1p%dist(k+v)=this%dos1p%dist(k+v)+deltaG(DBLE(v),DBLE(l))
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
            this%dos2p%dist(k+v)=this%dos2p%dist(k+v)+deltaG(DBLE(v),DBLE(l))*&
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
!          open(13,file='disp_modes.dat')
          dist=0.0d0
          do i=1,size(this%path)
           if(i.gt.1)then
            dist=dist+sqrt(sys%dist_rec(this%path(i)%k(:)-this%path(i-1)%k(:)))
           endif
           write(12,*) this%path(i)%k(:),dist,(this%path(i)%freq(s),s=1,size(this%path(i)%freq))
!           do s=1,sys%nats*3
!            write(13,*) 'K:',this%path(i)%k(:),'Mode:',s
!            do j=1,sys%nats*3
!             write(13,*) this%path(i)%k(:),dist,dble(this%path(i)%hess(j,s)),aimag(this%path(i)%hess(j,s))
!            enddo
!           enddo
          enddo
          close(12)
!          close(13)
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
        use units_parms
        use lattice_class
        use atoms_class
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: i,k,j,rest,l,s,t,ntemps
        double precision       :: coeff,mass1,Cv_tmp,temp,ratio
        double precision, allocatable :: Cv(:)

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

         if(mpi_id.eq.0)then 
          do i=1,1
           write(*,*) '         ','Gamma-point Phonons'
           do j=1,sys%nats*3
            write(*,*) '          ',j,this%list(i)%freq(j)
           enddo
          enddo
         endif

!         call this%smooth_bands()        
         this%list(1)%freq(1:3)=0.0d0           

         ! calc CV

         ntemps=300
         allocate(Cv(ntemps))
         Cv=0.0d0

         do t=1,ntemps

          temp=t

          j=this%k_start
          Cv_tmp=0.0d0

          do i=1,this%nloc
           do l=1,size(this%list(j)%freq)
            if(i.eq.1 .and. l.le.3) cycle
            ratio=this%list(j)%freq(l)/temp/kboltz
            Cv_tmp=Cv_tmp+kboltz*(ratio**2)*exp(ratio)/(exp(ratio)-1)**2
           enddo
           j=j+1
          enddo

          call mpi_allreduce(Cv_tmp,Cv(t),1,&
               mpi_double_precision,mpi_sum,mpi_comm_world,err)
 
         enddo

         Cv=Cv*0.0119627

         if(mpi_id.eq.0)then     
          open(133,file='Cv.dat')
          write(133,*) 0.0d0,0.0d0
          do i=1,size(Cv)
           write(133,*) i,Cv(i)  
          enddo
          close(133)
         endif

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

         ! broadcast calculated velocities

         do i=1,this%ntot
          do v=1,3
           call mpi_bcast(this%list(i)%vel(:,v),size(this%list(i)%vel,1),mpi_double_precision,this%mpi_conn(i),mpi_comm_world,err)
          enddo
         enddo
        
         if(mpi_id.eq.0)then 
           write(*,*) '     Group Velocities:'
          do i=1,this%ntot
            write(*,*) '     ',i,this%list(i)%k(:)
           do j=1,sys%nats*3
            write(*,*) '     ',j,this%list(i)%vel(j,:)
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

         call MPI_BCAST(this%effective_lt,1,mpi_logical,0,mpi_comm_world,err) 
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
        double precision              :: A,pi,B1(3),B2(3),B3(3),B4,mat(3,3)
        double complex                :: coeff
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

         if (sys%born_charges) then

          coeff=(0.0d0,0.0d0)
          do l=1,sys%nx*sys%ny*sys%nz      
           A=DOT_PRODUCT(sys%rcell(l,:),ph%k(:))
           coeff=coeff+exp(CMPLX(0.0d0,2*pi*A,8))
          enddo

          do i1=1,sys%nats
          do i2=1,sys%nats
           do s1=1,3
           do s2=1,3

            v1=(i1-1)*3+s1
            v2=(i2-1)*3+s2
             
             B1=0.0d0
             B2=0.0d0
             B3=0.0d0
             B4=0.0d0

             mat=sys%Zeff(sys%kind(i1),:,:)
             B1=matmul(mat,ph%k)
             mat=sys%Zeff(sys%kind(i2),:,:)
             B2=matmul(mat,ph%k)
             B3=matmul(sys%eps,ph%k)
             B4=dot_product(ph%k,B3)
             
                     ph%hess(v1,v2)=ph%hess(v1,v2)+&
                     4*pi*B1(s1)*B2(s2)/B4/sys%vol/sys%ntot&
                     *0.5291772109*0.5291772109*0.5291772109*coeff

           enddo
           enddo
          enddo
          enddo

         endif


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
!!!! 1822.89 turn the amu units of mass in a.u. (electron mass). Ener(i) alla fine è a.u.
!           do t=1,N  
!!!! the eigenvectors are transformed in cartesian displacements (angstrom) associated with the unit-less normal mode
!            hess(t,i)=hess(t,i)*0.5291772/dsqrt(mass2(t)*1822.89*ener(i))   
!!!! 1822.89 turn mass from amu to a.u. 
!!!! 0.5291 turn bohr in Ang
!           enddo

           ener(i)=sqrt(ener(i)/1822.89)    
           ener(i)=ener(i)*219474.6313702

!!!! 1/(c[a.u.]*bohr2cm*2pi) gives ener as frequency in cm-1
!!!! o in modo equivalente
!!!! 1/( 3E8*1E2*0.036752*4.13566E-15) = 1/(c(m/s)*(m2cm)*(eV2Eh)*h(eVs)) turn ener(i) (t^{-1} [a.u.]) in cm^{-1}

          else
                
           ener(i)=-sqrt(abs(ener(i))/1822.89)    
           ener(i)=ener(i)*219474.6313702

          endif
         enddo
       
        return
        end subroutine diaghess


        subroutine calc_IR(this,sys)
        use mpi
        use atoms_class
        implicit none
        class(brillouin)       :: this
        type(atoms_group)      :: sys
        integer                :: mpi_nproc,mpi_id,err
        integer                :: i,k,j,v,l
        double precision       :: mass1,coeff

         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)
         
         if(mpi_id.eq.0)then
          write(*,*) 'IR Intensities'
          write(*,*) 'warning! IR Intensities only implemented at Gamma',&
                     ' point'
         endif

         j=this%k_start

         do i=1,this%nloc

          allocate(this%list(j)%ir(3*sys%nats,3))
          this%list(j)%ir=0.0d0

          do k=1,size(this%list(j)%freq)

           if(j.eq.1 .and. k.le.3)cycle

           do v=1,sys%nats*3

            mass1=sys%mass(sys%kind((2+v)/3))
            coeff=0.5291772/dsqrt(mass1*1822.89*this%list(j)%freq(k)/219474.6313702)           

            this%list(j)%ir(k,1)=this%list(j)%ir(k,1)+&
                        coeff*this%list(j)%hess(v,k)*sys%dipole(v,1)
            this%list(j)%ir(k,2)=this%list(j)%ir(k,2)+&
                        coeff*this%list(j)%hess(v,k)*sys%dipole(v,2)
            this%list(j)%ir(k,3)=this%list(j)%ir(k,3)+&
                        coeff*this%list(j)%hess(v,k)*sys%dipole(v,3)

           enddo
          enddo 

          j=j+1
         enddo
        

         ! broadcast calculated bands 

         do i=1,this%ntot
          if(.not.allocated(this%list(i)%ir))then
           allocate(this%list(i)%ir(sys%nats*3,3))
          endif
          call mpi_bcast(this%list(i)%ir(:,1),size(this%list(i)%freq),mpi_double_precision,this%mpi_conn(i),MPI_COMM_WORLD,err)
          call mpi_bcast(this%list(i)%ir(:,2),size(this%list(i)%freq),mpi_double_precision,this%mpi_conn(i),MPI_COMM_WORLD,err)
          call mpi_bcast(this%list(i)%ir(:,3),size(this%list(i)%freq),mpi_double_precision,this%mpi_conn(i),MPI_COMM_WORLD,err)
         enddo

         if(mpi_id.eq.0)then 
          do i=1,this%ntot
            write(*,*) '      ',i,this%list(i)%k(:)
           do j=1,sys%nats*3
            write(*,*) '          ',j,this%list(i)%freq(j),this%list(i)%ir(j,:)
           enddo
          enddo
         endif
        
        return
        end subroutine calc_IR

        subroutine print_cart_disp(this,sys,kp,bn,step,nsteps)
        use mpi
        use atoms_class
        use units_parms
        implicit none
        class(brillouin)              :: this
        type(atoms_group)             :: sys
        integer                       :: kp,bn,j,jj,nsteps,i,s
        integer                       :: mpi_nproc,mpi_id,err
        double precision, allocatable :: disp(:,:)
        double precision              :: step

         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)         
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         if(mpi_id.eq.0)then

          allocate(disp(sys%nats,3))

          j=1
          do jj=1,sys%nats
           disp(jj,1)=bohr2ang/dsqrt(sys%mass(sys%kind(jj))*1822.89/219474.6313702)*&
                          dble(this%list(kp)%hess(j,bn))/sqrt(this%list(kp)%freq(bn))
           disp(jj,2)=bohr2ang/dsqrt(sys%mass(sys%kind(jj))*1822.89/219474.6313702)*&
                          dble(this%list(kp)%hess(j+1,bn))/sqrt(this%list(kp)%freq(bn))
           disp(jj,3)=bohr2ang/dsqrt(sys%mass(sys%kind(jj))*1822.89/219474.6313702)*&
                          dble(this%list(kp)%hess(j+2,bn))/sqrt(this%list(kp)%freq(bn))
                         
           j=j+3
          enddo

          write(*,*) "Writing Cartesian Displacements"
          open(13,file='Cart_disp.xyz')          

          do j=-nsteps,nsteps

           write(13,*) sys%nats
           write(13,*) 

           do i=1,sys%nats
            write(13,*) trim(sys%label(sys%kind(i))),(sys%x(i,s)+j*step*disp(i,s),s=1,3)
           enddo

          enddo

          do j=nsteps-1,-nsteps,-1

           write(13,*) sys%nats
           write(13,*) 

           do i=1,sys%nats
            write(13,*) trim(sys%label(sys%kind(i))),(sys%x(i,s)+j*step*disp(i,s),s=1,3)
           enddo

          enddo

          close(13)

         endif

        return
        end subroutine print_cart_disp


        end module phonons_class
