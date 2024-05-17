        module quantum_systems_class
        use mpi
        use blacs_utils
        use general_types_class
        use atoms_class
        use spins_dist_rs_class         
        use phonons_class
        use liuville_class

        implicit none
                
        type, extends(liuville_space) :: open_quantum_system
         ! phonon bath
         type(brillouin)                         :: phonons
         double precision                        :: temp
         double precision                        :: smear
         integer                                 :: type_smear=1 ! 1=Gaussian 0=Lorentzian
         type(dist_cmplx_mat)                    :: Vq
         ! lattice coupling
         type(atoms_group)                       :: lattice
         type(dist_cmplx_mat), allocatable       :: Vx(:)
         integer, allocatable                    :: map_s2a(:,:)
         integer                                 :: norder
         integer                                 :: nderiv         
         contains              
         procedure     ::  make_R21_lindbladian
         procedure     ::  make_R22_lindbladian
         procedure     ::  make_R41_lindbladian
         procedure     ::  X2Q
         procedure     ::  XY2QQ
        end type open_quantum_system

        type, extends(open_quantum_system) :: spin_quantum_system
         type(spins_group)                     :: spins
         double precision, allocatable         :: basis(:,:)
         integer, allocatable                  :: Hnodes(:,:)
         type(mat_cmplx), allocatable          :: SHrep(:,:)
         contains
         procedure   ::  make_spin_basis
         procedure   ::  make_SH_rep         
         procedure   ::  make_Hmat_nodes
         procedure   ::  get_Hij
         procedure   ::  make_Hmat
         procedure   ::  make_Smat
         procedure   ::  make_tinv
!         procedure   ::  make_Vx
!         procedure   ::  make_rot
         procedure   ::  set_dipolar
        end type spin_quantum_system

        contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   TRANSFORM LINEAR INTERACTION V FROM CARTESIAN TO PHONON BASIS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine X2Q(this,kk,bn)
        implicit none
        class(open_quantum_system)        :: this
        integer                           :: i,j,l,v,kk,bn
        double precision                  :: coeff_r,mass1
        double complex                    :: coeff
        double precision, allocatable     :: hess(:)

        allocate(hess(this%lattice%nats*3))

        this%Vq%mat=(0.0d0,0.0d0)

        l=1
        do i=1,this%lattice%nats
         mass1=this%lattice%mass(this%lattice%kind(i))*1822.89        
         coeff_r=0.5291772d0/sqrt(mass1*this%phonons%ntot*this%phonons%list(kk)%freq(bn)/219474.6313702d0)
         do j=1,3
          hess(l)=coeff_r*this%phonons%list(kk)%hess(l,bn)
          l=l+1
         enddo
        enddo

        do i=1,size(this%Vx)
         l=this%map_s2a(i,1)
         v=this%map_s2a(i,2)
         coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(this%phonons%list(kk)%k,this%lattice%rcell(v,:))
         this%Vq%mat=this%Vq%mat+this%Vx(i)%mat*hess(l)*exp(coeff)
        enddo

        deallocate(hess)

       return
       end subroutine X2Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   TRANSFORM QUADRATIC INTERACTION V FROM CARTESIAN TO PHONON BASIS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine XY2QQ(this,kk,bn,kk2,bn2)
        implicit none
        class(open_quantum_system)        :: this
        integer                           :: i,j,l,v,kk,bn,kk2,bn2,l2,v2
        double precision                  :: coeff_r,mass1
        double complex                    :: coeff,coeff2
        double precision, allocatable     :: hess(:),hess2(:)

        allocate(hess(this%lattice%nats*3))

        this%Vq%mat=(0.0d0,0.0d0)

        l=1
        do i=1,this%lattice%nats
         mass1=this%lattice%mass(this%lattice%kind(i))*1822.89        
         coeff_r=0.5291772d0/sqrt(mass1*this%phonons%ntot*this%phonons%list(kk)%freq(bn)/219474.6313702d0)
         do j=1,3
          hess(l)=coeff_r*this%phonons%list(kk)%hess(l,bn)
          l=l+1
         enddo
        enddo

        allocate(hess2(this%lattice%nats*3))

        l=1
        do i=1,this%lattice%nats
         mass1=this%lattice%mass(this%lattice%kind(i))*1822.89        
         coeff_r=0.5291772d0/sqrt(mass1*this%phonons%ntot*this%phonons%list(kk2)%freq(bn2)/219474.6313702d0)
         do j=1,3
          hess2(l)=coeff_r*this%phonons%list(kk2)%hess(l,bn2)
          l=l+1
         enddo
        enddo

        do i=1,size(this%Vx)
         l=this%map_s2a(i,1)
         v=this%map_s2a(i,2)
         l2=this%map_s2a(i,3)
         v2=this%map_s2a(i,4)
         coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(this%phonons%list(kk)%k,this%lattice%rcell(v,:))
         coeff2=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(this%phonons%list(kk2)%k,this%lattice%rcell(v2,:))
         if(l.eq.l2)then
          this%Vq%mat=this%Vq%mat+this%Vx(i)%mat*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)
         else
          this%Vq%mat=this%Vq%mat+this%Vx(i)%mat*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)
          this%Vq%mat=this%Vq%mat+this%Vx(i)%mat*hess(l2)*exp(coeff)*hess2(l)*exp(coeff2)
         endif
        enddo

        deallocate(hess)
        deallocate(hess2)

       return
       end subroutine XY2QQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN R21
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R21_lindbladian(this,min_ener,max_ener)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(open_quantum_system)        :: this
        integer                           :: t1,t2,rate,l,l2,i,j,v
        integer                           :: ph,phx,nstart,nloc,bn
        double precision                  :: min_ener,max_ener,freq
        double complex                    :: valc
        double complex, allocatable       :: Vmat(:,:)
        integer, allocatable              :: proc_grid(:)

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Building the second-order Lindbladian operator with linear coupling'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))
         if (.not.allocated(this%Vq%mat)) call this%Vq%set(this%Hdim,this%Hdim,NB,MB)

         ! run on k points and bn

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(this%phonons%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
          do bn=1,size(this%phonons%list(ph)%freq)

           if (ph.eq.1 .and. bn.le.3) cycle
           freq=this%phonons%list(ph)%freq(bn)
           if (freq.lt.min_ener) cycle
           if (freq.gt.max_ener) cycle           

           if(.not.allocated(this%phonons%list(ph)%hess) ) call this%phonons%list(ph)%diagD(this%lattice)

           if ( allocated(this%Vx) ) call this%X2Q(ph,bn)
         
           Vmat=(0.0d0,0.0d0)

           do l2=1,this%Hdim
            do l=1,this%Hdim
             call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
             Vmat(l2,l)=valc
            enddo
           enddo                              

         ! compute R21

           call this%make_R21(Vmat,this%temp,freq,this%smear,this%type_smear)

          enddo 
          ph=ph+1
         enddo

         do i=1,size(this%R21%mat,1)
          do j=1,size(this%R21%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(this%R21%mat(i,j),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           this%R21%mat(i,j)=valc
          enddo
         enddo

         this%R%mat=this%R%mat+this%R21%mat

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R21_lindbladian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN R22
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R22_lindbladian(this,min_ener,max_ener)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(open_quantum_system)        :: this
        integer                           :: t1,t2,rate,l,l2,i,j,v
        integer                           :: ph,phx,nstart,nloc,bn,ph2,bn2
        double precision                  :: min_ener,max_ener,freq,freq2
        double complex                    :: valc
        double complex, allocatable       :: Vmat(:,:)
        integer, allocatable              :: proc_grid(:)

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Building the second-order Lindbladian operator with quadratic coupling'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))
         if (.not.allocated(this%Vq%mat)) call this%Vq%set(this%Hdim,this%Hdim,NB,MB)

         ! run on k points and bn

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(this%phonons%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
          do bn=1,size(this%phonons%list(ph)%freq)

           if (ph.eq.1 .and. bn.le.3) cycle
           freq=this%phonons%list(ph)%freq(bn)
           if (freq.lt.min_ener) cycle
           if (freq.gt.max_ener) cycle           

           if(.not.allocated(this%phonons%list(ph)%hess) ) call this%phonons%list(ph)%diagD(this%lattice)

           do ph2=1,this%phonons%ntot
            do bn2=1,size(this%phonons%list(ph2)%freq)

             if (ph2.eq.1 .and. bn2.le.3) cycle
             freq2=this%phonons%list(ph2)%freq(bn2)
             if (freq2.lt.min_ener) cycle
             if (freq2.gt.max_ener) cycle           

             if(.not.allocated(this%phonons%list(ph2)%hess) ) call this%phonons%list(ph2)%diagD(this%lattice)

             if ( allocated(this%Vx) ) call XY2QQ(this,ph,bn,ph2,bn2)

             Vmat=(0.0d0,0.0d0)

             do l2=1,this%Hdim
              do l=1,this%Hdim
               call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
               Vmat(l2,l)=valc
              enddo
             enddo                   

             call this%make_R22(Vmat,this%temp,freq,freq2,this%smear,this%smear,this%type_smear)

            enddo ! bn2
           enddo ! ph2

          enddo 
          ph=ph+1
         enddo

         do i=1,size(this%R22%mat,1)
          do j=1,size(this%R22%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(this%R22%mat(i,j),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           this%R22%mat(i,j)=valc
          enddo
         enddo

         this%R%mat=this%R%mat+this%R22%mat

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R22_lindbladian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN R41
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R41_lindbladian(this,min_ener,max_ener,correction)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(open_quantum_system)        :: this
        integer                           :: t1,t2,rate,l,l2,i,j,v
        integer                           :: ph,phx,nstart,nloc,bn,ph2,bn2
        double precision                  :: min_ener,max_ener,freq,freq2
        double complex                    :: valc
        double complex, allocatable       :: Vmat(:,:),V2mat(:,:)
        integer, allocatable              :: proc_grid(:)
        logical                           :: correction

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Building the fourth-order Linbladian operator with linear coupling'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))
         allocate(V2mat(this%Hdim,this%Hdim))
         if (.not.allocated(this%Vq%mat)) call this%Vq%set(this%Hdim,this%Hdim,NB,MB)

         ! run on k points and bn

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(this%phonons%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
          do bn=1,size(this%phonons%list(ph)%freq)

           if (ph.eq.1 .and. bn.le.3) cycle
           freq=this%phonons%list(ph)%freq(bn)
           if (freq.lt.min_ener) cycle
           if (freq.gt.max_ener) cycle           

           if(.not.allocated(this%phonons%list(ph)%hess) ) call this%phonons%list(ph)%diagD(this%lattice)

           if ( allocated(this%Vx) ) call X2Q(this,ph,bn)

           Vmat=(0.0d0,0.0d0)

           do l2=1,this%Hdim
            do l=1,this%Hdim
             call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
             Vmat(l2,l)=valc
            enddo
           enddo                   

           do ph2=1,this%phonons%ntot
            do bn2=1,size(this%phonons%list(ph2)%freq)

             if ( (ph-1)*size(this%phonons%list(ph)%freq)+bn .ge. &
                  (ph2-1)*size(this%phonons%list(ph2)%freq)+bn2 ) cycle 

             if (ph2.eq.1 .and. bn2.le.3) cycle
             freq2=this%phonons%list(ph2)%freq(bn2)
             if (freq2.lt.min_ener) cycle
             if (freq2.gt.max_ener) cycle           

             if(.not.allocated(this%phonons%list(ph2)%hess) ) call this%phonons%list(ph2)%diagD(this%lattice)

             if ( allocated(this%Vx) ) call X2Q(this,ph2,bn2)
         
             V2mat=(0.0d0,0.0d0)

             do l2=1,this%Hdim
              do l=1,this%Hdim
               call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
               V2mat(l2,l)=valc
              enddo
             enddo                   

             call this%make_R41(Vmat,V2mat,this%temp,freq,freq2,this%smear,this%smear,this%type_smear,correction)

            enddo ! bn2
           enddo ! ph2

          enddo 
          ph=ph+1
         enddo

         do i=1,size(this%R41%mat,1)
          do j=1,size(this%R41%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(this%R41%mat(i,j),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           this%R41%mat(i,j)=valc
          enddo
         enddo

         this%R%mat=this%R%mat+this%R41%mat

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R41_lindbladian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE SPIN BASIS SET
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_spin_basis(this)
        implicit none
        class(spin_quantum_system)      :: this
        integer                         :: i,id,si
        integer                         :: t1,t2,rate

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Building Basis Set'
          flush(6)
         endif

         this%Hdim=1
         do i=1,this%spins%nspins
          this%Hdim=this%Hdim*(nint(2*this%spins%spin(this%spins%kind(i)))+1)
         enddo       
         allocate(this%basis(this%Hdim,this%spins%nspins))

         id=1
         si=1
         call build_spin_basis(this,id,si)
         this%basis=this%basis/2.0d0

         if(mpi_id.eq.0)then
          write(*,*) '     Total Hilbert space size: ',this%Hdim
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_spin_basis

        recursive subroutine build_spin_basis(this,id,si)
        implicit none
        class(spin_quantum_system)      :: this
        integer                         :: id,si,sj,i

 !        if ( id.eq.1 .and. si.eq.1 ) then
 !         do j=si+1,this%spins%nspins
 !          basis_tmp(j)=-nint(2*this%spins%spin(this%spins%kind(j)))
 !         enddo
 !        endif

         if(si.lt.this%spins%nspins)then

          do i=-nint(2*this%spins%spin(this%spins%kind(si))),nint(2*this%spins%spin(this%spins%kind(si))),2

           this%basis(id,si)=i

!           do j=si+1,this%spins%nspins
!            this%basis(id,j)=-nint(2*this%spins%spin(this%spins%kind(j)))
!           enddo

           sj=si+1  
           call build_spin_basis(this,id,sj)

          enddo

         else

          do i=-nint(2*this%spins%spin(this%spins%kind(si))),nint(2*this%spins%spin(this%spins%kind(si))),2
           this%basis(id,si)=i
           id=id+1
          enddo

         endif
       
        return
        end subroutine build_spin_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   SET UP THE DIPOLAR NETWORK AMONG THE SPINS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine set_dipolar(this,SH,ex_list)
        use spinham_class
        use mpi_utils
        implicit none
        class(spin_quantum_system)       :: this
        class(SpinHamiltonian)           :: SH
        integer                          :: i,j,s1,s2,ii,jj,v,l,m,ex,is
        integer, allocatable             :: ex_list(:,:)
        integer                          :: celli,cellj,Hdim1,Hdim2,Hdim12
        double precision                 :: dist0(3)
        double precision                 :: spin2(2),spin,psi(2),psi2(2,2)
        integer                          :: t1,t2,rate
        logical                          :: skip

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building Dipolar Network'
          flush(6)
         endif

         SH%nDdip=0

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%spins%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%spins%nspins_pr 
            do cellj=1,this%spins%ntot 

             s1=this%spins%nspins_pr*(celli-1)+ii
             s2=this%spins%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif

             if(s2.le.s1 .or. skip) cycle

             if(this%spins%dist(ii,celli,jj,cellj).le.SH%dipolar_thr)then

              if(this%spins%kind(s1).eq.SH%G(i)%kind .and. this%spins%kind(s2).eq.SH%G(j)%kind )then                 
               SH%nDdip=SH%nDdip+1
              endif

              if(this%spins%kind(s1).eq.SH%G(j)%kind .and. this%spins%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then
               SH%nDdip=SH%nDdip+1
              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         allocate(SH%Ddip(SH%nDdip))

         if(mpi_id.eq.0)  &
         write(*,*) '     Total Number of spin-spin dipolar interactions: ',SH%nDdip

         v=1

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%spins%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%spins%nspins_pr 
            do cellj=1,this%spins%ntot 

             s1=this%spins%nspins_pr*(celli-1)+ii
             s2=this%spins%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif


             if(s2.le.s1 .or. skip) cycle

             if(this%spins%dist(ii,celli,jj,cellj).le.SH%dipolar_thr)then

              dist0=this%spins%dist_vec_pbc(this%spins%x(s1,:),this%spins%x(s2,:))  

              if(this%spins%kind(s1).eq.SH%G(i)%kind .and. this%spins%kind(s2).eq.SH%G(j)%kind )then

               call SH%Ddip(v)%make_D(SH%G(i)%G,SH%G(j)%G, &
                    this%spins%bohr_mag(this%spins%kind(s1)),this%spins%bohr_mag(this%spins%kind(s2)),dist0,&
                    this%spins%dist(ii,celli,jj,cellj))
               SH%Ddip(v)%kind(1)=s1
               SH%Ddip(v)%kind(2)=s2

               v=v+1

              endif

              if(this%spins%kind(s1).eq.SH%G(j)%kind .and. this%spins%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then

               call SH%Ddip(v)%make_D(SH%G(j)%G,SH%G(i)%G,  &
                    this%spins%bohr_mag(this%spins%kind(s1)),this%spins%bohr_mag(this%spins%kind(s2)),dist0,&
                    this%spins%dist(ii,celli,jj,cellj))
               SH%Ddip(v)%kind(1)=s1
               SH%Ddip(v)%kind(2)=s2

               v=v+1

              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine set_dipolar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   COMPUTE MATRIX RAPRESENTATION OF EACH TERM OF THE SPIN HAMILTONIAN
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_SH_rep(this,SH,spin_id,spin_id2)
        use spinham_class
        implicit none
        class(spin_quantum_system) :: this
        class(SpinHamiltonian)     :: SH
        integer                    :: i,j,l2,l,s1,s2,l1
        integer                    :: ii,jj,ii1,ii2,jj1,jj2
        integer                    :: t1,t2,rate,a,b,c,d,spin_id,spin_id2
        double precision           :: spin2(2),psi2(2,2)
        double precision           :: spin1,psi1(2)

         if(allocated(this%SHrep))then
          do i=1,size(this%SHrep,1)
           do j=1,size(this%SHrep,2)
            call this%SHrep(i,j)%delete()
           enddo
          enddo
          deallocate(this%SHrep)
         endif

         ! allocate

         if(.not.allocated(this%SHrep))then
          allocate(this%SHrep(this%spins%nspins_pr,this%spins%nspins))
         endif

         do s1=1,this%spins%nspins_pr  
          if(allocated(this%SHrep(s1,s1)%mat)) deallocate(this%SHrep(s1,s1)%mat)
          ii=nint(2*this%spins%spin(this%spins%kind(s1)))+1
          allocate(this%SHrep(s1,s1)%mat(ii,ii))
          this%SHrep(s1,s1)%mat=(0.0d0,0.0d0)
         enddo

         do s1=1,this%spins%nspins_pr
          do s2=s1+1,this%spins%nspins
           if(allocated(this%SHrep(s1,s2)%mat)) deallocate(this%SHrep(s1,s2)%mat)
           ii=nint(2*this%spins%spin(this%spins%kind(s1)))+1
           ii=ii*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)
           allocate(this%SHrep(s1,s2)%mat(ii,ii))
           this%SHrep(s1,s2)%mat=(0.0d0,0.0d0)
          enddo
         enddo

         !

         do s1=1,this%spins%nspins_pr  
          if(spin_id.ne.spin_id2) cycle
          if(spin_id.ne.-1 .and. spin_id.ne.s1) cycle

          if(allocated(this%SHrep(s1,s1)%mat)) deallocate(this%SHrep(s1,s1)%mat)
          ii=nint(2*this%spins%spin(this%spins%kind(s1)))+1
          allocate(this%SHrep(s1,s1)%mat(ii,ii))
          this%SHrep(s1,s1)%mat=(0.0d0,0.0d0)

          do ii=1,nint(2*this%spins%spin(this%spins%kind(s1)))+1
           do jj=1,nint(2*this%spins%spin(this%spins%kind(s1)))+1

            psi1(1)=ii-this%spins%spin(this%spins%kind(s1))-1
            psi1(2)=jj-this%spins%spin(this%spins%kind(s1))-1
            spin1=this%spins%spin(this%spins%kind(s1))

        !  Stevens operators matrix elements
                   
            do i=1,SH%nO
             if(this%spins%kind(s1).eq.SH%O(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%O(i)%mat_elem(psi1,spin1)
             endif
            enddo

        !  Zeeman Stevens operators matrix elements

            do i=1,SH%nG
             if(this%spins%kind(s1).eq.SH%G(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%G(i)%mat_elem(psi1,spin1,this%spins%Bfield,this%spins%bohr_mag(this%spins%kind(s1)))
             endif
            enddo

        !  Single Ion Anisotropy operators matrix elements

            do i=1,SH%nDSI
             if(this%spins%kind(s1).eq.SH%DSI(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%DSI(i)%mat_elem(psi1,spin1)
             endif
            enddo


           enddo ! 1 spin basis
          enddo ! 1 spin basis

         enddo ! spins

         do s1=1,this%spins%nspins_pr
          if(spin_id.ne.-1 .and. spin_id.ne.s1) cycle
          do s2=s1+1,this%spins%nspins
           if(spin_id2.ne.-1 .and. spin_id2.ne.s2) cycle

           if(allocated(this%SHrep(s1,s2)%mat)) deallocate(this%SHrep(s1,s2)%mat)
           ii=nint(2*this%spins%spin(this%spins%kind(s1)))+1
           ii=ii*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)
           allocate(this%SHrep(s1,s2)%mat(ii,ii))
           this%SHrep(s1,s2)%mat=(0.0d0,0.0d0)

           ii=1
           do ii1=1,nint(2*this%spins%spin(this%spins%kind(s1)))+1
            do ii2=1,nint(2*this%spins%spin(this%spins%kind(s2)))+1

             jj=1
             do jj1=1,nint(2*this%spins%spin(this%spins%kind(s1)))+1
              do jj2=1,nint(2*this%spins%spin(this%spins%kind(s2)))+1

               psi2(1,1)=ii1-this%spins%spin(this%spins%kind(s1))-1
               psi2(1,2)=ii2-this%spins%spin(this%spins%kind(s2))-1
               psi2(2,1)=jj1-this%spins%spin(this%spins%kind(s1))-1
               psi2(2,2)=jj2-this%spins%spin(this%spins%kind(s2))-1
               spin2(1)=this%spins%spin(this%spins%kind(s1))
               spin2(2)=this%spins%spin(this%spins%kind(s2))

        !  Isotropic exchange operators matrix elements

               do i=1,SH%nJ
                if(this%spins%kind(s1).eq.SH%J(i)%kind(1) .and. this%spins%kind(s2).eq.SH%J(i)%kind(2) )then
                 if(this%spins%ntot.gt.1 .and. s2.gt.this%spins%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)
                 endif
                endif
                if(this%spins%kind(s1).eq.SH%J(i)%kind(2) .and. this%spins%kind(s2).eq.SH%J(i)%kind(1) .and. & 
                   SH%J(i)%kind(1).ne.SH%J(i)%kind(2) ) then
                 if(this%spins%ntot.gt.1 .and. s2.gt.this%spins%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)
                 endif
                endif
               enddo

        !  Exchange Anisotry operators matrix elements

               do i=1,SH%nD2S

                b=MOD((s2-1),this%spins%nspins_pr)+1
                d=INT((s2-1)/this%spins%nspins_pr)+1

                if(this%spins%kind(s1).eq.SH%D2S(i)%kind(1) .and. this%spins%kind(s2).eq.SH%D2S(i)%kind(2)) then
                 if(this%spins%ntot.gt.1 .and. s2.gt.this%spins%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)
                 endif
                endif
                if(this%spins%kind(s1).eq.SH%D2S(i)%kind(2) .and. this%spins%kind(s2).eq.SH%D2S(i)%kind(1) .and. & 
                   SH%D2S(i)%kind(1).ne.SH%D2S(i)%kind(2) ) then
                 if(this%spins%ntot.gt.1 .and. s2.gt.this%spins%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)
                 endif
                endif              

               enddo

        !  Dipolar Exchange operators matrix elements

               do i=1,SH%nDdip
                if( s1.eq.SH%Ddip(i)%kind(1) .and. s2.eq.SH%Ddip(i)%kind(2) ) then
                 if(this%spins%ntot.gt.1 .and. s2.gt.this%spins%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%Ddip(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%Ddip(i)%mat_elem(psi2,spin2)
                 endif
                endif
      
               enddo

               jj=jj+1
              enddo
             enddo

             ii=ii+1
            enddo
           enddo

          enddo ! s2           
         enddo ! s1

        return
        end subroutine make_SH_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   IDENTIFY POTENTIALLY NON-ZERO MATRIX ELEMENTS OF THE HAMILTONIAN
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_Hmat_nodes(this)
        use lists_class   
        use blacs_utils
        implicit none
        class(spin_quantum_system) :: this
        integer                    :: l,l2,ii,jj,s1,s2,indxl2g,i
        integer                    :: Nloc_row,Nloc_col,numroc
        integer                    :: spin(2),ndiff
        integer                    :: t1,t2,rate
        logical                    :: save_it
        type(list)                 :: nodes(4)


         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Calculation of non-zero elements of H0'
          flush(6)
         endif
         
         call nodes(1)%init()
         call nodes(2)%init()
         call nodes(3)%init()
         call nodes(4)%init()

         Nloc_row = NUMROC(this%Hdim,NB,myrow,0,nprow)
         Nloc_col = NUMROC(this%Hdim,MB,mycol,0,npcol)

         do ii=1,Nloc_row
          do jj=1,Nloc_col

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)
         
           ndiff=0
           save_it=.true.
           spin=0

           do s1=1,this%spins%nspins
            if(ABS(this%basis(l,s1)-this%basis(l2,s1)).gt.1.0d-8)then
             ndiff=ndiff+1
             if(ndiff.gt.2)then
              save_it=.false.
              exit
             endif
             spin(ndiff)=s1
            endif
           enddo

           if(save_it)then               
             
            s1=min(spin(1),spin(2))
            s2=max(spin(1),spin(2))

            if(s1.gt.this%spins%nspins_pr) goto 12

            call nodes(1)%add_node(ii)
            call nodes(2)%add_node(jj)
            call nodes(3)%add_node(s1)
            call nodes(4)%add_node(s2)

           endif

12         continue           

          enddo
         enddo

         allocate(this%Hnodes(nodes(1)%nelem,4))

         call nodes(1)%reboot()
         call nodes(2)%reboot()
         call nodes(3)%reboot()
         call nodes(4)%reboot()

         do ii=1,nodes(1)%nelem
          call nodes(1)%rd_val(this%Hnodes(ii,1))
          call nodes(2)%rd_val(this%Hnodes(ii,2))
          call nodes(3)%rd_val(this%Hnodes(ii,3))
          call nodes(4)%rd_val(this%Hnodes(ii,4))
          call nodes(1)%skip()
          call nodes(2)%skip()
          call nodes(3)%skip()
          call nodes(4)%skip()
         enddo 

         call nodes(1)%delete()
         call nodes(2)%delete()
         call nodes(3)%delete()
         call nodes(4)%delete()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_Hmat_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   COMPUTE A MATRIX ELEMENT OF THE HAMILTONIAN
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function get_Hij(this,node) result(val)
        implicit none
        class(spin_quantum_system) :: this
        integer                    :: s1,s2,l,l2,i,j
        integer                    :: ii,jj,indxl2g
        complex(8)                 :: val
        integer                    :: node(4)
        
         l=indxl2g(node(1),NB,myrow,0,nprow)
         l2=indxl2g(node(2),MB,mycol,0,npcol)

         val=(0.0d0,0.0d0)

         if(node(3).ne.0 .and. node(4).ne.0)then

          s1=node(3)
          s2=node(4)

          i=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
          j=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
          ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
          i=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
          j=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
          jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j

          val=val+this%SHrep(s1,s2)%mat(ii,jj)

         endif

         if(node(3).eq.0 .and. node(4).ne.0)then

          s2=node(4)
          if(s2.gt.this%spins%nspins_pr)then
           do s1=1,this%spins%nspins_pr
            if(s1.lt.s2)then
             i=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
            if(s1.gt.s2)then
             i=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
             j=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s1)))+1)+j
             i=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             j=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
             jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s1)))+1)+j
             val=val+this%SHrep(s2,s1)%mat(ii,jj)
            endif
            if(s1.eq.s2)then
             ii=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
           enddo
          else
           do s1=1,this%spins%nspins
            if(s1.lt.s2)then
             i=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
            if(s1.gt.s2)then
             i=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
             j=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s1)))+1)+j
             i=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             j=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
             jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s1)))+1)+j
             val=val+this%SHrep(s2,s1)%mat(ii,jj)
            endif
            if(s1.eq.s2)then
             ii=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
           enddo
          endif

         endif
        
         
         if(node(3).eq.0 .and. node(4).eq.0)then

          do s1=1,this%spins%nspins_pr
           do s2=s1,this%spins%nspins
            if(s1.ne.s2)then
             i=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spins%spin(this%spins%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spins%spin(this%spins%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spins%spin(this%spins%kind(s2)))+1)+j
            else
             ii=nint(this%basis(l,s1)+this%spins%spin(this%spins%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spins%spin(this%spins%kind(s2))+1)
            endif
            val=val+this%SHrep(s1,s2)%mat(ii,jj)
           enddo
          enddo

         endif

        return
        end function get_Hij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   COMPUTE ALL MATRIX ELEMENTS OF THE HAMILTONIAN
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_Hmat(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use spinham_class
        implicit none
        class(spin_quantum_system) :: this
        integer                    :: l,ii,jj
        integer                    :: t1,t2,rate

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building Hamiltonian Matrix'
          flush(6)
         endif
      
         call this%H0%set(this%Hdim,this%Hdim,NB,MB)  

         do l=1,size(this%Hnodes,1)
             
          ii=this%Hnodes(l,1)
          jj=this%Hnodes(l,2)
           
          this%H0%mat(ii,jj)=this%get_Hij(this%Hnodes(l,1:4))

         enddo ! l

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_Hmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   COMPUTE THE MATRIX ELEMENTS OF THE SPIN OPERATORS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_Smat(this,print_si)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spin_quantum_system)          :: this
        type(dist_cmplx_mat), allocatable   :: Sz(:),Sy(:),Sx(:)
        logical,allocatable                 :: a(:),skip
        integer, allocatable                :: print_si(:)
        integer                             :: i,j,v,k,s,t,ii,jj
        integer                             :: t1,t2,rate,s2print,indxl2g
        double precision                    :: aaa
        complex(8), allocatable             :: Mtmp(:)

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Calculation of the Spin operators matrix elements'
          flush(6)
         endif

         s2print=size(print_si)
         allocate(Sz(s2print))
         allocate(Sy(s2print))
         allocate(Sx(s2print))

         do i=1,s2print
          call Sx(i)%set(this%Hdim,this%Hdim,NB,MB)
          call Sy(i)%set(this%Hdim,this%Hdim,NB,MB)
          call Sz(i)%set(this%Hdim,this%Hdim,NB,MB)
          Sx(i)%mat=(0.0d0,0.0d0)
          Sy(i)%mat=(0.0d0,0.0d0)
          Sz(i)%mat=(0.0d0,0.0d0)
         enddo
                    
         allocate(a(this%spins%nspins))

         allocate(Mtmp(3))
         Mtmp=(0.0d0,0.0d0)

         do ii=1,size(Sz(1)%mat,1)
          do jj=1,size(Sz(1)%mat,2)
           s=indxl2g(ii,NB,myrow,0,nprow)
           t=indxl2g(jj,MB,mycol,0,npcol)
           if (s.eq.t) then                   
            do v=1,this%spins%nspins

             do k=1,s2print
              if(print_si(k).eq.v)then
               Sz(k)%mat(ii,jj)=this%basis(s,v)
              endif
              if(print_si(k).eq.-1)then
               Sz(k)%mat(ii,jj)=Sz(k)%mat(ii,jj)+this%basis(s,v)
              endif
             enddo
 
            enddo
           endif
          enddo
         enddo
          
         do ii=1,size(Sx(1)%mat,1)
          do jj=1,size(Sx(1)%mat,2)

           s=indxl2g(ii,NB,myrow,0,nprow)
           t=indxl2g(jj,MB,mycol,0,npcol)

           do v=1,this%spins%nspins

            skip=.true.
            do k=1,s2print
             if(print_si(k).eq.v) skip=.false.
             if(print_si(k).eq.-1) skip=.false.
            enddo

            if(skip) cycle

            a=.true.
            do k=1,this%spins%nspins
             if(k.ne.v)then
              if(ABS(this%basis(s,k)-this%basis(t,k)).gt.1.0d-8)then
               a(k)=.false.
               exit
              endif
             endif
            enddo

            if(ALL(a))then

! S+
             if(abs(this%basis(s,v)-this%basis(t,v)-1).lt.1.0E-06)then

              aaa=dsqrt((this%spins%spin(this%spins%kind(v))-this%basis(t,v))* &
                        (this%spins%spin(this%spins%kind(v))+this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)
              Mtmp(2)=cmplx(0.0d0,-1.0d0*aaa,8)

              do k=1,s2print
               if(print_si(k).eq.v)then
                Sx(k)%mat(ii,jj)=Mtmp(1)
                Sy(k)%mat(ii,jj)=Mtmp(2)
               endif
               if(print_si(k).eq.-1)then
                Sx(k)%mat(ii,jj)=Sx(k)%mat(ii,jj)+Mtmp(1)
                Sy(k)%mat(ii,jj)=Sy(k)%mat(ii,jj)+Mtmp(2)
               endif
              enddo

             endif

! S-
             if(abs(this%basis(s,v)-this%basis(t,v)+1).lt.1.0E-06)then

              aaa=dsqrt((this%spins%spin(this%spins%kind(v))+this%basis(t,v))* &
                        (this%spins%spin(this%spins%kind(v))-this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)
              Mtmp(2)=cmplx(0.0d0,aaa,8)

              do k=1,s2print
               if(print_si(k).eq.v)then
                Sx(k)%mat(ii,jj)=Mtmp(1)
                Sy(k)%mat(ii,jj)=Mtmp(2)
               endif
               if(print_si(k).eq.-1)then
                Sx(k)%mat(ii,jj)=Sx(k)%mat(ii,jj)+Mtmp(1)
                Sy(k)%mat(ii,jj)=Sy(k)%mat(ii,jj)+Mtmp(2)
               endif
              enddo

             endif

            endif

           enddo
          enddo
         enddo

         do s=1,s2print
          call this%to_eigenbasis(Sx(s))
          call this%to_eigenbasis(Sy(s))
          call this%to_eigenbasis(Sz(s))
         enddo

         if(allocated(this%QMOP))then
          do i=1,size(this%QMOP)
           deallocate(this%QMOP(i)%mat)           
          enddo
          deallocate(this%QMOP)           
         endif
         allocate(this%QMOP(3*s2print))

         do s=1,s2print
          call this%QMOP(3*(s-1)+1)%set(this%Hdim,this%Hdim,NB,MB)
          call this%QMOP(3*(s-1)+2)%set(this%Hdim,this%Hdim,NB,MB)
          call this%QMOP(3*(s-1)+3)%set(this%Hdim,this%Hdim,NB,MB)
          this%QMOP(3*(s-1)+1)%mat=Sx(s)%mat
          this%QMOP(3*(s-1)+2)%mat=Sy(s)%mat
          this%QMOP(3*(s-1)+3)%mat=Sz(s)%mat
          deallocate(Sx(s)%mat)
          deallocate(Sy(s)%mat)
          deallocate(Sz(s)%mat)
         enddo
                
         deallocate(Sx)
         deallocate(Sy)
         deallocate(Sz)
         deallocate(Mtmp)
         deallocate(a)
         

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return 
        end subroutine make_Smat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   CORRECT THE PHASE OF Hmat EIGENVECTORS 
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_tinv(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spin_quantum_system)          :: this
        type(dist_cmplx_mat)                :: Sx
        logical,allocatable                 :: a(:),skip
        integer                             :: i,j,v,k,s,t,ii,jj
        integer                             :: t1,t2,rate,s2print,indxl2g
        double precision                    :: aaa,phi
        complex(8), allocatable             :: Mtmp(:)

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Correcting the phase of Hmat eigenvectors'
          flush(6)
         endif


         call Sx%set(this%Hdim,this%Hdim,NB,MB)
         Sx%mat=(0.0d0,0.0d0)
                    
         allocate(a(this%spins%nspins))

         allocate(Mtmp(3))
         Mtmp=(0.0d0,0.0d0)
          
         do ii=1,size(Sx%mat,1)
          do jj=1,size(Sx%mat,2)

           s=indxl2g(ii,NB,myrow,0,nprow)
           t=indxl2g(jj,MB,mycol,0,npcol)

           do v=1,this%spins%nspins

            a=.true.
            do k=1,this%spins%nspins
             if(k.ne.v)then
              if(ABS(this%basis(s,k)-this%basis(t,k)).gt.1.0d-8)then
               a(k)=.false.
               exit
              endif
             endif
            enddo

            if(ALL(a))then

! S+
             if(abs(this%basis(s,v)-this%basis(t,v)-1).lt.1.0E-06)then

              aaa=dsqrt((this%spins%spin(this%spins%kind(v))-this%basis(t,v))* &
                        (this%spins%spin(this%spins%kind(v))+this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)

              Sx%mat(ii,jj)=Sx%mat(ii,jj)+Mtmp(1)

             endif

! S-
             if(abs(this%basis(s,v)-this%basis(t,v)+1).lt.1.0E-06)then

              aaa=dsqrt((this%spins%spin(this%spins%kind(v))+this%basis(t,v))* &
                        (this%spins%spin(this%spins%kind(v))-this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)

              Sx%mat(ii,jj)=Sx%mat(ii,jj)+Mtmp(1)

             endif

            endif

           enddo
          enddo
         enddo

         call this%to_eigenbasis(Sx)
        
         do ii=1,this%Hdim-1         
          phi=atan2(aimag(Sx%mat(i,i+1)),dble(Sx%mat(i,i+1)))
          this%H%mat(:,i+1)=this%H%mat(:,i+1)*exp(cmplx(0.0d0,-phi,8))
         enddo
                
         deallocate(Mtmp)
         deallocate(a)        

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return 
        end subroutine make_tinv

        end module quantum_systems_class
