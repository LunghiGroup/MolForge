        module quantum_systems_class
        use mpi
        use blacs_utils
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
         procedure     ::  make_R21_limbladian
         procedure     ::  make_R22_limbladian
         procedure     ::  make_R41_limbladian
         procedure     ::  X2Q
         procedure     ::  XY2QQ
        end type open_quantum_system

        type, extends(open_quantum_system) :: spin_quantum_system
         type(spins_group)     :: spins
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

        subroutine make_R21_limbladian(this,min_ener,max_ener)
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
          write(*,*) '     Building the second-order Limbladian operator with linear couplig'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))

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
         
           call this%to_eigenbasis(this%Vq)

           Vmat=(0.0d0,0.0d0)

           do l2=1,this%Hdim
            do l=1,this%Hdim
             call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
             Vmat(l2,l)=valc
            enddo
           enddo                   

         ! compute R21

           call make_R21(this,Vmat,this%temp,freq,this%smear,this%type_smear)

          enddo 
          ph=ph+1
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R21_limbladian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN R22
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R22_limbladian(this,min_ener,max_ener)
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
          write(*,*) '     Building the second-order Limbladian operator with quadratic couplig'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))

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
         
             call this%to_eigenbasis(this%Vq)

             Vmat=(0.0d0,0.0d0)

             do l2=1,this%Hdim
              do l=1,this%Hdim
               call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
               Vmat(l2,l)=valc
              enddo
             enddo                   

             call make_R22(this,Vmat,this%temp,freq,freq2,this%smear,this%smear,this%type_smear)

            enddo ! bn2
           enddo ! ph2

          enddo 
          ph=ph+1
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R22_limbladian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN R41
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R41_limbladian(this,min_ener,max_ener)
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

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Building the fourth-order Limbladian operator with linear couplig'
          flush(6)
         endif

         allocate(Vmat(this%Hdim,this%Hdim))
         allocate(V2mat(this%Hdim,this%Hdim))

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
         
           call this%to_eigenbasis(this%Vq)

           Vmat=(0.0d0,0.0d0)

           do l2=1,this%Hdim
            do l=1,this%Hdim
             call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
             Vmat(l2,l)=valc
            enddo
           enddo                   

           do ph2=1,this%phonons%ntot
            do bn2=1,size(this%phonons%list(ph2)%freq)

             if (ph2.eq.1 .and. bn2.le.3) cycle
             freq2=this%phonons%list(ph2)%freq(bn2)
             if (freq2.lt.min_ener) cycle
             if (freq2.gt.max_ener) cycle           

             if(.not.allocated(this%phonons%list(ph2)%hess) ) call this%phonons%list(ph2)%diagD(this%lattice)

             if ( allocated(this%Vx) ) call X2Q(this,ph2,bn2)
         
             call this%to_eigenbasis(this%Vq)

             V2mat=(0.0d0,0.0d0)

             do l2=1,this%Hdim
              do l=1,this%Hdim
               call pzelget('A',' ',valc,this%Vq%mat,l2,l,this%Vq%desc)
               V2mat(l2,l)=valc
              enddo
             enddo                   

             call make_R41(this,Vmat,V2mat,this%temp,freq,freq2,this%smear,this%smear,this%type_smear)

            enddo ! bn2
           enddo ! ph2

          enddo 
          ph=ph+1
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine make_R41_limbladian

        end module quantum_systems_class
