        program CrystalField
        use lapack_diag_simm
        implicit none
        interface
         subroutine read_orca_mat(N,Mat,filename)
         implicit none
         integer                       :: i,j,l,k,N,N1,N2
         character(len=100)            :: filename,skip
         double precision, allocatable :: MAT(:,:)
         end subroutine read_orca_mat
         subroutine projH(lmax,Nj,Ener,Jz,Bz,O,alpha,beta,gamma)
         use stevens_class
         use spinham_class 
         implicit none
         integer                        :: lmax,Odim,Nj
         integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
         double precision, allocatable  :: Ener(:)
         double precision               :: JMax,q1,q2,avg_ener,Bz,Bfield(3)
         double complex, allocatable    :: Jz(:,:),B(:),A(:,:),O(:),work(:)
         double complex                 :: mat_elem
         type(OSItensor), allocatable   :: Os(:)
         double precision, optional     :: alpha,beta,gamma
         end subroutine projH
         subroutine projdH(lmax,Nj,DHSOC,Jz,alpha,beta,gamma)
         use stevens_class
         use spinham_class 
         implicit none
         integer                        :: lmax,Odim,Nj
         integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
         double complex, allocatable    :: DHSOC(:,:)
         double precision               :: JMax,q1,q2,avg_ener
         double complex, allocatable    :: Jz(:,:),B(:),A(:,:),work(:)
         double complex                 :: mat_elem
         type(OSItensor), allocatable   :: Os(:)
         double precision, optional     :: alpha,beta,gamma
         end subroutine projdH
         subroutine diagH(Hdim,lmax,O,JMax)
         use lapack_diag_simm
         use stevens_class
         integer                        :: lmax,Hdim,i,j,k,l,q,s
         double precision               :: j1,j2,JMax
         double precision, allocatable  :: Ener(:)
         double complex, allocatable    :: Hmat(:,:),O(:)
         double complex                 :: val_tmp,val
         end subroutine diagH
        end interface
        double precision, allocatable :: SOCR(:,:),SOCI(:,:),EIG(:),Mat(:,:),Mj(:),EIG2(:),coeff(:,:) 
        double precision              :: norm,phi,alpha,beta,gamma,JMax,Bz
        double complex, allocatable   :: SOC(:,:),Lz(:,:),Sz(:,:),Jz(:,:),O(:),Sx(:,:),Lx(:,:),Jx(:,:),Jx2(:,:)
        double complex, allocatable   :: DHSOC(:,:)
        double precision,allocatable  :: valr(:),valc(:),phase(:)
        character(len=100)            :: filename,word
        integer                       :: i,j,k,N,Nj,k1,k2,lmax
        logical                       :: rotate_CF=.false.,spin_only=.false.,read_dH=.false.,add_field=.false.

        double precision, allocatable :: D(:),E(:)
        double complex, allocatable   :: Tau(:),work(:)
        integer                       :: info,lwork
        double precision              :: bohr_mag=-0.466867723

         if( iargc().eq.0)then
          write(*,*) 'CF Usage:'               
          write(*,*) '-JMult     : J Multiplicity of the ground state multiplet'
          write(*,*) '-CISIZE    : Size of CI basis set'               
          write(*,*) '-lmax      : Max Order of CF operators'               
          write(*,*) '-rot       : ZYZ Euler angoles (rad) for CF rotation'               
          write(*,*) '-spin_only : Assume <L>=0'               
          write(*,*) '-read_dH   : read the Vibronic coupling operator DH.dat'               
          write(*,*) '-Bz        : Apply a Magnetic Field along z'               
          stop
         endif

         do i=1,iargc()

          call getarg(i,word)
          
          select case (trim(word))

             case ('-JMult')
                 call getarg(i+1,word)
                 read(word,*) Nj
                 JMax=(Nj-1)/2.0d0

             case ('-CISIZE')
                 call getarg(i+1,word)
                 read(word,*) N

             case ('-lmax')
                 call getarg(i+1,word)
                 read(word,*) lmax

             case ('-Bz')
                 call getarg(i+1,word)
                 read(word,*) Bz
                 add_field=.true.

             case ('-rot')
                 call getarg(i+1,word)
                 read(word,*) alpha
                 call getarg(i+2,word)
                 read(word,*) beta
                 call getarg(i+3,word)
                 read(word,*) gamma
                 rotate_CF=.true.

             case ('-spin_only')
                 spin_only=.true.

             case ('-read_dH')
                 read_dH=.true.

          end select

         enddo

         allocate(SOC(N,N))
         allocate(phase(N))
         allocate(EIG(N))
         allocate(EIG2(N))
         allocate(Jz(Nj,Nj))
         allocate(Jx(Nj,Nj))
         allocate(Jx2(Nj,Nj))
         allocate(Mj(Nj))

         allocate(D(Nj))
         allocate(E(Nj))
         allocate(Tau(Nj))
         lwork=10000
         allocate(work(lwork))

         filename='PHASE.dat'
         open(11,file=filename)
         do i=1,N
          read(11,*) phase(i)
         enddo
         close(11)

         filename='SOCR.dat'
         call read_orca_mat(N,SOCR,filename)
         filename='SOCI.dat'
         call read_orca_mat(N,SOCI,filename)
         SOC=cmplx(SOCR,SOCI,8)

         filename='Lz.dat'
         call read_orca_mat(N,Mat,filename)
         Lz=cmplx(0.0d0,Mat,8)
         filename='Lx.dat'
         call read_orca_mat(N,Mat,filename)
         Lx=cmplx(0.0d0,Mat,8)

         filename='Sz.dat'
         call read_orca_mat(N,Mat,filename)
         Mat=Mat*0.5d0
         Sz=cmplx(Mat,0.0d0,8)

         filename='Sx.dat'
         call read_orca_mat(N,Mat,filename)
         Mat=Mat*0.5d0
         Sx=cmplx(Mat,0.0d0,8)

         if(add_field)then

          write(*,*) "Applying a Magnetic Field of ", Bz

          if(spin_only)then
           SOC=SOC+2.002319304362*Sz*2.12719108d-6*Bz
          else
           SOC=SOC+(Lz+2.002319304362*Sz)*2.12719108d-6*Bz
          endif
         endif

         do i=1,N
          do j=1,N
           SOC(i,j)=SOC(i,j)*phase(i)*phase(j)
          enddo
         enddo

         call new_diag(N,SOC,EIG)

         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'First Nj eigevanlues of the electronic H'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         open(11,file='SOC_EIGVEC_R.dat')
         open(12,file='SOC_EIGVEC_C.dat')
         open(13,file='SOC_EIGVAL.dat')

         do i=1,N
          write(11,*) (dble(SOC(i,j)),j=1,N)
          write(12,*) (aimag(SOC(i,j)),j=1,N)
          write(13,*) EIG(i)
         enddo

         close(11)
         close(12)
         close(13)

         do i=1,Nj
          EIG2(i)=(EIG(i)-EIG(1))*219474.63
          write(*,*) 'EIG: ',i,' ENER: ',EIG2(i)
         enddo

!         EIG=EIG2
        
         Jz=(0.0d0,0.0d0)
         Jx=(0.0d0,0.0d0)
         Jx2=(0.0d0,0.0d0)

         do i=1,Nj
          do j=1,Nj
           do k1=1,N
            do k2=1,N
               
             Jz(i,j)=Jz(i,j)&
                            +Sz(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j) !&

             if(.not. spin_only)then
             Jz(i,j)=Jz(i,j)&
                            +Lz(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j)
             endif

            enddo
           enddo
          enddo
         enddo

         do i=1,Nj
          do j=1,Nj
           do k1=1,N
            do k2=1,N
               
             Jx(i,j)=Jx(i,j)&
                            +Sx(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j) !&
!                            +Lx(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j)

             if(.not. spin_only)then
              Jx(i,j)=Jx(i,j)&
                            +Lx(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j)
             endif

            enddo
           enddo
          enddo
         enddo
        
         call new_diag(Nj,Jz,Mj)

         do i=1,Nj-1
          do k1=1,Nj
           do k2=1,Nj               
             Jx2(i,i+1)=Jx2(i,i+1)&
                            +Jx(k1,k2)*conjg(Jz(k1,i))*Jz(k2,i+1) 
           enddo
          enddo
          phi=atan2(aimag(Jx2(i,i+1)),dble(Jx2(i,i+1)))
          Jz(:,i+1)=Jz(:,i+1)*exp(cmplx(0.0d0,-phi,8))
         enddo
        
         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'Jz eigenvalues from electronic H NjxNj subspace'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         do i=1,Nj
          write(*,*) 'EIG: ',i,' Mj: ',Mj(i)
         enddo

         Jz=transpose(conjg(Jz))

!         open(12,file='eig.dat')
!         do i=1,Nj
!          read(12,*) EIG(i)
!         enddo
!         close(12)

!         open(12,file='Jz.dat')
!         allocate(coeff(4,2))

!         do i=1,Nj
!          read(12,*) coeff(1,1),coeff(1,2),coeff(2,1),coeff(2,2),coeff(3,1),coeff(3,2),coeff(4,1),coeff(4,2)
!          Jz(i,1)=cmplx(coeff(1,1),coeff(1,2),8)
!          Jz(i,2)=cmplx(coeff(2,1),coeff(2,2),8)
!          Jz(i,3)=cmplx(coeff(3,1),coeff(3,2),8)
!          Jz(i,4)=cmplx(coeff(4,1),coeff(4,2),8)
!         enddo
!         read(12,*)
!         do i=1,Nj
!          read(12,*) coeff(1,1),coeff(1,2),coeff(2,1),coeff(2,2),coeff(3,1),coeff(3,2),coeff(4,1),coeff(4,2)
!          Jz(i,5)=cmplx(coeff(1,1),coeff(1,2),8)
!          Jz(i,6)=cmplx(coeff(2,1),coeff(2,2),8)
!          Jz(i,7)=cmplx(coeff(3,1),coeff(3,2),8)
!          Jz(i,8)=cmplx(coeff(4,1),coeff(4,2),8)
!         enddo
!         read(12,*)
!         do i=1,Nj
!          read(12,*) coeff(1,1),coeff(1,2),coeff(2,1),coeff(2,2),coeff(3,1),coeff(3,2),coeff(4,1),coeff(4,2)
!          Jz(i,9)=cmplx(coeff(1,1),coeff(1,2),8)
!          Jz(i,10)=cmplx(coeff(2,1),coeff(2,2),8)
!          Jz(i,11)=cmplx(coeff(3,1),coeff(3,2),8)
!          Jz(i,12)=cmplx(coeff(4,1),coeff(4,2),8)
!         enddo
!         read(12,*)
!         do i=1,Nj
!          read(12,*) coeff(1,1),coeff(1,2),coeff(2,1),coeff(2,2),coeff(3,1),coeff(3,2),coeff(4,1),coeff(4,2)
!          Jz(i,13)=cmplx(coeff(1,1),coeff(1,2),8)
!          Jz(i,14)=cmplx(coeff(2,1),coeff(2,2),8)
!          Jz(i,15)=cmplx(coeff(3,1),coeff(3,2),8)
!          Jz(i,16)=cmplx(coeff(4,1),coeff(4,2),8)
!         enddo

         if(read_dH)then

          open(11,file='DH.dat')
          open(12,file='DHSOC.dat')

          allocate(DHSOC(N,N))
          allocate(valr(N))
          allocate(valc(N))
          DHSOC=(0.0d0,0.0d0)

          do i=1,N
           read(11,*) (valr(j),valc(j),j=1,N)
           do j=1,N
            DHSOC(i,j)=cmplx(valr(j),valc(j),8)
           enddo
           write(12,*) (dble(DHSOC(i,j)),aimag(DHSOC(i,j)),j=1,N)
          enddo

          close(11)
          close(12)
         
         endif

         if(rotate_CF)then                 
          call projH(lmax,Nj,EIG2,Jz,Bz,O,alpha,beta,gamma)
          if(read_dH) call projdH(lmax,Nj,DHSOC,Jz,alpha,beta,gamma)
         else
          call projH(lmax,Nj,EIG2,Jz,Bz,O)
          if(read_dH) call projdH(lmax,Nj,DHSOC,Jz)
         endif

         call diagH(Nj,lmax,O,JMax)

        return
        end program CrystalField

        subroutine diagH(Hdim,lmax,O,Jmax)
        use lapack_diag_simm
        use stevens_class
        integer                        :: lmax,Hdim,i,j,k,l,q,s
        double precision               :: j1,j2,JMax
        double precision, allocatable  :: Ener(:),Jzeig(:)
        double complex, allocatable    :: Hmat(:,:),O(:),Jz(:,:)
        double complex                 :: val_tmp,val

         allocate(Hmat(Hdim,Hdim))         
         allocate(Ener(Hdim))         

         j1=-JMax
         do i=1,Hdim
          j2=-JMax
          do j=1,Hdim

           val=(0.0d0,0.0d0)
           s=1
           do l=2,lmax,2
            do q=-l,l
             val_tmp=(0.0d0,0.0d0)
             call stevens_mat_elem(l,q,JMax,j1,JMax,j2,val_tmp)
             val=val+val_tmp*O(s)
             s=s+1
            enddo
           enddo

           Hmat(i,j)=val

           j2=j2+1.0d0
          enddo
          j1=j1+1.0d0
         enddo

         call new_diag(Hdim,Hmat,Ener)

         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'H eigenstates in the Jz eigenstates basis'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         do i=1,Hdim
          write(*,*) '#######',i,Ener(i)-Ener(1)
          do j=1,Hdim
           write(*,*) j,dble(Hmat(j,i)),aimag(Hmat(j,i)),dble(Hmat(j,i)*conjg(Hmat(j,i)))
          enddo
         enddo

!         allocate(Jz(Hdim,Hdim))
!         allocate(Jzeig(Hdim))

!         Jz=(0.0d0,0.0d0)
!         Jzeig=0.0d0

!         do i=1,Hdim
!          do j=1,Hdim
!           do l=1,Hdim            
!            Jz(i,j)=Jz(i,j)+conjg(Hmat(l,i))*Hmat(l,j)*(-JMax+(l-1))
!           enddo
!          enddo
!         enddo
          
!         call new_diag(Hdim,Jz,Jzeig)

!         write(*,*) '#################################################'
!         write(*,*) '#################################################'
!         write(*,*) 'Jz eigenstates in the H eigenstates basis'
!         write(*,*) '#################################################'
!         write(*,*) '#################################################'

!         do i=1,Hdim
!          write(*,*) '#######',i,Jzeig(i)
!          do j=1,Hdim
!           write(*,*) j,dble(Jz(j,i)),aimag(Jz(j,i)),dble(Jz(j,i)*conjg(Jz(j,i)))
!          enddo
!         enddo
         
        return
        end subroutine diagH

        subroutine projdH(lmax,Nj,DHSOC,Jz,alpha,beta,gamma)
        use stevens_class
        use spinham_class 
        implicit none
        integer                        :: lmax,Odim,Nj
        integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
        double complex, allocatable    :: DHSOC(:,:)
        double precision               :: JMax,q1,q2,avg_ener
        double complex, allocatable    :: Jz(:,:),B(:),A(:,:),work(:)
        double complex                 :: mat_elem
        type(OSItensor), allocatable   :: Os(:)
        double precision, optional     :: alpha,beta,gamma

         JMax=(Nj-1)/2.0d0

         Odim=0
         do i=2,lmax,2
          Odim=Odim+(2*i+1)
         enddo

         allocate(B(Nj*Nj))
         allocate(A(Nj*Nj,Odim))

         A=(0.0d0,0.0d0)        
         B=(0.0d0,0.0d0)        

         k=1
         do i1=1,Nj
          do i2=1,Nj
          
           do v=1,Nj
            do l=1,Nj
             B(k)=B(k)+DHSOC(l,v)*conjg(Jz(i2,v))*Jz(i1,l)
            enddo
           enddo

           q1=i1-((Nj-1)/2.0d0)-1
           q2=i2-((Nj-1)/2.0d0)-1

           s=1
           do l=2,lmax,2
            do j=-l,l
             call stevens_mat_elem(l,j,JMax,q1,JMax,q2,mat_elem)
             A(k,s)=mat_elem
             s=s+1
            enddo
           enddo

           k=k+1
          enddo
         enddo

         lwork=Nj*Nj+64*Nj*Nj+1000
         allocate(work(lwork))

         call zgels('N',Nj*Nj,Odim,1,A,Nj*Nj,B,Nj*Nj,WORK,LWORK,inf)
         if(inf.ne.0)then
          write(*,*) 'zgels failed',inf
          stop
         endif

         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'Derivatives of the Crystal Field Parameters:'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         s=1
         do l=2,lmax,2
          do j=-l,l
           write(*,*) l,j,dble(B(s))!,aimag(B(s))
           s=s+1
          enddo
         enddo

         if(present(alpha))then

          allocate(Os(lmax/2))

          write(*,*) '#################################################'
          write(*,*) '#################################################'
          write(*,*) 'Rotated Derivatives of the CF Parameters:'
          write(*,*) '#################################################'
          write(*,*) '#################################################'

          s=1
          ll=1
          do l=2,lmax,2
           Os(ll)%k=l
           allocate(Os(ll)%B(2*l+1))
           allocate(Os(ll)%q(2*l+1))
           v=1
           do j=-l,l
            Os(ll)%B(v)=B(s)
            Os(ll)%q(v)=j
            s=s+1
            v=v+1
           enddo
           call Os(ll)%rot(alpha,beta,gamma)
           v=1
           do j=-l,l
            write(*,*) Os(ll)%k,Os(ll)%q(v),dble(Os(ll)%B(v))
            v=v+1
           enddo
           ll=ll+1
          enddo

         endif

        return
        end subroutine projdH

        subroutine projH(lmax,Nj,Ener,Jz,Bz,O,alpha,beta,gamma)
        use stevens_class
        use spinham_class 
        implicit none
        integer                        :: lmax,Odim,Nj
        integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
        double precision, allocatable  :: Ener(:)
        double precision               :: JMax,q1,q2,avg_ener,Bz,Bfield(3)
        double complex, allocatable    :: Jz(:,:),B(:),A(:,:),O(:),work(:)
        double complex                 :: mat_elem,gtens_loc
        type(OSItensor), allocatable   :: Os(:)
        double precision, optional     :: alpha,beta,gamma

         Bfield=0.0d0
         Bfield(3)=Bz

         JMax=(Nj-1)/2.0d0

         Odim=0
         do i=2,lmax,2
          Odim=Odim+(2*i+1)
         enddo

         avg_ener=0.0d0
         do i=1,Nj
          avg_ener=avg_ener+Ener(i)
         enddo
         avg_ener=avg_ener/Nj

         allocate(O(Odim))
         allocate(B(Nj*Nj))
         allocate(A(Nj*Nj,Odim+1))

         A=(0.0d0,0.0d0)        
         B=(0.0d0,0.0d0)        
         O=(0.0d0,0.0d0)        

         k=1
         do i1=1,Nj
          do i2=1,Nj
          
           do v=1,Nj
            B(k)=B(k)+Ener(v)*conjg(Jz(i2,v))*Jz(i1,v)
           enddo

           q1=i1-((Nj-1)/2.0d0)-1
           q2=i2-((Nj-1)/2.0d0)-1

           s=1
           do l=2,lmax,2
            do j=-l,l
             call stevens_mat_elem(l,j,JMax,q1,JMax,q2,mat_elem)
             A(k,s)=mat_elem
             s=s+1
            enddo
           enddo

           do l=3,3
            do j=3,3
             mat_elem=gtens_loc(q1,q2,Jmax,Bfield,l,j)
             A(k,s)=mat_elem
             s=s+1
            enddo
           enddo

           k=k+1
          enddo
         enddo

         lwork=Nj*Nj+64*Nj*Nj+1000
         allocate(work(lwork))

         call zgels('N',Nj*Nj,Odim+1,1,A,Nj*Nj,B,Nj*Nj,WORK,LWORK,inf)
         if(inf.ne.0)then
          write(*,*) 'zgels failed',inf
          stop
         endif

         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'Crystal Field Parameters:'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         s=1
         do l=2,lmax,2
          do j=-l,l
           write(*,*) l,j,dble(B(s))!,aimag(B(s))
           O(s)=B(s)
           s=s+1
          enddo
         enddo
         
         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'G-tensor Parameters:'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         do l=3,3
          do j=3,3
           write(*,*) l,j,dble(B(s))!,aimag(B(s))
           s=s+1
          enddo
         enddo

         if(present(alpha))then

          allocate(Os(lmax/2))

          write(*,*) '#################################################'
          write(*,*) '#################################################'
          write(*,*) 'Rotated Crystal Field Parameters:'
          write(*,*) '#################################################'
          write(*,*) '#################################################'

          s=1
          ll=1
          do l=2,lmax,2
           Os(ll)%k=l
           allocate(Os(ll)%B(2*l+1))
           allocate(Os(ll)%q(2*l+1))
           v=1
           do j=-l,l
            Os(ll)%B(v)=B(s)
            Os(ll)%q(v)=j
            s=s+1
            v=v+1
           enddo
           call Os(ll)%rot(alpha,beta,gamma)
           v=1
           do j=-l,l
            write(*,*) Os(ll)%k,Os(ll)%q(v),dble(Os(ll)%B(v))
            v=v+1
           enddo
           ll=ll+1
          enddo

          s=1
          ll=1
          do l=2,lmax,2
           v=1
           do j=-l,l
            O(s)=Os(ll)%B(v)
            s=s+1
            v=v+1
           enddo
           ll=ll+1
          enddo

         endif

        return
        end subroutine projH

        subroutine read_orca_mat(N,Mat,filename)
        implicit none
        integer                       :: i,j,l,k,N,N1,N2
        character(len=100)            :: filename,skip
        double precision, allocatable :: MAT(:,:)

         open(11,file=trim(filename))
         
         N1=INT(N/6)
         N2=MOD(N,6)

         if(allocated(MAT)) deallocate(MAT)
         allocate(MAT(N,N))

         read(11,*)
         l=1
         do k=1,N1
          read(11,*)
          do i=1,N           
           read(11,*) skip,(MAT(i,l+j),j=0,5)
          enddo
          l=l+6
         enddo
         if(N2.gt.0)then
          read(11,*)
          do i=1,N
          read(11,*) skip,(MAT(i,l+j),j=0,N2-1)
          enddo
         endif

         close(11)

        return
        end subroutine read_orca_mat

        function gtens_loc(a1,b1,spin,B,i,j) result(val)
        implicit none
        double complex                          :: val,C
        integer                                 :: i,j
        double precision                        :: bohr_mag=-0.466867723
        double precision                        :: mapp(2),B(3),G(3,3),spin,a1,b1

         val=(0.0d0,0.0d0)
         G=0.0d0
         G(i,j)=1.0d0         
         mapp(1)=a1
         mapp(2)=b1

       !S-
         if(abs(mapp(1)-mapp(2)+1).lt.1.0E-06)then
          C=(0.0d0,0.0d0)
          C=0.5d0*G(1,1)*B(1)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*G(2,1)*B(2)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*G(3,1)*B(3)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          val=val+C
          C=(0.0d0,0.0d0)
          C=0.5d0*G(1,2)*B(1)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*G(2,2)*B(2)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*G(3,2)*B(3)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          val=val+C*CMPLX(0.0d0,1.0d0,8)
         endif
        !S+
         if(abs(mapp(1)-mapp(2)-1).lt.1.0E-06)then
          C=(0.0d0,0.0d0)
          C=0.5d0*G(1,1)*B(1)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*G(2,1)*B(2)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*G(3,1)*B(3)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          val=val+C!CMPLX(C,0,8)
          C=(0.0d0,0.0d0)
          C=0.5d0*G(1,2)*B(1)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*G(2,2)*B(2)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*G(3,2)*B(3)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          val=val+C*CMPLX(0.d0,-1.0d0,8)
         endif
        !Sz
         if(abs(mapp(1)-mapp(2)).lt.1.0E-06)then
           val=val+G(1,3)*B(1)*mapp(2)+G(2,3)*B(2)   &
                  *mapp(2)+G(3,3)*B(3)*mapp(2)
         endif

         val=-1.0d0*val*bohr_mag ! bohr magneton in cm-1/T

        return 
        end function gtens_loc
