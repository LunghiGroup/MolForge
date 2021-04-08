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
         subroutine projH(lmax,Nj,Ener,Jz,O,alpha,beta,gamma)
         use stevens_class
         use spinham_class 
         implicit none
         integer                        :: lmax,Odim,Nj
         integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
         double precision, allocatable  :: Ener(:)
         double precision               :: jj,q1,q2,avg_ener
         double complex, allocatable    :: Jz(:,:),B(:),A(:,:),O(:),work(:)
         double complex                 :: mat_elem
         type(OSItensor), allocatable   :: Os(:)
         double precision, optional     :: alpha,beta,gamma
         end subroutine projH
         subroutine diagH(Hdim,lmax,O)
         use lapack_diag_simm
         use stevens_class
         integer                        :: lmax,Hdim,i,j,k,l,q,s
         double precision               :: j1,j2,jj
         double precision, allocatable  :: Ener(:)
         double complex, allocatable    :: Hmat(:,:),O(:)
         double complex                 :: val_tmp,val
         end subroutine diagH
        end interface
        double precision, allocatable :: SOCR(:,:),SOCI(:,:),EIG(:),Mat(:,:),Mj(:),EIG2(:),coeff(:,:) 
        double precision              :: norm,phi,alpha,beta,gamma
        double complex, allocatable   :: SOC(:,:),Lz(:,:),Sz(:,:),Jz(:,:),O(:),Sx(:,:),Lx(:,:),Jx(:,:),Jx2(:,:)
        character(len=100)            :: filename,word
        integer                       :: i,j,k,N,Nj,k1,k2,lmax
        logical                       :: rotate_CF=.false.

        double precision, allocatable :: D(:),E(:)
        double complex, allocatable   :: Tau(:),work(:)
        integer                       :: info,lwork

         if( iargc().eq.0)then
          write(*,*) 'CF Usage:'               
          write(*,*) '-JMult   : J Multiplicity of the ground state multiplet'
          write(*,*) '-CISIZE  : Size of CI basis set'               
          write(*,*) '-lmax    : Max Order of CF operators'               
          write(*,*) '-rot     : ZYZ Euler angoles (rad) for CF rotation'               
          stop
         endif

         do i=1,iargc()

          call getarg(i,word)
          
          select case (trim(word))

             case ('-JMult')
                 call getarg(i+1,word)
                 read(word,*) Nj

             case ('-CISIZE')
                 call getarg(i+1,word)
                 read(word,*) N

             case ('-lmax')
                 call getarg(i+1,word)
                 read(word,*) lmax

             case ('-rot')
                 call getarg(i+1,word)
                 read(word,*) alpha
                 call getarg(i+2,word)
                 read(word,*) beta
                 call getarg(i+3,word)
                 read(word,*) gamma
                 rotate_CF=.true.

          end select

         enddo

         allocate(SOC(N,N))
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

         call new_diag(N,SOC,EIG)

         write(*,*) '#################################################'
         write(*,*) '#################################################'
         write(*,*) 'First Nj eigevanlues of the electronic H'
         write(*,*) '#################################################'
         write(*,*) '#################################################'

         do i=1,Nj
          EIG2(i)=(EIG(i)-EIG(1))*219474.63
          write(*,*) 'EIG: ',i,' ENER: ',EIG2(i)
         enddo

         EIG=EIG2
        
         Jz=(0.0d0,0.0d0)
         Jx=(0.0d0,0.0d0)
         Jx2=(0.0d0,0.0d0)

         do i=1,Nj
          do j=1,Nj
           do k1=1,N
            do k2=1,N
               
             Jz(i,j)=Jz(i,j)&
                            +Sz(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j) &
                            +Lz(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j)

            enddo
           enddo
          enddo
         enddo

         do i=1,Nj
          do j=1,Nj
           do k1=1,N
            do k2=1,N
               
             Jx(i,j)=Jx(i,j)&
                            +Sx(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j) &
                            +Lx(k1,k2)*conjg(SOC(k1,i))*SOC(k2,j)

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

         if(rotate_CF)then
          call projH(lmax,Nj,EIG,Jz,O,alpha,beta,gamma)
         else
          call projH(lmax,Nj,EIG,Jz,O)
         endif

         call diagH(Nj,lmax,O)

        return
        end program CrystalField

        subroutine diagH(Hdim,lmax,O)
        use lapack_diag_simm
        use stevens_class
        integer                        :: lmax,Hdim,i,j,k,l,q,s
        double precision               :: j1,j2,jj
        double precision, allocatable  :: Ener(:),Jzeig(:)
        double complex, allocatable    :: Hmat(:,:),O(:),Jz(:,:)
        double complex                 :: val_tmp,val

         allocate(Hmat(Hdim,Hdim))         
         allocate(Ener(Hdim))         
         jj=7.5d0

         j1=-7.5d0
         do i=1,Hdim
          j2=-7.5d0
          do j=1,Hdim

           val=(0.0d0,0.0d0)
           s=1
           do l=2,lmax,2
            do q=-l,l
             val_tmp=(0.0d0,0.0d0)
             call stevens_mat_elem(l,q,jj,j1,jj,j2,val_tmp)
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

         allocate(Jz(Hdim,Hdim))
         allocate(Jzeig(Hdim))

         Jz=(0.0d0,0.0d0)
         Jzeig=0.0d0

         do i=1,Hdim
          do j=1,Hdim
           do l=1,Hdim            
            Jz(i,j)=Jz(i,j)+conjg(Hmat(l,i))*Hmat(l,j)*(-7.5d0+(l-1))
           enddo
          enddo
         enddo
          
         call new_diag(Hdim,Jz,Jzeig)

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

        subroutine projH(lmax,Nj,Ener,Jz,O,alpha,beta,gamma)
        use stevens_class
        use spinham_class 
        implicit none
        integer                        :: lmax,Odim,Nj
        integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf,ll
        double precision, allocatable  :: Ener(:)
        double precision               :: jj,q1,q2,avg_ener
        double complex, allocatable    :: Jz(:,:),B(:),A(:,:),O(:),work(:)
        double complex                 :: mat_elem
        type(OSItensor), allocatable   :: Os(:)
        double precision, optional     :: alpha,beta,gamma

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
         allocate(A(Nj*Nj,Odim))

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
           jj=7.5d0
           s=1
           do l=2,lmax,2
            do j=-l,l
             call stevens_mat_elem(l,j,jj,q1,jj,q2,mat_elem)
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
