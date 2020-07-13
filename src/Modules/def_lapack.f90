        module lapack_inverse
        implicit none

        contains

        subroutine mat_inv(X,N)
        implicit none
        integer                                         :: info,lwork,N,ialloc
        integer, dimension(N)                           :: ipiv
        double precision,  dimension(N,N)               :: X_inv
        double complex,  dimension(N,N)                 :: cX_inv
        CLASS(*), intent(in), dimension(N,N)            :: X
        double precision, allocatable, dimension(:)     :: work


        select type(X)

        type is (complex(8))

        cX_inv=X

        lwork=(N)**2
        allocate(work(lwork),stat=ialloc)
        call zgetrf(N,N,cX_inv,N,IPIV,info)
        call zgetri(N,cX_inv,N,IPIV,work,lwork,info)

        X=cX_inv

        type is (double precision)

        X_inv=X

        lwork=(N)**2
        allocate(work(lwork),stat=ialloc)
        call dgetrf(N,N,X_inv,N,IPIV,info)
        call dgetri(N,X_inv,N,IPIV,work,lwork,info)

        X=X_inv

        end select

        return
        end subroutine mat_inv

        end module lapack_inverse

        module lapack_diag_simm
        implicit none

        contains

        subroutine new_diag(N,A,W)
        implicit none
        integer                                      :: i,j,INFO,s,k,N,ialloc
        integer                                      :: l,inf,infr
        CLASS(*),  dimension(N,N)                    :: A
        DOUBLE COMPLEX,    DIMENSION(N*(N+1)/2)      :: AP
        DOUBLE PRECISION,  DIMENSION(N)              :: W
        DOUBLE PRECISION, ALLOCATABLE,  DIMENSION(:) :: work
        DOUBLE COMPLEX,   ALLOCATABLE,  DIMENSION(:) :: cwork


        select type(A)

        type is (double precision)


        l=(3*N-1)
        allocate(work(l),stat=ialloc)
        call dsyev('V','U',N,A,N,W,work,l,infr)
        deallocate(work)
        if(infr.ne.0)then
        write(*,*) 'dsyev diagonalization failed'
        FLUSH(6)
        stop
        endif

        type is (complex(8))

        s=1
        do j=1,N
         do i=1,j
          AP(s)=A(i,j)
          s=s+1
         enddo
        enddo

        l=(2*N-1)+1000
        allocate(cwork(l),stat=ialloc)
        k=(3*N-2)+1000
        allocate(work(k),stat=ialloc)
        call zhpev('V','U',N,AP,W,A,N,cwork,work,inf)
        deallocate(cwork)
        deallocate(work)
        if(inf.ne.0)then
         write(*,*) 'zhpev diagonalization failed'
         stop
        endif

        end select

        return
        end subroutine new_diag

        end module lapack_diag_simm

        module lapack_diag_asimm
        implicit none

        contains

        subroutine new_diag2(N,A,W)
        implicit none
        integer                                      :: i,j,INFO,s,k,N,ialloc
        integer                                      :: l,inf,infr
        CLASS(*),  dimension(N,N)                    :: A
        DOUBLE COMPLEX,    DIMENSION(N*(N+1)/2)      :: AP
        DOUBLE COMPLEX,    DIMENSION(N,N)            :: VRI,VLI
        DOUBLE COMPLEX,    DIMENSION(N)              :: W
        DOUBLE PRECISION,  DIMENSION(N)              :: WR,WI
        DOUBLE PRECISION,  DIMENSION(N,N)            :: VL,VR
        DOUBLE PRECISION, ALLOCATABLE,  DIMENSION(:) :: work
        DOUBLE COMPLEX,   ALLOCATABLE,  DIMENSION(:) :: cwork


        select type(A)

        type is (double precision)


         l=4*N
         allocate(work(l),stat=ialloc)
         call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,l,infr)
         deallocate(work)
         if(infr.ne.0)then
          write(*,*) 'dgeev diagonalization failed'
          FLUSH(6)
          stop
         endif

         do j=1,N
          W(j)=CMPLX(WR(j),WI(j),8)
         enddo

         A=VR

        type is (complex(8))

         l=6*N
         allocate(cwork(l),stat=ialloc)
         allocate(work(l),stat=ialloc)
         call zgeev('V','V',N,A,N,W,VLI,N,VRI,N,cwork,l,work,infr)
         deallocate(cwork)
         deallocate(work)
         if(infr.ne.0)then
          write(*,*) 'zgeev diagonalization failed'
         endif
         A=VRI        

        end select

        return
        end subroutine new_diag2

        end module lapack_diag_asimm

