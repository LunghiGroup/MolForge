        module orthonorm_class
        implicit none

        contains


         subroutine schmidt_ort(N,Amat)
         implicit none
         double precision                     :: Amat(N,N),Bmat(N,N)
         double precision                     :: norm
         integer                              :: N,i,j,s,t,jj

          ! norm
          do j=1,N
           norm=0.0d0
           do i=1,N
            norm=norm+Amat(j,i)*Amat(j,i)         
           enddo
           Amat(j,:)=Amat(j,:)/sqrt(norm)
          enddo

          Bmat=0.0d0
          Bmat(1,:)=Amat(1,:)

          do j=2,N

           Bmat(j,:)=Amat(j,:)         
          ! project out previous directions
           do jj=1,j-1
            norm=0.0d0
            do s=1,N
             norm=norm+Bmat(j,s)*Bmat(jj,s)
            enddo
            do s=1,N
             Bmat(j,s)=Bmat(j,s)-norm*Bmat(jj,s)
            enddo
           enddo

          ! norm

           norm=0.0d0
           do i=1,N
            norm=norm+Bmat(j,i)*Bmat(j,i)         
           enddo
           Bmat(j,:)=Bmat(j,:)/sqrt(norm)

          enddo

          do i=1,N
           do j=1,N
            norm=0.0d0
            do s=1,N
             norm=norm+Bmat(i,s)*Bmat(j,s)
            enddo
            if(abs(norm).gt.1e-12)then
            write(*,*) i,j,norm
            endif
           enddo
          enddo
          do i=1,N
            write(*,*) i,Bmat(90,i)
          enddo

!         stop
  
          Amat=Bmat

         return
         end subroutine schmidt_ort


         subroutine OrthoNorm_Set(n_state,AP)
         use lapack_diag_simm
         implicit none
         double complex,  dimension(n_state,n_state)             :: AP
         integer                                                 :: n_state
         integer                                                 :: s,t,i,j,ialloc
         DOUBLE COMPLEX                                          :: sum
         DOUBLE COMPLEX, allocatable, dimension(:,:)             :: Smat,SmatInv,Shalfmat,ShalfmatInv,X
         DOUBLE PRECISION, allocatable, dimension(:)             :: Shalf
         double precision, allocatable, dimension(:)             :: work
!
         ALLOCATE(Smat(n_state,n_state),stat=ialloc)
         ALLOCATE(SmatInv(n_state,n_state),stat=ialloc)
         ALLOCATE(Shalf(n_state),stat=ialloc)
         ALLOCATE(Shalfmat(n_state,n_state),stat=ialloc)
         ALLOCATE(Shalfmatinv(n_state,n_state),stat=ialloc)
         ALLOCATE(X(n_state,n_state),stat=ialloc)

         Smat=(0.0d0,0.0d0)

         do t=1,n_state
          sum=(0.0d0,0.0d0)
          do s=1,n_state
           sum=sum+CONJG(AP(s,t))*AP(s,t)
          enddo
          do s=1,n_state
           AP(s,t)=AP(s,t)/dsqrt(DBLE(sum))
          enddo
         enddo

         do t=1,n_state
          do s=1,n_state
           do i=1,n_state
            Smat(s,t)=Smat(s,t)+CONJG(AP(i,s))*AP(i,t)
           enddo
          enddo
         enddo

         call new_diag(n_state,Smat,Shalf)

         do i=1,n_state
          do j=1,n_state
           if(i.eq.j)then
            if (Shalf(i).lt.1.0D-08)then
             write(*,*) 'Smat is singular, try another Model Space'
             FLUSH(6)
             stop
            else
             Shalfmat(i,j)=CMPLX(dsqrt(Shalf(i)),0.0d0,8)
             ShalfmatInv(i,j)=CMPLX(1.0d0/dsqrt(Shalf(i)),0.0d0,8)
            endif
           else
            Shalfmat(i,j)=(0.0d0,0.0d0)
            ShalfmatInv(i,j)=(0.0d0,0.0d0)
           endif
          enddo
         enddo

         X=MATMUL(Smat,ShalfmatInv)
         Smat=CONJG(Smat)
         Smat=TRANSPOSE(Smat)
         X=MATMUL(X,Smat)

         AP=MATMUL(X,AP)

         return
         end subroutine OrthoNorm_Set

        end module orthonorm_class
