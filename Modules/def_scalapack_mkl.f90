        module scalapack_diag_asimm_cmplx
        implicit none

        contains

        subroutine pzdiag2(N,A,W,VR,VL)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        integer                                      :: N,info,ilwork
        complex(8), allocatable                      :: tau(:)
        type(dist_cmplx_mat)                         :: Q,B,C,A
        type(dist_cmplx_mat)                         :: VL,VR,HC
        type(dist_cmplx_vec)                         :: V1,V2
        complex(8)                                   :: sum,val,W(N),WW(N)
        complex(8)                                   :: ZERO
        integer, allocatable                         :: iwork(:)
        integer                                      :: desca(9),descz(9),i,j,ILO,IHI
        logical, allocatable                         :: Sel(:)

         ILO=1 
         IHI=N

         !!!!!!!!!!!!!!!!!! PDGEBAL

!         call pdgebal('B',N,A%mat,A%desc,ILO,IHI,scale,info)
!         if(info.ne.0) write(*,*) 'pdgebal failed!'

         !!!!!!!!!!!!!!!!!! PDGEHRD

         allocate(tau(N))
         tau=(0.0d0,0.0d0)
         lwork=-1
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(1))

         call pzgehrd(N,ILO,IHI,A%mat,1,1,A%desc,tau,cwork,lwork,info)

         lwork=int(cwork(1))+100
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(lwork))

         call pzgehrd(N,ILO,IHI,A%mat,1,1,A%desc,tau,cwork,lwork,info)

         deallocate(cwork)
         if(info.ne.0) writE(*,*) 'pdgehrd failed'


         !!!!!!!!!!!!!!!!!! PDORMHR

         call Q%set(N,N,NB,MB)
         Q%mat=(0.0d0,0.0d0)
         do i=1,N                    
          call pzelset(Q%mat,i,i,Q%desc,(1.0d0,0.0d0))
         enddo

         lwork=-1
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(1))
         
         call pzunmhr('L','N',N,N,ILO,IHI,A%mat,1,1,A%desc,tau,Q%mat,1,1,Q%desc,cwork,lwork,info)

         lwork=int(cwork(1))+100
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(lwork))

         call pzunmhr('L','N',N,N,ILO,IHI,A%mat,1,1,A%desc,tau,Q%mat,1,1,Q%desc,cwork,lwork,info)

         deallocate(cwork)
         deallocate(tau)
         if(info.ne.0) writE(*,*) 'pdormhr failed'

         !!!!!!!!!!!!!!!!!! PDLASET

         ZERO=(0.0d0,0.0d0)
         CALL pzlaset( 'Lower triangular',N-2,N-2,ZERO,ZERO,A%mat,3,1,A%desc)

         !!!!!!!!!!!!!!!!!! PDHSEQR

!         ilwork=8000000
!         if(allocated(iwork)) deallocate(iwork)
!         allocate(iwork(ilwork))
!         lwork=8000000
!         if(allocated(work)) deallocate(work)
!         allocate(work(lwork))
         
!         call pdhseqr('S','V',N,ILO,IHI,A%mat,A%desc,WR,WI,Q%mat,Q%desc,work,lwork,iwork,ilwork,info)

!         deallocate(work)
!         deallocate(iwork)
!         if(info.ne.0) writE(*,*) 'pdhseqr failed'

         !!!!!!!!!!!!!!!!!! PZLAHQR

         call HC%set(N,N,NB,MB)
         HC%mat=A%mat
         VR%mat=Q%mat

         call A%dealloc()
         call Q%dealloc()

         lwork=8000000
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(lwork))

         ilwork=8000000
         if(allocated(iwork)) deallocate(iwork)
         allocate(iwork(ilwork))
         
         call pzlahqr(1,1,N,ILO,IHI,HC%mat,HC%desc,W,1,N,VR%mat,VR%desc,cwork,lwork,iwork,ilwork,info)

         deallocate(cwork)
         deallocate(iwork)
         if(info.ne.0) writE(*,*) 'pdlahqr failed'

         !!!!!!!!!!!!!!!!!! PZTREVC
         
         VL%mat=VR%mat

         lwork=8000000
         allocate(cwork(lwork))
         allocate(rwork(lwork))
         rwork=0.0d0
         cwork=(0.0d0,0.0d0)

         call pztrevc('B','B',Sel,N,HC%mat,HC%desc,VL%mat,VL%desc,VR%mat,VR%desc,N,N,cwork,rwork,info)
         if(info.ne.0) writE(*,*) 'pztrevc failed'

         call HC%dealloc()

         deallocate(cwork)
         deallocate(rwork)

         do i=1,N
          sum=(0.0d0,0.0d0)
          do j=1,N
           call pzelget('A',' ',val,VR%mat,j,i,VR%desc)
           sum=sum+val*conjg(val)
          enddo
          sum=sqrt(sum)
          do j=1,N
           call pzelget('A',' ',val,VR%mat,j,i,VR%desc)
           val=val/sum
           call pzelset(VR%mat,j,i,VR%desc,val)
          enddo
         enddo

         WW=(0.0d0,0.0d0)

         do i=1,N
          call pzdotc(N,WW(i),VR%mat,1,i,VR%desc,1,VL%mat,1,i,VL%desc,1)
         enddo

         do i=1,N
          do j=1,N
           call pzelget('A',' ',val,VL%mat,j,i,VL%desc)
           val=val/conjg(WW(i))
           call pzelset(VL%mat,j,i,VL%desc,val)           
          enddo
         enddo

         if(info.ne.0)then
          write(*,*) 'pzdiag2 failed!'
          stop
         endif

        return
        end subroutine pzdiag2

        end module scalapack_diag_asimm_cmplx

        module scalapack_diag_asimm
        implicit none

        contains

        subroutine pddiag2(N,A,W,VR,VL)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        integer                                      :: N,info,ilwork
        double precision, allocatable                :: tau(:)
        type(dist_dbl_mat)                           :: Q,A,B,C
        type(dist_cmplx_mat)                         :: VL,VR,HC
        type(dist_cmplx_vec)                         :: V1,V2
        complex(8)                                   :: sum,val,W(N),WW(N)
        double precision                             :: WR(N),WI(N),ZERO,SCALE(N)
        double precision                             :: VRL(N,N),VRR(N,N)
        integer, allocatable                         :: iwork(:)
        integer                                      :: desca(9),descz(9),i,j,ILO,IHI
        logical, allocatable                         :: Sel(:)

         ILO=1 
         IHI=N

         !!!!!!!!!!!!!!!!!! PDGEBAL

!         call pdgebal('B',N,A%mat,A%desc,ILO,IHI,scale,info)
!         if(info.ne.0) write(*,*) 'pdgebal failed!'

         !!!!!!!!!!!!!!!!!! PDGEHRD

         allocate(tau(N))
         tau=0.0d0
         lwork=-1
         if(allocated(work)) deallocate(work)
         allocate(work(1))

         call pdgehrd(N,ILO,IHI,A%mat,1,1,A%desc,tau,work,lwork,info)

         lwork=int(work(1))+100
         if(allocated(work)) deallocate(work)
         allocate(work(lwork))

         call pdgehrd(N,ILO,IHI,A%mat,1,1,A%desc,tau,work,lwork,info)

         deallocate(work)
         if(info.ne.0) writE(*,*) 'pdgehrd failed'


         !!!!!!!!!!!!!!!!!! PDORMHR

         call Q%set(N,N,NB,MB)
         Q%mat=0.0d0
         do i=1,N                    
          call pdelset(Q%mat,i,i,Q%desc,1.0d0)
         enddo

         lwork=-1
         if(allocated(work)) deallocate(work)
         allocate(work(1))
         
         call pdormhr('L','N',N,N,ILO,IHI,A%mat,1,1,A%desc,tau,Q%mat,1,1,Q%desc,work,lwork,info)

         lwork=int(work(1))+100
         if(allocated(work)) deallocate(work)
         allocate(work(lwork))

         call pdormhr('L','N',N,N,ILO,IHI,A%mat,1,1,A%desc,tau,Q%mat,1,1,Q%desc,work,lwork,info)

         deallocate(work)
         deallocate(tau)
         if(info.ne.0) writE(*,*) 'pdormhr failed'

         !!!!!!!!!!!!!!!!!! PDLASET

         ZERO=0.0d0
         CALL pdlaset( 'Lower triangular',N-2,N-2,ZERO,ZERO,A%mat,3,1,A%desc)


         !!!!!!!!!!!!!!!!!! PDHSEQR

!         ilwork=8000000
!         if(allocated(iwork)) deallocate(iwork)
!         allocate(iwork(ilwork))
!         lwork=8000000
!         if(allocated(work)) deallocate(work)
!         allocate(work(lwork))
         
!         call pdhseqr('S','V',N,ILO,IHI,A%mat,A%desc,WR,WI,Q%mat,Q%desc,work,lwork,iwork,ilwork,info)

!         deallocate(work)
!         deallocate(iwork)
!         if(info.ne.0) writE(*,*) 'pdhseqr failed'

         !!!!!!!!!!!!!!!!!! PZLAHQR

         call HC%set(N,N,NB,MB)
         HC%mat=A%mat
         VR%mat=Q%mat

         call A%dealloc()
         call Q%dealloc()

         lwork=8000000
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(lwork))

         ilwork=8000000
         if(allocated(iwork)) deallocate(iwork)
         allocate(iwork(ilwork))
         
         call pzlahqr(1,1,N,ILO,IHI,HC%mat,HC%desc,W,1,N,VR%mat,VR%desc,cwork,lwork,iwork,ilwork,info)

         deallocate(cwork)
         deallocate(iwork)
         if(info.ne.0) writE(*,*) 'pdlahqr failed'

         !!!!!!!!!!!!!!!!!! PZTREVC
         
         VL%mat=VR%mat

         lwork=8000000
         allocate(cwork(lwork))
         allocate(rwork(lwork))
         rwork=0.0d0
         cwork=(0.0d0,0.0d0)

         call pztrevc('B','B',Sel,N,HC%mat,HC%desc,VL%mat,VL%desc,VR%mat,VR%desc,N,N,cwork,rwork,info)
         if(info.ne.0) writE(*,*) 'pztrevc failed'

         call HC%dealloc()

         deallocate(cwork)
         deallocate(rwork)

         do i=1,N
          sum=(0.0d0,0.0d0)
          do j=1,N
           call pzelget('A',' ',val,VR%mat,j,i,VR%desc)
           sum=sum+val*conjg(val)
          enddo
          sum=sqrt(sum)
          do j=1,N
           call pzelget('A',' ',val,VR%mat,j,i,VR%desc)
           val=val/sum
           call pzelset(VR%mat,j,i,VR%desc,val)
          enddo
         enddo

         WW=(0.0d0,0.0d0)

         do i=1,N
          call pzdotc(N,WW(i),VR%mat,1,i,VR%desc,1,VL%mat,1,i,VL%desc,1)
         enddo

         do i=1,N
          do j=1,N
           call pzelget('A',' ',val,VL%mat,j,i,VL%desc)
           val=val/conjg(WW(i))
           call pzelset(VL%mat,j,i,VL%desc,val)           
          enddo
         enddo

         if(info.ne.0)then
          write(*,*) 'pddiag2 failed!'
          stop
         endif

        return
        end subroutine pddiag2

        end module scalapack_diag_asimm


        module scalapack_inv
        implicit none

        contains

        subroutine pzgeinv(N,A)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        integer                                      :: N,info
        type(dist_cmplx_mat)                         :: A
        integer, allocatable                         :: iwork(:),IPIV(:)


         allocate(IPIV(N))
         call pzgetrf(N,N,A%mat,1,1,A%desc,IPIV,info)
         if(info.ne.0) write(*,*) 'pzgetrf failed'

         if(allocated(cwork)) deallocate(cwork)
         if(allocated(iwork)) deallocate(iwork)
!         allocate(cwork(1))
!         allocate(iwork(1))
!         lwork=-1
!         lrwork=-1

!         call pzgetri(N,A%mat,1,1,A%desc,IPIV,cwork,lwork,iwork,lrwork,info)

!         if(allocated(cwork)) deallocate(cwork)
!         if(allocated(iwork)) deallocate(iwork)
         lwork=100*N
         lrwork=100*N
         allocate(cwork(lwork))
         allocate(iwork(lrwork))

         call pzgetri(N,A%mat,1,1,A%desc,IPIV,cwork,lwork,iwork,lrwork,info)

         if(info.ne.0) writE(*,*) 'pzgetri failed!'

         deallocate(IPIV)
         deallocate(cwork)
         deallocate(iwork)


        return
        end subroutine pzgeinv

        end module scalapack_inv



        module scalapack_diag_simm
        implicit none

        contains

        subroutine pzdiag(N,A,W,Z,desca,descz)
        implicit none
        integer                                      :: N,info,lwork,lrwork
        complex(8), allocatable                      :: A(:,:),Z(:,:)
        double precision                             :: W(N)
        double precision, allocatable                :: work(:)
        complex(8), allocatable                      :: cwork(:)!,crwork(:)
        double precision, allocatable                :: crwork(:)
        integer                                      :: desca(9),descz(9)
        integer                                      :: i,j
        complex(8)                                   :: val,sum

         lwork=-1
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(1))
         lrwork=-1
         if(allocated(crwork)) deallocate(crwork)
         allocate(crwork(1))

         call pzheev('V','U',N,A,1,1,desca,W,Z,1,1,descz,cwork,lwork,crwork,lrwork,info)

         lwork=2*int(cwork(1))+100
         if(allocated(cwork)) deallocate(cwork)
         allocate(cwork(lwork))
         lrwork=2*int(crwork(1))+100
         if(allocated(crwork)) deallocate(crwork)
         allocate(crwork(lrwork))

         call pzheev('V','U',N,A,1,1,desca,W,Z,1,1,descz,cwork,lwork,crwork,lrwork,info)
         
!         do i=1,N
!          sum=(0.0d0,0.0d0)
!          do j=1,N
!           call pzelget('A',' ',val,Z,j,i,descz)
!           sum=sum+val*conjg(val)
!          enddo
!          sum=sqrt(sum)
!          do j=1,N
!           call pzelget('A',' ',val,Z,j,i,descz)
!           val=val/sum
!           call pzelset(Z,j,i,descz,val)
!          enddo
!         enddo

         if(info.ne.0)then
          write(*,*) 'pzheev diag failed!'
          stop
         endif

        return
        end subroutine pzdiag

        end module scalapack_diag_simm
