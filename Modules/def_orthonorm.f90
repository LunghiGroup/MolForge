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

         subroutine lowdin_S(N,basis)
         use lapack_diag_simm
         implicit none
         double precision, allocatable :: basis(:,:)
         double precision, allocatable :: Smat(:,:),Smat_inv(:,:),Seig(:)
         integer                       :: N,i

         allocate(basis(N,N))
         allocate(Smat(N,N))
         allocate(Smat_inv(N,N))
         allocate(Seig(N))

         Smat=matmul(transpose(basis),basis)

         call new_diag(N,Smat,Seig)

         Smat_inv=0.0d0        
         do i=1,N
          Smat_inv(i,i)=1.0d0/sqrt(Seig(i))
         enddo

         Smat_inv=matmul(Smat,Smat_inv)
         Smat_inv=matmul(Smat_inv,transpose(Smat))        

         basis=matmul(basis,transpose(Smat_inv))

         return
         end subroutine lowdin_S

        end module orthonorm_class
