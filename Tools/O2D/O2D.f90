        program D2O
        double precision   :: D(3,3)
        double precision   :: O(5)         

         open(12,file='O2.dat')

         read(12,*) O(1)
         read(12,*) O(2)
         read(12,*) O(3)
         read(12,*) O(4)
         read(12,*) O(5)

         close(12)

         D=0.0d0

         D(3,3)=DBLE(O(3))*dsqrt(6d0)/3
         D(1,3)=DBLE(O(4))/dsqrt(2d0)     
         D(3,1)=D(1,3)
         D(2,3)=DBLE(O(2))/dsqrt(2d0)          
         D(3,2)=D(2,3)
         D(2,2)=-D(3,3)-DBLE(O(5))*dsqrt(2.d0)
         D(2,2)=D(2,2)/2.d0         
         D(1,2)=DBLE(O(1))/dsqrt(2.d0)
         D(2,1)=D(1,2)

         D(1,1)=-D(2,2)-D(3,3)

         write(*,*) D(1,1),D(1,2),D(1,3)
         write(*,*) D(2,1),D(2,2),D(2,3)
         write(*,*) D(3,1),D(3,2),D(3,3)


        return
        end program D2O


        subroutine traceless(N,A)
        implicit none
        integer                          :: N,s,t
        double precision                 :: trace
        double precision, DIMENSION(N,N) :: A

         trace=0.0d0
         do s=1,N
          trace=trace+A(s,s)
         enddo
         do s=1,N
          A(s,s)=A(s,s)-(trace/3.0d0)
         enddo

        return
        end subroutine traceless
        
        subroutine Dparam(A,D,E)
        implicit none
        integer                         :: s
        double precision                :: D,E
        double precision, DIMENSION(3)  :: A

         D=A(1)-0.5d0*(A(2)+A(3))
         E=abs(0.5d0*(A(2)-A(3)))
         if(abs(E/D).lt.0.333333)then
         return
         endif

         D=A(2)-0.5d0*(A(3)+A(1))
         E=abs(0.5d0*(A(3)-A(1)))
         if(abs(E/D).lt.0.333333)then
         return
         endif

         D=A(3)-0.5d0*(A(1)+A(2))
         E=abs(0.5d0*(A(1)-A(2)))
         if(abs(E/D).lt.0.333333)then
         return
         endif

         write(6,*) 'subroutine Dparam failed to calculate ', &
                    'D and E params'
         stop
        return
        end subroutine Dparam

