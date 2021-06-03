        module stevens_class
        implicit none


        contains

        subroutine stevens_mat_elem(k,q,j1,q1,j2,q2,mat_elem)
        use anglib
        implicit none
        integer             :: T1,M1,T,M,T2,M2,k,q
        double precision    :: j1,q1,j2,q2,a
        double complex      :: mat_elem

        mat_elem=(0.0,0.0)
        
        if(q.eq.0)then

        FLUSH(6)
         T2=NINT(2*j2)
         M2=NINT(2*q2)
         T=2*k
         M=2*q
         T1=NINT(2*j1)
         M1=NINT(2*q1)

         a=cleb(T2,M2,T,M,T1,M1)
         mat_elem=CMPLX(a,0,8)
         call spintens(k,j1,a)
         mat_elem=mat_elem*a/dsqrt(2*j1+1)

        endif


        if(q.gt.0)then

         T2=NINT(2*j2)
         M2=NINT(2*q2)
         T=2*k
         M=2*q
         T1=NINT(2*j1)
         M1=NINT(2*q1)

        if( (DABS(NINT(q/2.0d0)-q/2.0d0)) .lt. 1.0D-05)then
 
         a=cleb(T2,M2,T,M,T1,M1)
         M=-1*M
         a=a+cleb(T2,M2,T,M,T1,M1)
         mat_elem=CMPLX(a,0,8)
         call spintens(k,j1,a)
         mat_elem=mat_elem*a/dsqrt(4*j1+2)

        else

         a=-1*cleb(T2,M2,T,M,T1,M1)
         M=-1*M
         a=a+cleb(T2,M2,T,M,T1,M1)
         mat_elem=CMPLX(a,0,8)
         call spintens(k,j1,a)
         mat_elem=mat_elem*a/dsqrt(4*j1+2)

        endif

        endif


        if(q.lt.0)then

         T2=NINT(2*j2)
         M2=NINT(2*q2)
         T=2*k
         M=2*q
         T1=NINT(2*j1)
         M1=NINT(2*q1)

        if( (DABS(NINT(q/2.0d0)-q/2.0d0)) .lt. 1.0D-05)then

         a=cleb(T2,M2,T,M,T1,M1)
         M=-1*M
         a=a-cleb(T2,M2,T,M,T1,M1)
         mat_elem=CMPLX(0,a,8)
         call spintens(k,j1,a)
         mat_elem=mat_elem*a/dsqrt(4*j1+2)
        else

         a=cleb(T2,M2,T,M,T1,M1)
         M=-1*M
         a=a+cleb(T2,M2,T,M,T1,M1)
         mat_elem=CMPLX(0,a,8)
         call spintens(k,j1,a)
         mat_elem=mat_elem*a/dsqrt(4*j1+2)
        endif

        endif

        return
        end subroutine stevens_mat_elem


        subroutine spintens(K,S,B)
        use anglib
        implicit none
        integer          :: K,i
        DOUBLE PRECISION :: S,B

        B=1
       
        do i=-K,K
         B=B*(2*S+1+i)
        enddo
        
        B=B/((2**K)*factorial(2*K))

        B=dsqrt(B)
        B=B*factorial(K)

!        B=factorial(INT(2*S)+K+1)
!        B=B/factorial(2*K)
!        B=B/factorial(INT(2*S)-K)
!        B=dsqrt(B)
!        B=B*factorial(k)
!        B=B/(2**(K/2))


        return
        end

        end module stevens_class






