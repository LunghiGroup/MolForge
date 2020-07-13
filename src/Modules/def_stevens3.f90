        module stevens3_class
        implicit none

        contains

        subroutine spherical_mat_elem(k,q,j1,q1,j2,q2,mat_elem)
        use anglib
        integer             :: T1,M1,T,M,T2,M2,k,q
        double precision    :: j1,q1,j2,q2,redelem
        double complex      :: mat_elem

         mat_elem=(0.0d0,0.0d0)

         T2=NINT(2*j2)
         M2=NINT(2*q2)
         T=2*k
         M=2*q
         T1=NINT(2*j1)
         M1=NINT(2*q1)
          
         mat_elem=cleb(T2,M2,T,M,T1,M1)/sqrt(2*j1+1)
         call sto_redelem(k,j1,redelem)
         mat_elem=mat_elem*redelem

        return
        end subroutine spherical_mat_elem

        subroutine sto_redelem(k,j,matelem)
        use anglib
        implicit none
        double precision       :: j,matelem,val
        integer                :: k

         matelem=factorial(nint(2*j+k+1))
         matelem=matelem/factorial(2*k)
         matelem=matelem/factorial(nint(2*j-k))
         matelem=sqrt(matelem)
         matelem=(-1.0d0)**(k)*factorial(k)
         matelem=matelem*(-1.0d0)**(k)*sqrt(factorial(2*k))
         matelem=matelem/(2.0d0)**(k)
         matelem=matelem/factorial(k)

        return
        end subroutine sto_redelem

        subroutine stevens_mat_elem(k,q,j1,q1,j2,q2,mat_elem)
        use anglib
        implicit none
        integer             :: T1,M1,T,M,T2,M2,k,q
        double precision    :: j1,q1,j2,q2,a
        double complex      :: mat_elem,mat_elem_1,mat_elem_2

        mat_elem=(0.0d0,0.0d0)

        T2=NINT(2*j2)
        M2=NINT(2*q2)
        T=2*k
        M=2*q
        T1=NINT(2*j1)
        M1=NINT(2*q1)
        
        if(q.eq.0)then
         call spherical_mat_elem(k,q,j1,q1,j2,q2,mat_elem)
        endif

        if(q.gt.0)then          
         call spherical_mat_elem(k,q,j1,q1,j2,q2,mat_elem_1)
         call spherical_mat_elem(k,-1*q,j1,q1,j2,q2,mat_elem_2)
         mat_elem=((-1.0d0)**(q)*mat_elem_1+mat_elem_2)/sqrt(2.0d0)
        endif

        if(q.lt.0)then          
         call spherical_mat_elem(k,q,j1,q1,j2,q2,mat_elem_1)
         call spherical_mat_elem(k,-1*q,j1,q1,j2,q2,mat_elem_2)
         mat_elem=((-1.0d0)**(q+1)*mat_elem_1+mat_elem_2)/sqrt(2.0d0)
         mat_elem=cmplx(-1.0d0*aimag(mat_elem),dble(mat_elem),8)
        endif


        return
        end subroutine stevens_mat_elem

        end module stevens3_class






