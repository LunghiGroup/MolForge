        module rotations_class
        implicit none


        contains
        
        subroutine Rot_Wig(k,alpha,beta,gamma,rot)
        use anglib
        implicit none
        integer                                             :: v,t,s,mins,maxs,ialloc
        double precision                                    :: alpha,beta,gamma,a,c,k,djmn
        complex(8)                                          :: b
        complex(8), allocatable                             :: rot(:,:)


        if(allocated(rot)) deallocate(rot)
        allocate(rot(NINT(2*k+1),NINT(2*k+1)))
        rot=(0.0d0,0.0d0)

        do v=-NINT(2*k),NINT(2*k),2
         do t=-NINT(2*k),NINT(2*k),2

         mins=MAX( 0,NINT(t/2.-v/2.) )
         maxs=MIN( NINT(k+t/2.),NINT((k-v/2.)) )
 
         do s=mins,maxs
          a=factorial(NINT(k+t/2.0d0-s))*factorial(s)*factorial(NINT(v/2.-t/2.+s))*factorial(NINT(k-v/2.-s))
          a=((-1.0d0)**(v/2.0d0-t/2.0d0+s))/a
          c=NINT(2.0d0*k-v/2.0d0+t/2.0d0-2.0d0*s)
          a=a*((dcos(beta/2.0d0))**c)
          c=NINT(v/2.0d0-t/2.0d0+2.0d0*s)
          a=a*((dsin(beta/2.0d0))**c)
          rot(NINT(v/2.+k+1),NINT(t/2.0d0+k+1))=rot(NINT(v/2.0d0+k+1),NINT(t/2.0d0+k+1))+CMPLX(a,0.0d0,8)
         enddo

          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*& 
          CMPLX(dsqrt(factorial(NINT(k+v/2.))*factorial(NINT(k-v/2.))*factorial(NINT(k+t/2.))*factorial(NINT(k-t/2.))),0.0d0,8)

         enddo
        enddo

        do v=-NINT(2*k),NINT(2*k),2
         do t=-NINT(2*k),NINT(2*k),2
          b=CMPLX(0.0d0,(-1.0d0*(t/2.0d0)*gamma),8)
          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*EXP(b)
          b=CMPLX(0.0d0,(-1.0d0*(v/2.0d0)*alpha),8)
          rot(NINT(v/2.0d0+k+1.0d0),NINT(t/2.0d0+k+1.0d0))=EXP(b)*rot(NINT(v/2.+k+1),NINT(t/2.+k+1))
         enddo
        enddo

        rot=transpose(rot)

        return
        end subroutine Rot_Wig


        end module rotations_class

