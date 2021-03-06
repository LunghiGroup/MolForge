        subroutine Rot_Wig(k,alpha,beta,gamma,rot)
        use anglib
        implicit none
        integer                                             :: v,t,s,mins,maxs,ialloc
        double precision                                    :: alpha,beta,gamma,a,c,k,djmn
        double complex                                      :: b
        double complex, dimension(NINT(2*k+1),NINT(2*k+1))  :: rot

        do v=1,NINT(2*k+1)
         do t=1,NINT(2*k+1)
          rot(v,t)=(0.0d0,0.0d0)
         enddo
        enddo


        do v=-NINT(2*k),NINT(2*k),2
         do t=-NINT(2*k),NINT(2*k),2

         mins=MAX( 0,NINT(t/2.-v/2.) )
         maxs=MIN( NINT(k+t/2.),NINT((k-v/2.)) )
 
         do s=mins,maxs
          a=factorial(NINT(k+t/2.-s))*factorial(s)*factorial(NINT(v/2.-t/2.+s))*factorial(NINT(k-v/2.-s))
          a=((-1)**(v/2.-t/2.+s))/a
          c=NINT(2*k-v/2.+t/2.-2*s)
          a=a*((dcos(beta/2.0d0))**c)
          c=NINT(v/2.-t/2.+2*s)
          a=a*((dsin(beta/2.0d0))**c)
          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))+CMPLX(a,0.0,8)
         enddo

          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*& 
          CMPLX(dsqrt(factorial(NINT(k+v/2.))*factorial(NINT(k-v/2.))*factorial(NINT(k+t/2.))*factorial(NINT(k-t/2.))),0.0,8)

         enddo
        enddo

        do v=-NINT(2*k),NINT(2*k),2
         do t=-NINT(2*k),NINT(2*k),2
          b=CMPLX(0.0,(-1*(t/2.0d0)*gamma),8)
          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*EXP(b)
          b=CMPLX(0.0,(-1*(v/2.0d0)*alpha),8)
          rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=EXP(b)*rot(NINT(v/2.+k+1),NINT(t/2.+k+1))
         enddo
        enddo

        return
        end subroutine Rot_Wig



        subroutine Rot_O(k,O,alpha,beta,gamma)
        use variables, only                   : stev_oper
        implicit none
        integer                              :: k,l,t
        double complex                       :: rot(2*k+1,2*k+1)
        double precision                     :: alpha,beta,gamma
        TYPE(stev_oper), dimension(-k:k)     :: O,newO

        call Rot_Wig(DBLE(k),alpha,beta,gamma,rot)
        
        do l=-k,k

         newO(l)%k=O(l)%k
         newO(l)%q=O(l)%q

         if(newO(l)%q.eq.0)then
          newO(l)%val=O(l)%val
         endif

         if(newO(l)%q.gt.0)then
          newO(l)%val=O(l)%val-CMPLX(-DIMAG(O(-l)%val),DBLE(O(-l)%val),8)
          newO(l)%val=newO(l)%val*((-1)**(O(l)%q))
          newO(l)%val=newO(l)%val/dsqrt(2.0d0)
         endif

         if(newO(l)%q.lt.0)then
          newO(l)%val=O(-l)%val+CMPLX(-DIMAG(O(l)%val),DBLE(O(l)%val),8)
          newO(l)%val=newO(l)%val/dsqrt(2.0d0)
         endif

        enddo

        O=newO
        newO%val=(0.0,0.0)

        do l=-k,k
         do t=-k,k
          newO(l)%val=newO(l)%val+(CONJG(rot(l+k+1,t+k+1))*O(t)%val)
         enddo
        enddo

        do l=-k,k
         if(O(l)%q.eq.0)then
          O(l)%val=newO(l)%val
         endif
         if(O(l)%q.gt.0)then
          O(l)%val=((-1)**((O(l)%q)))*newO(l)%val+newO(-l)%val
          O(l)%val=O(l)%val*dsqrt(2.0d0)/2
         endif
         if(O(l)%q.lt.0)then
          O(l)%val=((-1)**(ABS(O(l)%q)))*newO(-l)%val-newO(l)%val
          O(l)%val=O(l)%val*dsqrt(2.0d0)/2
          O(l)%val=CMPLX(-DIMAG(O(l)%val),DBLE(O(l)%val),8)
         endif
        enddo

        return
        end subroutine Rot_O



