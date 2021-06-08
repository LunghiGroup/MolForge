      subroutine xyztoxray(cell,hinv)
      implicit none
      double precision      :: cell(6)
      double precision      :: deg2rad,pi,alpha,beta,gamma,a,b,c,V,hinv(3,3)

      
      pi=dacos(-1.0d0)
      deg2rad=pi/180.0d0

      a=cell(1)
      b=cell(2)
      c=cell(3)
      alpha=cell(4)*deg2rad
      beta=cell(5)*deg2rad
      gamma=cell(6)*deg2rad



      V=1-(cos(alpha)**2)
      V=V-(cos(beta)**2)
      V=V-(cos(gamma)**2)
      V=V+(2*cos(alpha)*cos(beta)*cos(gamma))
      V=dsqrt(V)

        write(*,*) alpha,beta,gamma,V
     
      hinv(1,1)=1/a
      hinv(1,2)=-cos(gamma)/(a*sin(gamma))
      hinv(1,3)=(cos(alpha)*cos(gamma)-cos(beta))/(a*V*sin(gamma))
      hinv(2,1)=0.0d0
      hinv(2,2)=1/(b*sin(gamma))
      hinv(2,3)=(cos(beta)*cos(gamma)-cos(alpha))/(b*V*sin(gamma))
      hinv(3,1)=0.0d0
      hinv(3,2)=0.0d0
      hinv(3,3)=sin(gamma)/(c*V)

      return
      end subroutine


      subroutine xraytoxyz(cell,h)
      implicit none
      double precision      :: cell(6)
      double precision      :: deg2rad,pi,alpha,beta,gamma,a,b,c,V,h(3,3)
      
      pi=acos(-1.0d0)
      deg2rad=pi/180.0d0

      a=cell(1)
      b=cell(2)
      c=cell(3)
      alpha=cell(4)*deg2rad
      beta=cell(5)*deg2rad
      gamma=cell(6)*deg2rad


      V=1-(cos(alpha))**2
      V=V-(cos(beta))**2
      V=V-(cos(gamma))**2
      V=V+(2*cos(alpha)*cos(beta)*cos(gamma))
      V=dsqrt(V)
      h(1,1)=a
      h(1,2)=b*cos(gamma)
      h(1,3)=c*cos(beta)
      h(2,1)=0
      h(2,2)=b*sin(gamma)
      h(2,3)=c*(cos(alpha)-(cos(beta)*cos(gamma)))/sin(gamma)
      h(3,1)=0
      h(3,2)=0
      h(3,3)=c*V/sin(gamma)


      return
      end subroutine

      subroutine dist_pbc(v1,v2,h,hinv,v,nrep)
      implicit none
      double precision                :: v1(3),v2(3),h(3,3),hinv(3,3)
      double precision                :: v(3),a(3),b(3)
      integer                         :: i,nrep


        forall(i=1:3)
         v(i)=0.0d0
        end forall

        a=v1
        b=v2
        a=MATMUL(hinv,a)
        b=MATMUL(hinv,b)

        v(1)=a(1)-b(1)
        v(2)=a(2)-b(2)
        v(3)=a(3)-b(3)

!        write(*,*) '--'
!        write(*,*) v

        do i =1,3
         v(i)=NINT(v(i)/nrep)
        enddo


        v(1)=a(1)-b(1)-v(1)*(nrep)
        v(2)=a(2)-b(2)-v(2)*(nrep)
        v(3)=a(3)-b(3)-v(3)*(nrep)


        v=MATMUL(h,v)

        return
        end subroutine dist_pbc


