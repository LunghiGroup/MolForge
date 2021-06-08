      program abc2cell
      implicit none
      integer           :: i
      double precision  :: a,b,c,alpha,beta,gamma,V,d(3,3),pi
      character(len=20) :: cellin,cellout


      if (iargc().eq.0)then
       write(*,*) "abc2cell Usage:"       
       write(*,*) "abc_file    : 1-line file with a,b,c,alpha,beta,gamma (Ang&Deg)"       
       write(*,*) "cell_file   : output file for cell"
      stop
      endif        

      call getarg(1,cellin)
      call getarg(2,cellout)

      open(11,file=cellin)
      open(12,file=cellout)

      read(11,*) a,b,c,alpha,beta,gamma

      pi=acos(-1.0d0)

      alpha=alpha*pi/180.0d0
      beta=beta*pi/180.0d0
      gamma=gamma*pi/180.0d0

      V=1-(cos(alpha)**2)
      V=V-(cos(beta)**2)
      V=V-(cos(gamma)**2)
      V=V+(2*cos(alpha)*cos(beta)*cos(gamma))
      V=dsqrt(V)

      d(1,1)=a
      d(1,2)=b*cos(gamma)
      d(1,3)=c*cos(beta)
      d(2,1)=0
      d(2,2)=b*sin(gamma)
      d(2,3)=c*(cos(alpha)-(cos(beta)*cos(gamma)))/sin(gamma)
      d(3,1)=0
      d(3,2)=0
      d(3,3)=c*V/sin(gamma)

      do i=1,3
       write(12,*) d(:,i)
      enddo

      return
      end program abc2cell
