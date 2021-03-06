        program cell2abc
        IMPLICIT NONE
        double precision :: at(3,3)
        double precision :: a,b,c, cosab, cosac, cosbc ,pi 
        double precision :: norm1, norm2, norm3
        character(len=100) :: cell_in,cell_out

         if (iargc().eq.0)then
          write(*,*) "cell2abc Usage:"       
          write(*,*) "cell_file   : input file with 3x3 cell"
          write(*,*) "abc_file    : 1-line file with a,b,c,alpha,beta,gamma (Ang&Deg)"       
          stop
         endif        

         call getarg(1,cell_in)
         call getarg(2,cell_out)

         open(10,file=trim(cell_in))
         open(11,file=trim(cell_out))

         read(10,*) at(1,1),at(1,2),at(1,3)
         read(10,*) at(2,1),at(2,2),at(2,3)
         read(10,*) at(3,1),at(3,2),at(3,3)
  
         norm1=sqrt(at(1,1)**2+at(1,2)**2+at(1,3)**2)
         norm2=sqrt(at(2,1)**2+at(2,2)**2+at(2,3)**2)
         norm3=sqrt(at(3,1)**2+at(3,2)**2+at(3,3)**2)
  
         a=norm1
         b=norm2
         c=norm3
  
         cosab=(at(1,1)*at(2,1)+at(1,2)*at(2,2)+at(1,3)*at(2,3))/norm1/norm2
         cosac=(at(1,1)*at(3,1)+at(1,2)*at(3,2)+at(1,3)*at(3,3))/norm1/norm3
         cosbc=(at(3,1)*at(2,1)+at(3,2)*at(2,2)+at(3,3)*at(2,3))/norm3/norm2

         pi=acos(-1.0d0)

         write(11,*) a,b,c,acos(cosbc)*180.0d0/pi,acos(cosac)*180.0d0/pi,acos(cosab)*180.0d0/pi

        return
        end program cell2abc
