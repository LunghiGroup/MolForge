        program mat2euler
        use lapack_diag_simm
        implicit none
        double precision :: mat(3,3),mat2(3,3),rot(3,3),eigenvals(3),prod(3),flip,newrot(3,3)
        double precision :: alpha,beta,gamma
        integer          :: x,y,z,s
        character(len=18) :: word,inp_file


         call getarg(1,inp_file)
         call getarg(2,word)

         open(13,file=inp_file)

         read(13,*) mat(1,:)
         read(13,*) mat(2,:)
         read(13,*) mat(3,:)

         close(13)

         rot=mat

         call new_diag(3,rot,eigenvals)
         
         x=1
         y=2
         z=3

         if(trim(word).eq.'min')then

          if( abs(eigenvals(1)-minval(eigenvals)).lt. 1.0e-6 )then
           z=1
           x=2
           y=3
          endif

          if( abs(eigenvals(2)-minval(eigenvals)).lt. 1.0e-6 )then
           z=2
           x=3
           y=1
          endif

          if( abs(eigenvals(3)-minval(eigenvals)).lt. 1.0e-6 )then
           z=3
           x=1
           y=2
          endif

         endif

         if(trim(word).eq.'max')then

          if( abs(eigenvals(1)-maxval(eigenvals)).lt. 1.0e-6 )then
           z=1
           x=2
           y=3
          endif

          if( abs(eigenvals(1)-maxval(eigenvals)).lt. 1.0e-6 )then
           z=2
           x=3
           y=1
          endif

          if( abs(eigenvals(1)-maxval(eigenvals)).lt. 1.0e-6 )then
           z=3
           x=1
           y=2
          endif

         endif

         newrot=rot
         rot(:,1)=newrot(:,x)
         rot(:,2)=newrot(:,y)
         rot(:,3)=newrot(:,z)

         prod(1)=rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2)
         prod(2)=rot(3,1)*rot(1,2)-rot(3,1)*rot(1,2)
         prod(3)=rot(1,1)*rot(2,2)-rot(2,1)*rot(1,2)
         flip=rot(1,3)*prod(1)+rot(2,3)*prod(2)+rot(3,3)*prod(3)

         if (flip.lt.0.d0)then
          write(6,*) 'Flipping the third eigenvector...'
          rot(:,3)=-1.0d0*rot(:,3)
         endif

         write(*,*) 'Eigenvalues:',eigenvals
         write(*,*) 'Eigenvector 1:',rot(:,1)
         write(*,*) 'Eigenvector 2:',rot(:,2)
         write(*,*) 'Eigenvector 3:',rot(:,3)

         open(13,file='rot.dat')

         write(13,*) rot(1,:)
         write(13,*) rot(2,:)
         write(13,*) rot(3,:)

         close(13)

         write(*,*) 'Rotation Matrix ---------'
         write(*,*) rot(1,:)
         write(*,*) rot(2,:)
         write(*,*) rot(3,:)
         write(*,*) '-------------------------'

         alpha=atan2(rot(2,3),rot(1,3))
         beta=acos(rot(3,3))
         gamma=atan2(-1.0d0*rot(3,2),rot(3,1))

         write(*,*) 'ZYZ Euler Angles:', alpha,beta,gamma
         
         mat2=0.0d0
         mat2=matmul(mat,rot)
         mat2=matmul(transpose(rot),mat2)

         open(13,file='mat2.dat')

         write(13,*) mat2(1,:)
         write(13,*) mat2(2,:)
         write(13,*) mat2(3,:)

         close(13)

         write(*,*) 'Rotated matrix ----------'
         write(*,*) mat2(1,:)
         write(*,*) mat2(2,:)
         write(*,*) mat2(3,:)
         write(*,*) '-------------------------'


        return
        end program mat2euler





