        program main
        use lapack_diag_simm
        implicit none
        double precision :: Dmat(3,3),Deig(3),O(5),newO(5),D,E








         call new_diag(3,Dmat,Deig)


         D=Eeig(3)-0.5*(Deig(1)+Deig(2))
         E=0.5*(Deig(1)-Deig(2))

        if( dabs(E/D).lt. 0.333333)then
         X=1
         Y=2
         Z=3
        else

         D=Deig(1)-0.5*(Dieg(2)+Deig(3))
         E=0.5*(Deig(2)-Deig(3))

         if( dabs(E/D).lt. 0.333333)then
          X=2
          Y=3
          Z=1
         else
          X=3
          Y=1
          Z=2
         endif

        endif

        alpha=atan2(Dmat(2,z),Dmat(1,Z))
        beta=acos(Dmat(3,z))
        gamma=atan2(-Dmat(3,y),Dmat(3,x))




        return
        end program main
