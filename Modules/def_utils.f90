        module units_parms
        implicit none

         double precision, parameter ::  pi=3.141592653589793 
         double precision, parameter ::  bohr2ang=0.52917721067 
         double precision, parameter ::  hartree2wavenum=4.5563D-6
         double precision, parameter ::  kboltz=0.6950347291 ! cm-1/K
         double precision, parameter ::  hplank=33.3571775619  ! cm-1*ps
         double precision, parameter ::  mag_vacuum=1.98643 ! T^2*A^3/cm-1 , = mu0/4pi = 10^-7 T^2*m^3/J

         contains

         subroutine order_array(A)
         implicit none
         double precision, allocatable    :: A(:),B(:)
         logical, allocatable             :: flag(:)
         integer                          :: i

          allocate(flag(size(A)))
          allocate(B(size(A)))

          flag=.true.
          B=0.0d0

          do i=1,size(A)          
           B(i)=minval(A,flag)
           flag(minloc(A,flag))=.false.
          enddo       

          A=B
          deallocate(B)

         return
         end subroutine order_array

         function bose(T,ener) result(stat)
         implicit none
         double precision :: T,ener,stat
          stat=dexp(ener/(T*kboltz))-1
          stat=1.0d0/stat
         return
         end function bose

         function pval(ener,N) result(val)
         implicit none
         double precision :: ener,N,val
          
          val=ener/(ener**2+N**2)

         return
         end function pval

         function d_pval(ener,N) result(val)
         implicit none
         double precision :: ener,N,val
          
          val=1/(ener**2+N**2)-2*ener**2/(ener**2+N**2)**2

         return
         end function d_pval

         function deltaL(ener,N) result(val)
         implicit none
         double precision :: ener,N,val

!!!      LORENTZIAN LINEWIDTH
          val=N/(ener**2+N**2)
          val=val/dacos(-1.0d0)

         return
         end function deltaL

         function d_deltaL(ener,N) result(val)
         implicit none
         double precision :: ener,N,val

!!!      LORENTZIAN LINEWIDTH
          val=-2*N*ener/(ener**2+N**2)**2
          val=val/dacos(-1.0d0)

         return
         end function d_deltaL

         function deltaG(ener,N) result(val)
         implicit none
         double precision :: ener,N,val

!!!      GAUSSIAN LINEWIDTH
          val=N*sqrt(pi)
          val=1.0d0/val
          val=val*dexp(-(ener**2/N**2))

         return
         end function deltaG

         function d_deltaG(ener,N) result(val)
         implicit none
         double precision :: ener,N,val

!!!      GAUSSIAN LINEWIDTH
          val=N*sqrt(pi)
          val=1.0d0/val
          val=val*dexp(-(ener**2/N**2))
          val=val*(-2/N**2)

         return
         end function d_deltaG

         function delta(Dtype,ener,N) result(val)
         implicit none
         double precision :: ener,N,val
         integer          :: Dtype
          
          if(Dtype.eq.1) val=deltaG(ener,N)
          if(Dtype.eq.0) val=deltaL(ener,N)

         return
         end function delta

         function d_delta(Dtype,ener,N) result(val)
         implicit none
         double precision :: ener,N,val
         integer          :: Dtype
          
          if(Dtype.eq.1) val=d_deltaG(ener,N)
          if(Dtype.eq.0) val=d_deltaL(ener,N) 

         return
         end function d_delta

         function deltaN(ener,N,NN) result(val)
         implicit none
         double precision :: ener,val,N,A,x
         integer          :: NN,i

          x=ener/(N*sqrt(2.0d0))
          val=0.0d0

          do i=0,NN
           A=((-1)**i)/(fact(i)*(4**i)*sqrt(pi))
           val=val+hermite(x,2*i)*exp(-x**2)*A
          enddo

          val=val/(N*sqrt(2.0d0))

         return
         end function deltaN

         function hermite(x,i) result(val)
         implicit none
         double precision        :: val,x
         integer                 :: NN,i,l

          val=0.0d0

          do l=0,i/2
           val=val+(((-1)**((i/2)-l) )/(fact(2*l)*fact((i/2)-l)))*((2*x)**(2*l))
          enddo

!          do i=1,NN,2
!           do l=0,i/2
!            val=val+(((-1)**(((i-1)/2)-l) )/(fact(2*l+1)*fact(((i-1)/2)-l)))*((2*x)**(2*l+1))    
!           enddo
!          enddo

          val=val*fact(i)

         return
         end function hermite

         function fact(i) result(val)
         implicit none
         integer :: i,l
         double precision :: val

          val=1.0d0

          if( i.lt.0)then
           write(*,*) 'Could not take the factorial of a negative value'
           stop
          endif

          if (i.eq.0)then
           val=1.0d0
          else
           do l=1,i
            val=val*l
           enddo
          endif

         return
         end function fact


        end module units_parms
