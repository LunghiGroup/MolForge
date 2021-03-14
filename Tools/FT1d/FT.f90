        program fourier
        implicit none
        integer                       :: N,M,omega
        integer                       :: i,l
        double precision              :: kconst,step
        double precision, allocatable :: corr(:)
        double complex, allocatable   :: ft(:)
         
         N=2000
         M=3000

         allocate(corr(N))
         allocate(ft(M))

         open(11,file='correlation.dat')         
         do i=1,N
          read(11,*) l,corr(i)
         enddo
         close(11)

         ft=(0.0d0,0.0d0)
         do omega=1,M
          kconst=2.0d0*acos(-1.0d0)*5.0d0/1000.0d0/33.357177561
          kconst=kconst*((omega-1)*1.0d0+4.0d0)
          do i=1,N
           if(i.eq.1 .or. i.eq.N)then
            ft(omega)=ft(omega)+0.5d0*corr(i)*exp(cmplx(0.0d0,kconst*(i-1),8))
           else
            ft(omega)=ft(omega)+corr(i)*exp(cmplx(0.0d0,kconst*(i-1),8))
           endif
          enddo
         enddo

         ft=ft*5.0d-3

         open(11,file='FT.dat')
         do i=1,M
          write(11,*) (i-1)*1d0+4.0d0,dble(ft(i)),aimag(ft(i))
         enddo
         close(11)         

        return
        end program fourier
