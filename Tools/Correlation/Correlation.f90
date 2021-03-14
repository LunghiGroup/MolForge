        program Correlation_Function
        implicit none
        integer                       :: nstep,bnstep,step,D,TMAX,MCOR,ialloc,death
        integer                       :: i,j,N,l
        double precision              :: sumA,sumB,norm,sigmaA,sigmaB
        double precision, allocatable :: Mid_Corr(:),Corr(:),A(:),B(:)
        character(len=20)             :: option,A_file_inp,B_file_inp
        logical                       :: do_norm_sigma=.false.,do_norm_unit=.false.,crosscorr=.false.

        if(iargc().eq.0)then
         write(*,*) '-A          : File containing property A'
         write(*,*) '-B          : File containing property B (optional for cross correlation'
         write(*,*) '-D          : Step of correlation function binning '
         write(*,*) '-T          : Length of the correlation function'
         write(*,*) '-step       : Number of steps between different origins for calculating average'
         write(*,*) '-norm_sigma : Normalize the function to <A(0)B(0)>=sigmaA*sigmaB '
         write(*,*) '-norm_unit  : Normalize the function to <A(0)B(0)>=1 '
         stop
        endif

        open(unit=12,file='Correlation.dat')

        do i=1,iargc()

         call  getarg(i,option)

         select case (option)

         case('-A')
          call getarg(i+1,option)
          read(option,*) A_file_inp
          open(unit=11,file=A_file_inp)

         case('-B')
          call getarg(i+1,option)
          read(option,*) B_file_inp
          open(unit=13,file=B_file_inp)
          crosscorr=.true.

         case('-D')
          call getarg(i+1,option)
          read(option,*) D

         case('-T')
          call getarg(i+1,option)
          read(option,*) TMAX

         case('-step')
          call getarg(i+1,option)
          read(option,*) step

         case('-norm_sigma')
          do_norm_sigma=.true.

         case('-norm_unit')
          do_norm_unit=.true.
   
         end select

        enddo

        NSTEP = 0
        BNSTEP=0

        do
         read(11,*,end=10)
         nstep = nstep + 1
        enddo
        10 continue
        rewind(11)

        if(crosscorr)then
         do
         read(13,*,end=9)
          bnstep = bnstep + 1
         enddo
         9 continue
         rewind(13)

         if(nstep.ne.bnstep)then
          write(*,*) 'Il numero di step di A e B è diverso!!!'
          stop
         endif
        endif

        mcor = tmax/D

        allocate(CORR(MCOR),stat=ialloc)
        allocate(MID_CORR(MCOR),stat=ialloc)
        allocate(A(NSTEP),stat=ialloc)

        sumA=0.0
        sigmaA=0.0
        do j = 1, NSTEP  
         read(11,*) A(j)
         sumA=sumA+A(j)
        enddo
        close(11)
        sumA=sumA/NSTEP

        do j=1,NSTEP
         A(j)=A(j)-sumA
        enddo

        do i=1,NSTEP
         sigmaA=sigmaA+A(i)**2
        enddo

        sigmaA=dsqrt(sigmaA/(NSTEP))

        write(*,*) 'Media Proprietà A          :',sumA
        write(*,*) 'Dev. Std.  Proprietà A     :',sigmaA

        if(crosscorr)then

         sumB=0.0
         sigmaB=0.0
         ALLOCATE (B(nstep),stat=ialloc)
         do j = 1,nstep
          read(13,*) B(j)
          sumB=sumB+B(j)
         enddo
         close(13)
         sumB=sumB/NSTEP
         do i=1,NSTEP
          B(i)=B(i)-sumB
         enddo
         do i=1,NSTEP
          sigmaB=sigmaB+B(i)**2
         enddo
         sigmaB=dsqrt(sigmaB/(NSTEP))
         write(*,*) 'Media Proprietà B          :',sumB
         write(*,*) 'Dev. Std.  Proprietà B     :',sigmaB

        endif

        N=0
        death=NSTEP-TMAX

        corr=0.0d0
        mid_corr=0.0d0

        if(crosscorr)then

         do j=1,death,step
          N=N+1 
          do i=1,TMAX,D
           MID_CORR(i)=A(j)*B(j+i-1)
           CORR(i)=(CORR(i)*(N-1)+MID_CORR(i))/N
          enddo
         enddo

         if(do_norm_sigma)then 
          norm=sigmaA*sigmaB
          corr=corr/norm
         endif

         if(do_norm_unit)then
          norm=corr(1)
          corr=corr/norm
         endif

        else

         do j=1,death,step
          N=N+1
          l=0
          do i=1,MCOR
           MID_CORR(i)=A(j)*A(j+l)
           CORR(i)=(CORR(i)*(N-1)+MID_CORR(i))/N
           l=l+D
          enddo
         enddo

         if(do_norm_sigma)then
          norm=sigmaA*sigmaA
          corr=corr/norm
         endif

         if(do_norm_unit)then
          norm=corr(1)
          corr=corr/norm
         endif

        endif

        l=0
        do i = 1, MCOR
         write(12,*) l,CORR(i)
         l=l+D
        enddo
        close(12)

        return
        end program Correlation_Function
