        program distribution
        use dist_class
        implicit none
        integer                       :: i,nsteps
        character(len=200)            :: option,filein,fileout
        type(dist1d)                  :: distfunc
        double precision              :: min_val,max_val,val,coeff
        logical                       :: donorm=.false.

        if(iargc().eq.0)then
         write(*,*) '-filein         : '
         write(*,*) '-fileout        : '
         write(*,*) '-sigma          : '
         write(*,*) '-step           : '
         write(*,*) '-minmax         : '
         write(*,*) '-norm           : '
         write(*,*) '-type_smear     : Lorentzian Smear = 0, Gaussian Smear = 1'
         stop
        endif

        do i=1,iargc()

         call  getarg(i,option)

         select case (option)

         case('-filein')
          call getarg(i+1,option)
          read(option,*) filein
          open(unit=11,file=filein)

         case('-fileout')
          call getarg(i+1,option)
          read(option,*) fileout
          open(unit=12,file=fileout)

         case('-sigma')
          call getarg(i+1,option)
          read(option,*) distfunc%sigma

         case('-step')
          call getarg(i+1,option)
          read(option,*) distfunc%step

         case('-minmax')
          call getarg(i+1,option)
          read(option,*) min_val
          call getarg(i+2,option)
          read(option,*) max_val

         case('-type_smear')
          call getarg(i+1,option)
          read(option,*) distfunc%type_smear

         case('-norm')
          donorm=.true.

         end select
        enddo

        distfunc%nsteps=nint((max_val-min_val)/distfunc%step)+1
        distfunc%shift=min_val
        call distfunc%alloc_dist()

        do
         read(11,*,end=10)
         nsteps = nsteps + 1
        enddo
        10 continue
        rewind(11)

        do i=1,nsteps
         read(11,*) val,coeff
         call distfunc%update_dist(val,coeff)
        enddo
        close(11)

        if(donorm) call distfunc%norm_dist()

        do i=1,distfunc%nsteps
         write(12,*) (i-1)*distfunc%step+distfunc%shift,distfunc%dist(i)
        enddo
        close(12)

        return
        end program distribution
