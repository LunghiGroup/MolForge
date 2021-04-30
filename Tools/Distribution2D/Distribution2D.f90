        program distribution2D
        use dist_class
        implicit none
        integer                       :: i,nsteps,j
        character(len=200)            :: option,filein,fileout
        type(dist2d)                  :: distfunc
        double precision              :: min_val1,max_val1,val1,coeff
        double precision              :: min_val2,max_val2,val2
        logical                       :: donorm=.false.

        if(iargc().eq.0)then
         write(*,*) '-filein         : '
         write(*,*) '-fileout        : '
         write(*,*) '-sigma1         : '
         write(*,*) '-step1          : '
         write(*,*) '-minmax1        : '
         write(*,*) '-sigma2         : '
         write(*,*) '-step2          : '
         write(*,*) '-minmax2        : '
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

         case('-sigma1')
          call getarg(i+1,option)
          read(option,*) distfunc%sigma1

         case('-step1')
          call getarg(i+1,option)
          read(option,*) distfunc%step1

         case('-minmax1')
          call getarg(i+1,option)
          read(option,*) min_val1
          call getarg(i+2,option)
          read(option,*) max_val1

         case('-sigma2')
          call getarg(i+1,option)
          read(option,*) distfunc%sigma2

         case('-step2')
          call getarg(i+1,option)
          read(option,*) distfunc%step2

         case('-minmax2')
          call getarg(i+1,option)
          read(option,*) min_val2
          call getarg(i+2,option)
          read(option,*) max_val2

         case('-type_smear')
          call getarg(i+1,option)
          read(option,*) distfunc%type_smear

         case('-norm')
          donorm=.true.

         end select
        enddo

        distfunc%nsteps1=nint((max_val1-min_val1)/distfunc%step1)+1
        distfunc%shift1=min_val1
        distfunc%nsteps2=nint((max_val2-min_val2)/distfunc%step2)+1
        distfunc%shift2=min_val2

        call distfunc%alloc_dist()

        do
         read(11,*,end=10)
         nsteps = nsteps + 1
        enddo
        10 continue
        rewind(11)

        do i=1,nsteps
         read(11,*) val1,val2,coeff
         call distfunc%update_dist(val1,val2,coeff)
        enddo
        close(11)

        if(donorm) call distfunc%norm_dist()

        do i=1,distfunc%nsteps1
         do j=1,distfunc%nsteps2
          write(12,*) (i-1)*distfunc%step1+distfunc%shift1,(j-1)*distfunc%step2+distfunc%shift2,distfunc%dist(i,j)
         enddo
         write(12,*)
        enddo
        close(12)

        return
        end program distribution2D
