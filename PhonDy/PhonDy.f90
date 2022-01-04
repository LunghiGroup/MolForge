        program Phonon_Dynamics
        use mpi
        use atoms_class
        use phonons_class
        use random_numbers_class
        implicit none

        type(atoms_group)              :: lattice
        type(brillouin)                :: kspace

        logical                        :: read_dipole=.false.
        character(len=50)              :: dipole_file,restart_file
        integer, allocatable           :: npts(:)
        double precision, allocatable  :: ki(:,:),kf(:,:)

        integer            :: i,t1,t2,err,mpi_nproc,mpi_id,npaths
        integer            :: nx,ny,nz,s
        integer            :: print_cart_kp,print_cart_bn,print_cart_nsteps=10
        integer            :: cart_dist_1,cart_dist_2
        double precision   :: print_cart_step=0.1
        logical            :: print_cart=.false.,print_amp=.false.,print_dist=.false.
        character(len=50)  :: arg
        double precision   :: rate,distmin=0.0d0,distmax=3500.0d0,step=1.0d0
        double precision   :: distmin2=0.0d0,distmax2=3500.0d0,step2=1.0d0
        logical            :: dodos2p=.false.,dodos1p=.false.

         call MPI_INIT(err)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)
         call system_clock(t1,rate)

         call init_random_seed()
 
         if(mpi_id.eq.0)then

          if(iargc().eq.0)then
           write(*,*) 'Usage:'
           write(*,*) '-restart         = File containing the force constants'
           write(*,*) '-smear           = Gaussian smear for dos in cm-1'
           write(*,*) '-nk              = <nk1> <nk2> <nk3>'
           write(*,*) '-cart_disp       = kp bn nsteps step'
           write(*,*) '-dos1p           = min max step'
           write(*,*) '-dos2p           = min max step'
           write(*,*) '-path            = file describing the BZ path'
           write(*,*) '-dipole          = file containing the dip. mom. derivatives'
           stop
          endif

          do i=1,iargc()

           call getarg(i,arg)

           if (arg.eq.'-path')then
            do_brillouin_path=.true.
            call getarg(i+1,arg)
            open(12,file=trim(arg))
            read(12,*) npaths
            allocate(npts(npaths))
            allocate(ki(3,npaths))
            allocate(kf(3,npaths))
            write(*,*) npaths
            do s=1,npaths
             read(12,*) npts(s),ki(1,s),ki(2,s),ki(3,s),kf(1,s),kf(2,s),kf(3,s)
            enddo
            close(12)
           endif

           if (arg.eq.'-nk')then
            do_brillouin_mesh=.true.
            call getarg(i+1,arg)
            read(arg,*) nx
            call getarg(i+2,arg)
            read(arg,*) ny
            call getarg(i+3,arg)
            read(arg,*) nz
           endif

           if (arg.eq.'-smear')then
            call getarg(i+1,arg)
            read(arg,*) smear
           endif

           if (arg.eq.'-cart_disp')then
            call getarg(i+1,arg)
            read(arg,*) print_cart_kp
            call getarg(i+2,arg)
            read(arg,*) print_cart_bn
            call getarg(i+3,arg)
            read(arg,*) print_cart_nsteps
            call getarg(i+4,arg)
            read(arg,*) print_cart_step
            print_cart=.true.
           endif

           if (arg.eq.'-dos1p')then
            call getarg(i+1,arg)
            read(arg,*) distmin
            call getarg(i+2,arg)
            read(arg,*) distmax
            call getarg(i+3,arg)
            read(arg,*) step
            dodos1p=.true.
           endif

           if (arg.eq.'-dos2p')then
            call getarg(i+1,arg)
            read(arg,*) distmin2
            call getarg(i+2,arg)
            read(arg,*) distmax2
            call getarg(i+3,arg)
            read(arg,*) step2
            dodos2p=.true.
           endif

           if (arg.eq.'-restart')then
            call getarg(i+1,arg)
            restart_file=arg
           endif

           if (arg.eq.'-dipole')then
            call getarg(i+1,arg)
            dipole_file=arg
            read_dipole=.true.
           endif

          enddo

         endif
         
          call mpi_bcast(do_brillouin_mesh,1,mpi_logical,0,MPI_COMM_WORLD,err)
          call mpi_bcast(do_brillouin_path,1,mpi_logical,0,MPI_COMM_WORLD,err)
          call mpi_bcast(read_dipole,1,mpi_logical,0,MPI_COMM_WORLD,err)

         if(do_brillouin_mesh)then
          call mpi_bcast(nx,1,mpi_integer,0,MPI_COMM_WORLD,err)
          call mpi_bcast(ny,1,mpi_integer,0,MPI_COMM_WORLD,err)
          call mpi_bcast(nz,1,mpi_integer,0,MPI_COMM_WORLD,err)
          call mpi_bcast(step,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          call mpi_bcast(distmin,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          call mpi_bcast(distmax,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          call mpi_bcast(dodos1p,1,mpi_logical,0,MPI_COMM_WORLD,err)
          call mpi_bcast(dodos2p,1,mpi_logical,0,MPI_COMM_WORLD,err)
          call mpi_bcast(step2,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          call mpi_bcast(distmin2,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          call mpi_bcast(distmax2,1,mpi_double_precision,0,MPI_COMM_WORLD,err)
          kspace%dos1p%step=step
          kspace%dos1p%sigma=smear
          kspace%dos1p%nsteps=1+INT((distmax-distmin)/step)
          kspace%dos2p%step=step2
          kspace%dos2p%sigma=smear
          kspace%dos2p%nsteps=1+INT((distmax2-distmin2)/step2)
         endif

         if(do_brillouin_path)then
          do i=1,size(npts)
           call mpi_bcast(npts(i),1,mpi_integer,0,MPI_COMM_WORLD,err)
           call mpi_bcast(ki(:,i),3,mpi_double_precision,0,MPI_COMM_WORLD,err)
           call mpi_bcast(kf(:,i),3,mpi_double_precision,0,MPI_COMM_WORLD,err)
          enddo
         endif

         call lattice%read_restart_file(restart_file)

         if(read_dipole)then
          call lattice%read_dipole_file(dipole_file)
         endif

         call lattice%atoms_bcast()
                
         if(do_brillouin_mesh)then
          call build_kspace(kspace,nx,ny,nz)
          call kspace%calc_bands(lattice)
          if(mpi_id.eq.0 .and. dodos1p) call kspace%calc_dos1p()      
          if(dodos2p) call kspace%calc_dos2p()
          if(read_dipole)then
           call kspace%calc_IR(lattice)
          endif
          if(print_cart)then
           call kspace%print_cart_disp(lattice,print_cart_kp,print_cart_bn,&
                        print_cart_step,print_cart_nsteps)
          endif
         endif

         if(do_brillouin_path)then
          call kspace%generate_path(npts,kf,ki)
          call kspace%calc_disps(lattice)
         endif

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) 'PhonDy total running time: ',real(t2-t1)/real(rate),'s'
          write(*,*) 'PhonDy finished correctly'
         endif

         call MPI_FINALIZE(err)

        return
        end program Phonon_Dynamics

        subroutine build_kspace(kspace,nx,ny,nz)
        use mpi
        use phonons_class
        use units_parms
        implicit none
        type(brillouin)    :: kspace
        integer            :: i,nx,ny,nz
        integer            :: err,mpi_nproc,mpi_id

         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         call kspace%generate_mesh(nx,ny,nz)

          boundary_scattering=.false.
          Length=1.0D5

          ntemps=25
          allocate(temp(ntemps))
          temp(1)=1.0d0
          temp(2)=2.0d0
          temp(3)=3.0d0
          temp(4)=4.0d0
          temp(5)=5.0d0
          temp(6)=6.0d0
          temp(7)=7.0d0
          temp(8)=8.0d0
          temp(9)=9.0d0
          temp(10)=10.0d0
          temp(11)=12.0d0
          temp(12)=14.0d0
          temp(13)=16.0d0
          temp(14)=18.0d0
          temp(15)=20.0d0
          temp(16)=25.0d0
          temp(17)=30.0d0
          temp(18)=35.0d0
          temp(19)=40.0d0
          temp(20)=45.0d0
          temp(21)=50.0d0
          temp(22)=55.0d0
          temp(23)=60.0d0
          temp(24)=65.0d0
          temp(25)=70.0d0
         

         if(mpi_id.eq.0)then

          write(*,*) '#################################'
          write(*,*) 'k-space mesh has been generate:',kspace%nx,kspace%ny,kspace%nz
          write(*,*) '#################################'
          write(*,*)
          write(*,*)
        
          write(*,*) '   Total number of k-space points:',kspace%ntot
          if (kspace%ntot.le.512)then
           do i=1,kspace%ntot 
            write(*,*) '      ',i,kspace%list(i)%k(1:3),kspace%list(i)%weight
           enddo
          endif
          write(*,*)
          write(*,*)

         endif

        return
        end subroutine build_kspace


