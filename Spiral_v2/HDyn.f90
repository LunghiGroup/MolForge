        program spiral2test
        use liuville_class
        use quantum_systems_class
        implicit none
        type(open_quantum_system)       :: sp
        double precision                :: step=10000.0d0,time=0.0d0,temp_rho=20.0d0,val(3),num
        integer                         :: start_step=0,dump_freq=1,nsteps=1000
        integer, allocatable            :: print_si(:)
        character(len=100)              :: type_rho0
        character(len=20)               :: rho_restart_file,lattice_restart
        double precision, allocatable   :: val_r(:),val_c(:)
        double complex, allocatable     :: Vx(:,:,:)
        double complex                  :: numc,coeff,trace

        character(len=20)                  :: option,input,output,arg
        integer                            :: i,j,l,key,mpi_nproc_spin,s,t,v,ncoeff
        integer                            :: dim_sizes(2),cart_world,coord(2),rank,aaa(2)
        integer, allocatable               :: map(:,:)
        logical                            :: reorder,wrap_around(2)
        type(dist_cmplx_mat)               :: AA,BB,CC
        integer                            :: N,M,NBl,MBl,numroc,info,ii,jj
        integer                            :: Nloc_row,Nloc_col,blacs_pnum

         call mpi_init(err)

         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)       
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         mpi_nproc_spin=mpi_nproc

         if(mpi_id.eq.0)then
          write(*,*) ''
          write(*,*) ''
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) '****                                                                        ****'
          write(*,*) '****                      SPIRAL@MolForge v1.0.0                            ****'
          write(*,*) '****                                                                        ****'
          write(*,*) '****               A first-principles spin dynamics software                ****'
          write(*,*) '****                                                                        ****'
          write(*,*) '****                      Author: Alessandro Lunghi                         ****'
          write(*,*) '****                        email: lunghia@tcd.ie                           ****'
          write(*,*) '****                                                                        ****'
          write(*,*) '****                                                                        ****'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) ''
          write(*,*) ''
          flush(6)
         endif

         if(mod(mpi_nproc,mpi_nproc_spin).ne.0)then
          write(*,*) 'Invalid CPUs decomposition'
          stop
         endif

         dim_sizes(2)=mpi_nproc_spin
         dim_sizes(1)=mpi_nproc/mpi_nproc_spin
         reorder=.false.
         wrap_around=.false.

         call mpi_cart_create(mpi_comm_world,2,dim_sizes,wrap_around,reorder,cart_world,err)
         call mpi_comm_rank(cart_world,rank,err)
         call mpi_cart_coords(cart_world,rank,2,coord,err)

         ! aggiungere genesi nuovo comm per blacs subset

         mpi_color=coord(1)
         key=mpi_id

         call mpi_comm_split(cart_world,mpi_color,key,mpi_blacs_world,err)

         call MPI_COMM_RANK(MPI_BLACS_WORLD,mpi_blacs_id,err)       
         call MPI_COMM_SIZE(MPI_BLACS_WORLD,mpi_blacs_nproc,err)

         mpi_color=coord(2)
         key=mpi_id

         call mpi_comm_split(cart_world,mpi_color,key,mpi_phonons_world,err)
         call MPI_COMM_RANK(MPI_PHONONS_WORLD,mpi_phonons_id,err)       
         call MPI_COMM_SIZE(MPI_PHONONS_WORLD,mpi_phonons_nproc,err)

!         call setup_blacs(mpi_nproc_spin)

         allocate(context(mpi_phonons_nproc))

         nprow=int(sqrt(dble(mpi_nproc_spin)))
         npcol=mpi_nproc_spin/nprow
         allocate(map(nprow,npcol))

         do i=1,mpi_phonons_nproc
          aaa(2)=0
          aaa(1)=i-1
          do jj=1,size(map,2)
           do ii=1,size(map,1)
            call mpi_cart_rank(cart_world,aaa,rank,err)
            map(ii,jj)=rank
            aaa(2)=aaa(2)+1
           enddo
          enddo
          call setup_multiblacs(mpi_nproc_spin,map,context(i))
         enddo
         call blacs_set_gridinfo()         

         if(myrow.ne.-1) mpi_color=1
         if(myrow.eq.-1) mpi_color=2

         key=mpi_blacs_id
         call mpi_comm_split(mpi_blacs_world,mpi_color,key,mpi_blacs_world,err)

         call MPI_COMM_RANK(MPI_BLACS_WORLD,mpi_blacs_id,err)
         call MPI_COMM_SIZE(MPI_BLACS_WORLD,mpi_blacs_nproc,err)

         key=mpi_phonons_id
         call mpi_comm_split(mpi_phonons_world,mpi_color,key,mpi_phonons_world,err)

         call MPI_COMM_RANK(MPI_PHONONS_WORLD,mpi_phonons_id,err)
         call MPI_COMM_SIZE(MPI_PHONONS_WORLD,mpi_phonons_nproc,err)

         if(myrow.eq.-1)goto 20

         if(mpi_id.eq.0)then

          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) 'Total Number of MPI processes: ',mpi_nproc
          write(*,*) 'Total Number of MPI processes dedicated to spin matrices: ',mpi_blacs_nproc
          write(*,*) 'Total Number of MPI processes dedicated to phonons matrices: ',mpi_phonons_nproc
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) ''
          flush(6)

         endif

         call getarg(1,arg)
         read(arg,*) sp%temp

! construct the hilbert space

         sp%Hdim=4
         call sp%H0%set(sp%Hdim,sp%Hdim,NB,MB)           
         sp%H0%mat=(0.0d0,0.0d0)

         open(12,file='eigen_Co.dat')
         do ii=1,sp%Hdim
          read(12,*) num
          numc=cmplx(num,0.0d0,8)
          call pzelset(sp%H0%mat,ii,ii,sp%H0%desc,numc)
         enddo
         close(12)

         call sp%diag_Hmat()

         allocate(val_r(sp%Hdim))
         allocate(val_c(sp%Hdim))

         open(12,file='J_SOC_R.dat')
         open(13,file='J_SOC_C.dat')

         allocate(sp%QMOP(1))
         call sp%QMOP(1)%set(sp%Hdim,sp%Hdim,NB,MB)

         do i=1,sp%Hdim
          read(12,*) val_r
          read(13,*) val_c
          sp%QMOP(1)%mat(i,:)=cmplx(val_r,val_c,8)
         enddo

         close(12)
         close(13)

         type_rho0='READ_FILE'
         rho_restart_file='rho_0.dat'
         call sp%make_rho(type_rho0,temp_rho,rho_restart_file)
         call sp%make_unitary_propagator(step)

! construct the liuville space

         call sp%make_basis_L()
         
         lattice_restart='FC2_Co'
         call sp%lattice%read_restart_file(lattice_restart)

         sp%phonons%nx=1
         sp%phonons%ny=1
         sp%phonons%nz=1
         sp%phonons%ntot=1
         call sp%phonons%generate_mesh(sp%phonons%nx,sp%phonons%ny,sp%phonons%nz)
         call sp%phonons%calc_bands(sp%lattice)

         sp%smear=10.0d0
         sp%norder=1
         sp%nderiv=213

         allocate(sp%map_s2a(sp%nderiv,sp%norder*2))
         allocate(sp%Vx(sp%nderiv))         

         allocate(Vx(sp%nderiv,sp%Hdim,sp%Hdim))

         open(12,file='grad_matrix_R_Co.dat')
         open(13,file='grad_matrix_C_Co.dat')
         read(12,*)  
         read(13,*) 

         sp%map_s2a(:,2)=1

         do i=1,sp%nderiv

          read(12,*) sp%map_s2a(i,1)
          read(13,*) sp%map_s2a(i,1)

          do j=1,sp%Hdim
           read(12,*) val_r
           read(13,*) val_c
           Vx(i,j,:)=cmplx(val_r,val_c,8)
          enddo
          do j=sp%Hdim+1,40
           read(12,*) val_r
           read(13,*) val_c
          enddo

!          trace=(0.0d0,0.0d0)

!          do ii=1,sp%Hdim
!           trace=trace+Vx(i,ii,ii)
!          enddo

!          do ii=1,sp%Hdim
!           Vx(i,ii,ii)=Vx(i,ii,ii)-trace/sp%Hdim
!          enddo

         enddo

         close(12)
         close(13)

         do s=1,sp%Hdim
          do t=1,sp%Hdim
           do l=0,2

            coeff=0.0d0
            ncoeff=0

            do i=1,sp%nderiv
             v=mod(sp%map_s2a(i,1)+2,3)
             if(v.ne.l) cycle
             coeff=coeff+Vx(i,s,t)
             ncoeff=ncoeff+1
            enddo
            do i=1,sp%nderiv
             v=mod(sp%map_s2a(i,1)+2,3)
             if(v.ne.l) cycle
             Vx(i,s,t)=Vx(i,s,t)-coeff/dble(ncoeff)
            enddo

           enddo
          enddo
         enddo
      
   ! Make Vx
        
         do i=1,sp%nderiv

          call sp%Vx(i)%set(sp%Hdim,sp%Hdim,NB,MB)           
          sp%Vx(i)%mat=(0.0d0,0.0d0)

          do ii=1,sp%Hdim
           do jj=1,sp%Hdim
            call pzelset(sp%Vx(i)%mat,ii,jj,sp%Vx(i)%desc,Vx(i,ii,jj))
           enddo
          enddo
  
          call sp%to_eigenbasis(sp%Vx(i))           
       
         enddo


!         call sp%make_R21_lindbladian(min_ener=5.0d0,max_ener=3000.0d0)
         call sp%make_R41_lindbladian(min_ener=5.0d0,max_ener=500.0d0)
         call sp%diag_limbladian(step)

! propagate

         call sp%propagate(start_step,time,nsteps,step,dump_freq)


20       continue
         call blacs_exit(-1)
         call mpi_finalize(err)
         stop

        return
        end program spiral2test

