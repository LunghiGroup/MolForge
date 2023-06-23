        program spiral2test
        use spin_phonon_class
        use liuville_class
        use quantum_systems_class
        implicit none
        type(spin_quantum_system)       :: sp
        double precision                :: step=10000.0d0,time=0.0d0,temp_rho=20.0d0,val(3),valr,euler(3)
        integer                         :: start_step=0,dump_freq=1,nsteps=1000
        integer, allocatable            :: print_si(:)
        character(len=100)              :: type_rho0
        character(len=20)               :: rho_restart_file,lattice_restart

        character(len=20)                  :: option,input,output,arg
        integer                            :: i,j,l,key,mpi_nproc_spin
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

         sp%spins%nspins=1
         sp%spins%nspins_pr=1
         sp%spins%nkinds=1
!        sp%spins%Bfield(3)=0.1d0
!         sp%spins%Bfield(1)=0.20860042062998443
!         sp%spins%Bfield(2)=2.8072557715195887e-09
!         sp%spins%Bfield(3)=-0.9780009532270372
!         sp%spins%Bfield=sp%spins%Bfield*0.01d0
         sp%spins%Bfield=0.0d0

!         euler(1)=3.1415810050546917       
!         euler(2)=0.29850593110274537       
!         euler(3)=-1.5708040439131328

         allocate(sp%spins%kind(1))
         allocate(sp%spins%spin(1))
         allocate(sp%spins%bohr_mag(1))

         sp%spins%kind(1)=1
         sp%spins%spin(1)=1.5d0
         sp%spins%bohr_mag(1)=-0.466867723d0
        
         sp%spins%SH%nG=1
         allocate(sp%spins%SH%G(1))
         sp%spins%SH%G(1)%kind=1

         sp%spins%SH%G(1)%G=0.0d0
         sp%spins%SH%G(1)%G(1,1)=2.0d0
         sp%spins%SH%G(1)%G(2,2)=2.0d0
         sp%spins%SH%G(1)%G(3,3)=2.0d0

!         sp%spins%SH%G(1)%G(1,1)=-0.791021    
!         sp%spins%SH%G(1)%G(1,2)=-0.189379    
!         sp%spins%SH%G(1)%G(1,3)=-1.907417
!         sp%spins%SH%G(1)%G(2,1)=0.054071   
!         sp%spins%SH%G(1)%G(2,2)=-0.225838   
!         sp%spins%SH%G(1)%G(2,3)=0.000025
!         sp%spins%SH%G(1)%G(3,1)=1.829858  
!         sp%spins%SH%G(1)%G(3,2)=0.438101    
!         sp%spins%SH%G(1)%G(3,3)=6.437077

         sp%spins%SH%nDSI=1
         allocate(sp%spins%SH%DSI(1))
         sp%spins%SH%DSI(1)%kind=1

!         sp%spins%SH%DSI(1)%D(1,1)=7.077304    
!         sp%spins%SH%DSI(1)%D(1,2)=-0.000004     
!         sp%spins%SH%DSI(1)%D(1,3)=3.769654
!         sp%spins%SH%DSI(1)%D(2,1)=-0.000004     
!         sp%spins%SH%DSI(1)%D(2,2)=7.316755    
!         sp%spins%SH%DSI(1)%D(2,3)=-0.000046
!         sp%spins%SH%DSI(1)%D(3,1)=3.769654    
!         sp%spins%SH%DSI(1)%D(3,2)=-0.000046    
!         sp%spins%SH%DSI(1)%D(3,3)=-4.013844 

         sp%spins%SH%nO=1
         allocate(sp%spins%SH%O(1))
         sp%spins%SH%O(1)%kind=1
         sp%spins%SH%O(1)%k=2
         
         allocate(sp%spins%SH%O(1)%B(5))
         allocate(sp%spins%SH%O(1)%q(5))
         sp%spins%SH%O(1)%q(1)=-2
         sp%spins%SH%O(1)%q(2)=-1
         sp%spins%SH%O(1)%q(3)=0
         sp%spins%SH%O(1)%q(4)=1
         sp%spins%SH%O(1)%q(5)=2
         sp%spins%SH%O(1)%B(1)=0.000000026083121591
         sp%spins%SH%O(1)%B(2)=0.000000041829917027
         sp%spins%SH%O(1)%B(3)=-6.6130296343394495
         sp%spins%SH%O(1)%B(4)=2.4106049909311422
         sp%spins%SH%O(1)%B(5)=-0.66640176704401499

!         call sp%spins%SH%rot(euler)

! construct the hilbert space

         call sp%make_spin_basis()                   
         call sp%make_SH_rep(sp%spins%SH,-1,-1)
         call sp%make_Hmat_nodes()
         call sp%make_Hmat()
         call sp%diag_Hmat()

         allocate(print_si(1))
         print_si(1)=1
         call sp%make_Smat(print_si)
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

         sp%temp=10.0d0
         sp%smear=10.0d0
         sp%norder=1
         sp%nderiv=213

!         sp%spins%SPH%nDSI=1 
!         allocate(sp%spins%SPH%DSI(1))
!         allocate(sp%spins%SPH%DSI_t(1))
!         sp%spins%SPH%DSI(1)%kind=1
!         sp%spins%SPH%DSI_t(1)%kind=1
!         allocate(sp%spins%SPH%DSI_t(1)%Dcart(sp%nderiv))
!         allocate(sp%map_s2a(sp%nderiv,sp%norder*2))

!         open(12,file='V_Cmat.dat')
!         read(12,*) 
!         do i=1,sp%nderiv
!          read(12,*) sp%map_s2a(i,:)
!          read(12,*) val
!          sp%spins%SPH%DSI_t(1)%Dcart(i)%D(1,:)=cmplx(val,0.0d0,8)
!          read(12,*) val
!          sp%spins%SPH%DSI_t(1)%Dcart(i)%D(2,:)=cmplx(val,0.0d0,8)
!          read(12,*) val
!          sp%spins%SPH%DSI_t(1)%Dcart(i)%D(3,:)=cmplx(val,0.0d0,8)
!         enddo
!         allocate(sp%spins%SPH%DSI_t(1)%map_s2a(sp%nderiv,sp%norder*2))
!         sp%spins%SPH%DSI_t(1)%map_s2a=sp%map_s2a
!         sp%spins%SPH%DSI_t(1)%norder=1
!         sp%spins%SPH%DSI_t(1)%nderiv=213
!         call sp%spins%SPH%DSI_t(1)%do_tinv()
      
         sp%spins%SPH%nO=1
         allocate(sp%spins%SPH%O(1))
         allocate(sp%spins%SPH%O_t(1))

         sp%spins%SPH%O(1)%kind=1
         sp%spins%SPH%O(1)%k=2
         sp%spins%SPH%O_t(1)%kind=1
         sp%spins%SPH%O_t(1)%k=2
         allocate(sp%spins%SPH%O_t(1)%Ocart(sp%nderiv))

         allocate(sp%map_s2a(sp%nderiv,sp%norder*2))

         allocate(sp%spins%SPH%O(1)%B(5))
         allocate(sp%spins%SPH%O(1)%q(5))

         open(12,file='B_values.dat')
         read(12,*) 
         do i=1,sp%nderiv
          allocate(sp%spins%SPH%O_t(1)%Ocart(i)%B(5))
          allocate(sp%spins%SPH%O_t(1)%Ocart(i)%q(5))
          sp%spins%SPH%O_t(1)%Ocart(i)%kind=1
          sp%spins%SPH%O_t(1)%Ocart(i)%k=2
          read(12,*) sp%map_s2a(i,:)
          do j=1,5
           read(12,*) sp%spins%SPH%O_t(1)%Ocart(i)%q(j),valr
           sp%spins%SPH%O_t(1)%Ocart(i)%B(j)=valr
          enddo
         enddo
         allocate(sp%spins%SPH%O_t(1)%map_s2a(sp%nderiv,sp%norder*2))
         sp%spins%SPH%O_t(1)%map_s2a=sp%map_s2a
         sp%spins%SPH%O_t(1)%norder=1
         sp%spins%SPH%O_t(1)%nderiv=sp%nderiv
         call sp%spins%SPH%O_t(1)%do_tinv()

   ! Make Vx

         allocate(sp%Vx(sp%nderiv))         

         do i=1,sp%nderiv

          call sp%Vx(i)%set(sp%Hdim,sp%Hdim,NB,MB)           
          sp%Vx(i)%mat=(0.0d0,0.0d0)

          sp%spins%SPH%O(1)%B=sp%spins%SPH%O_t(1)%Ocart(i)%B
          sp%spins%SPH%O(1)%q=sp%spins%SPH%O_t(1)%Ocart(i)%q
!          sp%spins%SPH%DSI(1)%D=sp%spins%SPH%DSI_t(1)%Dcart(i)%D

!          call sp%spins%SPH%rot(euler)

          call sp%make_SH_rep(sp%spins%SPH,1,1) 
 
          do l=1,size(sp%Hnodes,1) 
 
           ii=sp%Hnodes(l,1) 
           jj=sp%Hnodes(l,2) 
             
           sp%Vx(i)%mat(ii,jj)=sp%get_Hij(sp%Hnodes(l,1:4))

          enddo ! l            
 
          call sp%to_eigenbasis(sp%Vx(i)) 

         enddo

!         call sp%make_R21_limbladian(min_ener=5.0d0,max_ener=3000.0d0)
         call sp%make_R41_limbladian(min_ener=5.0d0,max_ener=500.0d0)
         call sp%diag_limbladian(step)

! propagate

         call sp%propagate(start_step,time,nsteps,step,dump_freq)


20       continue
         call blacs_exit(-1)
         call mpi_finalize(err)
         stop

        return
        end program spiral2test

