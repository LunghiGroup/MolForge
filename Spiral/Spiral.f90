        program Spiral_Main
        use mpi
        use mpi_utils
        use blacs_utils
        use control_variables
        use spins_dist_rs_class
        use hilbert_dist_class
        use atoms_class
        use phonons_class
        use spindy_parser
        implicit none
        type(spins_hilbert)                :: spindy,gsh
        type(atoms_group)                  :: lattice
        type(brillouin)                    :: phondy
        integer                            :: t1,t2,rate
        character(len=20)                  :: option,input,output,arg
        integer                            :: i,j,l,key,mpi_nproc_spin

        integer :: dim_sizes(2),cart_world,coord(2),rank,aaa(2)
        integer, allocatable :: map(:,:)
        logical :: reorder,wrap_around(2)

        type(dist_cmplx_mat)          :: AA,BB,CC
        integer                       :: N,M,NBl,MBl,numroc,info,ii,jj
        integer                       :: Nloc_row,Nloc_col,blacs_pnum

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
          write(*,*) '****                      SPIRAL@MolForge v1.0_beta.0                       ****'
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

         if ( iargc().eq.0 )then
          if(mpi_id.eq.0)then
           write(*,*) '-i                  : input file'
           write(*,*) '-o                  : output file'
           write(*,*) '-spin_cpu           : number of MPI processes allocated to Spin matrices'
          endif
          call blacs_exit(-1)
          call mpi_finalize(err)
          stop
         endif

         do i=1,iargc()

          call getarg(i,arg)

          if (arg.eq.'-spin_cpu')then
           call getarg(i+1,arg)
           read(arg,*) mpi_nproc_spin
          endif

         enddo

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

         if(mpi_id.eq.0)then

          call system_clock(t1,rate)         

          do l=1,iargc()

           call getarg(l,option)

           select case (trim(option))

            case('-i')
             call getarg(l+1,input)

            case('-o')

            call getarg(l+1,output)
            open(6,file=output,recl=1056)

           end select

          enddo


          open(10,file=input)
          rewind(10)
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) 'Parsing input file: ',trim(input)
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) ''
          flush(6)

         endif

         do 

          if(mpi_id.eq.0) call parse_input(input,spindy,phondy,lattice,gsh)
         
          call mpi_bcast(operation,100,mpi_character,0,mpi_comm_world,err)

          if(operation.eq.'END') exit

          call do_action(spindy,phondy,lattice,gsh)

         enddo

         if(mpi_id.eq.0)then

          call system_clock(t2)
          write(*,*) ''
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'
          write(*,*) 'SPIRAL correctly ended '
          write(*,*) 'Total running time: ',real(t2-t1)/real(rate),'s'
          write(*,*) '********************************************************************************'
          write(*,*) '********************************************************************************'

          close(10)        
          close(6)

         endif

20       continue
         call blacs_exit(-1)
         call mpi_finalize(err)
         stop

        return
        end program Spiral_Main
     

