        program get_mols
        use mpi
        use mpi_utils
        use atoms_class
        implicit none
        integer                        :: i,j,frames=1,ii,step=1,skip=0,centre_atom_id=1
        logical                        :: new_type,centre_atom=.false.
        double precision               :: cell(3,3)
        character(len=100)             :: word,inp_file,cell_file
        character(len=10), allocatable :: label(:)
        type(atoms_group)              :: sys
        
         call MPI_INIT(err)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         if( iargc().eq.0)then
          write(*,*) 'Findmols Usage:'
          write(*,*) '-xyz          : Name of the xyz file to read from'
          write(*,*) '-cell         : 3x3 Matrix of the Cell'
          write(*,*) '-frames       : Number of frames in the xyz file'
          write(*,*) '-step         : Sampling Frequency from xyz file'
          write(*,*) '-skip_steps   : Skip N steps before Sampling the xyz file'
          write(*,*) '-centre       : Set the centre of the cell around one atom'
          call MPI_FINALIZE(err)
          stop
         endif

         do i=1,iargc()
          call getarg(i,word)
          
          select case (trim(word))

             case ('-frames')
                 call getarg(i+1,word)
                 read(word,*) frames

             case ('-step')
                 call getarg(i+1,word)
                 read(word,*) step

             case ('-xyz')
                 call getarg(i+1,word)
                 read(word,*) inp_file

             case ('-cell')
                 call getarg(i+1,word)
                 read(word,*) cell_file

             case ('-centre')
                 call getarg(i+1,word)
                 read(word,*) centre_atom_id
                 centre_atom=.true.

          end select

         enddo

         open(11,file=trim(cell_file))
         open(10,file=inp_file)

         do ii=1,frames
                    
          cell=0.0d0
          read(11,*) cell(1,1),cell(1,2),cell(1,3)
          read(11,*) cell(2,1),cell(2,2),cell(2,3)
          read(11,*) cell(3,1),cell(3,2),cell(3,3)

          call sys%read_cell(cell)
          sys%nx=1
          sys%ny=1
          sys%nz=1
          sys%ntot=1
          call sys%do_supercell()

          read(10,*) sys%nats
          read(10,*)

          allocate(sys%x(sys%nats,3))
          if(allocated(label)) deallocate(label)
          allocate(label(sys%nats))

          do i=1,sys%nats
           read(10,*) label(i),sys%x(i,:)
          enddo

          if (ii.gt.skip .and. (ii.eq.(skip+1) .or. mod(ii+skip,step).eq.0))then      

           allocate(sys%kind(sys%nats))
           sys%kind(1)=1
           sys%nkinds=1

           do i=2,sys%nats
            new_type=.true.
            do j=1,i-1

             if(label(i).eq.label(j))then
              sys%kind(i)=sys%kind(j)
              new_type=.false.
              exit
             endif

            enddo

            if(new_type)then
             sys%nkinds=sys%nkinds+1
             sys%kind(i)=sys%nkinds
            endif
           enddo

           allocate(sys%label(sys%nkinds))

           do i=1,sys%nkinds
            do j=1,sys%nats
             if(i.eq.sys%kind(j))then
              sys%label(i)=label(j)
              exit
             endif
            enddo
           enddo

           allocate(sys%mass(sys%nkinds))
           sys%mass=0.0d0
           do i=1,sys%nkinds
            if(trim(sys%label(i)).eq.'C') sys%mass(i)=12.010700225830078
            if(trim(sys%label(i)).eq.'Co') sys%mass(i)=58.9331950000000
            if(trim(sys%label(i)).eq.'P') sys%mass(i)=30.97376200000000
            if(trim(sys%label(i)).eq.'F') sys%mass(i)=18.98840300000000
            if(trim(sys%label(i)).eq.'Dy') sys%mass(i)=162.5000000000000
            if(trim(sys%label(i)).eq.'H') sys%mass(i)=1.0078999996185303
            if(trim(sys%label(i)).eq.'V') sys%mass(i)=50.941501617431641
            if(trim(sys%label(i)).eq.'S') sys%mass(i)=32.064998626708984
            if(trim(sys%label(i)).eq.'P') sys%mass(i)=30.973760000000000
            if(trim(sys%label(i)).eq.'Se') sys%mass(i)=78.96000000000000
            if(trim(sys%label(i)).eq.'N') sys%mass(i)=14.006699562072754
            if(trim(sys%label(i)).eq.'O') sys%mass(i)=15.999899864196777
           enddo

           if(centre_atom) call sys%wrap_geo(centre_atom_id)
           call sys%find_mols()

          endif
        
          call sys%delete()

         enddo ! frames

         close(10)
         close(11)

         call MPI_FINALIZE(err)

        return
        end program get_mols
