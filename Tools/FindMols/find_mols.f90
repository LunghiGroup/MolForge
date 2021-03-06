        program get_mols
        use mpi
        use mpi_utils
        use atoms_class
        implicit none
        integer                        :: i,j,frames=1,ii,step=1,skip=0,centre_atom_id=1
        logical                        :: new_type,centre_atom=.false.,remap_mols=.false.,reorder_mols=.false.
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
          write(*,*) '-remap_mols   : Find Molecules and remap them around pbc'
          write(*,*) '-reorder_mols : Find Molecules and order list '
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

             case ('-remap_mols') 
                 remap_mols=.true.

             case ('-reorder_mols') 
                 reorder_mols=.true.

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
            call get_mass(sys%label(i),sys%mass(i))
           enddo

           if(centre_atom) call sys%wrap_geo(centre_atom_id)
           if(remap_mols .or. reorder_mols)  call sys%find_mols(reorder_mols,remap_mols)

           write(*,*) sys%nats
           write(*,*) sys%molid
           do i=1,sys%nats
             write(*,"(a2,2x,3(f10.6,2x))") sys%label(sys%kind(i)),sys%x(i,:)
           enddo

          endif
        
          call sys%delete()

         enddo ! frames

         close(10)
         close(11)

         call MPI_FINALIZE(err)

        return
        end program get_mols
