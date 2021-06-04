        program get_bis
        use mpi
        use mpi_utils
        use atoms_class
        implicit none
        integer                        :: atom,Jmax,i
        double precision               :: rcut
        double precision               :: cell(3,3)
        character(len=100)             :: inp_file,cell_file,word
        type(atoms_group)              :: sys
        
         call MPI_INIT(err)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         if( iargc().eq.0)then
          write(*,*) 'Get_Bispectrum Usage:'
          write(*,*) '-xyz          : Name of the xyz file to read from'
          write(*,*) '-cell         : 3x3 Matrix of the Cell'
          write(*,*) '-atom         : Number of frames in the xyz file'
          write(*,*) '-Jmax         : Sampling Frequency from xyz file'
          write(*,*) '-rcut         : Skip N steps before Sampling the xyz file'
          call MPI_FINALIZE(err)
          stop
         endif

         do i=1,iargc()
          call getarg(i,word)
          
          select case (trim(word))

             case ('-xyz')
                 call getarg(i+1,word)
                 read(word,*) inp_file

             case ('-cell')
                 call getarg(i+1,word)
                 read(word,*) cell_file

             case ('-atom')
                 call getarg(i+1,word)
                 read(word,*) atom

             case ('-Jmax')
                 call getarg(i+1,word)
                 read(word,*) Jmax

             case ('-rcut')
                 call getarg(i+1,word)
                 read(word,*) rcut

          end select

         enddo

         open(11,file=trim(cell_file))

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

         close(11)

         call sys%read_xyz(11,inp_file)
         call sys%build_neighbour_list(rcut)
         call sys%build_descriptors('BIS',Jmax,rcut)

         do i=1,size(sys%at_desc(atom)%desc)
          write(*,*) sys%at_desc(atom)%desc(i)
         enddo

         call sys%delete()
         call MPI_FINALIZE(err)

        return
        end program get_bis
