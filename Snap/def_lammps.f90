        module lammps_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use atoms_class
        implicit none

        type,extends(atoms_group)         :: lammps_obj
         type(C_ptr)                      :: lmp
         contains
         !procedure     :: atoms_to_lammps_obj        
         procedure     :: setup_lammps_lattice
         !procedure     :: get_bis
         !procedure     :: get_der_bis
        end type lammps_obj

        contains

        subroutine setup_lammps_lattice(this,nkinds)
        use atoms_class
        implicit none
        class(lammps_obj),intent(in)          :: this
        integer,intent(in)                    :: nkinds
        character(len=100)                    :: a1,a2,a3
        character(len=200)                    :: region_string
        character(len=250)                    :: create_atoms_string,mass_string,create_box_string
        integer                               :: j

         write(a1,*) 'a1',this%cell(1,1),this%cell(1,2),this%cell(1,3)
         write(a2,*) 'a2',this%cell(2,1),this%cell(2,2),this%cell(2,3)
         write(a3,*) 'a3',this%cell(3,1),this%cell(3,2),this%cell(3,3)
         write(region_string,*) -this%origin(1),this%cell(1,1)-this%origin(1),-this%origin(2), &
         this%cell(2,2)-this%origin(2),-this%origin(3),this%cell(3,3)-this%origin(3),this%cell(2,1),&
          this%cell(3,1),this%cell(3,2)

         call lammps_command(this%lmp,"units real")
         call lammps_command(this%lmp,"dimension 3")
         call lammps_command(this%lmp,"boundary p p p")
         call lammps_command(this%lmp,"atom_style charge")
         call lammps_command(this%lmp,"neighbor 0.3 bin")
         call lammps_command(this%lmp,'lattice custom 1.0'&
          //trim(a1)//trim(a2)//trim(a3)//' basis 0 0 0')
         call lammps_command(this%lmp,'region id_1 prism'//trim(region_string)//' units box')
         call lammps_command(this%lmp,"atom_modify map yes")
         
         write(create_box_string,*) 'create_box',nkinds,'id_1'
        call lammps_command(this%lmp,trim(create_box_string))

        do j=1,this%nats

         create_atoms_string=""

         write(create_atoms_string,*)'create_atoms',this%kind(j),'single',(this%x(j,1)),(this%x(j,2)),(this%x(j,3)),&
           'remap yes units box'
         call lammps_command(this%lmp,trim(create_atoms_string))

        end do

        do j=1,this%nkinds

         write(mass_string,*)"mass",j,this%mass(j)
         call lammps_command(this%lmp,mass_string)

        end do

        if (nkinds.ne.this%nkinds) then
         do j=this%nkinds+1,nkinds
          write(mass_string,*)"mass",j,'1'
          call lammps_command(this%lmp,mass_string)
         end do
        end if


       end subroutine setup_lammps_lattice
        
       subroutine atoms_to_lammps_obj(set,nconfig,file_input,len_file_inp)
       type(lammps_obj), allocatable       :: set(:)
       integer                              :: nconfig, i,j,nats,ntypes
       character(len=100),allocatable       :: tmp(:,:)
       character(len=100),dimension(10)     :: tmp_cell_nkinds
       character(len=80)                    :: err_string
       character(len=50)                    :: filename
       character(len=*)                     :: file_input
       integer,intent(in)                   :: len_file_inp
       character(len=150)                   :: file_inp
       integer                              :: ierror    
       integer                              :: status
       character(len=80)                    :: err_msg

        file_inp=file_input(1:len_file_inp)
        allocate(set(nconfig))
        open(unit=1,file=trim(file_inp),status="old",action="read", iostat=ierror, iomsg=err_string)
        
        if (ierror /= 0) then
          write(*,*) err_string
        end if

        do i=1,nconfig
         read(1,*,iostat=status,iomsg=err_msg) nats
         
         if(status /= 0) then
          write(*,*) err_msg
        end if

         read(1,*)tmp_cell_nkinds(:)
         allocate(tmp(nats,6))

         do j=1,nats
          read(1,*) tmp(j,:)
         end do

         filename="file_temporaneo"
         open(unit=2, file = trim(filename), status="replace",action="write",iostat=ierror,iomsg=err_string)
         
         !if (ierror /= 0) then
          !write(*,*) err_string
         !end if         
         
          write(2,*)nats
          write(2,*)trim(tmp_cell_nkinds(1)),' ',trim(tmp_cell_nkinds(2)),' ',trim(tmp_cell_nkinds(3)),' '&
           ,trim(tmp_cell_nkinds(4)),' ',trim(tmp_cell_nkinds(5)),' ',trim(tmp_cell_nkinds(6)),' ',trim(tmp_cell_nkinds(7))&
           ,' ',trim(tmp_cell_nkinds(8)),' ',trim(tmp_cell_nkinds(9)),' ',trim(tmp_cell_nkinds(10))
          do j=1,nats
           write(2,*)trim(tmp(j,1)),' ',trim(tmp(j,2)),' ',trim(tmp(j,3)),' ',trim(tmp(j,4)),' ',trim(tmp(j,5)),&
            ' ',trim(tmp(j,6))
          end do

         deallocate(tmp)
         close(2)
         call set(i)%read_extended_xyz(3,trim(filename))
      
        end do
        close(1)

        end subroutine atoms_to_lammps_obj
        end module
