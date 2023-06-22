        module lammps_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use atoms_class
        implicit none

        type,extends(atoms_group)         :: lammps_obj
         type(C_ptr)                      :: lmp
         real(kind=dbl)                   :: cutoff
         integer                          :: twojmax
         contains
         procedure     :: initialize => init_lammps_obj
         procedure     :: finalize => finalize_lammps_obj
         procedure     :: setup => setup_lammps_obj
         procedure     :: get_desc => get_bis
         procedure     :: get_der_desc => get_der_bis
        end type lammps_obj

        contains
        
        subroutine init_lammps_obj(this)
        implicit none
        class(lammps_obj),intent(in)         :: this
                 
        call lammps_open_no_mpi("lmp -screen none -log log.simple",this%lmp)
        
        end subroutine init_lammps_obj

        subroutine finalize_lammps_obj(this)
        implicit none
        class(lammps_obj),intent(in)         :: this

        call lammps_close(this%lmp)

        end subroutine finalize_lammps_obj


        subroutine setup_lammps_obj(this,nkinds)
        implicit none
        class(lammps_obj),intent(in)          :: this
        integer,intent(in)                    :: nkinds
        character(len=100)                    :: a1,a2,a3
        character(len=300)                    :: region_string
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

       end subroutine setup_lammps_obj
 

       subroutine get_bis(this)
       implicit none
       class(lammps_obj),intent(inout)                       :: this
       integer                                               :: i,j,k,pos,m
       integer (C_int), dimension(:),pointer                 :: id
       integer                                               :: nlocal
       integer                                               :: components
       integer,allocatable                                   :: store_kind(:)
       real (C_double), dimension(:,:), pointer              :: bispec => NULL()
       real (C_double), dimension(:,:), pointer              :: der_bis => NULL()
       character(len=150)                                    :: cutoff_string,der_bis_string
       
       call lammps_command(this%lmp,"pair_style zero 30")
       call lammps_command(this%lmp,"pair_coeff * *")
       write(cutoff_string,*)'compute bispec all sna/atom',this%cutoff,'1',this%twojmax
       
       do i=1,this%nkinds
        cutoff_string=trim(cutoff_string)//' 0.5'
       end do
       
       do i=1,this%nkinds
        cutoff_string=trim(cutoff_string)//' 1'
       end do

       call lammps_command(this%lmp,trim(cutoff_string))
       call lammps_command(this%lmp,"run 0")
       call lammps_extract_atom(id,this%lmp,"id")
       call number_bispec(this%twojmax,components)
       allocate(this%at_desc(this%nats))
       
       call lammps_extract_compute(bispec,this%lmp,'bispec',LMP_STYLE_ATOM,LMP_TYPE_ARRAY)
       
       do k=1,this%nats
         
         allocate(this%at_desc(k)%desc(components))
         pos=FINDLOC(id,k,1)
         this%at_desc(k)%desc(1)=1.0
         
         do m=1,components-1
          this%at_desc(k)%desc(m+1)=bispec(m,pos)
         end do
       
       end do
       call lammps_command(this%lmp,"uncompute bispec")
       
       end subroutine get_bis
        
       subroutine get_der_bis(this)
       implicit none
       class(lammps_obj),intent(inout)                       :: this
       integer                                               :: i,j,k,pos,m
       integer (C_int), dimension(:),pointer                 :: id
       integer                                               :: nlocal
       integer                                               :: components
       integer,allocatable                                   :: store_kind(:)
       real (C_double), dimension(:,:), pointer              :: x,bis_comp,bispec,der_bis
       character(len=800)                                    :: cutoff_string,der_bis_string, &
                                                                mass_string

       call lammps_command(this%lmp,"pair_style zero 30")
       call lammps_command(this%lmp,"pair_coeff * *")
       write(der_bis_string,*)'compute der_bis all snad/atom',this%cutoff,'1',this%twojmax

       do i=1,this%nkinds
        der_bis_string=trim(der_bis_string)//' 0.5'
       end do

       do i=1,this%nkinds
        der_bis_string=trim(der_bis_string)//' 1'
       end do
       
       call lammps_command(this%lmp,trim(der_bis_string))
       call lammps_command(this%lmp,"run 0")
       call lammps_extract_compute(der_bis,this%lmp,'der_bis',LMP_STYLE_ATOM,LMP_TYPE_ARRAY)
       call lammps_extract_atom(id,this%lmp,"id")
             
       call number_bispec(this%twojmax,components)
       allocate(this%der_at_desc(this%nats))
 
       do k=1,this%nats

        allocate(this%der_at_desc(k)%desc(3*this%nkinds*(components-1)))
        pos=FINDLOC(id,k,1)
        this%der_at_desc(k)%desc(:)=der_bis(:,pos)
      
       end do

       call lammps_command(this%lmp,"uncompute der_bis")

       end subroutine get_der_bis

       subroutine import_lammps_obj_list(set,nconfig,file_input,len_file_inp)
       implicit none
       type(lammps_obj), allocatable        :: set(:)
       integer                              :: nconfig, i,j,nats
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
        
        end subroutine import_lammps_obj_list

        subroutine number_bispec(twojmax,components)
        implicit none
        integer,intent(in)            :: twojmax
        integer,intent(out)           :: components
        real(kind=dbl)                :: order_coeff

        if (modulo(twojmax,2)==0) then
         order_coeff=(twojmax/2.0)+1
         components=order_coeff*(order_coeff+1)*(2*order_coeff+1)/6.0
        else
         order_coeff=(twojmax+1)/2.0
         components=order_coeff*(order_coeff+1)*(order_coeff+2)/3.0
        end if

        components=components + 1

        end subroutine number_bispec

        end module
