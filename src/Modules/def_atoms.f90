        module atoms_class         
        use lattice_class
        use general_types_class 
        use descriptors_class
        use ffs_class
        use mlmodels_class
        implicit none

        type fc_branch
         integer                       :: order
         integer                       :: nfcs
         double precision, allocatable :: val(:)
         integer, allocatable          :: nat(:,:)
         integer, allocatable          :: cell(:,:)
        end type fc_branch

        type, extends(bravais_lattice) :: atoms_group
         integer                         :: nats 
         integer                         :: nkinds
         integer, pointer                :: kind(:)         
         integer, allocatable            :: molid(:)
         logical, allocatable            :: good_ids(:)
         double precision, allocatable   :: x(:,:)         
         double precision, allocatable   :: v(:,:)         
         double precision, allocatable   :: mass(:)
         double precision, allocatable   :: charge(:)
         double precision, allocatable   :: dipole(:,:)
         double precision, allocatable   :: dist(:,:,:,:)
         double precision, allocatable   :: rij(:,:)
         type(vector_int), allocatable   :: neigh(:) 
         double precision                :: neigh_cutoff=3.7d0
         character(len=10), allocatable  :: label(:)
         logical                         :: born_charges=.false. 
         double precision, allocatable   :: fcs2(:,:,:,:,:)
         type(fc_branch), allocatable    :: fcs3
         type(descriptor), pointer       :: at_desc(:)
         type(force_field), pointer      :: FF
         type(mlmodel), pointer          :: ML
         contains
         procedure        :: delete => delete_atoms_group
         procedure        :: build_descriptors
         procedure        :: build_neighbour_list
         procedure        :: read_restart_file
         procedure        :: read_dipole_file
         procedure        :: read_structure_file
         procedure        :: read_extended_xyz
         procedure        :: read_xyz
         procedure        :: write_restart_file
         procedure        :: atoms_bcast
         procedure        :: dist_ij
         procedure        :: find_mols
         procedure        :: remove_overlap
         procedure        :: wrap_geo
         procedure        :: cart2frac
         procedure        :: frac2cart
         procedure        :: init_fcs2
         procedure        :: init_fcs3
        end type atoms_group

        contains

        subroutine remove_overlap(this)
        use lists_class
        implicit none
        class(atoms_group)               :: this
        double precision, allocatable    :: new_geo(:,:)
        double precision                 :: rad,bond_thr
        character(len=10), allocatable   :: new_labels(:)
        integer                          :: i,j

         call this%dist_ij()                   
         if(.not. allocated(this%good_ids)) allocate(this%good_ids(this%nats))
         this%good_ids=.true.  

         do i=1,this%nats
          if(.not. this%good_ids(i)) cycle
          do j=i+1,this%nats
           if(.not. this%good_ids(j)) cycle

           call get_cov_radius(this%label(this%kind(i)),rad)
           bond_thr=rad
           call get_cov_radius(this%label(this%kind(j)),rad)
           bond_thr=bond_thr+rad
           if(this%dist(i,1,j,1).lt.bond_thr-0.30d0)then
            this%good_ids(j)=.false.
           endif

          enddo
         enddo
           
!         allocate(new_geo(count(good_ids),3))
!         allocate(new_labels(count(good_ids)))
        
!         do i=1,this%nats
!          if(good_ids(i))then
!           write(*,*) this%label(this%kind(i)),this%x(i,:) 
!           new_geo(i,:)=this%x(i,:) 
!           new_labels(i)=this%label(this%kind(i)) 
!          endif
!         enddo

        return
        end subroutine remove_overlap

        subroutine read_dipole_file(this,dipole_file)
        use mpi
        use units_parms
        implicit none
        class(atoms_group)            :: this
        character(len=*)              :: dipole_file
        integer                       :: i,j
        integer                       :: err,mpi_nproc,mpi_id
           
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,err)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_id,err)

         if(mpi_id.eq.0)then

          write(*,*) 'Reading Dipole Derivarives:',dipole_file
          write(*,*)
          write(*,*)

          open(unit=12,file=dipole_file)
          allocate(this%dipole(3*this%nats,3))

          do i=1,3*this%nats
           read(12,*) (this%dipole(i,j),j=1,3)
          enddo

          close(12)

         endif

        return
        end subroutine read_dipole_file

        subroutine delete_atoms_group(this)
        implicit none
        class(atoms_group)               :: this
 
         if(associated(this%kind)) deallocate(this%kind)
         this%kind=>null()
         if(allocated(this%x)) deallocate(this%x)
         if(allocated(this%v)) deallocate(this%v)
         if(allocated(this%mass)) deallocate(this%mass)
         if(allocated(this%charge)) deallocate(this%charge)
         if(allocated(this%dist)) deallocate(this%dist)
         if(allocated(this%label)) deallocate(this%label)

          ! add removal of descriptors, neighbours, and FF

        return
        end subroutine delete_atoms_group

        subroutine build_neighbour_list(this,cutoff)
        use lists_class
        implicit none
        class(atoms_group)               :: this
        integer                          :: i,j
        double precision, optional       :: cutoff
        type(list)                       :: listid
        
         call this%dist_ij()

         if(allocated(this%neigh))then
          do i=1,size(this%neigh)
           if(allocated(this%neigh(i)%v)) deallocate(this%neigh(i)%v)
          enddo
          deallocate(this%neigh)
         endif
         allocate(this%neigh(this%nats))         

         if (present(cutoff)) this%neigh_cutoff=cutoff
         
         do i=1,this%nats
          call listid%init()
          do j=1,this%nats
           if(this%dist(i,1,j,1).lt.this%neigh_cutoff .and. i.ne.j) call listid%add_node(j)
          enddo
          allocate(this%neigh(i)%v(listid%nelem))
          call listid%reboot()
          do j=1,listid%nelem
           call listid%rd_val(this%neigh(i)%v(j))
           call listid%skip()
          enddo
          call listid%delete()
         enddo

        return
        end subroutine build_neighbour_list

        subroutine build_descriptors(this,kind_desc,Jmax,r0)
        implicit none
        class(atoms_group)               :: this
        integer                          :: i,j,s,id
        double precision, optional       :: r0
        integer, optional                :: Jmax
        double precision, allocatable    :: neighbours(:,:)
        character(len=*)                 :: kind_desc
        type(bispectrum)                 :: bis
        type(cartesian)                  :: cart

         if(trim(kind_desc).eq.'CART')then
          allocate(this%at_desc(1))
          call cart%get_desc(this%x)
          this%at_desc(1)%desc=cart%desc
          this%at_desc(1)%size_desc=cart%size_desc
          deallocate(cart%desc)            
         endif

         if(trim(kind_desc).eq.'BIS')then

          bis%BI%max_order=Jmax
          call bis%BI%setup()

          allocate(this%at_desc(this%nats))

          do i=1,this%nats
           allocate(neighbours(size(this%neigh(i)%v),3))
           do j=1,size(this%neigh(i)%v)
            id=this%neigh(i)%v(j)
            neighbours(j,:)=this%x(i,:)-this%x(id,:)
           enddo
           call bis%get_desc(neighbours,r0)
           this%at_desc(i)%desc=bis%desc
           this%at_desc(i)%size_desc=bis%size_desc
           deallocate(bis%desc)
           call bis%BI%reset()
           deallocate(neighbours)
          enddo

          call bis%BI%delete()

         endif

        return
        end subroutine build_descriptors

        subroutine read_extended_xyz(this,IOid,filename)
        implicit none
        class(atoms_group)            :: this
        integer                       :: i,j,IOid
        double precision, allocatable :: mass(:)
        character(len=100), optional  :: filename
        character(len=5), allocatable :: label(:)

         if(present(filename)) open(IOid,file=trim(filename))

          read(IOid,*) this%nats
          read(IOid,*) this%cell(1,:),this%cell(2,:),this%cell(3,:),&
                     this%nkinds
          allocate(this%x(this%nats,3)) 
          allocate(label(this%nats)) 
          allocate(mass(this%nats)) 
          allocate(this%label(this%nkinds)) 
          allocate(this%mass(this%nkinds)) 
          allocate(this%kind(this%nats)) 
         
         do i=1,this%nats
          read(IOid,*) label(i),this%x(i,:),this%kind(i),mass(i)
         enddo
         
         do i=1,this%nkinds
          do j=1,this%nats
           if(this%kind(j).eq.i)then
            this%mass(i)=mass(j)
            this%label(i)=trim(label(j))
            exit
           endif
          enddo
         enddo

         deallocate(mass)
         deallocate(label)
         
         if(present(filename)) close(IOid)

         this%touched=.true.
         call this%cell2abc()
         call this%init()
         call this%do_supercell()
         call this%dist_ij()
         call this%build_neighbour_list()

        return
        end subroutine read_extended_xyz        

        subroutine read_xyz(this,IOid,filename)
        implicit none
        class(atoms_group)            :: this
        integer                       :: i,j,IOid
        double precision, allocatable :: mass(:)
        character(len=100), optional  :: filename
        character(len=5), allocatable :: label(:)
        logical                       :: new_type

         if(present(filename)) open(IOid,file=trim(filename))

         read(IOid,*) this%nats
         read(IOid,*) 
         allocate(this%x(this%nats,3)) 
         allocate(label(this%nats)) 
         allocate(this%kind(this%nats)) 
         
         do i=1,this%nats
          read(IOid,*) label(i),this%x(i,:)
         enddo

         if(present(filename)) close(IOid)

         this%kind(1)=1
         this%nkinds=1

         do i=2,this%nats
          new_type=.true.
          do j=1,i-1

           if(label(i).eq.label(j))then
            this%kind(i)=this%kind(j)
            new_type=.false.
            exit
           endif

          enddo            
          if(new_type)then
           this%nkinds=this%nkinds+1
           this%kind(i)=this%nkinds
          endif
         enddo

         allocate(this%label(this%nkinds))

         do i=1,this%nkinds
          do j=1,this%nats 
           if(i.eq.this%kind(j))then
            this%label(i)=label(j)
            exit
           endif
          enddo
         enddo
         
         deallocate(label)

         allocate(this%mass(this%nkinds))
         this%mass=0.0d0
         do i=1,this%nkinds
          call get_mass(this%label(i),this%mass(i))
         enddo                

        return
        end subroutine read_xyz        

        subroutine wrap_geo(this,j)
        implicit none
        class(atoms_group)      :: this
        integer                 :: i,j,celli,cellj,v1,v2,v,nat0
        double precision        :: c(3),a(3),b(3)

         call this%cart2frac()
         
         do i=2,this%nats
          a=this%x(j,:)
          b=this%x(i,:)
          c(1)=a(1)-b(1)
          c(1)=nint(c(1)/dble(this%nx))*this%nx
          c(2)=a(2)-b(2)
          c(2)=nint(c(2)/dble(this%ny))*this%ny
          c(3)=a(3)-b(3)
          c(3)=nint(c(3)/dble(this%nz))*this%nz
          this%x(i,:)=this%x(i,:)+c(:)
         enddo

         call this%frac2cart()

        return
        end subroutine wrap_geo

        subroutine find_mols(this,reorder,remap)
        use lists_class
        use sparse_class
        implicit none
        class(atoms_group)            :: this
        integer                       :: i,j,celli,cellj,v1,v2,v,N,jj,l
        logical, allocatable          :: check(:)
        integer, allocatable          :: mapp(:),blc(:),NearNeigh(:),new_kind(:)
        double precision              :: rad,bond_thr,a(3),b(3),c(3),c2(3)
        double precision, allocatable :: new_geo(:,:)
        type(csr_mat_int)             :: CN
        type(list)                    :: AI,AJ,Aval
        logical                       :: print_flag,reorder,remap
        type(list)                    :: r,r2,queue
        integer                       :: dim_block,pos
        integer                       :: k,m,x
        class(*), pointer             :: arrow

         call this%dist_ij()

         N=this%ntot*this%nats

         allocate(NearNeigh(N))
         NearNeigh=0

         allocate(CN%AI(N+1))
         CN%AI(1)=0

         call AJ%init()
         call Aval%init()

         do v1=1,this%nats
         do celli=1,this%ntot
          i=(celli-1)*this%nats+v1

          CN%AI(i+1)=CN%AI(i)

          do v2=1,this%nats
          do cellj=1,this%ntot
           j=(cellj-1)*this%nats+v2

           if(j.eq.i)cycle
           call get_cov_radius(this%label(this%kind(v1)),rad)
           bond_thr=rad
           call get_cov_radius(this%label(this%kind(v2)),rad)
           bond_thr=bond_thr+rad
           if(this%dist(v1,celli,v2,cellj).lt.bond_thr+0.30d0)then
            CN%AI(i+1)=CN%AI(i+1)+1
            NearNeigh(i)=NearNeigh(i)+1
            call AJ%add_node(j)
            call Aval%add_node(1)
           endif

          enddo
          enddo
         enddo
         enddo

         call AJ%reboot()
         call Aval%reboot()
         CN%nzel=AJ%nelem
         allocate(CN%AJ(AJ%nelem))
         allocate(CN%A(AJ%nelem))

         do i=1,AJ%nelem
          call AJ%rd_val(CN%AJ(i))
          call Aval%rd_val(CN%A(i))
          call AJ%skip()
          call Aval%skip()
         enddo

         call AJ%delete()
         call Aval%delete()

         call this%cart2frac()

         call r%init()      ! elementi del blocco
         call r2%init()     ! dimensione blocco
         call queue%init()  
                   
         N=size(CN%AI)-1
         allocate(check(N))
         check=.false.
       
        ! k scorre su gli elementi non nulli

         k=1

         do while ( .not. all(check) .or. k.le.N  )      
        
          if( check(k) ) then

           k=k+1

          else

           check(k)=.true.
           call queue%add_node(k)
           call r%add_node(k)
           dim_block=1
                                 
           do while ( queue%nelem .gt. 0  ) 

!        prendi il primo in lista e cerca i vicini

            select type (arrow=>queue%head%key) 
            type is (integer)
             v=arrow
            end select
 
            do i=1,CN%AI(v+1)-CN%AI(v)

             pos=CN%AJ(i+CN%AI(v))
                
             if(check(pos))then

             else
              check(pos)=.true.
              call queue%add_node(pos)
              call r%add_node(pos)
              dim_block=dim_block+1          

              if( remap )then
               a=this%x(v,:)
               b=this%x(pos,:)
               c(1)=a(1)-b(1)
               c(2)=a(2)-b(2)
               c(3)=a(3)-b(3)
               c(1)=nint(c(1)/dble(this%nx))*this%nx
               c(2)=nint(c(2)/dble(this%ny))*this%ny
               c(3)=nint(c(3)/dble(this%nz))*this%nz
               this%x(pos,:)=this%x(pos,:)+c(:)
              endif

             endif
            enddo

            call queue%rm
           enddo ! while queue .ne. 0

           k=k+1
        
           call r2%add_node(dim_block)

          endif
         enddo ! while check

!!!     compatta liste e crea matrice permutazioni e determina
!!!     dimensioni blocchi

         check=.true.

         allocate(mapp(N))
         allocate(blc(r2%nelem+1))

         blc(1)=0

         call r2%reboot
         call r%reboot       

         m=1

         do v=1,r2%nelem
  
          blc(v+1)=blc(v)

          select type (arrow=>r2%node%key)
          type is (integer)
           l=arrow
          end select

          call r2%skip        

          do k=1,l
                         
           select type (arrow=>r%node%key)
           type is (integer)
            x=arrow
           end select

           call r%skip 

           if ( check(x) ) then
            mapp(m)=x
            check(x)=.false.
            m=m+1
            blc(v+1)=blc(v+1)+1
           endif

          enddo ! ciclo su blocco
         enddo ! ciclo su graphs
       
         if( allocated(CN%AC)) deallocate(CN%AC)
         call r%delete()
         call r2%delete()
         call queue%delete()

         call this%frac2cart()


         if ( reorder ) then

          allocate(new_geo(this%nats,3))
          allocate(new_kind(this%nats))

          if(allocated(this%molid))deallocate(this%molid)
          allocate(this%molid(size(blc)))
          this%molid=blc

          v=1
          do j=1,size(blc)-1
           do i=1+blc(j),blc(j+1)
            new_geo(v,:)=this%x(mapp(i),:)
            new_kind(v)=this%kind(mapp(i))
 !           this%molid(v)=j
            v=v+1
           enddo
          enddo

          this%x=new_geo
          this%kind=new_kind
          call this%dist_ij()

         endif

        return
        end subroutine find_mols

        subroutine dist_ij(this)
        implicit none
        class(atoms_group)      :: this
        integer                 :: i,j,celli,cellj,v1,v2,v
        double precision        :: c(3),a(3),b(3)

         if(allocated(this%dist)) deallocate(this%dist)
         allocate(this%dist(this%nats,this%ntot,this%nats,this%ntot))

         call this%cart2frac()
         do i=1,this%nats
          do j=1,this%nats
           do celli=1,this%ntot
            do cellj=1,this%ntot
             v1=(celli-1)*this%nats+i
             v2=(cellj-1)*this%nats+j
             do v=1,3
              a(v)=this%x(i,v)+this%rcell(celli,v)
              b(v)=this%x(j,v)+this%rcell(cellj,v)
             enddo
             c(1)=a(1)-b(1)
             c(1)=c(1)-nint(c(1)/dble(this%nx))*this%nx
             c(2)=a(2)-b(2)
             c(2)=c(2)-nint(c(2)/dble(this%ny))*this%ny
             c(3)=a(3)-b(3)
             c(3)=c(3)-nint(c(3)/dble(this%nz))*this%nz
             this%dist(i,celli,j,cellj)=this%dist_dir(c)
            enddo
           enddo
          enddo
         enddo
         call this%frac2cart()

        return
        end subroutine dist_ij

        subroutine atoms_bcast(this)
        use MPI
        use mpi_utils
        implicit none
        class(atoms_group)   :: this
        integer              :: ncell,i,i1,i2,s1,s2
        integer              :: nfcs
        logical              :: read_fc3=.false.


         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)       

         if(mpi_nproc.gt.1)then

          call this%lattice_bcast()        

          call MPI_BCAST(this%nats,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%nkinds,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%born_charges,1,mpi_logical,0,mpi_comm_world,err)

          if(.not.allocated(this%x))then
           allocate(this%x(this%nats,3))
          endif

          do i=1,this%nats
           call MPI_BCAST(this%x(i,:),3,mpi_double_precision,0,mpi_comm_world,err)
          enddo

          if(.not.associated(this%kind))then
           allocate(this%kind(this%nats))
          endif

          call MPI_BCAST(this%kind,this%nats,mpi_integer,0,mpi_comm_world,err)

          if(.not.allocated(this%label))then
           allocate(this%label(this%nkinds))
          endif

          call MPI_BCAST(this%label,10,mpi_char,0,mpi_comm_world,err)

          if(.not.allocated(this%mass))then
           allocate(this%mass(this%nkinds))
          endif

          call MPI_BCAST(this%mass,this%nkinds,mpi_double_precision,0,mpi_comm_world,err)
          
          if(.not.allocated(this%fcs2))then
           call this%init_fcs2()
          endif
        
          ncell=this%nx*this%ny*this%nz

          do i1=1,this%nats
           do i2=1,this%nats
            do s1=1,3
             do s2=1,3
              call MPI_BCAST(this%fcs2(:,i1,s1,i2,s2),ncell,mpi_double_precision,0,mpi_comm_world,err)
             enddo
            enddo
           enddo
          enddo

          if(mpi_id.eq.0)then
           if(allocated(this%fcs3))then
            read_fc3=.true.           
            nfcs=this%fcs3%nfcs
           endif
          endif

          call mpi_bcast(read_fc3,1,mpi_logical,0,mpi_comm_world,err)  
           
          if(read_fc3)then

           call mpi_bcast(nfcs,1,mpi_integer,0,mpi_comm_world,err)      

           if(.not.allocated(this%fcs3))then
            call this%init_fcs3(nfcs)
           endif

           do i=1,this%fcs3%nfcs
            call MPI_BCAST(this%fcs3%val,this%fcs3%nfcs,mpi_double_precision,0,mpi_comm_world,err) 
           enddo

           call MPI_BCAST(this%fcs3%nat(:,1),this%fcs3%nfcs,mpi_integer,0,mpi_comm_world,err) 
           call MPI_BCAST(this%fcs3%nat(:,2),this%fcs3%nfcs,mpi_integer,0,mpi_comm_world,err) 
           call MPI_BCAST(this%fcs3%nat(:,3),this%fcs3%nfcs,mpi_integer,0,mpi_comm_world,err) 
           call MPI_BCAST(this%fcs3%cell(:,1),this%fcs3%nfcs,mpi_integer,0,mpi_comm_world,err) 
           call MPI_BCAST(this%fcs3%cell(:,2),this%fcs3%nfcs,mpi_integer,0,mpi_comm_world,err) 

          endif

         endif

        return       
        end subroutine atoms_bcast


        subroutine  init_fcs3(this,nfcs)
        implicit none
        class(atoms_group) :: this
        integer            :: nfcs
         allocate(this%fcs3)
         this%fcs3%nfcs=nfcs
         this%fcs3%order=3
         allocate(this%fcs3%val(nfcs))
         allocate(this%fcs3%cell(nfcs,2))
         allocate(this%fcs3%nat(nfcs,3))
         this%fcs3%val=0.0d0
        return
        end subroutine init_fcs3


        subroutine  init_fcs2(this)
        implicit none
        class(atoms_group) :: this
        integer            :: ncells
         ncells=this%nx*this%ny*this%nz
         allocate(this%fcs2(ncells,this%nats,3,this%nats,3))
         this%fcs2=0.0d0
        return
        end subroutine init_fcs2

        subroutine frac2cart(this)
        implicit none
        class(atoms_group) :: this
        integer            :: i
         do i=1,this%nats
          this%x(i,:)=matmul(this%J,this%x(i,:))
         enddo
        return
        end subroutine frac2cart

        subroutine cart2frac(this)
        implicit none
        class(atoms_group) :: this
        integer            :: i
         do i=1,this%nats
          this%x(i,:)=matmul(this%Jinv,this%x(i,:))
         enddo
        return
        end subroutine cart2frac

        subroutine read_structure_file(this,structure_file)
        use mpi
        use mpi_utils
        use units_parms
        implicit none
        class(atoms_group)            :: this
        character(len=50)             :: skip
        character(len=*)              :: structure_file
        integer                       :: i,i1,i2,i3,s1,s2,s3,v1,v2,v3,counter
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)


         if(mpi_id.eq.0)then

         write(*,*) 'Reading Structure File:',structure_file
         write(*,*)
         write(*,*)

         open(unit=11,file=structure_file)

         read(11,*) this%nkinds,this%nats
         write(*,*) '   Detected ',this%nkinds,' kinds and ',&
                    this%nats,' atoms'
         write(*,*) 
         write(*,*)

         allocate(this%x(this%nats,3))
         allocate(this%kind(this%nats))

         allocate(this%label(this%nkinds))
         allocate(this%mass(this%nkinds))

         read(11,*) (this%cell(1,i),i=1,3)
         read(11,*) (this%cell(2,i),i=1,3)
         read(11,*) (this%cell(3,i),i=1,3)

         this%cell=this%cell*bohr2ang

         do i=1,this%nkinds
          read(11,*)  skip,this%label(i),this%mass(i)
         enddo
         
         do i=1,this%nats
          read(11,*)  skip,this%kind(i),this%x(i,1),this%x(i,2),this%x(i,3)
         enddo

         this%touched=.true.
         call this%cell2abc()

         write(*,*) '#################################'
         write(*,*) 'Bravais Lattice Properties:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)
        
         write(*,*) '   Unit Cell Vectors (Angstrom):'
         write(*,*) '      ',this%cell(1,:)
         write(*,*) '      ',this%cell(2,:)
         write(*,*) '      ',this%cell(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Unit Cell Parameters (Angstrom, Deg):'
         write(*,*) '        a=',this%a
         write(*,*) '        b=',this%b
         write(*,*) '        c=',this%c
         write(*,*) '    alpha=',this%alpha*180.0d0/acos(-1.0d0)
         write(*,*) '     beta=',this%beta*180.0d0/acos(-1.0d0)
         write(*,*) '    gamma=',this%gamma*180.0d0/acos(-1.0d0)
         write(*,*)
         write(*,*)
        
         call this%init()

         write(*,*) '   Jacobian Tensor:'
         write(*,*) '      ',this%J(1,:)
         write(*,*) '      ',this%J(2,:)
         write(*,*) '      ',this%J(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Inverse Jacobian Tensor:'
         write(*,*) '      ',this%Jinv(1,:)
         write(*,*) '      ',this%Jinv(2,:)
         write(*,*) '      ',this%Jinv(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Covariant Metric Tensor:'
         write(*,*) '      ',this%Gcov(1,:)
         write(*,*) '      ',this%Gcov(2,:)
         write(*,*) '      ',this%Gcov(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Controvariant Metric Tensor:'
         write(*,*) '      ',this%Gcon(1,:)
         write(*,*) '      ',this%Gcon(2,:)
         write(*,*) '      ',this%Gcon(3,:)
         write(*,*)
         write(*,*)

         call this%frac2cart()

         write(*,*) '#################################'
         write(*,*) 'Unit Cell Atomic Kinds:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)

         write(*,*)
         do i=1,this%nkinds
          write(*,*) '      ','kind=',this%label(i),'mass=',this%mass(i)
         enddo
         write(*,*)
         write(*,*)
       
         write(*,*) '#################################'
         write(*,*) 'Unit Cell Cartesian Coordinates:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)

         write(*,*)
         do i=1,this%nats   
          write(*,*) '      ',this%label(this%kind(i)),this%x(i,1:3)
         enddo
         write(*,*)
         write(*,*)

         read(11,*) this%born_charges

         if(this%born_charges)then
          write(*,*) 'Input with Born_Charges not implemented yet'
          stop
         endif

         read(11,*) this%nx,this%ny,this%nz
         write(*,*) '   Detected Real-Space Supercell:',this%nx,this%ny,this%nz 
         write(*,*) 
         write(*,*)

         call this%do_supercell()

         close(11)

         endif

        return
        end subroutine read_structure_file

        subroutine read_restart_file(this,restart_file,anharmonic_file)
        use mpi
        use mpi_utils
        use units_parms
        implicit none
        class(atoms_group)            :: this
        character(len=50)             :: skip
        character(len=*)              :: restart_file
        character(len=*),optional     :: anharmonic_file
        integer                       :: i,i1,i2,i3,s1,s2,s3,v1,v2,v3,counter
        logical                       :: read_fc3=.false.
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)


         if(mpi_id.eq.0)then

         if(present(anharmonic_file)) read_fc3=.true.

         write(*,*) 'Reading Structure File:',restart_file
         write(*,*)
         write(*,*)

         open(unit=11,file=restart_file)

         read(11,*) this%nkinds,this%nats
         write(*,*) '   Detected ',this%nkinds,' kinds and ',&
                    this%nats,' atoms'
         write(*,*) 
         write(*,*)

         allocate(this%x(this%nats,3))
         allocate(this%kind(this%nats))

         allocate(this%label(this%nkinds))
         allocate(this%mass(this%nkinds))

         read(11,*) (this%cell(1,i),i=1,3)
         read(11,*) (this%cell(2,i),i=1,3)
         read(11,*) (this%cell(3,i),i=1,3)

         this%cell=this%cell*bohr2ang

         do i=1,this%nkinds
          read(11,*)  skip,this%label(i),this%mass(i)
         enddo
         
         do i=1,this%nats
          read(11,*)  skip,this%kind(i),this%x(i,1),this%x(i,2),this%x(i,3)
         enddo


         this%touched=.true.
         call this%cell2abc()

         write(*,*) '#################################'
         write(*,*) 'Bravais Lattice Properties:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)
        
         write(*,*) '   Unit Cell Vectors (Angstrom):'
         write(*,*) '      ',this%cell(1,:)
         write(*,*) '      ',this%cell(2,:)
         write(*,*) '      ',this%cell(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Unit Cell Parameters (Angstrom, Deg):'
         write(*,*) '        a=',this%a
         write(*,*) '        b=',this%b
         write(*,*) '        c=',this%c
         write(*,*) '    alpha=',this%alpha*180.0d0/acos(-1.0d0)
         write(*,*) '     beta=',this%beta*180.0d0/acos(-1.0d0)
         write(*,*) '    gamma=',this%gamma*180.0d0/acos(-1.0d0)
         write(*,*)
         write(*,*)
        
         call this%init()

         write(*,*) '   Jacobian Tensor:'
         write(*,*) '      ',this%J(1,:)
         write(*,*) '      ',this%J(2,:)
         write(*,*) '      ',this%J(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Inverse Jacobian Tensor:'
         write(*,*) '      ',this%Jinv(1,:)
         write(*,*) '      ',this%Jinv(2,:)
         write(*,*) '      ',this%Jinv(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Covariant Metric Tensor:'
         write(*,*) '      ',this%Gcov(1,:)
         write(*,*) '      ',this%Gcov(2,:)
         write(*,*) '      ',this%Gcov(3,:)
         write(*,*)
         write(*,*)

         write(*,*) '   Controvariant Metric Tensor:'
         write(*,*) '      ',this%Gcon(1,:)
         write(*,*) '      ',this%Gcon(2,:)
         write(*,*) '      ',this%Gcon(3,:)
         write(*,*)
         write(*,*)

         call this%frac2cart()

         write(*,*) '#################################'
         write(*,*) 'Unit Cell Atomic Kinds:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)

         write(*,*)
         do i=1,this%nkinds
          write(*,*) '      ','kind=',this%label(i),'mass=',this%mass(i)
         enddo
         write(*,*)
         write(*,*)
       
         write(*,*) '#################################'
         write(*,*) 'Unit Cell Cartesian Coordinates:'
         write(*,*) '#################################'
         write(*,*)
         write(*,*)

         write(*,*)
         do i=1,this%nats   
          write(*,*) '      ',this%label(this%kind(i)),this%x(i,1:3)
         enddo
         write(*,*)
         write(*,*)

         read(11,*) this%born_charges

         if(this%born_charges)then
          write(*,*) 'Input with Born_Charges not implemented yet'
          stop
         endif

         read(11,*) this%nx,this%ny,this%nz
         write(*,*) '   Detected Real-Space Supercell:',this%nx,this%ny,this%nz 
         write(*,*) 
         write(*,*)

         call this%do_supercell()
         call this%init_fcs2()         

         do s1=1,3
         do s2=1,3

          do v1=1,this%nats
          do v2=1,this%nats

           read(11,*) 
           counter=1

           do i1=1,this%nz
           do i2=1,this%ny
           do i3=1,this%nx
            read(11,*) skip,skip,skip,this%fcs2(counter,v1,s1,v2,s2)
            this%fcs2(counter,v1,s1,v2,s2)=this%fcs2(counter,v1,s1,v2,s2)/2.0d0
            counter=counter+1
           enddo
           enddo
           enddo

          enddo
          enddo

         enddo
         enddo

         close(11)

         write(*,*) '   2nd Order Force Constants read correctly.'
         write(*,*) 
         write(*,*)

         if(read_fc3)then

        
          write(*,*) 'Reading 3rd order FCs from File:',anharmonic_file
          write(*,*)
          write(*,*)

          open(unit=11,file=anharmonic_file)

          read(11,*) counter
          call this%init_fcs3(counter)         

          do i=1,counter
         
           read(11,*) this%fcs3%cell(i,1), &
                      this%fcs3%cell(i,2), &
                      this%fcs3%nat(i,1),  &
                      this%fcs3%nat(i,2),  &
                      this%fcs3%nat(i,3),  &
                      this%fcs3%val(i)   

          enddo

          close(11)

          write(*,*) '   3nd Order Force Constants read correctly.'
          write(*,*) 
          write(*,*)

         end if ! on read fc3

         endif


        return 
        end subroutine read_restart_file
       
        subroutine write_restart_file(this,restart_file)
        use mpi
        use mpi_utils
        use units_parms
        implicit none
        class(atoms_group)            :: this
        character(len=50)             :: skip
        character(len=*)              :: restart_file
        integer                       :: i,i1,i2,i3,s1,s2,s3,v1,v2,v3,counter
           
         call MPI_COMM_SIZE(mpi_comm_world,mpi_nproc,err)
         call MPI_COMM_RANK(mpi_comm_world,mpi_id,err)

         if(mpi_id.eq.0)then

         open(unit=11,file=restart_file)

         write(11,*) this%nkinds,this%nats

         write(11,*) (this%cell(1,i)/bohr2ang,i=1,3)
         write(11,*) (this%cell(2,i)/bohr2ang,i=1,3)
         write(11,*) (this%cell(3,i)/bohr2ang,i=1,3)

         do i=1,this%nkinds
          write(11,*)  i,this%label(i),this%mass(i)
         enddo
         
         call this%cart2frac()

         do i=1,this%nats
          write(11,*)  i,this%kind(i),this%x(i,1),this%x(i,2),this%x(i,3)
         enddo

         write(11,*) this%born_charges

         write(11,*) this%nx,this%ny,this%nz


         do s1=1,3
         do s2=1,3

          do v1=1,this%nats
          do v2=1,this%nats

           write(11,*) s1,s2,v1,v2 
           counter=1

           do i1=1,this%nz
           do i2=1,this%ny
           do i3=1,this%nx
            write(11,*) i1,i2,i3,this%fcs2(counter,v1,s1,v2,s2)*2.0d0
            counter=counter+1
           enddo
           enddo
           enddo

          enddo
          enddo

         enddo
         enddo

         close(11)

         write(*,*) '   Restart File correctly wrote.'
         write(*,*) 
         write(*,*)

         endif

        return 
        end subroutine write_restart_file

        subroutine get_mass(label,mass)
        implicit none
        character(len=2)      :: label
        double precision      :: mass

         if(trim(label).eq.'H') mass=1.00784d0
         if(trim(label).eq.'He') mass=4.002602d0
         if(trim(label).eq.'Li') mass=6.938d0
         if(trim(label).eq.'Be') mass=9.0121831d0
         if(trim(label).eq.'B') mass=10.806d0
         if(trim(label).eq.'C') mass=12.0096d0
         if(trim(label).eq.'N') mass=14.00643d0
         if(trim(label).eq.'O') mass=15.99903d0
         if(trim(label).eq.'F') mass=18.998403163d0
         if(trim(label).eq.'Ne') mass=20.1797d0
         if(trim(label).eq.'Na') mass=22.98976928d0
         if(trim(label).eq.'Mg') mass=24.304d0
         if(trim(label).eq.'Al') mass=26.9815385d0
         if(trim(label).eq.'Si') mass=28.084d0
         if(trim(label).eq.'P') mass=30.973761998d0
         if(trim(label).eq.'S') mass=32.059d0
         if(trim(label).eq.'Cl') mass=35.446d0
         if(trim(label).eq.'Ar') mass=39.948d0
         if(trim(label).eq.'K') mass=39.0983d0
         if(trim(label).eq.'Ca') mass=40.078d0
         if(trim(label).eq.'Sc') mass=44.955908d0
         if(trim(label).eq.'Ti') mass=47.867d0
         if(trim(label).eq.'V') mass=50.9415d0
         if(trim(label).eq.'Cr') mass=51.9961d0
         if(trim(label).eq.'Mn') mass=54.938044d0
         if(trim(label).eq.'Fe') mass=55.845d0
         if(trim(label).eq.'Co') mass=58.933194d0
         if(trim(label).eq.'Ni') mass=58.6934d0
         if(trim(label).eq.'Cu') mass=63.546d0
         if(trim(label).eq.'Zn') mass=65.38d0
         if(trim(label).eq.'Ga') mass=69.723d0
         if(trim(label).eq.'Ge') mass=72.630d0
         if(trim(label).eq.'As') mass=74.921595d0
         if(trim(label).eq.'Se') mass=78.971d0
         if(trim(label).eq.'Br') mass=79.901d0
         if(trim(label).eq.'Kr') mass=83.798d0
         if(trim(label).eq.'Rb') mass=85.4678d0
         if(trim(label).eq.'Sr') mass=87.62d0
         if(trim(label).eq.'Y') mass=88.90584d0
         if(trim(label).eq.'Zr') mass=91.224d0
         if(trim(label).eq.'Nb') mass=92.90637d0
         if(trim(label).eq.'Mo') mass=95.95d0
         if(trim(label).eq.'Tc') mass=98.0d0
         if(trim(label).eq.'Ru') mass=101.07d0
         if(trim(label).eq.'Rh') mass=102.90550d0
         if(trim(label).eq.'Pd') mass=106.42d0
         if(trim(label).eq.'Ag') mass=107.8682d0
         if(trim(label).eq.'Cd') mass=112.414d0
         if(trim(label).eq.'In') mass=114.818d0
         if(trim(label).eq.'Sn') mass=118.710d0
         if(trim(label).eq.'Sb') mass=121.760d0
         if(trim(label).eq.'Te') mass=127.60d0
         if(trim(label).eq.'I') mass=126.90447d0
         if(trim(label).eq.'Xe') mass=131.293d0
         if(trim(label).eq.'Cs') mass=132.90545196d0
         if(trim(label).eq.'Ba') mass=137.327d0
         if(trim(label).eq.'La') mass=138.90547d0
         if(trim(label).eq.'Ce') mass=140.116d0
         if(trim(label).eq.'Pr') mass=140.90766d0
         if(trim(label).eq.'Nd') mass=144.242d0
         if(trim(label).eq.'Pm') mass=145.0d0
         if(trim(label).eq.'Sm') mass=150.36d0
         if(trim(label).eq.'Eu') mass=151.964d0
         if(trim(label).eq.'Gd') mass=157.25d0
         if(trim(label).eq.'Tb') mass=158.92535d0
         if(trim(label).eq.'Dy') mass=162.500d0
         if(trim(label).eq.'Ho') mass=164.93033d0
         if(trim(label).eq.'Er') mass=167.259d0
         if(trim(label).eq.'Tm') mass=168.93422d0
         if(trim(label).eq.'Yb') mass=173.054d0
         if(trim(label).eq.'Lu') mass=174.9668d0
         if(trim(label).eq.'Hf') mass=178.49d0
         if(trim(label).eq.'Ta') mass=180.94788d0
         if(trim(label).eq.'W') mass=183.84d0
         if(trim(label).eq.'Re') mass=186.207d0
         if(trim(label).eq.'Os') mass=190.23d0
         if(trim(label).eq.'Ir') mass=192.217d0
         if(trim(label).eq.'Pt') mass=195.084d0
         if(trim(label).eq.'Au') mass=196.966569d0
         if(trim(label).eq.'Hg') mass=200.592d0
         if(trim(label).eq.'Tl') mass=204.382d0
         if(trim(label).eq.'Pb') mass=207.2d0
         if(trim(label).eq.'Bi') mass=208.98040d0
         if(trim(label).eq.'Po') mass=209.0d0
         if(trim(label).eq.'Rn') mass=222.0d0
         if(trim(label).eq.'Fr') mass=223.0d0
         if(trim(label).eq.'Ra') mass=226.0d0
         if(trim(label).eq.'Th') mass=232.0d0
         if(trim(label).eq.'Pa') mass=231.03588d0
         if(trim(label).eq.'U') mass=238.02891d0

        return
        end subroutine get_mass

        subroutine get_cov_radius(label,radius)
        implicit none
        character(len=2)      :: label
        double precision      :: radius

         if(trim(label).eq.'H') radius=0.31d0
         if(trim(label).eq.'He') radius=0.28d0
         if(trim(label).eq.'Li') radius=1.28d0
         if(trim(label).eq.'Be') radius=0.96d0
         if(trim(label).eq.'B') radius=0.84d0
         if(trim(label).eq.'C') radius=0.76d0
         if(trim(label).eq.'N') radius=0.71d0
         if(trim(label).eq.'O') radius=0.66d0
         if(trim(label).eq.'F') radius=0.57d0
         if(trim(label).eq.'Ne') radius=0.58d0
         if(trim(label).eq.'Na') radius=1.66d0
         if(trim(label).eq.'Mg') radius=1.41d0
         if(trim(label).eq.'Al') radius=1.21d0
         if(trim(label).eq.'Si') radius=1.11d0
         if(trim(label).eq.'P') radius=1.07d0
         if(trim(label).eq.'S') radius=1.05d0
         if(trim(label).eq.'Cl') radius=1.02d0
         if(trim(label).eq.'Ar') radius=1.06d0
         if(trim(label).eq.'K') radius=2.03d0
         if(trim(label).eq.'Ca') radius=1.76d0
         if(trim(label).eq.'Sc') radius=1.70d0
         if(trim(label).eq.'Ti') radius=1.60d0
         if(trim(label).eq.'V') radius=1.53d0
         if(trim(label).eq.'Cr') radius=1.39d0
         if(trim(label).eq.'Mn') radius=1.61d0
         if(trim(label).eq.'Fe') radius=1.52d0
         if(trim(label).eq.'Co') radius=1.50d0
         if(trim(label).eq.'Ni') radius=1.24d0
         if(trim(label).eq.'Cu') radius=1.32d0
         if(trim(label).eq.'Zn') radius=1.22d0
         if(trim(label).eq.'Ga') radius=1.22d0
         if(trim(label).eq.'Ge') radius=1.20d0
         if(trim(label).eq.'As') radius=1.19d0
         if(trim(label).eq.'Se') radius=1.20d0
         if(trim(label).eq.'Br') radius=1.20d0
         if(trim(label).eq.'Kr') radius=1.16d0
         if(trim(label).eq.'Rb') radius=2.20d0
         if(trim(label).eq.'Sr') radius=1.95d0
         if(trim(label).eq.'Y') radius=1.90d0
         if(trim(label).eq.'Zr') radius=1.75d0
         if(trim(label).eq.'Nb') radius=1.64d0
         if(trim(label).eq.'Mo') radius=1.54d0
         if(trim(label).eq.'Tc') radius=1.47d0
         if(trim(label).eq.'Ru') radius=1.46d0
         if(trim(label).eq.'Rh') radius=1.42d0
         if(trim(label).eq.'Pd') radius=1.39d0
         if(trim(label).eq.'Ag') radius=1.45d0
         if(trim(label).eq.'Cd') radius=1.44d0
         if(trim(label).eq.'In') radius=1.42d0
         if(trim(label).eq.'Sn') radius=1.39d0
         if(trim(label).eq.'Sb') radius=1.39d0
         if(trim(label).eq.'Te') radius=1.38d0
         if(trim(label).eq.'I') radius=1.39d0
         if(trim(label).eq.'Xe') radius=1.40d0
         if(trim(label).eq.'Cs') radius=2.44d0
         if(trim(label).eq.'Ba') radius=2.15d0
         if(trim(label).eq.'La') radius=2.07d0
         if(trim(label).eq.'Ce') radius=2.04d0
         if(trim(label).eq.'Pr') radius=2.03d0
         if(trim(label).eq.'Nd') radius=2.01d0
         if(trim(label).eq.'Pm') radius=1.99d0
         if(trim(label).eq.'Sm') radius=1.98d0
         if(trim(label).eq.'Eu') radius=1.98d0
         if(trim(label).eq.'Gd') radius=1.96d0
         if(trim(label).eq.'Tb') radius=1.94d0
         if(trim(label).eq.'Dy') radius=1.92d0
         if(trim(label).eq.'Ho') radius=1.92d0
         if(trim(label).eq.'Er') radius=1.89d0
         if(trim(label).eq.'Tm') radius=1.90d0
         if(trim(label).eq.'Yb') radius=1.87d0
         if(trim(label).eq.'Lu') radius=1.87d0
         if(trim(label).eq.'Hf') radius=1.75d0
         if(trim(label).eq.'Ta') radius=1.70d0
         if(trim(label).eq.'W') radius=1.62d0
         if(trim(label).eq.'Re') radius=1.51d0
         if(trim(label).eq.'Os') radius=1.44d0
         if(trim(label).eq.'Ir') radius=1.41d0
         if(trim(label).eq.'Pt') radius=1.36d0
         if(trim(label).eq.'Au') radius=1.36d0
         if(trim(label).eq.'Hg') radius=1.32d0
         if(trim(label).eq.'Tl') radius=1.45d0
         if(trim(label).eq.'Pb') radius=1.46d0
         if(trim(label).eq.'Bi') radius=1.48d0
         if(trim(label).eq.'Po') radius=1.40d0
         if(trim(label).eq.'At') radius=1.50d0
         if(trim(label).eq.'Rn') radius=1.50d0
         if(trim(label).eq.'Fr') radius=2.60d0
         if(trim(label).eq.'Ra') radius=2.21d0
         if(trim(label).eq.'Ac') radius=2.15d0
         if(trim(label).eq.'Th') radius=2.06d0
         if(trim(label).eq.'Pa') radius=2.00d0
         if(trim(label).eq.'U') radius=1.96d0
         if(trim(label).eq.'Np') radius=1.90d0
         if(trim(label).eq.'Pu') radius=1.87d0
         if(trim(label).eq.'Am') radius=1.80d0
         if(trim(label).eq.'Cm') radius=1.69d0
         
        return
        end subroutine get_cov_radius

        end module atoms_class
