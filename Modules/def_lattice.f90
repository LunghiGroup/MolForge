        module lattice_class         
        implicit none

        type bravais_lattice
         logical          :: touched=.false.
         double precision :: a=0.0d0
         double precision :: b=0.0d0
         double precision :: c=0.0d0
         double precision :: alpha=0.0d0
         double precision :: beta=0.0d0
         double precision :: gamma=0.0d0
         double precision :: vol=0.0d0
         double precision :: cell(3,3)=0.0d0
         double precision, allocatable :: origin (:)
         integer          :: nx=1
         integer          :: ny=1
         integer          :: nz=1
         integer          :: ntot=1
         double precision, allocatable :: rcell(:,:)
         double precision :: J(3,3)=0.0d0
         double precision :: Jinv(3,3)=0.0d0
         double precision :: Gcon(3,3)=0.0d0
         double precision :: Gcov(3,3)=0.0d0
         contains
         procedure        :: init => init_lattice
         procedure        :: read_cell
         procedure        :: lattice_bcast
         procedure        :: abc2cell
         procedure        :: cell2abc
         procedure        :: calc_Js
         procedure        :: calc_Gs
         procedure        :: calc_vol
         procedure        :: do_supercell
         procedure        :: dist_rec
         procedure        :: dist_dir
         procedure        :: dist_pbc
         procedure        :: dist_vec_pbc
        end type bravais_lattice

        contains

        subroutine read_cell(this,cell)
        implicit none
        class(bravais_lattice)   :: this
        double precision         :: cell(3,3)

         this%cell=cell
         this%touched=.true.
         call this%init()
         call this%cell2abc

        return
        end subroutine read_cell

        subroutine lattice_bcast(this)
        use MPI
        use mpi_utils
        implicit none
        class(bravais_lattice)   :: this

         if(mpi_nproc.gt.1)then

          call MPI_BCAST(this%touched,1,mpi_logical,0,mpi_comm_world,err) 
          call MPI_BCAST(this%a,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%b,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%c,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%alpha,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%beta,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%gamma,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%vol,1,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%cell(1,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%cell(2,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%cell(3,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%nx,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%ny,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%nz,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%ntot,1,mpi_integer,0,mpi_comm_world,err) 
          call MPI_BCAST(this%J(1,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%J(2,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%J(3,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Jinv(1,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Jinv(2,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Jinv(3,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcon(1,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcon(2,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcon(3,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcov(1,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcov(2,:),3,mpi_double_precision,0,mpi_comm_world,err) 
          call MPI_BCAST(this%Gcov(3,:),3,mpi_double_precision,0,mpi_comm_world,err) 

          if(.not.allocated(this%rcell))then
           call this%do_supercell()
          endif

         endif

        return       
        end subroutine lattice_bcast

        function dist_vec_pbc(this,v1,v2) result (dist)
        implicit none
        class(bravais_lattice)  :: this
        double precision        :: v1(3),v2(3),a(3),b(3),c(3),dist(3)
        integer                 :: i
                
         a=matmul(this%Jinv,v1)
         b=matmul(this%Jinv,v2)
         c(1)=a(1)-b(1)
         c(1)=c(1)-nint(c(1)/dble(this%nx))*this%nx
         c(2)=a(2)-b(2)
         c(2)=c(2)-nint(c(2)/dble(this%ny))*this%ny
         c(3)=a(3)-b(3)
         c(3)=c(3)-nint(c(3)/dble(this%nz))*this%nz

         dist=matmul(this%J,c)

        return
        end function dist_vec_pbc

        function dist_pbc(this,v1,v2) result (dist)
        implicit none
        class(bravais_lattice)  :: this
        double precision        :: v1(3),v2(3),a(3),b(3),c(3),dist
        integer                 :: i
                
         a=matmul(this%Jinv,v1)
         b=matmul(this%Jinv,v2)
         c(1)=a(1)-b(1)
         c(1)=c(1)-nint(c(1)/dble(this%nx))*this%nx
         c(2)=a(2)-b(2)
         c(2)=c(2)-nint(c(2)/dble(this%ny))*this%ny
         c(3)=a(3)-b(3)
         c(3)=c(3)-nint(c(3)/dble(this%nz))*this%nz
         dist=this%dist_dir(c)

        return
        end function dist_pbc

        function dist_dir(lat,q)  result(abc)
        implicit none
        class(bravais_lattice)   :: lat
        double precision         :: abc,q(3)
        integer                  :: i,j
          abc=0.0d0
          do i=1,3
           do j=1,3
            abc=abc+lat%Gcov(i,j)*q(i)*q(j)
           enddo
          enddo
          abc=sqrt(abc)
        return 
        end function dist_dir

        function dist_rec(lat,q)  result(abc)
        implicit none
        class(bravais_lattice)  :: lat
        double precision        :: abc,q(3)
        integer                 :: i,j
          abc=0.0d0
          do i=1,3
           do j=1,3
            abc=abc+lat%Gcon(i,j)*q(i)*q(j)
           enddo
          enddo
          abc=sqrt(abc)
        return 
        end function dist_rec

        subroutine init_lattice(this)
        implicit none
        class(bravais_lattice) :: this        
         if(this%touched)then
          call this%calc_Js()
          call this%calc_Gs()
          call this%calc_vol()
         else
          write(*,*) 'Error, no bravais_lattice has been defined'
          stop
         endif
        return
        end subroutine init_lattice

        subroutine calc_Js(this)
        use lapack_inverse
        implicit none
        class(bravais_lattice) :: this
         this%J=transpose(this%cell)
         this%Jinv=this%J
         call mat_inv(this%Jinv,3)                                     
        return
        end subroutine calc_Js

        subroutine calc_Gs(this)
        use lapack_inverse
        implicit none
        class(bravais_lattice) :: this
         this%Gcov=matmul(transpose(this%J),this%J)
         this%Gcon=this%Gcov
         call mat_inv(this%Gcon,3)                                     
        return
        end subroutine calc_Gs

        subroutine calc_vol(this)
        implicit none
        class(bravais_lattice) :: this
         this%vol=1-(cos(this%alpha))**2
         this%vol=this%vol-(cos(this%beta))**2
         this%vol=this%vol-(cos(this%gamma))**2
         this%vol=this%vol+(2*cos(this%alpha)*cos(this%beta)*cos(this%gamma))
         this%vol=sqrt(this%vol)
         this%vol=this%vol*this%a*this%b*this%c
        return
        end subroutine calc_vol

        subroutine do_supercell(this)
        implicit none
        class(bravais_lattice) :: this
        integer        :: v,i,s,k,ialloc

         v=1
         if(.not. allocated(this%rcell))then
          allocate(this%rcell((this%nx*this%ny*this%nz),3),stat=ialloc)
         endif

         do k=0,this%nz-1
          do s=0,this%ny-1
           do i=0,this%nx-1
            this%rcell(v,1)=i-nint(DBLE(i)/DBLE(this%nx))*this%nx
            this%rcell(v,2)=s-nint(DBLE(s)/DBLE(this%ny))*this%ny
            this%rcell(v,3)=k-nint(DBLE(k)/DBLE(this%nz))*this%nz
            v=v+1  
           enddo
          enddo
         enddo

         this%ntot=this%nx*this%ny*this%nz

        return
        end subroutine do_supercell

        subroutine abc2cell(this)
        implicit none
        class(bravais_lattice) :: this              
        return
        end subroutine abc2cell


        subroutine cell2abc(this)
        implicit none
        class(bravais_lattice) :: this
        double precision :: norm1,norm2,norm3,cosab,cosac,cosbc
         norm1=SQRT( this%cell(1,1)**2 + this%cell(1,2)**2 + this%cell(1,3)**2 )
         norm2=SQRT( this%cell(2,1)**2 + this%cell(2,2)**2 + this%cell(2,3)**2 )
         norm3=SQRT( this%cell(3,1)**2 + this%cell(3,2)**2 + this%cell(3,3)**2 )

         this%a=norm1
         this%b=norm2
         this%c=norm3

         cosab=(this%cell(1,1)*this%cell(2,1) + &
                this%cell(1,2)*this%cell(2,2) + &
                this%cell(1,3)*this%cell(2,3))/norm1/norm2
         cosac=(this%cell(1,1)*this%cell(3,1) + &
                this%cell(1,2)*this%cell(3,2) + &
                this%cell(1,3)*this%cell(3,3))/norm1/norm3
         cosbc=(this%cell(3,1)*this%cell(2,1) + &
                this%cell(3,2)*this%cell(2,2) + &
                this%cell(3,3)*this%cell(2,3))/norm3/norm2

         this%alpha=acos(cosbc)
         this%beta=acos(cosac)
         this%gamma=acos(cosab)        
        return
        end subroutine cell2abc


        end module lattice_class
