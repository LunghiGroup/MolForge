        module blacs_utils
        use mpi
        use mpi_utils
        implicit none

         integer              :: nprow,npcol
         integer              :: myrow,mycol
         integer              :: mpi_color
         integer              :: NB=6,MB=6
         integer              :: lwork,lrwork,icntx
         integer,allocatable  :: context(:)
         double precision, allocatable :: work(:),rwork(:)
         complex(8), allocatable       :: cwork(:),crwork(:)

         type :: dist_mat
          integer                 :: desc(9)
          integer                 :: N,M
          logical                 :: alloc=.false.
         end type dist_mat

         type,extends(dist_mat) :: dist_dbl_mat
          double precision, allocatable :: mat(:,:)
          contains
          procedure               :: set =>  setup_blacs_dbl_matrix
          procedure               :: dealloc =>  dealloc_dbl
          procedure               :: get_nze => get_nze_dbl
         end type  dist_dbl_mat

         type ,extends(dist_mat):: dist_cmplx_mat
          complex(8), allocatable :: mat(:,:)
          complex(8)              :: csr_mat_cmplx 
          contains
          procedure               :: set =>  setup_blacs_cmplx_matrix
          procedure               :: dealloc =>  dealloc_cmplx
          procedure               :: get_nze => get_nze_cmplx
          procedure               :: raze => raze_cmplx
          procedure               :: print_mat => print_mat_cmplx
          procedure               :: set_val => set_val_cmplx
          procedure               :: read_mat => read_mat_cmplx
         end type  dist_cmplx_mat

         type ,extends(dist_mat):: dist_cmplx_vec
          complex(8), allocatable :: vec(:)
          contains
          procedure               :: set =>  setup_blacs_cmplx_vector
          procedure               :: dealloc =>  dealloc_cmplx_vector
         end type  dist_cmplx_vec


         type ,extends(dist_mat):: dist_int_mat
          integer, allocatable :: mat(:,:)
          contains
          procedure               :: set =>  setup_blacs_int_matrix
          procedure               :: dealloc =>  dealloc_int
         end type  dist_int_mat

         type ,extends(dist_mat):: dist_logic_mat
          logical, allocatable    :: mat(:,:)
          contains
          procedure               :: set =>  setup_blacs_logical_matrix
          procedure               :: dealloc =>  dealloc_logical
         end type  dist_logic_mat


        contains

         subroutine set_val_cmplx(this,i,j,val)
         use mpi
         use mpi_utils
         implicit none
         class(dist_cmplx_mat) :: this
         integer               :: i,j
         double complex        :: val

          call pzelset(this%mat,i,j,this%desc,val)

         return
         end subroutine set_val_cmplx

         subroutine read_mat_cmplx(this,filename)
         use mpi
         use mpi_utils
         implicit none
         class(dist_cmplx_mat)   :: this
         integer                 :: l,i,j,Nelem
         double complex          :: val
         double precision        :: valr,valc
         character(len=30)       :: filename

          if(mpi_id.eq.0) open(15,file=trim(filename))
          if(mpi_id.eq.0) read(15,*) Nelem
          
          this%mat=(0.0d0,0.0d0)

          do l=1,Nelem
           if(mpi_id.eq.0) read(15,*) i,j,valr,valc
           if(mpi_id.eq.0) val=cmplx(valr,valc,8)          
           call mpi_bcast(val,1,mpi_double_complex,0,mpi_blacs_world,err)
           call mpi_bcast(i,1,mpi_integer,0,mpi_blacs_world,err)
           call mpi_bcast(j,1,mpi_integer,0,mpi_blacs_world,err)
           call pzelset(this%mat,i,j,this%desc,val)
          enddo

          if(mpi_id.eq.0) close(15)

         return
         end subroutine read_mat_cmplx

         subroutine print_mat_cmplx(this,filename)
         use mpi
         use mpi_utils
         implicit none
         class(dist_cmplx_mat)   :: this
         integer                 :: l,j
         double complex          :: val
         character(len=30)       :: filename

          if(mpi_id.eq.0) open(15,file=trim(filename))
          if(mpi_id.eq.0) write(15,*) this%N*this%M

          do l=1,this%N
           do j=1,this%M
            call pzelget('A',' ',val,this%mat,l,j,this%desc)
            if(mpi_id.eq.0) write(15,*) l,j,dble(val),aimag(val)
           enddo
          enddo

          if(mpi_id.eq.0) flush(15)
          if(mpi_id.eq.0) close(15)

         return
         end subroutine print_mat_cmplx

         function get_nze_dbl(this,thr) result(nze)
         use mpi
         use mpi_utils
         implicit none
         class(dist_dbl_mat) :: this
         integer             :: nze,i,j
         double precision    :: thr

          nze=0

          do i=1,size(this%mat,1)
           do j=1,size(this%mat,2)
            if(abs(this%mat(i,j)).gt.thr) nze=nze+1
           enddo
          enddo

          call mpi_allreduce(nze,nze,1,mpi_integer,mpi_sum,mpi_blacs_world,err)

         return
         end function get_nze_dbl

         subroutine raze_cmplx(this,thr)
         use mpi
         use mpi_utils
         implicit none        
         class(dist_cmplx_mat) :: this
         integer               :: i,j
         double precision      :: thr

          do i=1,size(this%mat,1)
           do j=1,size(this%mat,2)
            if( abs(dble(this%mat(i,j))).lt.thr .and. & 
               abs(aimag(this%mat(i,j))).lt.thr  ) this%mat(i,j)=(0.0d0,0.0d0)
           enddo
          enddo

         return
         end subroutine raze_cmplx

         function get_nze_cmplx(this,thr) result(nze)
         use mpi
         use mpi_utils
         implicit none        
         class(dist_cmplx_mat) :: this
         integer               :: nze,i,j
         double precision      :: thr

          nze=0

          do i=1,size(this%mat,1)
           do j=1,size(this%mat,2)
            if( abs(dble(this%mat(i,j))).gt.thr .or. & 
               abs(aimag(this%mat(i,j))).gt.thr  ) nze=nze+1
           enddo
          enddo

          call mpi_allreduce(nze,nze,1,mpi_integer,mpi_sum,mpi_blacs_world,err)

         return
         end function get_nze_cmplx

         subroutine dealloc_cmplx_vector(this)
         implicit none
         class(dist_cmplx_vec) :: this
          deallocate(this%vec)
         return
         end subroutine dealloc_cmplx_vector

         subroutine dealloc_dbl(this)
         implicit none
         class(dist_dbl_mat) :: this
          deallocate(this%mat)
         return
         end subroutine dealloc_dbl

         subroutine dealloc_cmplx(this)
         implicit none
         class(dist_cmplx_mat) :: this
          deallocate(this%mat)
         return
         end subroutine dealloc_cmplx

         subroutine dealloc_int(this)
         implicit none
         class(dist_int_mat) :: this
          deallocate(this%mat)
         return
         end subroutine dealloc_int

         subroutine dealloc_logical(this)
         implicit none
         class(dist_logic_mat) :: this
          deallocate(this%mat)
         return
         end subroutine dealloc_logical

         subroutine setup_blacs(mpi_nproc_loc)
         use mpi
         use mpi_utils
         implicit none
         integer                :: mpi_nproc_loc

          nprow=int(sqrt(dble(mpi_nproc_loc)))
          npcol=mpi_nproc_loc/nprow

          if(.not.allocated(context)) allocate(context(1))
          call blacs_get( -1, 0, context(1) )
          call blacs_gridinit(context(1),'R',nprow,npcol)
          call blacs_gridinfo(context(1),nprow,npcol,myrow,mycol)

         return
         end subroutine setup_blacs

         subroutine setup_multiblacs(mpi_nproc_loc,map,cntx)
         use mpi
         use mpi_utils
         implicit none
         integer                :: mpi_nproc_loc,cntx
         integer, allocatable   :: map(:,:)

          nprow=int(sqrt(dble(mpi_nproc_loc)))
          npcol=mpi_nproc_loc/nprow

          call blacs_get( -1, 0, cntx )
          call blacs_gridmap(cntx,map,nprow,nprow,npcol)
          call blacs_gridinfo(cntx,nprow,npcol,myrow,mycol)

         return
         end subroutine setup_multiblacs

         subroutine dismiss_blacs()
         implicit none

          if(myrow.ne.-1) call blacs_gridexit(context(1))
          call blacs_exit(-1)

         return
         end subroutine dismiss_blacs

         subroutine dismiss_multiblacs(cntx)
         implicit none
         integer :: cntx

          if(myrow.ne.-1) call blacs_gridexit(cntx)
          call blacs_exit(-1)

         return
         end subroutine dismiss_multiblacs

         subroutine blacs_set_gridinfo()
         implicit none
         integer                :: i

          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow,npcol,myrow,mycol)
           if(myrow.ne.-1) exit
          enddo

         return
         end subroutine blacs_set_gridinfo

        !
        ! matrices set-ups
        !

         subroutine setup_blacs_int_matrix(this,N,M,NBl,MBl)
         implicit none
         class(dist_int_mat)            :: this
         integer                        :: N,M,NBl,MBl,numroc,info
         integer                        :: Nloc_row,Nloc_col,i
         integer                        :: myrow_loc,mycol_loc
         integer                        :: nprow_loc,npcol_loc

          this%N=N
          this%M=M
          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow_loc,npcol_loc,myrow_loc,mycol_loc)
           if(myrow_loc.ne.-1)then
            Nloc_row = NUMROC(N,NBl,myrow_loc,0,nprow_loc)
            Nloc_col = NUMROC(M,MBl,mycol_loc,0,npcol_loc)
            allocate(this%mat(Nloc_row,Nloc_col))
            call descinit(this%desc,N,M,NBl,MBl,0,0,context(i),Nloc_row,info)
           endif
          enddo

          this%mat=0
          this%alloc=.true.

         return
         end subroutine setup_blacs_int_matrix

         subroutine setup_blacs_logical_matrix(this,N,M,NBl,MBl)
         implicit none
         class(dist_logic_mat)           :: this
         integer                        :: N,M,NBl,MBl,numroc,info
         integer                        :: Nloc_row,Nloc_col,i
         integer                        :: myrow_loc,mycol_loc
         integer                        :: nprow_loc,npcol_loc
         
          this%N=N
          this%M=M
          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow_loc,npcol_loc,myrow_loc,mycol_loc)
           if(myrow_loc.ne.-1)then
            Nloc_row = NUMROC(N,NBl,myrow_loc,0,nprow_loc)
            Nloc_col = NUMROC(M,MBl,mycol_loc,0,npcol_loc)
            allocate(this%mat(Nloc_row,Nloc_col))
            call descinit(this%desc,N,M,NBl,MBl,0,0,context(i),Nloc_row,info)
           endif
          enddo

          this%mat=.false.
          this%alloc=.true.

         return
         end subroutine setup_blacs_logical_matrix

         subroutine setup_blacs_cmplx_matrix(this,N,M,NBl,MBl)
         implicit none
         class(dist_cmplx_mat)          :: this
         integer                        :: N,M,NBl,MBl,numroc,info
         integer                        :: Nloc_row,Nloc_col,i
         integer                        :: myrow_loc,mycol_loc
         integer                        :: nprow_loc,npcol_loc

          this%N=N
          this%M=M
          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow_loc,npcol_loc,myrow_loc,mycol_loc)
           if(myrow_loc.ne.-1)then
            Nloc_row = NUMROC(N,NBl,myrow_loc,0,nprow_loc)
            Nloc_col = NUMROC(M,MBl,mycol_loc,0,npcol_loc)
            allocate(this%mat(Nloc_row,Nloc_col))
            call descinit(this%desc,N,M,NBl,MBl,0,0,context(i),Nloc_row,info)
           endif
          enddo

          this%mat=(0.0d0,0.0d0)
          this%alloc=.true.

         return
         end subroutine setup_blacs_cmplx_matrix

         subroutine setup_blacs_cmplx_vector(this,N,NBl)
         implicit none
         class(dist_cmplx_vec)          :: this
         integer                        :: N,NBl,numroc,info
         integer                        :: Nloc_row,i
         integer                        :: myrow_loc,mycol_loc
         integer                        :: nprow_loc,npcol_loc

          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow_loc,npcol_loc,myrow_loc,mycol_loc)
           if(myrow_loc.ne.-1)then
            Nloc_row = NUMROC(N,NBl,myrow_loc,0,nprow_loc)
            allocate(this%vec(Nloc_row))
            call descinit(this%desc,N,1,NBl,1,0,0,context(i),Nloc_row,info)
           endif
          enddo

         return
         end subroutine setup_blacs_cmplx_vector

         subroutine setup_blacs_dbl_matrix(this,N,M,NBl,MBl)
         implicit none
         class(dist_dbl_mat)            :: this
         integer                        :: N,M,NBl,MBl,numroc,info
         integer                        :: Nloc_row,Nloc_col,i
         integer                        :: myrow_loc,mycol_loc
         integer                        :: nprow_loc,npcol_loc

          this%N=N
          this%M=M
          do i=1,size(context)
           call blacs_gridinfo(context(i),nprow_loc,npcol_loc,myrow_loc,mycol_loc)
           if(myrow_loc.ne.-1)then
            Nloc_row = NUMROC(N,NBl,myrow_loc,0,nprow_loc)
            Nloc_col = NUMROC(M,MBl,mycol_loc,0,npcol_loc)
            allocate(this%mat(Nloc_row,Nloc_col))
            call descinit(this%desc,N,M,NBl,MBl,0,0,context(i),Nloc_row,info)
           endif
          enddo

          this%mat=0.0d0
          this%alloc=.true.

         return
         end subroutine setup_blacs_dbl_matrix

        end module blacs_utils
