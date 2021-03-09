        module sparse_class
        use lists_class
        implicit none

        type csr_mat
         integer                   :: nzel
         integer, allocatable      :: AI(:)
         integer, allocatable      :: AJ(:)
         integer, allocatable      :: AC(:)
         contains
         procedure     :: init => init_0
         procedure     :: delete => delete_0 
         procedure     :: check => check_ij
         procedure     :: do_coord
         procedure     :: get_pos
         procedure     :: read_ij
         procedure     :: mult_dense
         generic       :: mult => mult_dense
         generic       :: rd => read_ij
         procedure     :: block => mat2blocks
        end type csr_mat

        type, extends(csr_mat) :: csr_mat_cmplx 
         complex(8), allocatable     :: A(:)
         contains
         procedure     :: read_ij_cmplx
         generic       :: rd => read_ij_cmplx
         procedure     :: init => init_cmplx
         procedure     :: delete => delete_cmplx
         procedure     :: todense => todense_cmplx
         procedure     :: tosparse => tosparse_cmplx
         procedure     :: trans => transpose_cmplx
         procedure     :: mult_dense_cmplx
         generic       :: mult => mult_dense_cmplx
        end type

        type, extends(csr_mat) :: csr_mat_dbl 
         double precision, allocatable     :: A(:)
         contains
         procedure     :: read_ij_dbl
         generic       :: rd => read_ij_dbl
         procedure     :: init => init_dbl
         procedure     :: mult_dense_dbl
         generic       :: mult => mult_dense_dbl
        end type

        type, extends(csr_mat) :: csr_mat_int
         integer, allocatable     :: A(:)
         contains
         procedure     :: read_ij_int
         generic       :: rd => read_ij_int
         procedure     :: init => init_int
         procedure     :: mult_dense_int
         generic       :: mult => mult_dense_int
        end type

        contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function check_ij(A,i,j) result (got_it)
        implicit none
        class(csr_mat)          :: A
        integer                 :: i,j        
        integer                 :: l,n
        logical                 :: got_it
       
        if(i.gt.(size(A%AI)-1) .or. i.le.0) then
         write(*,*) i,j
         stop
        endif

        if(j.gt.(size(A%AI)-1) .or. j.le.0) then
         write(*,*) i,j
         stop
        endif

        got_it=.false.

        do l=1,A%AI(i+1)-A%AI(i)         
         if(A%AJ(A%AI(i)+l).eq.j)then             
            got_it=.true.
            exit
         endif
        enddo         

        return
        end function check_ij


        function get_pos(A,i,j) result (n)
        implicit none
        class(csr_mat)      :: A
        integer                 :: i,j        
        integer                 :: l,n
       
        if(i.gt.(size(A%AI)-1) .or. i.le.0) stop
        if(j.gt.(size(A%AI)-1) .or. j.le.0) stop

        n=0

        do l=1,A%AI(i+1)-A%AI(i)         
         if(A%AJ(A%AI(i)+l).eq.j)then             
            n=A%AI(i)+l
            exit
         endif
        enddo
         
        return
        end function get_pos


        subroutine do_coord(A)
        implicit none
        class(csr_mat)          :: A
        integer                 :: i
        integer                 :: l,v
       
        v=1

        if(.not.allocated(A%AC))then
         allocate(A%AC(size(A%AJ)))

         do i=1,size(A%AI)-1           
          do l=1,A%AI(i+1)-A%AI(i)         
           A%AC(v)=i
           v=v+1
          enddo
         enddo
 
        endif

        return
        end subroutine do_coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine tosparse_cmplx(this,N,A)
        use mpi
        use mpi_utils
        use blacs_utils
        use lists_class
        implicit none
        class(csr_mat_cmplx)    :: this
        type(dist_cmplx_mat)    :: A
        type(list)              :: AI,AJ,Aval
        integer                 :: N,i,j,ii,jj
        complex(8)              :: val

         allocate(this%AI(N+1))
         this%AI(1)=0

         call AJ%init()
         call Aval%init()

         do i=1,N
          this%AI(i+1)=this%AI(i)
          do j=1,N             
           
           call pzelget('A',' ',val,A%mat,i,j,A%desc)

           if(abs(dble(val)).gt.1.0d-10.or.abs(aimag(val)).gt.1.0d-10)then
            this%AI(i+1)=this%AI(i+1)+1  
            call AJ%add_node(j)
            call Aval%add_node(val)
           endif

          enddo
         enddo

         call AJ%reboot()
         call Aval%reboot()

         allocate(this%AJ(AJ%nelem))
         allocate(this%A(AJ%nelem))

         do i=1,AJ%nelem         
          call AJ%rd_val(this%AJ(i)) 
          call Aval%rd_val(this%A(i)) 
          call AJ%skip()
          call Aval%skip()
         enddo

         call AJ%delete()
         call Aval%delete()
   
        return
        end subroutine tosparse_cmplx 

        subroutine todense_cmplx(this,A)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(csr_mat_cmplx)    :: this
        type(dist_cmplx_mat)    :: A
        integer                 :: i,j,ii,jj

         call A%set(size(this%AI)-1,size(this%AI)-1,NB,MB)
         A%mat=(0.0d0,0.0d0)

         do i=1,size(this%AI)-1
          do j=1,this%AI(i+1)-this%AI(i)
           ii=i
           jj=this%AJ(this%AI(i)+j)
           call pzelset(A%mat,ii,jj,A%desc, &
                         this%A(this%AI(i)+j))
          enddo
         enddo       

        return
        end subroutine todense_cmplx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine read_ij(A,i,j)
        implicit none   
        class(csr_mat) :: A
        integer        :: i,j
        return
        end subroutine read_ij

        subroutine read_ij_dbl(A,i,j,val)
        implicit none
        class(csr_mat_dbl)      :: A
        integer                 :: i,j        
        double precision        :: val
        integer                 :: l
       
        if(i.gt.(size(A%AI)-1) .or. i.le.0) stop
        if(j.gt.(size(A%AI)-1) .or. j.le.0) stop

        val=0.0d0

        do l=1,A%AI(i+1)-A%AI(i)         
         if(A%AJ(A%AI(i)+l).eq.j)then
              val=A%A(A%AI(i)+l)
              exit
         endif
        enddo

        return
        end subroutine read_ij_dbl
        
        subroutine read_ij_cmplx(A,i,j,val)
        implicit none
        class(csr_mat_cmplx)    :: A
        integer                 :: i,j        
        complex(8)              :: val
        integer                 :: l
        integer                 :: t1,t2,rate

        if(i.gt.(size(A%AI)-1) .or. i.le.0) then
         write(*,*) 'i',i
         flush(6)
         stop
        endif
        if(j.gt.(size(A%AI)-1) .or. j.le.0) then
         write(*,*) 'j',j
         stop
        endif

        val=(0.0d0,0.0d0)

        do l=1,A%AI(i+1)-A%AI(i)         
         if(A%AJ(A%AI(i)+l).eq.j)then
              val=A%A(A%AI(i)+l)
              exit
         endif
        enddo

        return
        end subroutine read_ij_cmplx

        subroutine read_ij_int(A,i,j,val)
        implicit none
        class(csr_mat_int)      :: A
        integer                 :: i,j        
        integer                 :: val
        integer                 :: l
       
        if(i.gt.(size(A%AI)-1) .or. i.le.0) stop
        if(j.gt.(size(A%AI)-1) .or. j.le.0) stop

        val=0

        do l=1,A%AI(i+1)-A%AI(i)         
         if(A%AJ(A%AI(i)+l).eq.j)then
              val=A%A(A%AI(i)+l)
              exit
         endif
        enddo

        return
        end subroutine read_ij_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine init_0(this,NI,NJ)
        implicit none
        class(csr_mat)  :: this
        integer         :: NI,NJ        
         this%nzel=NI
         allocate(this%AI(NI+1))
         allocate(this%AJ(NJ))
        return
        end subroutine init_0

        subroutine init_int(this,NI,NJ)
        implicit none
        class(csr_mat_int)  :: this
        integer               :: NI,NJ
         allocate(this%A(NJ))       
         allocate(this%AI(NI+1))
         allocate(this%AJ(NJ))
        return
        end subroutine init_int

        subroutine init_dbl(this,NI,NJ)
        implicit none
        class(csr_mat_dbl)    :: this
        integer               :: NI,NJ
         allocate(this%A(NJ))       
         allocate(this%AI(NI+1))
         allocate(this%AJ(NJ))
        return
        end subroutine init_dbl

        subroutine init_cmplx(this,NI,NJ)
        implicit none
        class(csr_mat_cmplx)  :: this
        integer               :: NI,NJ
         allocate(this%A(NJ))       
         allocate(this%AI(NI+1))
         allocate(this%AJ(NJ))
        return
        end subroutine init_cmplx

        subroutine delete_0(this)
        implicit none
        class(csr_mat)  :: this                
         this%nzel=0
         if(allocated(this%AI)) deallocate(this%AI)
         if(allocated(this%AJ)) deallocate(this%AJ)
         if(allocated(this%AC)) deallocate(this%AC)
        return
        end subroutine delete_0

        subroutine delete_cmplx(this)
        implicit none
        class(csr_mat_cmplx)  :: this                
         this%nzel=0
         if(allocated(this%AI)) deallocate(this%AI)
         if(allocated(this%AJ)) deallocate(this%AJ)
         if(allocated(this%A)) deallocate(this%A)
         if(allocated(this%AC)) deallocate(this%AC)
        return
        end subroutine delete_cmplx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine transpose_cmplx(this)
        implicit none
        class(csr_mat_cmplx) :: this
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        integer              :: i,j,k,counter
       
         call mat%init()
         call matI%init()
         call matJ%init()

         call matI%add_node(0)
         counter=0

         do k=1,size(this%AI)-1
          do i=1,size(this%AJ)
           
           if(this%AJ(i).eq.k)then
            
            counter=counter+1
            call mat%add_node(conjg(this%A(i)))

            do j=1,size(this%AI)
             if(this%AI(j).ge.i)then
              call matJ%add_node(j-1)
              exit
             endif
            enddo

           endif

          enddo
          call matI%add_node(counter)
         enddo

         call this%delete()
         call this%init(matI%nelem-1,matJ%nelem)
         call mat%reboot()
         call matI%reboot()
         call matJ%reboot()

         do i=1,matI%nelem
          call matI%rd_val(this%AI(i))
          call matI%skip
         enddo

         do i=1,matJ%nelem
          call matJ%rd_val(this%AJ(i))
          call mat%rd_val(this%A(i))
          call matJ%skip
          call mat%skip
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete
           
        return
        end subroutine transpose_cmplx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine mult_dense(A,C)
        implicit none 
        class(csr_mat) :: A,C         
        return
        end subroutine mult_dense

        subroutine mult_dense_dbl(A,B,C)
        implicit none 
        class(csr_mat_dbl) :: A,C         
        double precision, pointer :: B(:,:)
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        double precision     :: Cij
        integer              :: i,v,j,counter

         call mat%init()
         call matI%init()
         call matJ%init()

         counter=0
         call matI%add_node(counter)

         do i=1,size(A%AI)-1
          do v=1,size(B,1)
           Cij=0
           do j=1,A%AI(i+1)-A%AI(i)
            Cij=Cij+A%A(A%AI(i)+j)*B(j,v)
           enddo
           if(Cij.ne.0)then
            call mat%add_node(Cij)
            call matJ%add_node(v)
            counter=counter+1
           endif
          enddo
          call matI%add_node(counter)
         enddo       

         call C%init(matI%nelem-1,matJ%nelem)
         call mat%reboot
         call matI%reboot
         call matJ%reboot

         i=1
         do while (matI%nelem.gt.0)
          select type ( arrow => matI%node%key )
           type is (integer)
           C%AI(i)=arrow
          end select
          call matI%rm
          i=i+1
         enddo

         i=1
         do while (matJ%nelem.gt.0)
          select type ( arrow => matJ%node%key )
           type is (integer)
           C%AJ(i)=arrow
          end select
          select type ( arrow => mat%node%key )
           type is (double precision)
            C%A(i)=arrow
          end select
          call matJ%rm
          call mat%rm
          i=i+1
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete

        return
        end subroutine mult_dense_dbl

        subroutine mult_dense_cmplx(A,B,C)
        implicit none 
        class(csr_mat_cmplx) :: A,C         
        complex(8), pointer :: B(:,:)
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        complex(8)           :: Cij
        integer              :: i,v,j,counter

         call mat%init()
         call matI%init()
         call matJ%init()

         counter=0
         call matI%add_node(counter)

         do i=1,size(A%AI)-1
          do v=1,size(B,1)
           Cij=0
           do j=1,A%AI(i+1)-A%AI(i)
            Cij=Cij+A%A(A%AI(i)+j)*B(j,v)
           enddo
           if(Cij.ne.0)then
            call mat%add_node(Cij)
            call matJ%add_node(v)
            counter=counter+1
           endif
          enddo
          call matI%add_node(counter)
         enddo       

         call C%init(matI%nelem-1,matJ%nelem)
         call mat%reboot
         call matI%reboot
         call matJ%reboot

         i=1
         do while (matI%nelem.gt.0)
          select type ( arrow => matI%node%key )
           type is (integer)
           C%AI(i)=arrow
          end select
          call matI%rm
          i=i+1
         enddo

         i=1
         do while (matJ%nelem.gt.0)
          select type ( arrow => matJ%node%key )
           type is (integer)
           C%AJ(i)=arrow
          end select
          select type ( arrow => mat%node%key )
           type is (complex(8))
            C%A(i)=arrow
          end select
          call matJ%rm
          call mat%rm
          i=i+1
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete

        return
        end subroutine mult_dense_cmplx

        subroutine mult_dense_int(A,B,C)
        implicit none 
        class(csr_mat_int) :: A,C         
        integer, pointer   :: B(:,:)
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        integer              :: Cij,i,v,j,counter

         call mat%init()
         call matI%init()
         call matJ%init()

         counter=0
         call matI%add_node(counter)

         do i=1,size(A%AI)-1
          do v=1,size(B,1)
           Cij=0
           do j=1,A%AI(i+1)-A%AI(i)
            Cij=Cij+A%A(A%AI(i)+j)*B(j,v)
           enddo
           if(Cij.ne.0)then
            call mat%add_node(Cij)
            call matJ%add_node(v)
            counter=counter+1
           endif
          enddo
          call matI%add_node(counter)
         enddo       

         call C%init(matI%nelem-1,matJ%nelem)
         call mat%reboot
         call matI%reboot
         call matJ%reboot

         i=1
         do while (matI%nelem.gt.0)
          select type ( arrow => matI%node%key )
           type is (integer)
           C%AI(i)=arrow
          end select
          call matI%rm
          i=i+1
         enddo

         i=1
         do while (matJ%nelem.gt.0)
          select type ( arrow => matJ%node%key )
           type is (integer)
           C%AJ(i)=arrow
          end select
          select type ( arrow => mat%node%key )
           type is (integer)
            C%A(i)=arrow
          end select
          call matJ%rm
          call mat%rm
          i=i+1
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete

        return
        end subroutine mult_dense_int




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  MAT2BLOCKS
!!!
!!!  Given a sparse matrix A it computes all the disconnected graph of
!!!  its elements and save the permutation of basis vector in mapp(:) and
!!!  the size of blocks in blc(:)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine mat2blocks(A,blc,mapp)
        implicit none
        class(csr_mat)         :: A
        integer, allocatable   :: mapp(:),blc(:)
        type(list)             :: r,r2,queue
        integer                :: dim_block,N,pos
        integer                :: i,j,l,k,m,x,v
        logical, allocatable   :: check(:)
        class(*), pointer      :: arrow
        
          
        call r%init()      ! elementi del blocco
        call r2%init()     ! dimensione blocco
        call queue%init()  

                   
        N=size(A%AI)-1
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
 
        do i=1,A%AI(v+1)-A%AI(v)

!!       !!!! ATT MODIFICA FATTA 6/6/2019
!         pos=A%AJ(v)
         pos=A%AJ(i+A%AI(v))
                
         if(check(pos))then

         else
          check(pos)=.true.
          call queue%add_node(pos)
          call r%add_node(pos)
          dim_block=dim_block+1          
         endif

        enddo

        call queue%rm
        enddo

        k=k+1
        
        call r2%add_node(dim_block)

        endif

        enddo

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
       

        if( allocated(A%AC)) deallocate(A%AC)
        call r%delete()
        call r2%delete()
        call queue%delete()

        return
        end subroutine mat2blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  MULT_SPARSE MULT_SPARSE_CMPLX MULT_SPARSE_DBL MULT_SPARSE_INT
!!!
!!!  Specialised routines for sparse-sparse matrices multiplication
!!!  mult_sparse redirect to the type specific mult_sparse_XXX routines
!!! !according to the type of A and B, which must be the same
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine mult_sparse_cmplx(cross,A,B,C)
        implicit none
        type(csr_mat_cmplx), intent(in)   :: A
        type(csr_mat_cmplx), intent(in)   :: B
        type(csr_mat_cmplx), intent(out)  :: C   
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        integer              :: i,v,j,counter
        integer              :: t1,t2,rate
        complex(8)           :: val,Cij
        logical              :: cross

         call system_clock(t1,rate)        

         if(size(A%AI).ne.size(B%AI))then
          write(*,*) size(A%AI),size(B%AI)
          flush(6)
          stop
         endif

         if (cross) call A%trans()

         call mat%init()
         call matI%init()
         call matJ%init()

         call matI%add_node(0)
         counter=0

         do i=1,size(A%AI)-1
          do v=1,size(B%AI)-1
           Cij=(0.0d0,0.0d0)
           do j=1,A%AI(i+1)-A%AI(i)
            call B%rd(A%AJ(A%AI(i)+j),v,val)
            Cij=Cij+val*A%A(A%AI(i)+j)
            val=(0.0d0,0.0d0)
           enddo
           if(abs(dble(Cij)).gt.1.0e-9.or.abs(aimag(Cij)).gt.1.0e-9)then
            call mat%add_node(Cij)
            call matJ%add_node(v)
            counter=counter+1
           endif
          enddo
          call matI%add_node(counter)
         enddo       

         call C%init(matI%nelem-1,matJ%nelem)
         call mat%reboot()
         call matI%reboot()
         call matJ%reboot()

         do i=1,matI%nelem
          call matI%rd_val(C%AI(i))
          call matI%skip
         enddo

         do i=1,matJ%nelem
          call matJ%rd_val(C%AJ(i))
          call mat%rd_val(C%A(i))
          call matJ%skip
          call mat%skip
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete

         call system_clock(t2)
         write(*,*) 'Sparse mult done in ',real(t2-t1)/real(rate),'s'
         flush(6)
        
        return
        end subroutine mult_sparse_cmplx

        subroutine mult2_sparse_cmplx(cross,A2,B2,C)
        use mpi
        use mpi_utils
        implicit none
        type(csr_mat_cmplx), intent(in)   :: A2
        type(csr_mat_cmplx), intent(in)   :: B2
        type(csr_mat_cmplx), intent(out)  :: C   
        type(csr_mat_cmplx)               :: A,B
        type(list)           :: mat
        type(list)           :: matI
        type(list)           :: matJ
        integer              :: i,v,j,counter,k,l,val2
        integer              :: t1,t2,rate
        complex(8)           :: val,Cij
        logical              :: cross,check
        integer              :: nloc,nstart,nnz
        integer,allocatable  :: proc_grid(:)

         if(mpi_blacs_id.eq.0)then
          call system_clock(t1,rate)        
         endif

         if(size(A2%AI).ne.size(B2%AI))then
          write(*,*) size(A2%AI),size(B2%AI)
          flush(6)
          stop
         endif

         call mpi_dist_nprocess(size(A2%AI)-1,nloc,nstart,proc_grid,mpi_blacs_world)

         call A%init(size(A2%AI)-1,size(A2%AJ))
         call B%init(size(B2%AI)-1,size(B2%AJ))

         A%AI=A2%AI
         A%AJ=A2%AJ
         A%A=A2%A

         B%AI=B2%AI
         B%AJ=B2%AJ
         B%A=B2%A

         if (cross) call A%trans()
         call B%trans()

         call mat%init()
         call matI%init()
         call matJ%init()

         counter=0

         i=nstart
         do l=1,nloc
          do v=1,size(B%AI)-1
           j=1
           k=1
           Cij=(0.0d0,0.0d0)
           check=.false.
           do while (j.le.(A%AI(i+1)-A%AI(i)) .and. &
                     k.le.(B%AI(v+1)-B%AI(v)) )

            if(A%AJ(A%AI(i)+j).lt.B%AJ(B%AI(v)+k))then
             j=j+1
             else if(A%AJ(A%AI(i)+j).gt.B%AJ(B%AI(v)+k))then
             k=k+1
            else 
             Cij=Cij+A%A(A%AI(i)+j)*conjg(B%A(B%AI(v)+k))
             j=j+1
             k=k+1
             check=.true.
            endif

           enddo
           if( abs(dble(Cij)).gt.1.0e-9.or.abs(aimag(Cij)).gt.1.0e-9 )then
            call mat%add_node(Cij)
            call matJ%add_node(v)
            counter=counter+1
           endif
          enddo
          call matI%add_node(counter)
          i=i+1
         enddo       
         
         ! allgather pieces of C%AI(i)

         call mpi_allreduce(matJ%nelem,nnz,1, &
                  mpi_integer,MPI_SUM,mpi_blacs_world,err)    

         call C%init(size(A2%AI)-1,nnz)
         call mat%reboot()
         call matI%reboot()
         call matJ%reboot()

         C%AI(1)=0
         counter=0
         do i=1,size(C%AI)-1
          if(mpi_blacs_id.eq.proc_grid(i))then
           call matI%rd_val(val2)
           call matI%skip
           C%AI(i+1)=val2+counter
          endif
          call mpi_bcast(C%AI(i+1),1,mpi_integer,proc_grid(i),mpi_blacs_world,err)
          if(i.lt.size(C%AI)-1)then
           if(proc_grid(i+1).ne.proc_grid(i)) counter=C%AI(i+1)
          endif
          call mpi_bcast(counter,1,mpi_integer,proc_grid(i),mpi_blacs_world,err)
         enddo

         counter=1
         do j=0,mpi_blacs_nproc-1
          if(mpi_blacs_id.eq.j) nnz=mat%nelem
          call mpi_bcast(nnz,1,mpi_integer,j,mpi_blacs_world,err)
          do i=1,nnz           
           if(mpi_blacs_id.eq.j)then
            call matJ%rd_val(val2)
            call matJ%skip
            call mat%rd_val(val)
            call mat%skip
            C%AJ(counter)=val2
            C%A(counter)=val
           endif
           call mpi_bcast(C%AJ(counter),1,mpi_integer,j,mpi_blacs_world,err)
           call mpi_bcast(C%A(counter),1,mpi_double_complex,j,mpi_blacs_world,err)
           counter=counter+1
          enddo
         enddo

         call mat%delete
         call matI%delete
         call matJ%delete

         if(mpi_blacs_id.eq.0)then
          call system_clock(t2)
          write(*,*) 'Sparse mult done in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif
        
        return
        end subroutine mult2_sparse_cmplx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        
        end module sparse_class


!!!!    add block decomposition
!!!!    add diagonalization matrices

