        module descriptors_class                
        use bispectrum_class                 
        implicit none

         type descriptor 
          double precision, allocatable :: desc(:)
          double precision              :: rot(3,3)
          integer                       :: size_desc
         end type descriptor

         type, extends(descriptor) :: coul_mat
          contains 
          procedure                     :: get_desc => get_Cmat
         end type coul_mat

         type, extends(descriptor) :: cartesian
          contains 
          procedure                     :: get_desc => get_cart
         end type cartesian

         type, extends(descriptor) :: bispectrum
          type(bispectrum_components)           :: BI
          contains 
          procedure                             :: get_desc => get_bispectrum
         end type bispectrum

        contains
         
         subroutine get_bispectrum(this,geo,r0)
         use bispectrum_class                 
         implicit none
         class(bispectrum)              :: this
         double precision, allocatable  :: geo(:,:)
         double precision               :: vec(3),vec2(3),weight,r0
         integer                        :: i,j,v,l,Jmax
         integer                        :: t1,t2
         double precision               :: rate
 
          weight=1.0d0

          do l=1,this%BI%max_order+1
           do i=1,size(geo,1)
            vec=geo(i,:)
            call this%BI%SH(l)%get_val(vec,r0,weight)
           enddo
          enddo

          call this%BI%get_bi()

          allocate(this%desc(size(this%BI%B)))
          this%desc=dble(this%BI%B)

         return
         end subroutine get_bispectrum

         subroutine get_cart(this,geo)
         implicit none
         class(cartesian)              :: this
         double precision, allocatable :: geo(:,:)
         integer                       :: i,j,v
                                       
          this%size_desc=size(geo,1)*3

          allocate(this%desc(this%size_desc))
          this%desc=0.0d0

          v=1
          do i=1,size(geo,1)
           do j=1,3
            this%desc(v)=geo(i,j)
            v=v+1
           enddo
          enddo

         return
         end subroutine get_cart

         subroutine get_Cmat(this,geo,Nat)
         use lapack_diag_simm
         implicit none
         class(coul_mat)               :: this
         double precision, allocatable :: geo(:,:),Cmat(:,:)
         integer, allocatable          :: Nat(:)
         double precision              :: r2
         integer                       :: i,j,v
                                       
          this%size_desc=size(geo,1)

          allocate(Cmat(size(geo,1),size(geo,1)))
          allocate(this%desc(this%size_desc))

          Cmat=0.0d0
          this%desc=0.0d0

          do i=1,size(geo,1)
           do j=1,size(geo,1)

            if(i.ne.j)then
             r2=(geo(i,1)-geo(j,1))**2
             r2=r2+(geo(i,2)-geo(j,2))**2
             r2=r2+(geo(i,3)-geo(j,3))**2
             Cmat(i,j)=nat(i)*nat(j)/sqrt(r2)
            else
             Cmat(i,j)=0.5d0*nat(i)**(2.4d0)
            endif

           enddo
          enddo

          call new_diag(this%size_desc,Cmat,this%desc)

          deallocate(Cmat)

         return
         end subroutine get_Cmat


        end module descriptors_class

