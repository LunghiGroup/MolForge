        module spinham_class         
        use lists_class
        implicit none

        double precision                         :: ZERO=1.0E-18
        double complex                           :: ZEROC=(1.0E-8,1.0E-8)

        type, extends(list) :: Jiso_list
        contains
        procedure                 :: add_node => add_Jnode        
        procedure                 :: rd_node => rd_Jnode
        end type Jiso_list 

        type, extends(list) :: DSItensor_list
        contains
        procedure                 :: add_node => add_DSInode        
        procedure                 :: rd_node => rd_DSInode
        end type DSItensor_list 

        type, extends(list) :: D2Stensor_list
        contains
        procedure                 :: add_node => add_D2Snode        
        procedure                 :: rd_node => rd_D2Snode
        end type D2Stensor_list 

        type, extends(list) :: Gtensor_list
        contains
        procedure                 :: add_node => add_Gnode        
        procedure                 :: rd_node => rd_Gnode
        end type Gtensor_list 
        
        type, extends(list) :: OSItensor_list
        contains
        procedure                 :: add_node => add_Onode        
        procedure                 :: rd_node => rd_Onode
        end type OSItensor_list 

        type :: cart_tensor
         complex(8)               :: D(3,3)
         complex(8), allocatable  :: mat(:,:)
         double precision         :: cutoff
         contains
         procedure            :: dump => printD
         procedure            :: traceless 
         procedure            :: get_norm => get_D_norm
         procedure            :: cart2stev
         procedure            :: stev2cart
         procedure            :: rot => rot_cart_tensor
        end type cart_tensor

        type, extends(cart_tensor) :: DSItensor
         integer             :: kind
         contains
         procedure           :: mat_elem => DSI_matelem
        end type DSItensor

        type, extends(cart_tensor) :: D2Stensor
         integer             :: kind(2)
         contains
         procedure           :: mat_elem => D2S_matelem
        end type D2Stensor

        type, extends( D2Stensor) :: Ddipolar
         contains        
         procedure           :: make_D   => make_D_dipolar
        end type Ddipolar
        
        type :: Hmat
         double complex, allocatable  :: H(:,:)
         integer                      :: Hdim
        end type Hmat

        type :: OSItensor
         integer                  :: k
         integer,allocatable      :: q(:)
         complex(8), allocatable  :: B(:)
         complex(8), allocatable  :: mat(:,:)
         integer                  :: kind
         contains
         procedure            :: mat_elem => OSI_matelem
         procedure            :: rot      => rot_O
         procedure            :: get_norm => get_O_norm
         procedure            :: stev2spher
         procedure            :: spher2stev
        end type OSItensor

        type :: Jiso
         complex(8)               :: J
         integer                  :: kind(2)
         complex(8), allocatable  :: mat(:,:)
         contains
         procedure            :: mat_elem => Jiso_matelem
         procedure            :: get_norm => get_J_norm
        end type Jiso

        type :: Gtensor
         complex(8)               :: G(3,3)
         complex(8), allocatable  :: mat(:,:)
         integer                  :: kind
         contains
         procedure           :: mat_elem => G_matelem
         procedure           :: get_norm => get_G_norm
         procedure           :: rotg     => rot_g_tensor
        end type Gtensor

        type :: SpinHamiltonian
         integer                      :: nH=0
         type(Hmat),allocatable       :: H(:)
         integer                      :: nO=0
         type(OSItensor),allocatable  :: O(:)
         integer                      :: nJ=0
         type(Jiso),allocatable       :: J(:)
         integer                      :: nG=0
         type(Gtensor),allocatable    :: G(:)
         integer                      :: nDSI=0
         type(DSItensor),allocatable  :: DSI(:)
         integer                      :: nD2S=0
         type(D2Stensor),allocatable  :: D2S(:)
         integer                      :: nDdip=0
         type(Ddipolar),allocatable   :: Ddip(:)
         double precision             :: dipolar_thr=10.0D6
         logical                      :: make_dipolar=.false.
         contains
         procedure                    :: spinham_bcast
         procedure                    :: rot => rotate_sh
        end type SpinHamiltonian

        contains

        subroutine rotate_sh(this,euler)
        implicit none
        class(SpinHamiltonian)  :: this
        double precision        :: euler(3)
        integer                 :: v

         do v=1,this%nO
          call this%O(v)%rot(euler(1),euler(2),euler(3))
         enddo        

         do v=1,this%nG
          call this%G(v)%rotg(euler(1),euler(2),euler(3))
         enddo        

         do v=1,this%nDSI
          call this%DSI(v)%rot(euler(1),euler(2),euler(3))
         enddo        

         do v=1,this%nD2S
          call this%D2S(v)%rot(euler(1),euler(2),euler(3))
         enddo        

        return
        end subroutine rotate_sh

        subroutine spinham_bcast(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(SpinHamiltonian)  :: this
        integer                 :: i,l,v

         call mpi_bcast(this%nH,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nO,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nJ,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nG,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nDSI,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nD2S,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%make_dipolar,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(this%dipolar_thr,1,mpi_double_precision,0,mpi_comm_world,err)
         
         if(this%nH.gt.0)then
          if(.not.allocated(this%H)) allocate(this%H(this%nH))
          do i=1,this%nH
           call mpi_bcast(this%H(i)%Hdim,1,mpi_integer,0,mpi_comm_world,err)
           do l=1,this%H(i)%Hdim
            do v=1,this%H(i)%Hdim
             call mpi_bcast(this%H(i)%H(v,l),1,mpi_complex,0,mpi_comm_world,err)
            enddo
           enddo
          enddo
         endif

         if(this%nJ.gt.0)then
          if(.not.allocated(this%J)) allocate(this%J(this%nJ))
          do i=1,this%nJ
           call mpi_bcast(this%J(i)%J,1,mpi_complex,0,mpi_comm_world,err)
           call mpi_bcast(this%J(i)%kind,2,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif

         if(this%nG.gt.0)then
          if(.not.allocated(this%G)) allocate(this%G(this%nG))
          do i=1,this%nG
           do l=1,3
            do v=1,3
             call mpi_bcast(this%G(i)%G(v,l),1,mpi_complex,0,mpi_comm_world,err)
            enddo
           enddo
           call mpi_bcast(this%G(i)%kind,1,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif

         if(this%nDSI.gt.0)then
          if(.not.allocated(this%DSI)) allocate(this%DSI(this%nDSI))
          do i=1,this%nDSI
           do l=1,3
            do v=1,3
             call mpi_bcast(this%DSI(i)%D(l,v),1,mpi_complex,0,mpi_comm_world,err)
            enddo
           enddo
           call mpi_bcast(this%DSI(i)%kind,1,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif

         if(this%nD2S.gt.0)then
          if(.not.allocated(this%D2S)) allocate(this%D2S(this%nD2S))
          do i=1,this%nD2S
           do l=1,3
            do v=1,3            
             call mpi_bcast(this%D2S(i)%D(v,l),1,mpi_complex,0,mpi_comm_world,err)
            enddo
           enddo
           call mpi_bcast(this%D2S(i)%kind,2,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif

         if(this%nO.gt.0)then
          if(.not.allocated(this%O)) allocate(this%O(this%nO))
          do i=1,this%nO           
           call mpi_bcast(this%O(i)%k,1,mpi_integer,0,mpi_comm_world,err)
           if(.not. allocated(this%O(i)%B) ) allocate(this%O(i)%B(2*this%O(i)%k+1))
           if(.not. allocated(this%O(i)%q) ) allocate(this%O(i)%q(2*this%O(i)%k+1))
           call mpi_bcast(this%O(i)%B,2*this%O(i)%k+1,mpi_complex,0,mpi_comm_world,err)
           call mpi_bcast(this%O(i)%q,2*this%O(i)%k+1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%O(i)%kind,1,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif


        return
        end subroutine spinham_bcast


        subroutine make_D_dipolar (this,G1,G2,bohr_mag1,bohr_mag2,x1,dist)
        use stevens_class
        use units_parms 
        implicit none
        class(Ddipolar)                          :: this
        complex(8)                               :: G1(3,3),G2(3,3)
        double precision                         :: x1(3),x2(3),dist
        double precision                         :: bohr_mag1,bohr_mag2
        integer                                  :: s,t,k,l

         this%D=(0.0d0,0.0d0)
         
         do s=1,3
          do t=1,3
           do k=1,3
            this%D(s,t)=this%D(s,t)+G1(s,k)*G2(k,t)
           enddo
          enddo
         enddo

         do s=1,3
          do t=1,3
           do k=1,3
            do l=1,3
             this%D(s,t)=this%D(s,t)-3.0d0*G1(s,k)*x1(k)*x1(l)*G2(l,t)/dist**2
            enddo
           enddo
          enddo
         enddo

         this%D=(this%D*mag_vacuum*bohr_mag1*bohr_mag2)/dist**3
         call this%traceless()
        
        return
        end subroutine make_D_dipolar

        subroutine get_D_norm(this,val)
        implicit none
        class(cart_tensor)                       :: this
        double precision                         :: val
        integer                                  :: i,j

         val=0.0d0
         do i=1,3
          do j=1,3
           val=val+this%D(i,j)*conjg(this%D(i,j))
          enddo
         enddo

        return
        end subroutine get_D_norm

        subroutine get_G_norm(this,val)
        implicit none
        class(Gtensor)                           :: this
        double precision                         :: val
        integer                                  :: i,j

         val=0.0d0
         do i=1,3
          do j=1,3
           val=val+this%G(i,j)*conjg(this%G(i,j))
          enddo
         enddo

        return
        end subroutine get_G_norm

        subroutine get_O_norm(this,val)
        implicit none
        class(OSItensor)                         :: this
        double precision                         :: val
        integer                                  :: i

         val=0.0d0
         do i=1,size(this%B)
          val=val+this%B(i)*conjg(this%B(i))
         enddo

        return
        end subroutine get_O_norm

        subroutine get_J_norm(this,val)
        implicit none
        class(Jiso)                              :: this
        double precision                         :: val

         val=0.0d0
         val=val+this%J*conjg(this%J)

        return
        end subroutine get_J_norm

        function G_matelem(this,mapp,spin,B,bohr_mag) result(val)
        implicit none
        class(Gtensor)                           :: this
        complex(8)                               :: val,C
        double precision                         :: mapp(2)
        double precision                         :: spin,B(3),bohr_mag
        integer                                  :: i

         val=(0.0d0,0.0d0)

       !S-
         if(abs(mapp(1)-mapp(2)+1).lt.1.0E-06)then
          C=(0.0d0,0.0d0)
          C=0.5d0*this%G(1,1)*B(1)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*this%G(2,1)*B(2)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*this%G(3,1)*B(3)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          val=val+C!CMPLX(C,0,8)
          C=(0.0d0,0.0d0)
          C=0.5d0*this%G(1,2)*B(1)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*this%G(2,2)*B(2)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          C=C+0.5d0*this%G(3,2)*B(3)*dsqrt((spin+mapp(2))*(spin-mapp(2)+1))
          val=val+C*CMPLX(0.0d0,1.0d0,8)
         endif
        !S+
         if(abs(mapp(1)-mapp(2)-1).lt.1.0E-06)then
          C=(0.0d0,0.0d0)
          C=0.5d0*this%G(1,1)*B(1)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*this%G(2,1)*B(2)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*this%G(3,1)*B(3)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          val=val+C!CMPLX(C,0,8)
          C=(0.0d0,0.0d0)
          C=0.5d0*this%G(1,2)*B(1)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*this%G(2,2)*B(2)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          C=C+0.5d0*this%G(3,2)*B(3)*dsqrt((spin-mapp(2))*(spin+mapp(2)+1))
          val=val+C*CMPLX(0.d0,-1.0d0,8)
         endif
        !Sz
         if(abs(mapp(1)-mapp(2)).lt.1.0E-06)then
!           val=val+CMPLX(this%G(1,3)*B(1)*mapp(2)+this%G(2,3)*B(2)   &
!                  *mapp(2)+this%G(3,3)*B(3)*mapp(2),0.0d0,8)
           val=val+this%G(1,3)*B(1)*mapp(2)+this%G(2,3)*B(2)   &
                  *mapp(2)+this%G(3,3)*B(3)*mapp(2)
         endif

         val=-1.0d0*val*bohr_mag ! bohr magneton in cm-1/T
        
        return
        end function G_matelem


        function OSI_matelem(this,mapp,spin) result(val)
        use stevens_class
        implicit none
        class(OSItensor)                         :: this
        complex(8)                               :: val,val_tmp
        double precision                         :: mapp(2),R,C
        double precision                         :: spin
        integer                                  :: i

         val=(0.0d0,0.0d0)
         
         do i=1,2*this%k+1
          val_tmp=(0.0d0,0.0d0)          
          call stevens_mat_elem(this%k,this%q(i),spin,mapp(1),spin,mapp(2),val_tmp)
          val=val+val_tmp*this%B(i)
         enddo

        return
        end function OSI_matelem


        function DSI_matelem(this,mapp,spin) result(val)
        implicit none
        class(DSItensor)                         :: this
        complex(8)                               :: val,R,C
        double precision                         :: mapp(2)
        double precision                         :: spin        

         val=(0.0d0,0.0d0)

!	SzSz
         if(abs(mapp(1)-mapp(2)).lt.ZERO)then
          R=this%D(3,3)*mapp(1)*mapp(2)
          val=val+R!CMPLX(R,0.0d0,8)
         endif
        
!	S+S+
         if(abs(mapp(1)-mapp(2)-2.0d0).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)-0.25d0*this%D(2,2))
          R=R*dsqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))
          R=R*dsqrt((spin-mapp(2)-1.0d0)*(spin+mapp(2)+2.0d0))
          C=(-0.25d0*this%D(1,2)-0.25d0*this%D(2,1))
          C=C*dsqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))
          C=C*dsqrt((spin-mapp(2)-1.0d0)*(spin+mapp(2)+2.0d0))
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	S-S-
         if(abs(mapp(1)-mapp(2)+2.0d0).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)-0.25d0*this%D(2,2))
          R=R*dsqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))
          R=R*dsqrt((spin+mapp(2)-1.0d0)*(spin-mapp(2)+2.0d0))
          C=(0.25d0*this%D(1,2)+0.25d0*this%D(2,1))
          C=C*dsqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))
          C=C*dsqrt((spin+mapp(2)-1.0d0)*(spin-mapp(2)+2.0d0))
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	S+S-
         if(abs(mapp(1)-mapp(2)).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)+0.25d0*this%D(2,2))
          R=R*(spin+mapp(2))*(spin-mapp(2)+1.0d0)
          C=(0.25d0*this%D(1,2)-0.25d0*this%D(2,1))
          C=C*(spin+mapp(2))*(spin-mapp(2)+1.0d0)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	S-S+
         if(abs(mapp(1)-mapp(2)).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)+0.25d0*this%D(2,2))*(spin-mapp(2))
          R=R*(spin+mapp(2)+1.0d0)
          C=(-0.25d0*this%D(1,2)+0.25d0*this%D(2,1))*(spin-mapp(2))
          C=C*(spin+mapp(2)+1.0d0)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	S+Sz
         if(abs(mapp(1)-mapp(2)-1.0d0).lt.ZERO)THEN
          R=0.5d0*this%D(1,3)
          R=R*sqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))*mapp(2)
          C=-0.5d0*this%D(2,3)
          C=C*sqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))*mapp(2)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif
        
!	S-Sz
         if(abs(mapp(1)-mapp(2)+1.0d0).lt.ZERO)THEN
          R=0.5d0*this%D(1,3)
          R=R*sqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))*mapp(2)
          C=0.5d0*this%D(2,3)
          C=C*sqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))*mapp(2)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	SzS+
         if(abs(mapp(1)-mapp(2)-1.0d0).lt.ZERO)THEN
          R=0.5d0*this%D(3,1)
          R=R*sqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))*(mapp(2)+1.0d0)
          C=-0.5d0*this%D(3,2)
          C=C*sqrt((spin-mapp(2))*(spin+mapp(2)+1.0d0))*(mapp(2)+1.0d0)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

!	SzS-
         if(abs(mapp(1)-mapp(2)+1.0d0).lt.ZERO)THEN
           R=0.5d0*this%D(3,1)
           R=R*sqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))*(mapp(2)-1.0d0)
           C=0.5d0*this%D(3,2)
           C=C*sqrt((spin+mapp(2))*(spin-mapp(2)+1.0d0))*(mapp(2)-1.0d0)
          val=val+R+C*CMPLX(0.0d0,1.0d0,8)
         endif

        return
        end function DSI_matelem


        function D2S_matelem(this,mapp,spin) result(val)
        implicit none
        class(D2Stensor)                         :: this
        complex(8)                               :: val
        double precision                         :: mapp(2,2),R,C
        double precision                         :: spin(2)

         val=(0.0d0,0.0d0) 

!	SzSz
         if(abs(mapp(1,1)-mapp(2,1)).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)).lt.ZERO)THEN
           R=this%D(3,3)*mapp(2,1)*mapp(2,2)
           val=val+R !CMPLX(R,0.0d0,8)
         endif
        
!	S+S+
         if(abs(mapp(1,1)-mapp(2,1)-1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)-1).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)-0.25d0*this%D(2,2))
          R=R*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))
          R=R*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))
          C=(-0.25d0*this%D(1,2)-0.25d0*this%D(2,1))
          C=C*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))
          C=C*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))
          val=val+CMPLX(R,C,8)
         endif

!	S-S-
         if(abs(mapp(1,1)-mapp(2,1)+1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)+1).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)-0.25d0*this%D(2,2))
          R=R*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))
          R=R*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))
          C=(0.25d0*this%D(1,2)+0.25d0*this%D(2,1))
          C=C*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))
          C=C*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))
          val=val+CMPLX(R,C,8)
         endif

!	S+S-
         if(abs(mapp(1,1)-mapp(2,1)-1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)+1).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)+0.25d0*this%D(2,2))
          R=R*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))
          R=R*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))
          C=(0.25d0*this%D(1,2)-0.25d0*this%D(2,1))
          C=C*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))
          C=C*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))
          val=val+CMPLX(R,C,8)
         endif

!	S-S+
         if(abs(mapp(1,1)-mapp(2,1)+1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)-1).lt.ZERO)THEN
          R=(0.25d0*this%D(1,1)+0.25d0*this%D(2,2))
          R=R*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))
          R=R*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))
          C=(-0.25d0*this%D(1,2)+0.25d0*this%D(2,1))
          C=C*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1)) 
          C=C*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))
          val=val+CMPLX(R,C,8)
         endif

!	S+Sz
         if(abs(mapp(1,1)-mapp(2,1)-1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)).lt.ZERO)THEN
          R=0.5d0*this%D(1,3)
          R=R*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))*mapp(2,2)
          C=-0.5d0*this%D(2,3)
          C=C*dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))*mapp(2,2)
          val=val+CMPLX(R,C,8)
         endif
        
!	S-Sz
         if(abs(mapp(1,1)-mapp(2,1)+1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)).lt.ZERO)THEN
          R=0.5d0*this%D(1,3)
          R=R*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))*mapp(2,2)
          C=0.5d0*this%D(2,3)
          C=C*dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))*mapp(2,2)
          val=val+CMPLX(R,C,8)
         endif

!	SzS+
         if(abs(mapp(1,1)-mapp(2,1)).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)-1).lt.ZERO)THEN
          R=0.5d0*this%D(3,1)
          R=R*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))*mapp(2,1)
          C=-0.5d0*this%D(3,2)
          C=C*dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))*mapp(2,1)
          val=val+CMPLX(R,C,8)
         endif

!	SzS-
         if(abs(mapp(1,1)-mapp(2,1)).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)+1).lt.ZERO)THEN
          R=0.5d0*this%D(3,1)
          R=R*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))*mapp(2,1)
          C=0.5d0*this%D(3,2)
          C=C*dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))*mapp(2,1)
          val=val+CMPLX(R,C,8)
         endif

        return
        end function D2S_matelem


        function Jiso_matelem(this,mapp,spin) result(val)
        implicit none
        class(Jiso)                              :: this
        complex(8)                               :: val
        double precision                         :: mapp(2,2),R,C
        double precision                         :: spin(2)

         val=(0.0d0,0.0d0)
        
         if(abs(mapp(1,1)-mapp(2,1)).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)).lt.ZERO)THEN
          val=this%J*mapp(1,1)*mapp(1,2)
         endif
                
         if(abs(mapp(1,1)-mapp(2,1)-1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)+1).lt.ZERO)THEN
          val=this%J*                                           &
          dsqrt((spin(1)-mapp(2,1))*(spin(1)+mapp(2,1)+1))*     & 
          dsqrt((spin(2)+mapp(2,2))*(spin(2)-mapp(2,2)+1))*0.5
         endif

         if(abs(mapp(1,1)-mapp(2,1)+1).lt.ZERO.and.abs(mapp(1,2)-mapp(2,2)-1).lt.ZERO)THEN
          val=this%J*                                           & 
          dsqrt((spin(1)+mapp(2,1))*(spin(1)-mapp(2,1)+1))*     & 
          dsqrt((spin(2)-mapp(2,2))*(spin(2)+mapp(2,2)+1))*0.5
         endif

        return
        end function Jiso_matelem

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!
        !!!!! Rotation functions
        !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine rot_g_tensor(this,alpha,beta,gamma)
        implicit none
        class(Gtensor)                         :: this
        double complex                         :: rotG(3,3)
        double complex                         :: rotmat(3,3)
        double precision                       :: alpha,beta,gamma

         rotmat(1,1)=cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) 
         rotmat(2,1)=cos(alpha)*sin(gamma)+cos(beta)*cos(gamma)*sin(alpha)
         rotmat(3,1)=-cos(gamma)*sin(beta)

         rotmat(1,2)=-cos(gamma)*sin(alpha)-cos(alpha)*cos(beta)*sin(gamma)
         rotmat(2,2)=cos(alpha)*cos(gamma)-cos(beta)*sin(gamma)*sin(alpha)
         rotmat(3,2)=sin(beta)*sin(gamma)

         rotmat(1,3)=cos(alpha)*sin(beta)
         rotmat(2,3)=sin(alpha)*sin(beta)
         rotmat(3,3)=cos(beta)

         rotG=matmul(transpose(rotmat),this%G)
         rotG=matmul(rotG,rotmat)

         this%G=rotG

        return
        end subroutine rot_g_tensor

        subroutine rot_cart_tensor(this,alpha,beta,gamma)
        implicit none
        class(cart_tensor)                     :: this
        double complex                         :: rotD(3,3)
        double complex                         :: rotmat(3,3)
        double precision                       :: alpha,beta,gamma
                                 

         rotmat(1,1)=cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) 
         rotmat(2,1)=cos(alpha)*sin(gamma)+cos(beta)*cos(gamma)*sin(alpha)
         rotmat(3,1)=-cos(gamma)*sin(beta)

         rotmat(1,2)=-cos(gamma)*sin(alpha)-cos(alpha)*cos(beta)*sin(gamma)
         rotmat(2,2)=cos(alpha)*cos(gamma)-cos(beta)*sin(gamma)*sin(alpha)
         rotmat(3,2)=sin(beta)*sin(gamma)

         rotmat(1,3)=cos(alpha)*sin(beta)
         rotmat(2,3)=sin(alpha)*sin(beta)
         rotmat(3,3)=cos(beta)

         rotD=matmul(transpose(rotmat),this%D)
         rotD=matmul(rotD,rotmat)
         
         this%D=rotD

        return
        end subroutine rot_cart_tensor

        subroutine rot_O(O,alpha,beta,gamma)
        use stevens_class
        use rotations_class
        implicit none
        class(OSItensor)                     :: O
        type(OSItensor)                      :: newO
        integer                              :: k,l
        double complex, allocatable          :: rot(:,:)
        double precision                     :: alpha,beta,gamma

         call Rot_Wig(DBLE(O%k),alpha,beta,gamma,rot)

         newO%k=O%k
         allocate(newO%q(2*newO%k+1))
         allocate(newO%B(2*newO%k+1))
         newO%q=O%q
         newO%B=(0.0d0,0.0d0)

         do l=1,2*O%k+1

          if(newO%q(l).eq.0)then
           newO%B(l)=O%B(l)
          endif

          if(newO%q(l).gt.0)then
           k=newO%k-newO%q(l)+1
           newO%B(l)=O%B(l)+CMPLX(-AIMAG(O%B(k)),DBLE(O%B(k)),8)
           newO%B(l)=newO%B(l)*((-1)**(O%q(k)))
           newO%B(l)=newO%B(l)/dsqrt(2.0d0)
          endif

          if(newO%q(l).lt.0)then
           k=newO%k-newO%q(l)+1
           newO%B(l)=O%B(k)-CMPLX(-AIMAG(O%B(l)),DBLE(O%B(l)),8)
           newO%B(l)=newO%B(l)/dsqrt(2.0d0)
          endif

        enddo

        O%B=newO%B
        newO%B=(0.0d0,0.0d0)
        newO%B=matmul(rot,O%B)
        O%B=(0.0d0,0.0d0)

        do l=1,2*O%k+1
         if(O%q(l).eq.0)then
          O%B(l)=newO%B(l)
         endif
         if(O%q(l).gt.0)then
          k=newO%k-newO%q(l)+1
          O%B(l)=(-1.0d0)**(O%q(l))*newO%B(l)+newO%B(k)
          O%B(l)=O%B(l)*dsqrt(2.0d0)/2.0d0
         endif
         if(O%q(l).lt.0)then
          k=newO%k-newO%q(l)+1
          O%B(l)=(-1.0d0)**(O%q(k)+1)*newO%B(k)+newO%B(l)
          O%B(l)=O%B(l)*dsqrt(2.0d0)/2.0d0
          O%B(l)=CMPLX(-AIMAG(O%B(l)),DBLE(O%B(l)),8)
         endif
        enddo

        deallocate(rot)

        return
        end subroutine rot_O

!        subroutine rot_O(O,alpha,beta,gamma)
!        use stevens_class
!        use rotations_class
!        implicit none
!        class(OSItensor)                     :: O
!        type(OSItensor)                      :: newO
!        integer                              :: k,l
!        double complex, allocatable          :: rot(:,:)
!        double precision                     :: alpha,beta,gamma

!         call Rot_Wig(DBLE(k),alpha,beta,gamma,rot)
        
!         newO%k=O%k
!         allocate(newO%q(2*newO%q))
!         do l=1,2*O%k+1

!          newO%q(l)=O%q(l)

!          if(newO%q(l).eq.0)then
!           newO%B(l)=O%B(l)
!          endif

!          if(newO%q(l).gt.0)then
!           k=newO%k-newO%q(l)+1
!           newO%B(l)=O%B(l)-CMPLX(-AIMAG(O%B(k)),DBLE(O%B(k)),8)
!           newO%B(l)=newO%B(l)*((-1)**(O%q(l)))
!           newO%B(l)=newO%B(l)/dsqrt(2.0d0)
!          endif

!          if(newO%q(l).lt.0)then
!           k=newO%k-newO%q(l)+1
!           newO%B(l)=O%B(k)+CMPLX(-AIMAG(O%B(l)),DBLE(O%B(l)),8)
!           newO%B(l)=newO%B(l)/dsqrt(2.0d0)
!          endif

!        enddo

!        O%q=newO%q
!        newO%B=(0.0d0,0.0d0)
        
!        newO%B=matmul(rot,O%B)

!        do l=1,2*O%k+1
!         if(O%q(l).eq.0)then
!          O%B(l)=newO%B(l)
!         endif
!         if(O%q(l).gt.0)then
!          k=newO%k-newO%q(l)+1
!          O%B(l)=((-1.0d0)**((O%q(l))))*newO%B(l)+newO%B(k)
!          O%B(l)=O%B(l)*dsqrt(2.0d0)/2.0d0
!         endif
!         if(O%q(l).lt.0)then
!          k=newO%k-newO%q(l)+1
!          O%B(l)=((-1.0d0)**(ABS(O%q(l))))*newO%B(k)-newO%B(l)
!          O%B(l)=O%B(l)*dsqrt(2.0d0)/2.0d0
!          O%B(l)=CMPLX(-AIMAG(O%B(l)),DBLE(O%B(l)),8)
!         endif
!        enddo

!        return
!        end subroutine rot_O

        subroutine spher2stev(this,O)
        implicit none
        class(OSItensor)        :: this
        type(OSItensor)         :: O
        integer                 :: l

         O%k=this%k
         if(allocated(O%q)) deallocate(O%q)
         if(allocated(O%B)) deallocate(O%B)

         do l=1,2*this%k+1
          O%q(l)=this%q(l)
          O%B(l)=(0.0d0,0.0d0)
         enddo

         do l=1,2*this%k+1

          if(this%q(l).eq.0)then
           O%B(l)=this%B(l)
          endif
          if(this%q(l).gt.0)then
           O%B(l)=((-1)**((O%q(l))))*this%B(l)+this%B(2*this%k+2-l)
           O%B(l)=O%B(l)*dsqrt(2.0d0)/2
          endif
          if(this%q(l).lt.0)then
           O%B(l)=((-1)**(ABS(O%q(2*this%k+2-l)+1)))*this%B(2*this%k+2-l)+this%B(l)
           O%B(l)=O%B(l)*dsqrt(2.0d0)/2
           O%B(l)=CMPLX(-DIMAG(O%B(l)),DBLE(O%B(l)),8)
          endif

         enddo

        return
        end subroutine spher2stev

        subroutine stev2spher(this,O)
        implicit none
        class(OSItensor)        :: this
        type(OSItensor)         :: O
        integer                 :: l

         O%k=this%k
         if(allocated(O%q)) deallocate(O%q)
         if(allocated(O%B)) deallocate(O%B)

         do l=1,2*this%k+1
          O%q(l)=this%q(l)
          O%B(l)=(0.0d0,0.0d0)
         enddo

         do l=1,2*this%k+1

          if(O%q(l).eq.0)then
           O%B(l)=this%B(l)
          endif           

          if(O%q(l).gt.0)then
           O%B(l)=this%B(l)-CMPLX(-DIMAG(this%B(2*this%k+1-l+1)),DBLE(this%B(l-this%k-1)),8)
           O%B(l)=O%B(l)*((-1)**(this%B(2*this%k+1-l+1)))
           O%B(l)=O%B(l)/dsqrt(2.0d0)
          endif

          if(O%q(l).lt.0)then
           O%B(l)=this%B(2*this%k+1-l+1)+CMPLX(-DIMAG(this%B(l)),DBLE(this%B(l)),8)
           O%B(l)=O%B(l)/dsqrt(2.0d0)
          endif

         enddo


        return
        end subroutine stev2spher
  
        !!!!!!!!!!!!! Cartesian tensors functions

        subroutine stev2cart(this,O)
        implicit none
        class(cart_tensor)      :: this
        type(OSItensor)         :: O
        integer                 :: i
      
         this%D(3,3)=DBLE(O%B(3))*sqrt(6.0d0)/3.0d0

         this%D(1,3)=dble((O%B(2)-O%B(4)))/2.0d0
         this%D(3,1)=this%D(1,3)

         this%D(2,3)=-1.0d0*aimag(O%B(2)+O%B(4))/2.0d0
         this%D(3,2)=this%D(2,3)

         this%D(1,2)=aimag(-O%B(1)+O%B(5))/2.0d0
         this%D(2,1)=this%D(1,2)

         this%D(1,1)=(dble(O%B(1)+O%B(5))-this%D(3,3))/2.0d0
         this%D(2,2)=-this%D(3,3)-this%D(1,1)
  
        return
        end subroutine stev2cart

        subroutine cart2stev(this,O)
        implicit none
        class(cart_tensor)      :: this
        type(OSItensor)         :: O
        integer                 :: i

          O%k=2
          if(allocated(O%q)) deallocate(O%q)
          if(allocated(O%B)) deallocate(O%B)
          allocate(O%q(2*O%k+1))
          allocate(O%B(2*O%k+1))
          do i=-O%k,O%k
           O%q(i+O%k+1)=i
          enddo

          O%B(3)=3*this%D(3,3)-(this%D(1,1)+this%D(2,2)+this%D(3,3))/sqrt(6.0d0)
          O%B(2)=+0.5d0*((this%D(1,3)+this%D(3,1))+cmplx(0.0d0,1.0d0,8)*(-1.0d0*(this%D(2,3)+this%D(3,2))))
          O%B(4)=-0.5d0*((this%D(1,3)+this%D(3,1))+cmplx(0.0d0,1.0d0,8)*(this%D(2,3)+this%D(3,2)))
          O%B(1)=0.5d0*((this%D(1,1)-this%D(2,2))+cmplx(0.0d0,-1.0d0,8)*(this%D(1,2)+this%D(2,1)))
          O%B(5)=0.5d0*((this%D(1,1)-this%D(2,2))+cmplx(0.0d0,1.0d0,8)*(this%D(1,2)+this%D(2,1)))

        return
        end subroutine cart2stev

        subroutine printD(this)
        use lapack_diag_simm
        implicit none
        class(cart_tensor)      :: this
        integer                 :: s,l
        double precision        :: autoval(3),projD(3,3)
        double precision        :: axialD,rhombicE
        
         projD=this%D

         write(6,*) '############################'
         write(6,*) 'D tensor'
        do s=1,3
         write(6,*) (projD(s,l),l=1,3)
        enddo
         write(6,*) '###########################'

        call new_diag(3,projD,autoval)

         write(6,*) '############################'
         write(6,*) 'Projected D Autoket and Autoval'
         write(6,*) projD(1,1),autoval(1)
         write(6,*) projD(2,1)
         write(6,*) projD(3,1)
         write(6,*) projD(1,2),autoval(2)
         write(6,*) projD(2,2)
         write(6,*) projD(3,2)
         write(6,*) projD(1,3),autoval(3)
         write(6,*) projD(2,3)
         write(6,*) projD(3,3)
         write(6,*) '###########################'

        call Dparam(autoval,axialD,rhombicE)
         write(6,*) 'Axial parameter   D: ',axialD
         write(6,*) 'Rhombic parameter E: ',rhombicE
         write(6,*) 'Ratio parameter E/D: ',abs(rhombicE/axialD)
       

        return
        end subroutine printD


        subroutine Dparam(A,D,E)
        implicit none
        integer               :: s
        real*8                :: D,E
        real*8, DIMENSION(3)  :: A

         D=A(1)-0.5*(A(2)+A(3))
         E=abs(0.5*(A(2)-A(3)))
         if(abs(E/D).lt.0.333333)then
         return
         endif


         D=A(2)-0.5*(A(3)+A(1))
         E=abs(0.5*(A(3)-A(1)))
         if(abs(E/D).lt.0.333333)then
         return
         endif


         D=A(3)-0.5*(A(1)+A(2))
         E=abs(0.5*(A(1)-A(2)))
         if(abs(E/D).lt.0.333333)then
         return
         endif

         write(6,*) 'subroutine Dparam failed to calculate ', &
                    'D and E params'
        stop
        end subroutine Dparam

        subroutine traceless(this)
        implicit none
        class(cart_tensor)     :: this
        integer                :: s
        double precision       :: trace

         trace=0.0d0
         do s=1,3
          trace=trace+this%D(s,s)
         enddo
         do s=1,3
          this%D(s,s)=this%D(s,s)-(trace/3.0d0)
         enddo

        return
        end subroutine traceless


        !!!!!!!!!!! List functions

        subroutine add_DSInode(this_list,val)
        implicit none
        class(DSItensor_list)          :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (DSItensor)
            allocate(DSItensor::this_list%head%key)
            allocate(DSItensor::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (DSItensor)
            allocate(DSItensor::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (DSItensor)
           select type (arrow=>tmp_node%key)
            type is (DSItensor)
            arrow%D=val%D
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_DSInode


        subroutine add_D2Snode(this_list,val)
        implicit none
        class(D2Stensor_list)          :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (D2Stensor)
            allocate(D2Stensor::this_list%head%key)
            allocate(D2Stensor::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (D2Stensor)
            allocate(D2Stensor::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (D2Stensor)
           select type (arrow=>tmp_node%key)
            type is (D2Stensor)
            arrow%D=val%D
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_D2Snode

        subroutine add_Gnode(this_list,val)
        implicit none
        class(Gtensor_list)          :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (Gtensor)
            allocate(Gtensor::this_list%head%key)
            allocate(Gtensor::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (Gtensor)
            allocate(Gtensor::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (Gtensor)
           select type (arrow=>tmp_node%key)
            type is (Gtensor)
            arrow%G=val%G
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Gnode


        subroutine add_Jnode(this_list,val)
        implicit none
        class(Jiso_list)            :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (Jiso)
            allocate(Jiso::this_list%head%key)
            allocate(Jiso::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (Jiso)
            allocate(Jiso::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (Jiso)
           select type (arrow=>tmp_node%key)
            type is (Jiso)
            arrow%J=val%J
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Jnode

        subroutine add_Onode(this_list,val)
        implicit none
        class(OSItensor_list)            :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (OSItensor)
            allocate(OSItensor::this_list%head%key)
            allocate(OSItensor::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (OSItensor)
            allocate(OSItensor::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (OSItensor)
           select type (arrow=>tmp_node%key)
            type is (OSItensor)
            arrow%k=val%k
            arrow%q=val%q
            arrow%B=val%B
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Onode

        subroutine rd_Jnode(this,J)
        implicit none
        class(Jiso_list) :: this
        class(Jiso)      :: J

         select type (bho=>this%node%key)
          type is (Jiso)
           J%J=bho%J
           J%kind=bho%kind
         end select

        return
        end subroutine rd_Jnode

        subroutine rd_DSInode(this,DSI)
        implicit none
        class(DSItensor_list)  :: this
        class(DSItensor)      :: DSI
        class(*),pointer       :: bho

         select type (bho=>this%node%key)
          type is (DSItensor)
           DSI%D=bho%D
           DSI%kind=bho%kind
         end select

        return
        end subroutine rd_DSInode

        subroutine rd_D2Snode(this,D2S)
        implicit none
        class(D2Stensor_list)  :: this
        class(D2Stensor)      :: D2S
        class(*),pointer       :: bho

         select type (bho=>this%node%key)
          type is (D2Stensor)
           D2S%D=bho%D
           D2S%kind=bho%kind
         end select

        return
        end subroutine rd_D2Snode

        subroutine rd_Gnode(this,G)
        implicit none
        class(Gtensor_list)  :: this
        class(Gtensor)      :: G
        class(*),pointer       :: bho

         select type (bho=>this%node%key)
          type is (Gtensor)
           G%G=bho%G
           G%kind=bho%kind
         end select

        return
        end subroutine rd_Gnode

        subroutine rd_Onode(this,O)
        implicit none
        class(OSItensor_list)  :: this
        class(OSItensor)      :: O
        class(*),pointer       :: bho

         select type (bho=>this%node%key)
          type is (OSItensor)
           O%k=bho%k
           O%q=bho%q
           O%B=bho%B
           O%kind=bho%kind
         end select

        return
        end subroutine rd_Onode

        end module spinham_class

