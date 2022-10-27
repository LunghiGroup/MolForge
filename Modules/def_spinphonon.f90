        module spin_phonon_class
        use spinham_class
        use lists_class
        use dist_class
        implicit none

        type, extends(list) :: Jt_list
        contains
        procedure                 :: add_node => add_Jtnode        
        procedure                 :: rd_node => rd_Jtnode
        end type Jt_list 

        type, extends(list) :: Gt_list
        contains
        procedure                 :: add_node => add_Gtnode        
        procedure                 :: rd_node => rd_Gtnode
        end type Gt_list 

        type, extends(list) :: DSIt_list
        contains
        procedure                 :: add_node => add_DSItnode        
        procedure                 :: rd_node => rd_DSItnode
        end type DSIt_list 

        type, extends(list) :: D2St_list
        contains
        procedure                 :: add_node => add_D2Stnode        
        procedure                 :: rd_node => rd_D2Stnode
        end type D2St_list 

        type, extends(list) :: Ot_list
        contains
        procedure                 :: add_node => add_Otnode        
        procedure                 :: rd_node => rd_Otnode
        end type Ot_list 

        type :: sph_thermos
         integer, allocatable           :: map_s2a(:,:)
         integer                        :: nderiv
         integer                        :: norder
        end type sph_thermos

        type, extends(sph_thermos) :: Jthermos
         type(Jiso), allocatable  :: Jcart(:)
         integer                  :: kind(2)
         contains
         procedure     :: cart2brill => cart2brill_J
         procedure     :: cart2brill2 => cart2brill2_J
        end type Jthermos

        type, extends(sph_thermos) :: D2Sthermos
         type(D2Stensor), allocatable  :: Dcart(:)
         integer                       :: kind(2)
         contains
         procedure     :: cart2brill => cart2brill_D2S
         procedure     :: cart2brill2 => cart2brill2_D2S
         procedure     :: do_tinv => do_tinv_D2S
        end type D2Sthermos

        type, extends(sph_thermos) :: DSIthermos
         type(DSItensor), allocatable  :: Dcart(:)
         integer                       :: kind
         contains
         procedure     :: cart2brill => cart2brill_DSI
         procedure     :: cart2brill2 => cart2brill2_DSI
         procedure     :: do_tinv => do_tinv_DSI
        end type DSIthermos

        type, extends(sph_thermos) :: Gthermos
         type(Gtensor), allocatable    :: Gcart(:)
         integer                       :: kind
         contains
         procedure     :: cart2brill => cart2brill_G
         procedure     :: cart2brill2 => cart2brill2_G
         procedure     :: do_tinv => do_tinv_G
        end type Gthermos

        type, extends(sph_thermos) :: Othermos
         type(OSItensor), allocatable  :: Ocart(:)
         integer                       :: kind
         integer                       :: k
         contains
         procedure     :: cart2brill => cart2brill_OSI
         procedure     :: cart2brill2 => cart2brill2_OSI
         procedure     :: do_tinv => do_tinv_O
        end type Othermos

        type, extends(sph_thermos)  :: dDthermos
         type(Ddipolar)               :: dDcart(6)
         integer                      :: kind(2)
         contains        
         procedure                :: cart2brill => cart2brill_dD
         procedure                :: cart2brill2 => cart2brill2_dD
         procedure                :: make_dD   => make_dD_dipolar
         procedure                :: make_ddD   => make_ddD_dipolar
        end type dDthermos

        type, extends(SpinHamiltonian) :: SpinPhononHamiltonian
         type(Othermos),allocatable     :: O_t(:)
         type(Jthermos),allocatable     :: J_t(:)
         type(Gthermos),allocatable     :: G_t(:)
         type(DSIthermos),allocatable   :: DSI_t(:)
         type(D2Sthermos),allocatable   :: D2S_t(:)
         type(dDthermos),allocatable    :: Ddip_t(:)
         type(dist1D),allocatable       :: O_dist(:)
         type(dist1D),allocatable       :: J_dist(:)
         type(dist1D),allocatable       :: G_dist(:)
         type(dist1D),allocatable       :: DSI_dist(:)
         type(dist1D),allocatable       :: D2S_dist(:)
         type(dist1D),allocatable       :: Ddip_dist(:)
         integer, allocatable           :: mapp(:)
         contains
         procedure                      :: cart2brill
         procedure                      :: cart2brill2
         procedure                      :: spinphonon_bcast
         procedure                      :: get_dists
         procedure                      :: alloc_dists
         procedure                      :: dump_dists
         procedure                      :: merge_dists
         procedure                      :: remove_dists
        end type SpinPhononHamiltonian

        contains

        subroutine remove_dists(this)
        use dist_class
        implicit none
        class(SpinPhononHamiltonian) :: this          
        integer                      :: l

         if(this%nO.gt.0)then
          do l=1,this%nO
           call this%O_dist(l)%delete_dist()
          enddo
         endif
         if(this%nJ.gt.0)then
          do l=1,this%nJ
           call this%J_dist(l)%delete_dist()
          enddo
         endif
         if(this%nG.gt.0)then
          do l=1,this%nG
           call this%G_dist(l)%delete_dist()
          enddo
         endif
         if(this%nDSI.gt.0)then
          do l=1,this%nDSI
           call this%DSI_dist(l)%delete_dist()
          enddo
         endif
         if(this%nD2S.gt.0)then
          do l=1,this%nD2S
           call this%D2S_dist(l)%delete_dist()
          enddo
         endif
         if(this%nDdip.gt.0)then
          do l=1,this%nDdip
           call this%Ddip_dist(l)%delete_dist()
          enddo
         endif

        return
        end subroutine remove_dists

        subroutine merge_dists(this)
        use dist_class
        use mpi_utils
        use blacs_utils
        implicit none
        class(SpinPhononHamiltonian) :: this          
        integer                      :: l,s
        double precision             :: val


         if(this%nO.gt.0)then
          do l=1,this%nO
           do s=1,this%O_dist(l)%nsteps
             call mpi_allreduce(this%O_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%O_dist(l)%dist(s)=val
           enddo
          enddo
         endif

         if(this%nJ.gt.0)then
          do l=1,this%nJ
           do s=1,this%J_dist(l)%nsteps
             call mpi_allreduce(this%J_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%J_dist(l)%dist(s)=val
           enddo
          enddo
         endif

         if(this%nG.gt.0)then
          do l=1,this%nG
           do s=1,this%G_dist(l)%nsteps
             call mpi_allreduce(this%G_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%G_dist(l)%dist(s)=val
           enddo
          enddo
         endif

         if(this%nDSI.gt.0)then
          do l=1,this%nDSI
           do s=1,this%DSI_dist(l)%nsteps
             call mpi_allreduce(this%DSI_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%DSI_dist(l)%dist(s)=val
           enddo
          enddo
         endif

         if(this%nD2S.gt.0)then
          do l=1,this%nD2S
           do s=1,this%D2S_dist(l)%nsteps
             call mpi_allreduce(this%D2S_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%D2S_dist(l)%dist(s)=val
           enddo
          enddo
         endif

         if(this%nDdip.gt.0)then
          do l=1,this%nDdip
           do s=1,this%Ddip_dist(l)%nsteps
             call mpi_allreduce(this%Ddip_dist(l)%dist(s),val,1,&
                   mpi_double_precision,mpi_sum,mpi_phonons_world,err)
            this%Ddip_dist(l)%dist(s)=val
           enddo
          enddo
         endif

        return
        end subroutine merge_dists

        subroutine dump_dists(this,file_name)
        use dist_class
        implicit none
        class(SpinPhononHamiltonian) :: this          
        integer                      :: l,v,s
        character(len=*)             :: file_name
        character(len=1)             :: id

         v=1

         if(this%nO.gt.0)then
          do l=1,this%nO
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%O_dist(l)%nsteps
            write(111,*) (s-1)*this%O_dist(l)%step,this%O_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif
         if(this%nJ.gt.0)then
          do l=1,this%nJ
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%J_dist(l)%nsteps
            write(111,*) (s-1)*this%J_dist(l)%step,this%J_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif
         if(this%nG.gt.0)then
          do l=1,this%nG
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%G_dist(l)%nsteps
            write(111,*) (s-1)*this%G_dist(l)%step,this%G_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif
         if(this%nDSI.gt.0)then
          do l=1,this%nDSI
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%DSI_dist(l)%nsteps
            write(111,*) (s-1)*this%DSI_dist(l)%step,this%DSI_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif
         if(this%nD2S.gt.0)then
          do l=1,this%nD2S
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%D2S_dist(l)%nsteps
            write(111,*) (s-1)*this%D2S_dist(l)%step,this%D2S_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif
         if(this%nDdip.gt.0)then
          do l=1,this%nDdip
           write(id,"(I1)") v
           open(111,file=trim(file_name)//trim(id)//'.dat')
           do s=1,this%Ddip_dist(l)%nsteps
            write(111,*) (s-1)*this%Ddip_dist(l)%step,this%Ddip_dist(l)%dist(s)
           enddo
           close(111)
           v=v+1
          enddo
         endif

        return
        end subroutine dump_dists

        subroutine alloc_dists(this,sigma,step,nsteps)
        use dist_class
        implicit none
        class(SpinPhononHamiltonian) :: this          
        integer                      :: l,nsteps
        double precision             :: sigma,step

         if(this%nO.gt.0)then
          allocate(this%O_dist(this%nO))
          do l=1,this%nO
           this%O_dist(l)%sigma=sigma
           this%O_dist(l)%step=step
           this%O_dist(l)%nsteps=nsteps
           call this%O_dist(l)%alloc_dist()
          enddo
         endif
         if(this%nJ.gt.0)then
          allocate(this%J_dist(this%nJ))
          do l=1,this%nJ
           this%J_dist(l)%sigma=sigma
           this%J_dist(l)%step=step
           this%J_dist(l)%nsteps=nsteps
           call this%J_dist(l)%alloc_dist()
          enddo
         endif
         if(this%nG.gt.0)then
          allocate(this%G_dist(this%nG))
          do l=1,this%nG
           this%G_dist(l)%sigma=sigma
           this%G_dist(l)%step=step
           this%G_dist(l)%nsteps=nsteps
           call this%G_dist(l)%alloc_dist()
          enddo
         endif
         if(this%nDSI.gt.0)then
          allocate(this%DSI_dist(this%nDSI))
          do l=1,this%nDSI
           this%DSI_dist(l)%sigma=sigma
           this%DSI_dist(l)%step=step
           this%DSI_dist(l)%nsteps=nsteps
           call this%DSI_dist(l)%alloc_dist()
          enddo
         endif
         if(this%nD2S.gt.0)then
          allocate(this%D2S_dist(this%nD2S))
          do l=1,this%nD2S
           this%D2S_dist(l)%sigma=sigma
           this%D2S_dist(l)%step=step
           this%D2S_dist(l)%nsteps=nsteps
           call this%D2S_dist(l)%alloc_dist()
          enddo
         endif
         if(this%nDdip.gt.0)then
          allocate(this%Ddip_dist(this%nDdip))  
          do l=1,this%nDdip
           this%Ddip_dist(l)%sigma=sigma
           this%Ddip_dist(l)%step=step
           this%Ddip_dist(l)%nsteps=nsteps
           call this%Ddip_dist(l)%alloc_dist()
          enddo
         endif

        return
        end subroutine alloc_dists

        subroutine get_dists(this,freq)
        use dist_class
        implicit none
        class(SpinPhononHamiltonian) :: this          
        double precision             :: freq,val
        integer                      :: l

         do l=1,this%nJ
          call this%J(l)%get_norm(val)
          call this%J_dist(l)%update_dist(freq,val)
         enddo

         do l=1,this%nO
          call this%O(l)%get_norm(val)
          call this%O_dist(l)%update_dist(freq,val)
         enddo
         
         do l=1,this%nG
          call this%G(l)%get_norm(val)
          call this%G_dist(l)%update_dist(freq,val)
         enddo

         do l=1,this%nDSI
          call this%DSI(l)%get_norm(val)
          call this%DSI_dist(l)%update_dist(freq,val)
         enddo

         do l=1,this%nD2S
          call this%D2S(l)%get_norm(val)          
          call this%D2S_dist(l)%update_dist(freq,val)
         enddo

         do l=1,this%nDdip
          call this%Ddip(l)%get_norm(val)
          call this%Ddip_dist(l)%update_dist(freq,val)
         enddo

        return
        end subroutine get_dists

        subroutine do_tinv_O(this)
        implicit none
        class(Othermos)                          :: this
        integer                                  :: i,l,t,s,v,ncoeff
        double precision                         :: coeff
                         
        if(this%norder.eq.1)then

         do s=1,this%k*2+1
          do l=0,2
           coeff=0.0d0
           ncoeff=0
           do i=1,size(this%Ocart)
            v=mod(this%map_s2a(i,1)+2,3)
            if(v.ne.l) cycle
            coeff=coeff+this%Ocart(i)%B(s)
            ncoeff=ncoeff+1
           enddo
           do i=1,size(this%Ocart)
            v=mod(this%map_s2a(i,1)+2,3)
            if(v.ne.l) cycle
            this%Ocart(i)%B(s)=this%Ocart(i)%B(s)-coeff/dble(ncoeff)
           enddo
          enddo
         enddo

        else
         write(*,*) '2nd Order Translational Invariance for Stevens',&
                      ' Operators not yet Implemented'
        endif

        return
        end subroutine do_tinv_O

        subroutine do_tinv_D2S(this)
        implicit none
        class(D2Sthermos)                        :: this
        integer                                  :: i,l,t,s,v,ncoeff,vv
        double precision                         :: coeff
                         
         if(this%norder.eq.1)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              coeff=coeff+this%Dcart(i)%D(s,t)
              ncoeff=ncoeff+1
             enddo
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              this%Dcart(i)%D(s,t)=this%Dcart(i)%D(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo

         endif

         if(this%norder.eq.2)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              coeff=coeff+this%Dcart(i)%D(s,t)
              ncoeff=ncoeff+1
             enddo
!             write(*,*) l+1,coeff/dble(ncoeff)
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              this%Dcart(i)%D(s,t)=this%Dcart(i)%D(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo

         endif

        return
        end subroutine do_tinv_D2S

        subroutine do_tinv_DSI(this)
        implicit none
        class(DSIthermos)                        :: this
        integer                                  :: i,l,t,s,v,ncoeff,vv
        double precision                         :: coeff
          
         
         if(this%norder.eq.1)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              coeff=coeff+this%Dcart(i)%D(s,t)
              ncoeff=ncoeff+1
             enddo
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              this%Dcart(i)%D(s,t)=this%Dcart(i)%D(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo
          return
         endif
       
         if(this%norder.eq.2)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              coeff=coeff+this%Dcart(i)%D(s,t)
              ncoeff=ncoeff+1
             enddo
!             write(*,*) l+1,coeff/dble(ncoeff)
             do i=1,size(this%Dcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              this%Dcart(i)%D(s,t)=this%Dcart(i)%D(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo

         endif

        return
        end subroutine do_tinv_DSI

        subroutine do_tinv_G(this)
        implicit none
        class(Gthermos)                          :: this
        integer                                  :: i,l,t,s,v,ncoeff,vv
        double precision                         :: coeff
           
         if(this%norder.eq.1)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Gcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              coeff=coeff+this%Gcart(i)%G(s,t)
              ncoeff=ncoeff+1
             enddo
             do i=1,size(this%Gcart)
              v=mod(this%map_s2a(i,1)+2,3)
              if(v.ne.l) cycle
              this%Gcart(i)%G(s,t)=this%Gcart(i)%G(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo

         endif
       
         if(this%norder.eq.2)then

          do s=1,3
           do t=1,3
            do l=0,2
             coeff=0.0d0
             ncoeff=0
             do i=1,size(this%Gcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              coeff=coeff+this%Gcart(i)%G(s,t)
              ncoeff=ncoeff+1
             enddo
!             write(*,*) l+1,coeff/dble(ncoeff)
             do i=1,size(this%Gcart)
              v=mod(this%map_s2a(i,1)+2,3)
              vv=mod(this%map_s2a(i,3)+2,3)
              if(v.ne.l .or. vv.ne.l) cycle
              this%Gcart(i)%G(s,t)=this%Gcart(i)%G(s,t)-coeff/dble(ncoeff)
             enddo
            enddo
           enddo
          enddo

         endif

        return
        end subroutine do_tinv_G

        subroutine make_dD_dipolar (this,G1,G2,bohr_mag1,bohr_mag2,x1,dist)
        use stevens_class
        use units_parms 
        implicit none
        class(dDthermos)                        :: this
        complex(8)                               :: G1(3,3),G2(3,3)
        double precision                         :: x1(3),x2(3),dist
        double precision                         :: bohr_mag1,bohr_mag2
        integer                                  :: s,t,k,l,is,js,ii

         ii=1

         do js=1,-1,-2
          do is=1,3
           
           this%dDcart(ii)%D=0.0d0
         
           do s=1,3
            do t=1,3
             do k=1,3
              this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)-3.0d0*G1(s,k)*G2(k,t)*x1(is)*js/dist
             enddo
            enddo
           enddo

           do s=1,3
            do t=1,3
             do l=1,3
              this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)-3.0d0*G1(s,is)*x1(l)*G2(l,t)*js/dist &
                                    -3.0d0*G1(s,l)*x1(l)*G2(is,t)*js/dist
             enddo
            enddo
           enddo

           do s=1,3
            do t=1,3
             do k=1,3
              do l=1,3
               this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)+15.0d0*G1(s,k)*x1(k)*x1(l)*G2(l,t)*x1(is)*js/dist**3
              enddo
             enddo
            enddo
           enddo

           this%dDcart(ii)%D=(this%dDcart(ii)%D*mag_vacuum*bohr_mag1*bohr_mag2)/dist**4
           call this%dDcart(ii)%traceless()

           ii=ii+1

          enddo
         enddo

        return
        end subroutine make_dD_dipolar

        subroutine make_ddD_dipolar (this,G1,G2,bohr_mag1,bohr_mag2,x1,dist)
        use stevens_class
        use units_parms 
        implicit none
        class(dDthermos)                         :: this
        complex(8)                               :: G1(3,3),G2(3,3)
        double precision                         :: x1(3),x2(3),dist
        double precision                         :: bohr_mag1,bohr_mag2
        integer                                  :: s,t,k,l,is,js,ii,it,jt

         ii=1

         do js=1,-1,-2
         do is=1,3
          do jt=1,-1,-2
          do it=1,3
           
           this%dDcart(ii)%D=0.0d0
         
           do s=1,3
            do t=1,3
             do k=1,3
              this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)&
                        +15.0d0*G1(s,k)*G2(k,t)*x1(is)*js*x1(it)*jt/dist**2
             enddo
            enddo
           enddo

           if(it.eq.is)then
            do s=1,3
             do t=1,3
              do k=1,3
               this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)&
                         -3.0d0*G1(s,k)*G2(k,t)*js*jt
              enddo
             enddo
            enddo
           endif

           do s=1,3
            do t=1,3
              this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t) &
                -3.0d0*G1(s,is)*G2(it,t)*js*jt   &
                -3.0d0*G1(s,it)*G2(is,t)*js*jt   
            enddo
           enddo

           do s=1,3
            do t=1,3
             do l=1,3
              this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t) &
                +15.0d0*G1(s,l)*x1(l)*G2(is,t)*js*x1(it)*jt/dist**2   &
                +15.0d0*G1(s,is)*x1(l)*G2(l,t)*js*x1(it)*jt/dist**2   &
                +15.0d0*G1(s,it)*x1(l)*G2(l,t)*js*x1(is)*jt/dist**2   &
                +15.0d0*G1(s,l)*x1(l)*G2(it,t)*js*x1(is)*jt/dist**2   &
                +15.0d0*G1(s,l)*x1(it)*G2(l,t)*js*x1(is)*jt/dist**2   
             enddo
            enddo
           enddo

           do s=1,3
            do t=1,3
             do k=1,3
              do l=1,3
               this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)&
                -105.0d0*G1(s,k)*x1(k)*x1(l)*G2(l,t)*x1(is)*js*x1(it)*jt/dist**4
              enddo
             enddo
            enddo
           enddo

           if(is.eq.it)then
            do s=1,3
             do t=1,3
              do k=1,3
               do l=1,3
                this%dDcart(ii)%D(s,t)=this%dDcart(ii)%D(s,t)&
                 +15.0d0*G1(s,k)*x1(k)*x1(l)*G2(l,t)*js*jt/dist**2
               enddo
              enddo
             enddo
            enddo
           endif

           this%dDcart(ii)%D=(this%dDcart(ii)%D*mag_vacuum*bohr_mag1*bohr_mag2)/dist**5
           call this%dDcart(ii)%traceless()

           ii=ii+1

          enddo
          enddo
         enddo
         enddo
        
        return
        end subroutine make_ddD_dipolar

        subroutine spinphonon_bcast(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(SpinPhononHamiltonian)  :: this
        integer                       :: i,l,j,nmapp

         call mpi_bcast(this%nO,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nJ,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nG,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nDSI,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%nD2S,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(this%make_dipolar,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(this%dipolar_thr,1,mpi_double_precision,0,mpi_comm_world,err)
         if(this%make_dipolar .and. mpi_id.eq.0) nmapp=size(this%mapp)
         if(this%make_dipolar) &
          call mpi_bcast(nmapp,1,mpi_integer,0,mpi_comm_world,err)
         if(this%make_dipolar .and. .not.allocated(this%mapp)) allocate(this%mapp(nmapp)) 
         if(this%make_dipolar) &
         call mpi_bcast(this%mapp,nmapp,mpi_integer,0,mpi_comm_world,err)
         
         if(this%nJ.gt.0)then
          if(.not.allocated(this%J)) allocate(this%J(this%nJ))
          if(.not.allocated(this%J_t)) allocate(this%J_t(this%nJ))
          do i=1,this%nJ
           call mpi_bcast(this%J_t(i)%nderiv,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%J_t(i)%norder,1,mpi_integer,0,mpi_comm_world,err)
           if(.not.allocated(this%J_t(i)%Jcart)) allocate(this%J_t(i)%Jcart(this%J_t(i)%nderiv))
           if(.not.allocated(this%J_t(i)%map_s2a)) allocate(this%J_t(i)%map_s2a(this%J_t(i)%nderiv,2*this%J_t(i)%norder))
           do j=1,this%J_t(i)%nderiv
            call mpi_bcast(this%J_t(i)%map_s2a(j,:),this%J_t(i)%norder*2,mpi_integer,0,mpi_comm_world,err)
            call mpi_bcast(this%J_t(i)%Jcart(j)%J,1,mpi_double_precision,0,mpi_comm_world,err)
            call mpi_bcast(this%J_t(i)%Jcart(j)%kind,2,mpi_integer,0,mpi_comm_world,err)
           enddo
          enddo
         endif

         if(this%nG.gt.0)then

          if(.not.allocated(this%G)) allocate(this%G(this%nG))
          if(.not.allocated(this%G_t)) allocate(this%G_t(this%nG))
          do i=1,this%nG
           call mpi_bcast(this%G_t(i)%nderiv,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%G_t(i)%norder,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%G_t(i)%kind,1,mpi_integer,0,mpi_comm_world,err)
           if(.not.allocated(this%G_t(i)%Gcart)) allocate(this%G_t(i)%Gcart(this%G_t(i)%nderiv))
           if(.not.allocated(this%G_t(i)%map_s2a)) allocate(this%G_t(i)%map_s2a(this%G_t(i)%nderiv,2*this%G_t(i)%norder))
           do j=1,this%G_t(i)%nderiv
            call mpi_bcast(this%G_t(i)%map_s2a(j,:),this%G_t(i)%norder*2,mpi_integer,0,mpi_comm_world,err)
            do l=1,3
             call mpi_bcast(this%G_t(i)%Gcart(j)%G(l,:),3,mpi_double_complex,0,mpi_comm_world,err)
            enddo
            call mpi_bcast(this%G_t(i)%Gcart(j)%kind,1,mpi_integer,0,mpi_comm_world,err)
           enddo
          enddo

          do i=1,size(this%G_t)
           call this%G_t(i)%do_tinv()
          enddo

         endif

         if(this%nDSI.gt.0)then
          if(.not.allocated(this%DSI)) allocate(this%DSI(this%nDSI))
          if(.not.allocated(this%DSI_t)) allocate(this%DSI_t(this%nDSI))
          do i=1,this%nDSI
           call mpi_bcast(this%DSI_t(i)%nderiv,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%DSI_t(i)%norder,1,mpi_integer,0,mpi_comm_world,err)
           if(.not.allocated(this%DSI_t(i)%Dcart)) allocate(this%DSI_t(i)%Dcart(this%DSI_t(i)%nderiv))
           if(.not.allocated(this%DSI_t(i)%map_s2a)) allocate(this%DSI_t(i)%map_s2a(this%DSI_t(i)%nderiv,2*this%DSI_t(i)%norder))
           do j=1,this%DSI_t(i)%nderiv
            call mpi_bcast(this%DSI_t(i)%map_s2a(j,:),this%DSI_t(i)%norder*2,mpi_integer,0,mpi_comm_world,err)
            do l=1,3
             call mpi_bcast(this%DSI_t(i)%Dcart(j)%D(l,:),3,mpi_double_complex,0,mpi_comm_world,err)
            enddo
            call mpi_bcast(this%DSI_t(i)%Dcart(j)%kind,2,mpi_integer,0,mpi_comm_world,err)
           enddo
          enddo
          do i=1,size(this%DSI_t)
           call this%DSI_t(i)%do_tinv()
          enddo
         endif


         if(this%nD2S.gt.0)then
          if(.not.allocated(this%D2S)) allocate(this%D2S(this%nD2S))
          if(.not.allocated(this%D2S_t)) allocate(this%D2S_t(this%nD2S))
          do i=1,this%nD2S
           call mpi_bcast(this%D2S_t(i)%nderiv,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%D2S_t(i)%norder,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%D2S_t(i)%kind,2,mpi_integer,0,mpi_comm_world,err)
           if(.not.allocated(this%D2S_t(i)%Dcart)) allocate(this%D2S_t(i)%Dcart(this%D2S_t(i)%nderiv))
           if(.not.allocated(this%D2S_t(i)%map_s2a)) allocate(this%D2S_t(i)%map_s2a(this%D2S_t(i)%nderiv,2*this%D2S_t(i)%norder))
           do j=1,this%D2S_t(i)%nderiv
            call mpi_bcast(this%D2S_t(i)%map_s2a(j,:),this%D2S_t(i)%norder*2,mpi_integer,0,mpi_comm_world,err)
            do l=1,3
             call mpi_bcast(this%D2S_t(i)%Dcart(j)%D(l,:),3,mpi_double_complex,0,mpi_comm_world,err)
            enddo
            call mpi_bcast(this%D2S_t(i)%Dcart(j)%kind,2,mpi_integer,0,mpi_comm_world,err)
           enddo
          enddo
          do i=1,size(this%D2S_t)
           call this%D2S_t(i)%do_tinv()
          enddo
         endif

         if(this%nO.gt.0)then
          if(.not.allocated(this%O)) allocate(this%O(this%nO))
          if(.not.allocated(this%O_t)) allocate(this%O_t(this%nO))
          do i=1,this%nO           
           call mpi_bcast(this%O_t(i)%nderiv,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%O_t(i)%norder,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%O_t(i)%k,1,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(this%O_t(i)%kind,1,mpi_integer,0,mpi_comm_world,err)
           if(.not.allocated(this%O_t(i)%Ocart)) allocate(this%O_t(i)%Ocart(this%O_t(i)%nderiv))
           if(.not.allocated(this%O_t(i)%map_s2a)) allocate(this%O_t(i)%map_s2a(this%O_t(i)%nderiv,2*this%O_t(i)%norder))
           do j=1,this%O_t(i)%nderiv
            call mpi_bcast(this%O_t(i)%map_s2a(j,:),this%O_t(i)%norder*2,mpi_integer,0,mpi_comm_world,err)
            call mpi_bcast(this%O_t(i)%Ocart(j)%k,1,mpi_integer,0,mpi_comm_world,err)
            if(.not. allocated(this%O_t(i)%Ocart(j)%B) ) &
                allocate(this%O_t(i)%Ocart(j)%B(2*this%O_t(i)%Ocart(j)%k+1))
            if(.not. allocated(this%O_t(i)%Ocart(j)%q) ) &
                allocate(this%O_t(i)%Ocart(j)%q(2*this%O_t(i)%Ocart(j)%k+1))
            do l=1,2*this%O_t(i)%Ocart(j)%k+1
             call mpi_bcast(this%O_t(i)%Ocart(j)%B(l),1,mpi_double_precision,0,mpi_comm_world,err)
             call mpi_bcast(this%O_t(i)%Ocart(j)%q(l),1,mpi_integer,0,mpi_comm_world,err)
            enddo
            call mpi_bcast(this%O_t(i)%Ocart(j)%kind,1,mpi_integer,0,mpi_comm_world,err)
           enddo
          enddo
          do i=1,size(this%O_t)
           call this%O_t(i)%do_tinv()
          enddo
         endif

        return
        end subroutine spinphonon_bcast

        subroutine cart2brill(this,rcell,sys,phondy,ki,kj)
        use phonons_class
        use atoms_class
        use proj_disp_class
        use units_parms
        implicit none
        class(SpinPhononHamiltonian)    :: this
        class(brillouin)                :: phondy
        class(atoms_group)              :: sys
        type(molecule)                  :: mol
        integer                         :: ki,kj,l,i,j
        complex(8), allocatable         :: hess(:)
        double precision                :: k(3),mass1,coeff
        double precision, allocatable   :: rcell(:,:)
        integer                         :: loc_nats,i1,s1
        double precision, allocatable   :: mass(:),geo(:,:)
        logical                         :: project_intern=.false.

         if(allocated(hess)) deallocate(hess)
         allocate(hess(size(phondy%list(ki)%hess,1)))
         hess=(0.0d0,0.0d0)

         if(project_intern)then

          loc_nats=30

          allocate(mass(loc_nats))
          allocate(geo(loc_nats,3))
          geo=sys%x(1:loc_nats,1:3)

          call mol%def_mol(geo,mass)
          geo=0.0d0
          do j=1,loc_nats*3
           mass1=sys%mass(sys%kind((2+j)/3))
           coeff=bohr2ang/dsqrt(mass1*1822.89/219474.6313702)
           i1=(2+j)/3
           s1=mod(j-1,3)+1
           geo(i1,s1)=coeff*dble(phondy%list(ki)%hess(j,kj))/sqrt(phondy%list(ki)%freq(kj))
          enddo

          geo=geo+sys%x(1:loc_nats,:)

          call mol%def_mol_dist(geo)
          call mol%proj_disp()
          
          l=1
          do i=1,loc_nats
           do j=1,3
            hess(l)=cmplx(mol%cart_int(i,j),0.0d0,8)
!            mass1=sys%mass(sys%kind((2+j)/3))
!            coeff=bohr2ang/dsqrt(mass1*1822.89/219474.6313702)
!            write(*,*) 'Mode',ki,kj,hess(l),coeff*dble(phondy%list(ki)%hess(l,kj))/sqrt(phondy%list(ki)%freq(kj))
            l=l+1
           enddo
          enddo

         else
         
!          writE(*,*) '####',ki,kj,phondy%liest(ki)%hess(l,kj)

          l=1
          do i=1,sys%nats
           mass1=sys%mass(sys%kind(i))*1822.89        
           coeff=0.5291772d0/sqrt(mass1*phondy%ntot*phondy%list(ki)%freq(kj)/219474.6313702d0)
           do j=1,3
            hess(l)=coeff*phondy%list(ki)%hess(l,kj)
            l=l+1
           enddo
          enddo

         endif

         do l=1,this%nJ
          call this%J_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%J(l))
         enddo

         do l=1,this%nG
          call this%G_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%G(l))
         enddo

         do l=1,this%nO
          call this%O_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%O(l))
         enddo

         do l=1,this%nDSI
          call this%DSI_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%DSI(l))
         enddo

         do l=1,this%nD2S
          call this%D2S_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%D2S(l))
         enddo
         
         do l=1,this%nDdip
          call this%Ddip_t(l)%cart2brill(hess,phondy%list(ki)%k,rcell,this%Ddip(l))
         enddo

        return
        end subroutine cart2brill

        subroutine cart2brill2(this,rcell,sys,phondy,ki,kj,ki2,kj2)
        use phonons_class
        use atoms_class
        use proj_disp_class
        use units_parms
        implicit none
        class(SpinPhononHamiltonian)    :: this
        class(brillouin)                :: phondy
        class(atoms_group)              :: sys
        type(molecule)                  :: mol
        integer                         :: ki,kj,l,i,j,ki2,kj2
        complex(8), allocatable         :: hess(:),hess2(:)
        double precision                :: k(3),mass1,coeff
        double precision, allocatable   :: rcell(:,:)
        integer                         :: loc_nats,i1,s1
        double precision, allocatable   :: mass(:),geo(:,:)
        logical                         :: project_intern=.false.

         if(allocated(hess)) deallocate(hess)
         allocate(hess(size(phondy%list(ki)%hess,1)))
         hess=(0.0d0,0.0d0)
         if(allocated(hess2)) deallocate(hess2)
         allocate(hess2(size(phondy%list(ki2)%hess,1)))
         hess2=(0.0d0,0.0d0)
         
         l=1
         do i=1,sys%nats
          mass1=sys%mass(sys%kind(i))*1822.89
          coeff=0.5291772d0/sqrt(mass1*phondy%ntot*phondy%list(ki)%freq(kj)/219474.6313702d0)
          do j=1,3
           hess(l)=coeff*phondy%list(ki)%hess(l,kj)
           l=l+1
          enddo
         enddo
         l=1
         do i=1,sys%nats
          mass1=sys%mass(sys%kind(i))*1822.89
          coeff=0.5291772d0/sqrt(mass1*phondy%ntot*phondy%list(ki2)%freq(kj2)/219474.6313702d0)
          do j=1,3
           hess2(l)=coeff*phondy%list(ki2)%hess(l,kj2)
           l=l+1
          enddo
         enddo

         do l=1,this%nJ
          call this%J_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%J(l))
         enddo

         do l=1,this%nG
          call this%G_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%G(l))
         enddo

         do l=1,this%nO
          call this%O_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%O(l))
         enddo

         do l=1,this%nDSI
          call this%DSI_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%DSI(l))
         enddo

         do l=1,this%nD2S
          call this%D2S_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%D2S(l))
         enddo
         
         do l=1,this%nDdip
          call this%Ddip_t(l)%cart2brill2(hess,phondy%list(ki)%k,hess2,phondy%list(ki2)%k,rcell,this%Ddip(l))
         enddo

        return
        end subroutine cart2brill2

        subroutine cart2brill_J(this,hess,k,rcell,Jtmp)
        use spinham_class
        implicit none
        class(Jthermos)                          :: this
        type(Jiso)                               :: Jtmp
        integer                                  :: i,l,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                         
         Jtmp%kind=this%Jcart(1)%kind
         Jtmp%J=0.0d0

         do i=1,size(this%Jcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          Jtmp%J=Jtmp%J+this%Jcart(i)%J*hess(l)*exp(coeff)
         enddo

        return
        end subroutine cart2brill_J

        subroutine cart2brill2_J(this,hess,k,hess2,k2,rcell,Jtmp)
        use spinham_class
        implicit none
        class(Jthermos)                          :: this
        type(Jiso)                               :: Jtmp
        integer                                  :: i,l,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         Jtmp%kind=this%Jcart(1)%kind
         Jtmp%J=0.0d0

         do i=1,size(this%Jcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)
          coeff2=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           Jtmp%J=Jtmp%J+this%Jcart(i)%J*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)
          else
           Jtmp%J=Jtmp%J+this%Jcart(i)%J*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)+&
                         this%Jcart(i)%J*hess(l2)*exp(coeff)*hess2(l)*exp(coeff2)
          endif
         enddo

        return
        end subroutine cart2brill2_J

        subroutine cart2brill_dD(this,hess,k,rcell,dDtmp)
        use spinham_class
        implicit none
        class(dDthermos)                         :: this
        type(Ddipolar)                           :: dDtmp
        integer                                  :: i,l,t,s,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                         
         dDtmp%kind(1)=this%kind(1)
         dDtmp%kind(2)=this%kind(2)
         dDtmp%D=0.0d0

         do i=1,size(this%dDcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          do s=1,3
           do t=1,3           
            dDtmp%D(s,t)=dDtmp%D(s,t)+this%dDcart(i)%D(s,t)*hess(l)*exp(coeff)
           enddo
          enddo
         enddo

        return
        end subroutine cart2brill_dD

        subroutine cart2brill2_dD(this,hess,k,hess2,k2,rcell,dDtmp)
        use spinham_class
        implicit none
        class(dDthermos)                         :: this
        type(Ddipolar)                           :: dDtmp
        integer                                  :: i,l,t,s,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         dDtmp%kind(1)=this%kind(1)
         dDtmp%kind(2)=this%kind(2)
         dDtmp%D=(0.0d0,0.0d0)

         do i=1,size(this%dDcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)
          coeff2=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           do s=1,3
            do t=1,3
             dDtmp%D(s,t)=dDtmp%D(s,t)+this%dDcart(i)%D(s,t)*hess(l)*exp(coeff)*&
                          hess2(l2)*exp(coeff2)
            enddo
           enddo
          else
           do s=1,3
            do t=1,3
             dDtmp%D(s,t)=dDtmp%D(s,t)+this%dDcart(i)%D(s,t)*hess(l)*exp(coeff)* &
                          hess2(l2)*exp(coeff2)+ &
                          this%dDcart(i)%D(s,t)*hess(l2)*exp(coeff)* &
                          hess2(l)*exp(coeff2)
            enddo
           enddo
          endif
         enddo

        return
        end subroutine cart2brill2_dD

        subroutine cart2brill_D2S(this,hess,k,rcell,D2Stmp)
        use spinham_class
        implicit none
        class(D2Sthermos)                        :: this
        type(D2Stensor)                          :: D2Stmp
        integer                                  :: i,l,t,s,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                         
         D2Stmp%kind=this%Dcart(1)%kind
         D2Stmp%D=(0.0d0,0.0d0)

         do i=1,size(this%Dcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          do s=1,3
           do t=1,3
            D2Stmp%D(s,t)=D2Stmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)
           enddo
          enddo
         enddo

        return
        end subroutine cart2brill_D2S

        subroutine cart2brill2_D2S(this,hess,k,hess2,k2,rcell,D2Stmp)
        use spinham_class
        implicit none
        class(D2Sthermos)                        :: this
        type(D2Stensor)                          :: D2Stmp
        integer                                  :: i,l,t,s,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         D2Stmp%kind=this%Dcart(1)%kind
         D2Stmp%D=0.0d0

         do i=1,size(this%Dcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)
          coeff2=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           do s=1,3
            do t=1,3
             D2Stmp%D(s,t)=D2Stmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)*&
                           hess2(l2)*exp(coeff2)
            enddo
           enddo
          else
           do s=1,3
            do t=1,3
             D2Stmp%D(s,t)=D2Stmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)*&
                           hess2(l2)*exp(coeff2)+&
                           this%Dcart(i)%D(s,t)*hess(l2)*exp(coeff)*&
                           hess2(l)*exp(coeff2)
            enddo
           enddo
          endif
         enddo

        return
        end subroutine cart2brill2_D2S

        subroutine cart2brill_DSI(this,hess,k,rcell,DSItmp)
        use spinham_class
        implicit none
        class(DSIthermos)                        :: this
        type(DSItensor)                          :: DSItmp
        integer                                  :: i,l,t,s,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                         
         DSItmp%kind=this%Dcart(1)%kind
         DSItmp%D=0.0d0

         do i=1,size(this%Dcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          do s=1,3
           do t=1,3
            DSItmp%D(s,t)=DSItmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)
           enddo
          enddo
         enddo

        return
        end subroutine cart2brill_DSI

        subroutine cart2brill2_DSI(this,hess,k,hess2,k2,rcell,DSItmp)
        use spinham_class
        implicit none
        class(DSIthermos)                        :: this
        type(DSItensor)                          :: DSItmp
        integer                                  :: i,l,t,s,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         DSItmp%kind=this%Dcart(1)%kind
         DSItmp%D=0.0d0

         do i=1,size(this%Dcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           do s=1,3
            do t=1,3
             DSItmp%D(s,t)=DSItmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)*&
                                         hess2(l2)*exp(coeff2)
            enddo
           enddo
          else
           do s=1,3
            do t=1,3
             DSItmp%D(s,t)=DSItmp%D(s,t)+this%Dcart(i)%D(s,t)*hess(l)*exp(coeff)*&
                                         hess2(l2)*exp(coeff2)+&
                                         this%Dcart(i)%D(s,t)*hess(l2)*exp(coeff)*&
                                         hess2(l)*exp(coeff2)
            enddo
           enddo
          endif
         enddo

        return
        end subroutine cart2brill2_DSI

        subroutine cart2brill_G(this,hess,k,rcell,Gtmp)
        use spinham_class
        implicit none
        class(Gthermos)                          :: this
        type(Gtensor)                            :: Gtmp
        integer                                  :: i,l,t,s,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                         
         Gtmp%kind=this%Gcart(1)%kind
         Gtmp%G=(0.0d0,0.0d0)

         do i=1,size(this%Gcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)          
          coeff=cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          do s=1,3
           do t=1,3
            Gtmp%G(s,t)=Gtmp%G(s,t)+&
              this%Gcart(i)%G(s,t)*hess(l)*exp(coeff)
           enddo
          enddo
         enddo

        return
        end subroutine cart2brill_G

        subroutine cart2brill2_G(this,hess,k,hess2,k2,rcell,Gtmp)
        use spinham_class
        implicit none
        class(Gthermos)                          :: this
        type(Gtensor)                            :: Gtmp
        integer                                  :: i,l,t,s,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         Gtmp%kind=this%Gcart(1)%kind
         Gtmp%G=(0.0d0,0.0d0)

         do i=1,size(this%Gcart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)          
          coeff=cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)          
          coeff2=cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           do s=1,3
            do t=1,3
             gtmp%g(s,t)=gtmp%g(s,t)+&
               this%gcart(i)%g(s,t)*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)
            enddo
           enddo
          else
           do s=1,3
            do t=1,3
             gtmp%g(s,t)=gtmp%g(s,t)+&
               this%gcart(i)%g(s,t)*hess(l)*exp(coeff)*hess2(l2)*exp(coeff2)+&
               this%gcart(i)%g(s,t)*hess(l2)*exp(coeff)*hess2(l)*exp(coeff2)
            enddo
           enddo
          endif

         enddo

        return
        end subroutine cart2brill2_G

        subroutine cart2brill_OSI(this,hess,k,rcell,Otmp)
        use spinham_class
        implicit none
        class(Othermos)                          :: this
        type(OSItensor)                          :: Otmp
        integer                                  :: i,l,t,v
        complex(8), allocatable                  :: hess(:)
        complex(8)                               :: coeff
        double precision                         :: k(3)
        double precision, allocatable            :: rcell(:,:)
                        
         Otmp%k=this%Ocart(1)%k
         Otmp%kind=this%Ocart(1)%kind
         if(.not. allocated(Otmp%q)) allocate(Otmp%q(Otmp%k*2+1))
         if(.not. allocated(Otmp%B)) allocate(Otmp%B(Otmp%k*2+1))
         Otmp%B=(0.0d0,0.0d0)
         Otmp%q=this%Ocart(1)%q

         do i=1,size(this%Ocart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          do t=1,2*Otmp%k+1
           Otmp%B(t)=Otmp%B(t)+this%Ocart(i)%B(t)*hess(l)*exp(coeff)   
          enddo
         enddo

        return
        end subroutine cart2brill_OSI

        subroutine cart2brill2_OSI(this,hess,k,hess2,k2,rcell,Otmp)
        use spinham_class
        implicit none
        class(Othermos)                          :: this
        type(OSItensor)                          :: Otmp
        integer                                  :: i,l,t,v,l2,v2
        complex(8), allocatable                  :: hess(:),hess2(:)
        complex(8)                               :: coeff,coeff2
        double precision                         :: k(3),k2(3)
        double precision, allocatable            :: rcell(:,:)
                         
         Otmp%k=this%Ocart(1)%k
         Otmp%kind=this%Ocart(1)%kind

         if(allocated(Otmp%q)) deallocate(Otmp%q)
         if(allocated(Otmp%B)) deallocate(Otmp%B)
         allocate(Otmp%q(Otmp%k*2+1))
         allocate(Otmp%B(Otmp%k*2+1))
         Otmp%B=(0.0d0,0.0d0)
         l=1
         do i=-2*Otmp%k,2*Otmp%k
          Otmp%q(l)=i
          l=l+1
         enddo

         do i=1,size(this%Ocart)
          l=this%map_s2a(i,1)
          v=this%map_s2a(i,2)
          coeff=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k,rcell(v,:))
          l2=this%map_s2a(i,3)
          v2=this%map_s2a(i,4)
          coeff2=cmplx(0.0d0,1.0d0,8)*2*acos(-1.0d0)*DOT_PRODUCT(k2,rcell(v2,:))
          if(l.eq.l2)then
           do t=1,2*Otmp%k+1
            Otmp%B(t)=Otmp%B(t)+this%Ocart(i)%B(t)*hess(l)*exp(coeff)*&
                         hess2(l2)*exp(coeff2)
           enddo
          else
           do t=1,2*Otmp%k+1
            Otmp%B(t)=Otmp%B(t)+this%Ocart(i)%B(t)*hess(l)*exp(coeff)*&
                         hess2(l2)*exp(coeff2)+&
                         this%Ocart(i)%B(t)*hess(l2)*exp(coeff)*&
                         hess2(l)*exp(coeff2)
           enddo
          endif
         enddo

        return
        end subroutine cart2brill2_OSI

        subroutine add_Jtnode(this_list,val)
        implicit none
        class(Jt_list)                 :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: i

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (Jthermos)
            allocate(Jthermos::this_list%head%key)
            allocate(Jthermos::this_list%node%key)
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
           type is (Jthermos)
            allocate(Jthermos::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (Jthermos)
           select type (arrow=>tmp_node%key)
            type is (Jthermos)
            arrow%nderiv=val%nderiv
            arrow%norder=val%norder
            allocate(arrow%Jcart(arrow%nderiv))
            allocate(arrow%map_s2a(arrow%nderiv,2*arrow%norder))
            do i=1,arrow%nderiv
             arrow%Jcart(i)%J=val%Jcart(i)%J
             arrow%Jcart(i)%kind=val%Jcart(i)%kind
            enddo
            arrow%map_s2a=val%map_s2a
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Jtnode

        subroutine rd_Jtnode(this,J)
        implicit none
        class(Jt_list)     :: this
        class(Jthermos)    :: J
        integer            :: i

         select type (bho=>this%node%key)
          type is (Jthermos)
           J%nderiv=bho%nderiv
           J%norder=bho%norder
           allocate(J%Jcart(J%nderiv))
           allocate(J%map_s2a(J%nderiv,J%norder*2))
           do i=1,J%nderiv
            J%Jcart(i)%J=bho%Jcart(i)%J
            J%Jcart(i)%kind=bho%Jcart(i)%kind
           enddo
           J%map_s2a=bho%map_s2a
           J%kind=bho%kind
         end select

        return
        end subroutine rd_Jtnode

        subroutine add_Gtnode(this_list,val)
        implicit none
        class(Gt_list)                 :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: i

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (Gthermos)
            allocate(Gthermos::this_list%head%key)
            allocate(Gthermos::this_list%node%key)
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
           type is (Gthermos)
            allocate(Gthermos::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (Gthermos)
           select type (arrow=>tmp_node%key)
            type is (Gthermos)
            arrow%nderiv=val%nderiv
            arrow%norder=val%norder
            allocate(arrow%Gcart(arrow%nderiv))
            allocate(arrow%map_s2a(arrow%nderiv,2*arrow%norder))
            do i=1,arrow%nderiv
             arrow%Gcart(i)%G=val%Gcart(i)%G
             arrow%Gcart(i)%kind=val%Gcart(i)%kind
            enddo
            arrow%map_s2a=val%map_s2a
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Gtnode

        subroutine rd_Gtnode(this,G)
        implicit none
        class(Gt_list)     :: this
        class(Gthermos)    :: G
        integer            :: i

         select type (bho=>this%node%key)
          type is (Gthermos)
           G%nderiv=bho%nderiv
           G%norder=bho%norder
           allocate(G%Gcart(G%nderiv))
           allocate(G%map_s2a(G%nderiv,G%norder*2))
           do i=1,G%nderiv
            G%Gcart(i)%G=bho%Gcart(i)%G
            G%Gcart(i)%kind=bho%Gcart(i)%kind
           enddo
           G%map_s2a=bho%map_s2a
           G%kind=bho%kind
         end select

        return
        end subroutine rd_Gtnode

        subroutine add_DSItnode(this_list,val)
        implicit none
        class(DSIt_list)               :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: i

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (DSIthermos)
            allocate(DSIthermos::this_list%head%key)
            allocate(DSIthermos::this_list%node%key)
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
           type is (DSIthermos)
            allocate(DSIthermos::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (DSIthermos)
           select type (arrow=>tmp_node%key)
            type is (DSIthermos)
            arrow%nderiv=val%nderiv
            arrow%norder=val%norder
            allocate(arrow%Dcart(arrow%nderiv))
            allocate(arrow%map_s2a(arrow%nderiv,2*arrow%norder))
            do i=1,arrow%nderiv
             arrow%Dcart(i)%D=val%Dcart(i)%D
             arrow%Dcart(i)%kind=val%Dcart(i)%kind
            enddo
            arrow%map_s2a=val%map_s2a
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_DSItnode

        subroutine rd_DSItnode(this,D)
        implicit none
        class(DSIt_list)   :: this
        class(DSIthermos)  :: D
        integer            :: i

         select type (bho=>this%node%key)
          type is (DSIthermos)
           D%nderiv=bho%nderiv
           D%norder=bho%norder
           allocate(D%Dcart(D%nderiv))
           allocate(D%map_s2a(D%nderiv,D%norder*2))
           do i=1,D%nderiv
            D%Dcart(i)%D=bho%Dcart(i)%D
            D%Dcart(i)%kind=bho%Dcart(i)%kind
           enddo
           D%map_s2a=bho%map_s2a
           D%kind=bho%kind
         end select

        return
        end subroutine rd_DSItnode

        subroutine add_D2Stnode(this_list,val)
        implicit none
        class(D2St_list)               :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: i

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (D2Sthermos)
            allocate(D2Sthermos::this_list%head%key)
            allocate(D2Sthermos::this_list%node%key)
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
           type is (D2Sthermos)
            allocate(D2Sthermos::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (D2Sthermos)
           select type (arrow=>tmp_node%key)
            type is (D2Sthermos)
            arrow%nderiv=val%nderiv
            arrow%norder=val%norder
            allocate(arrow%Dcart(arrow%nderiv))
            allocate(arrow%map_s2a(arrow%nderiv,2*arrow%norder))
            do i=1,arrow%nderiv
             arrow%Dcart(i)%D=val%Dcart(i)%D
             arrow%Dcart(i)%kind=val%Dcart(i)%kind
            enddo
            arrow%map_s2a=val%map_s2a
            arrow%kind=val%kind
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_D2Stnode

        subroutine rd_D2Stnode(this,D)
        implicit none
        class(D2St_list)   :: this
        class(D2Sthermos)  :: D
        integer            :: i

         select type (bho=>this%node%key)
          type is (D2Sthermos)
           D%nderiv=bho%nderiv
           D%norder=bho%norder
           allocate(D%Dcart(D%nderiv))
           allocate(D%map_s2a(D%nderiv,D%norder*2))
           do i=1,D%nderiv
            D%Dcart(i)%D=bho%Dcart(i)%D
            D%Dcart(i)%kind=bho%Dcart(i)%kind
           enddo
           D%map_s2a=bho%map_s2a
           D%kind=bho%kind
         end select

        return
        end subroutine rd_D2Stnode

        subroutine add_Otnode(this_list,val)
        implicit none
        class(Ot_list)                 :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node
        integer                        :: i

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (Othermos)
            allocate(Othermos::this_list%head%key)
            allocate(Othermos::this_list%node%key)
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
           type is (Othermos)
            allocate(Othermos::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif
        
         if ( present(val) ) then
          select type (val)
                  
          type is (Othermos)
           select type (arrow=>tmp_node%key)
            type is (Othermos)
            arrow%nderiv=val%nderiv
            arrow%norder=val%norder
            allocate(arrow%Ocart(arrow%nderiv))
            allocate(arrow%map_s2a(arrow%nderiv,arrow%norder*2))
            do i=1,arrow%nderiv
             arrow%Ocart(i)%k=val%Ocart(i)%k
             allocate(arrow%Ocart(i)%B(2*arrow%Ocart(i)%k+1))
             allocate(arrow%Ocart(i)%q(2*arrow%Ocart(i)%k+1))
             arrow%Ocart(i)%B=val%Ocart(i)%B
             arrow%Ocart(i)%q=val%Ocart(i)%q
             arrow%Ocart(i)%kind=val%Ocart(i)%kind
            enddo
            arrow%map_s2a=val%map_s2a
            arrow%kind=val%kind
            arrow%k=val%k
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_Otnode

        subroutine rd_Otnode(this,O)
        implicit none
        class(Ot_list)   :: this
        class(Othermos)  :: O
        integer          :: i

         select type (bho=>this%node%key)
          type is (Othermos)
           O%nderiv=bho%nderiv
           O%norder=bho%norder
           allocate(O%Ocart(O%nderiv))
           allocate(O%map_s2a(O%nderiv,O%norder*2))
           do i=1,O%nderiv
            O%Ocart(i)%k=bho%Ocart(i)%k
            allocate(O%Ocart(i)%B(2*O%Ocart(i)%k+1))
            allocate(O%Ocart(i)%q(2*O%Ocart(i)%k+1))
            O%Ocart(i)%B=bho%Ocart(i)%B
            O%Ocart(i)%q=bho%Ocart(i)%q
            O%Ocart(i)%kind=bho%Ocart(i)%kind
           enddo
           O%map_s2a=bho%map_s2a
           O%kind=bho%kind
           O%k=bho%k
         end select

        return
        end subroutine rd_Otnode


        end module spin_phonon_class
