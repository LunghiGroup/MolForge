        module proj_disp_class
        implicit none


         type molecule        
          double precision, allocatable    :: cart(:,:)
          double precision, allocatable    :: cart_eq(:,:)
          double precision, allocatable    :: cart_rot(:,:)
          double precision, allocatable    :: cart_int(:,:) 
          double precision                 :: tr(3)
          double precision, allocatable    :: mass(:) 
          double precision                 :: Mtot
          double precision                 :: lamb(3,3)
          double precision                 :: rot(3,3)
          double precision                 :: theta_cov(3)
          double precision                 :: theta_con(3)
          double precision                 :: psi
          double precision                 :: cmass(3)
          double precision                 :: cmass_eq(3)
          double precision                 :: Mpsi(3,3)
          double precision                 :: Mpsi2(3,3)
          double precision, allocatable    :: jaco(:,:,:)
          integer                          :: nat
          double precision                 :: inertia(3,3)
          double precision                 :: inertia_inv(3,3)
          double precision                 :: M(3,3,3)
          contains          
          procedure                        :: def_mol
          procedure                        :: def_mol_dist
          procedure                        :: proj_disp
          procedure                        :: get_trasl
          procedure                        :: get_M
          procedure                        :: get_inertia
          procedure                        :: get_rot
          procedure                        :: get_internals
          procedure                        :: get_jaco
          procedure                        :: get_norms
         end type molecule

        contains

         subroutine def_mol(this,geo,mass)
         implicit none
         class(molecule)                   :: this
         double precision, allocatable     :: geo(:,:),mass(:)

          this%nat=size(mass)

          allocate(this%mass(this%nat))
          this%mass=mass
          
          allocate(this%cart_eq(this%nat,3))
          this%cart_eq=geo

         return
         end subroutine def_mol


         subroutine def_mol_dist(this,geo)
         implicit none
         class(molecule)                   :: this
         double precision, allocatable     :: geo(:,:)
         
          if(.not.allocated(this%cart)) allocate(this%cart(this%nat,3))
          this%cart=geo

         return
         end subroutine def_mol_dist


         subroutine proj_disp(this)
         implicit none
         class(molecule)                   :: this
         integer                           :: i

          if(.not.allocated(this%cart_rot)) allocate(this%cart_rot(this%nat,3))
          if(.not.allocated(this%cart_int)) allocate(this%cart_int(this%nat,3))

          call this%get_trasl          ! calcola lo spostamente del baricentro t
          call this%get_M              ! definisce le matrici di rotazione infinitesima
          call this%get_inertia        ! calcola il tensore di inerazia e la sua inversa
          call this%get_rot            ! calcola la matrice rot relativa alle sole rotazioni
          call this%get_internals      ! calcola i soli spostamente interni come int=tot-rot-trasl

         return
         end subroutine proj_disp

         subroutine get_M(this)
         implicit none
         class(molecule)                   :: this

          this%M(1,1,1)=0.0d0
          this%M(1,1,2)=0.0d0
          this%M(1,1,3)=0.0d0
          this%M(1,2,1)=0.0d0
          this%M(1,2,2)=0.0d0
          this%M(1,2,3)=-1.0d0
          this%M(1,3,1)=0.0d0
          this%M(1,3,2)=1.0d0
          this%M(1,3,3)=0.0d0

          this%M(2,1,1)=0.0d0
          this%M(2,1,2)=0.0d0
          this%M(2,1,3)=1.0d0
          this%M(2,2,1)=0.0d0
          this%M(2,2,2)=0.0d0
          this%M(2,2,3)=0.0d0
          this%M(2,3,1)=-1.0d0
          this%M(2,3,2)=0.0d0
          this%M(2,3,3)=0.0d0

          this%M(3,1,1)=0.0d0
          this%M(3,1,2)=-1.0d0
          this%M(3,1,3)=0.0d0
          this%M(3,2,1)=1.0d0
          this%M(3,2,2)=0.0d0
          this%M(3,2,3)=0.0d0
          this%M(3,3,1)=0.0d0
          this%M(3,3,2)=0.0d0
          this%M(3,3,3)=0.0d0

         return
         end subroutine get_M

         subroutine get_inertia(this)
         use lapack_inverse
         implicit none
         class(molecule)                   :: this
         integer          :: i

          this%inertia=0.0d0
          this%inertia_inv=0.0d0

          do i=1,this%nat
           this%inertia(1,1)=this%inertia(1,1)+this%mass(i)*&
                        ((this%cart_eq(i,2)-this%cmass_eq(2))**2+(this%cart_eq(i,3)-this%cmass_eq(3))**2)
           this%inertia(1,2)=this%inertia(1,2)-this%mass(i)*&
                        (this%cart_eq(i,1)-this%cmass_eq(1))*(this%cart_eq(i,2)-this%cmass_eq(2))
           this%inertia(1,3)=this%inertia(1,3)-this%mass(i)*&
                        (this%cart_eq(i,1)-this%cmass_eq(1))*(this%cart_eq(i,3)-this%cmass_eq(3))
           this%inertia(2,2)=this%inertia(2,2)+this%mass(i)*&
                        ((this%cart_eq(i,1)-this%cmass_eq(1))**2+(this%cart_eq(i,3)-this%cmass_eq(3))**2)
           this%inertia(2,1)=this%inertia(2,1)-this%mass(i)*&
                        (this%cart_eq(i,1)-this%cmass_eq(1))*(this%cart_eq(i,2)-this%cmass_eq(2))
           this%inertia(2,3)=this%inertia(2,3)-this%mass(i)*&
                        (this%cart_eq(i,2)-this%cmass_eq(2))*(this%cart_eq(i,3)-this%cmass_eq(3))
           this%inertia(3,3)=this%inertia(3,3)+this%mass(i)*&
                        ((this%cart_eq(i,1)-this%cmass_eq(1))**2+(this%cart_eq(i,2)-this%cmass_eq(2))**2)
           this%inertia(3,1)=this%inertia(3,1)-this%mass(i)*&
                        (this%cart_eq(i,3)-this%cmass_eq(3))*(this%cart_eq(i,1)-this%cmass_eq(1))
           this%inertia(3,2)=this%inertia(3,2)-this%mass(i)*&
                        (this%cart_eq(i,3)-this%cmass_eq(3))*(this%cart_eq(i,2)-this%cmass_eq(2))
          enddo

          this%inertia_inv(1:3,1:3)=this%inertia(1:3,1:3)

          call mat_inv(this%inertia_inv,3)

         return
         end subroutine get_inertia


         subroutine get_trasl(this)
         implicit none
         class(molecule)                   :: this
         integer          :: i,s

          this%tr=0.0d0
          this%Mtot=0.0d0
          this%cmass=0.0d0
          this%cmass_eq=0.0d0

          do i=1,this%nat
           this%Mtot=this%Mtot+this%mass(i)
          enddo

          do s=1,3
           do i=1,this%nat
            this%cmass(s)=this%cmass(s)+this%mass(i)*this%cart(i,s)
           enddo
          enddo

          do s=1,3
           do i=1,this%nat
            this%cmass_eq(s)=this%cmass_eq(s)+this%mass(i)*this%cart_eq(i,s)
           enddo
          enddo

          do s=1,3
           do i=1,this%nat
            this%tr(s)=this%tr(s)+this%mass(i)*(this%cart(i,s)-this%cart_eq(i,s))
           enddo
          enddo

          this%cmass_eq=this%cmass_eq/this%Mtot
          this%cmass=this%cmass/this%Mtot
          this%tr=this%tr/this%Mtot

         return
         end subroutine get_trasl


         subroutine get_jaco(this)
         implicit none
         class(molecule)                   :: this
         integer         :: s,i,t,v,ialloc

          if(.not. allocated(this%jaco)) allocate(this%jaco(3,this%nat,3))
          this%jaco=0.0d0

          do s=1,3
           do i=1,this%nat
            do t=1,3
             do v=1,3
              this%jaco(s,i,t)=this%jaco(s,i,t)+this%M(s,t,v)*(this%cart_eq(i,v)-this%cmass_eq(v))
             enddo
            enddo
           enddo
          enddo

         return
         end subroutine get_jaco

         subroutine get_rot(this)
         implicit none
         class(molecule)                   :: this
         integer                           :: s,t,v,i,iter
         double precision                  :: aa(this%nat,3),psiv(3)

          this%psi=0.0d0
          this%rot=0.0d0
          this%lamb=0.0d0
          this%Mpsi=0.0d0
          this%Mpsi2=0.0d0
          iter=0

          do s=1,3
           this%lamb(s,s)=1.0d0
           this%rot(s,s)=1.0d0
          enddo

          call this%get_jaco   ! calcola la jacobiana covariante cart -> rot_cov

         ! calcolo iterativo della matrice di rotazione relativa allo spostamento cartesiano

          do while ( (abs(this%theta_cov(1)).gt.1.0d-8) .or. (abs(this%theta_cov(2)).gt.1.0d-8) .or. &
                     (abs(this%theta_cov(3)).gt.1.0d-8)  .or. (iter.lt.5) )

          aa=0.0d0
          this%theta_cov=0.0d0
          this%psi=0.0d0
          psiv=0.0d0
          this%lamb=0.0d0
          this%Mpsi=0.0d0
          this%Mpsi2=0.0d0

          do i=1,this%nat
           do s=1,3
            do v=1,3
             aa(i,s)=aa(i,s)+this%rot(s,v)*(this%cart_eq(i,v)-this%cmass_eq(v))
            enddo
           enddo
          enddo

          do s=1,3
           do i=1,this%nat
            do t=1,3
              this%theta_cov(s)=this%theta_cov(s)+this%mass(i)*this%jaco(s,i,t)*(this%cart(i,t)-aa(i,t)-this%cmass(t))
            enddo
           enddo
          enddo

          this%theta_con=MATMUL(this%inertia_inv,this%theta_cov)

          do s=1,3
           this%psi=this%psi+this%theta_con(s)*this%theta_con(s)
          enddo

          this%psi=dsqrt(this%psi)

          do i=1,3
           psiv(i)=this%theta_con(i)/this%psi
          enddo

          this%lamb(1,1)=psiv(1)*psiv(1)*(1-dcos(this%psi))+dcos(this%psi)
          this%lamb(1,2)=psiv(1)*psiv(2)*(1-dcos(this%psi))-psiv(3)*dsin(this%psi)
          this%lamb(1,3)=psiv(1)*psiv(3)*(1-dcos(this%psi))+psiv(2)*dsin(this%psi)
          this%lamb(2,1)=psiv(2)*psiv(1)*(1-dcos(this%psi))+psiv(3)*dsin(this%psi)
          this%lamb(2,2)=psiv(2)*psiv(2)*(1-dcos(this%psi))+dcos(this%psi)
          this%lamb(2,3)=psiv(2)*psiv(3)*(1-dcos(this%psi))-psiv(1)*dsin(this%psi)
          this%lamb(3,1)=psiv(3)*psiv(1)*(1-dcos(this%psi))-psiv(2)*dsin(this%psi)
          this%lamb(3,2)=psiv(3)*psiv(2)*(1-dcos(this%psi))+psiv(1)*dsin(this%psi)
          this%lamb(3,3)=psiv(3)*psiv(3)*(1-dcos(this%psi))+dcos(this%psi)

          this%rot=MATMUL(this%lamb,this%rot)
          iter=iter+1

          enddo

          this%cart_rot=0.0d0

          do i=1,this%nat
           do s=1,3
            do t=1,3
             this%cart_rot(i,s)=this%cart_rot(i,s)+this%rot(s,t)*(this%cart_eq(i,t)-this%cmass_eq(t))
            enddo
           enddo
!           this%cart_rot(i,1:3)=this%cart_rot(i,1:3)+this%cmass_eq(1:3)
          enddo

         return
         end subroutine get_rot


         subroutine get_internals(this)
         implicit none
         class(molecule)                   :: this
         integer        :: i,t,v,ialloc

          if(.not.allocated(this%cart_int)) allocate(this%cart_int(this%nat,3))
          this%cart_int=0.0d0

          do i=1,this%nat
           do t=1,3
            do v=1,3
             this%cart_int(i,t)=this%cart_int(i,t)+this%rot(v,t)*(this%cart(i,v)-this%cmass(v))
            enddo
           enddo
          enddo

          do i=1,this%nat
           do t=1,3
            this%cart_int(i,t)=this%cart_int(i,t)-this%cart_eq(i,t)+this%cmass_eq(t)
           enddo
          enddo

!          do i=1,this%nat
!           do t=1,3
!            this%cart_int(i,t)=this%cart(i,t)-this%cart_eq(i,t)-this%tr(t)
!           enddo
!           this%cart_int(i,:)=matmul(transpose(this%rot),this%cart_int(i,:))
!          enddo

         return
         end subroutine get_internals

         subroutine get_norms(this,norm)
         implicit none
         class(molecule)                   :: this
         integer                           :: i,s,t
         double precision                  :: norm(5)
         double precision, allocatable     :: inter_rot(:,:)

          norm=0.0d0

          allocate(inter_rot(this%nat,3))
          inter_rot=0.0d0

          do i=1,this%nat
           do s=1,3
            do t=1,3
             inter_rot(i,s)=inter_rot(i,s)+this%rot(s,t)*this%cart_int(i,t)
            enddo
           enddo
          enddo


          do i=1,this%nat
           do s=1,3
            norm(1)=norm(1)+(this%cart(i,s)-this%cart_eq(i,s))**2
            norm(2)=norm(2)+(this%tr(s))**2
            norm(3)=norm(3)+(this%cart_rot(i,s)-this%cart_eq(i,s)+this%cmass_eq(s))**2
            norm(4)=norm(4)+(inter_rot(i,s))**2
            norm(5)=norm(5)+(this%cart(i,s)-this%cart_rot(i,s) &
                             -inter_rot(i,s)-this%cmass(s))**2
           enddo
          enddo

!          norm=sqrt(norm)         

         return
         end subroutine get_norms

         
         end module proj_disp_class
