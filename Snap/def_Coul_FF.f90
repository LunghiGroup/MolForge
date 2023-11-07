        module Coul_FF_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use parameters_class
        use potential_class
        implicit none

        type,extends(potential)                      :: COUL_FF
        real(kind=dbl), allocatable                  :: beta_dip(:)
        real(kind=dbl), allocatable                  :: beta_forces(:)
        integer                                      :: tot_kinds
        integer                                      :: num_bisp
        real(kind=dbl)                               :: r_screen

        contains

        procedure                        :: import_coeff => import_Coul_coeff
        procedure                        :: get_fval  => get_Coul_energy
        procedure                        :: get_fgrad => get_Coul_force
        procedure                        :: import_coeff_forces
        procedure                        :: set_charges

        end type COUL_FF

        contains
        
        subroutine import_Coul_coeff(this)
        implicit none
        class(COUL_FF)                                :: this
        integer                                       :: i

        allocate(this%beta_dip(this%num_bisp*this%tot_kinds))
        open(222,file='snapcoeff_dipoles',action='read')
         do i=1,this%num_bisp*this%tot_kinds
          read(222,*) this%beta_dip(i)
         end do
        close(222)

        end subroutine import_Coul_coeff        

        subroutine import_coeff_forces(this)
        implicit none
        class(COUL_FF)                                :: this
        integer                                       :: i,k

        allocate(this%beta_forces((this%num_bisp-1)*this%tot_kinds))
        open(222,file='snapcoeff_dipoles',action='read')
        do i=1,this%tot_kinds
         read(222,*)
         do k=1,this%num_bisp-1
          read(222,*) this%beta_forces((i-1)*(this%num_bisp-1)+k)
         end do
        end do
        close(222)

        end subroutine import_coeff_forces
        
        subroutine set_charges(this)
        implicit none
        class(COUL_FF)                                :: this
        integer                                       :: j,m
        
        call this%frame%initialize()
        call this%frame%setup(this%frame%nkinds)
        call this%frame%get_desc("DIPOLE")
        call this%frame%finalize()

        allocate(this%frame%charge(this%frame%nats))
        this%frame%charge=0.0
        

        do j=1,this%frame%nats
         do m=1,this%num_bisp
        
          this%frame%charge(j)=this%frame%charge(j)+this%frame%at_desc_dip(j)%desc(m)*this%beta_dip&
                  ((this%frame%kind(j)-1)*this%num_bisp+m)

         end do
        end do
        
        do j=1,this%frame%nats
         deallocate(this%frame%at_desc_dip(j)%desc)
        end do
        
        deallocate(this%frame%at_desc_dip)
        end subroutine set_charges

        subroutine get_Coul_energy(this,vec,val)
        implicit none
        class(COUL_FF)                                          :: this
        real(kind=dbl)                                          :: val
        real(kind=dbl), allocatable                             :: vec(:)        
        integer                                                 :: i,j,k
        real(kind=dbl)                                          :: screen
        real(kind=dbl),dimension(:,:),allocatable               :: rel_dist

        ! this subroutine assumes that the charges have already been assigned with the set_charges subroutine
        val=0.0

        allocate(rel_dist(this%frame%nats,this%frame%nats))
        rel_dist=0.0
        call this%frame%dist_ij

        do j=1,this%frame%nats-1
         do k=j+1,this%frame%nats
        
           rel_dist(j,k)=A_to_B*this%frame%dist(j,1,k,1)
           if (rel_dist(j,k) <= this%r_screen) then
            screen=0.5*(1-dcos(PI*rel_dist(j,k)/this%r_screen))
           else
            screen=1.0
           end if

          val=val+screen*(this%frame%charge(j)*this%frame%charge(k)/rel_dist(j,k))
        end do
       end do
       deallocate(rel_dist)


       val=val*Har_to_Kc

       end subroutine get_Coul_energy

       subroutine get_Coul_force(this,vec,val,grad)
       implicit none
       class(COUL_FF)                                  :: this
       real(kind=dbl)                                  :: val
       real(kind=dbl), allocatable                     :: vec(:),grad(:)
       real(kind=dbl)                                  :: dump
       integer                                         :: i,j,k,l,m
       character(len=150)                              :: atom_string
       real(kind=dbl),allocatable                      :: rel_dist(:,:),coord(:,:)
       
       !beware that this subroutine needs to be called after import_coeff_forces
       !at every stage of the subroutine we are acutally calculating the forces
       !just at the end we will perform a change in sign to obtain the gradients

       call this%frame%initialize()
       call this%frame%setup(this%frame%nats)
       do i=1,this%frame%nats
        write(atom_string,*) 'set atom',i,'type',i
        call lammps_command(this%frame%lmp,trim(atom_string))
       end do
       
       call this%frame%get_der_desc("DIPOLE",this%frame%nats)
       call this%frame%finalize()
       
       if (.not.allocated(rel_dist)) allocate(rel_dist(this%frame%nats,this%frame%nats))
       if (.not.allocated(coord)) allocate(coord(this%frame%nats,3))
       if (.not.allocated(grad)) allocate(grad(this%frame%nats*3))
       grad=0.0

       call this%frame%dist_ij

       coord=A_to_B*this%frame%x
        
       do i=1,this%frame%nats
        do j=1,this%frame%nats
         if (i.ne.j) then 
          do m=1,3
         
           rel_dist(i,j)=A_to_B*this%frame%dist(i,1,j,1)
          
           if (rel_dist(i,j) <= this%r_screen) then
            dump=0.5*(1-dcos(PI*rel_dist(i,j)/this%r_screen))
           else
            dump=1.0
           end if
           !
           !here we implement a simple Coulomb interaction with the additional dumping. 

           grad((i-1)*3+m)=grad((i-1)*3+m)+dump*(this%frame%charge(i)*this%frame%charge(j)&
           *(coord(i,m)-coord(j,m)))/(rel_dist(i,j)**3)
           
           !
           !derivation of the screening part of the fictitious Coulomb interaction
          
           if (rel_dist(i,j) <= this%r_screen) then

            grad((i-1)*3+m)=grad((i-1)*3+m)-0.5*this%frame%charge(i)*this%frame%charge(j)*&
            (dsin(PI*rel_dist(i,j)/this%r_screen))*PI*(coord(i,m)-coord(j,m))/(this%r_screen*(rel_dist(i,j)**2))

           end if
          
          end do                
         end if
        end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=1,this%frame%nats
        this%frame%der_at_desc_dip(i)%desc=this%frame%der_at_desc_dip(i)%desc/A_to_B
        do l=1,this%frame%nats
         do j=1,this%frame%nats
         
          if (l.ne.j) then

           rel_dist(l,j)=A_to_B*this%frame%dist(l,1,j,1)
           if (rel_dist(l,j) <= this%r_screen) then
            dump=0.5*(1-dcos(PI*rel_dist(l,j)/this%r_screen))
           else
            dump=1.0
           end if
           
           do m=1,3
            do k=1,this%num_bisp-1
             grad((i-1)*3+m)=grad((i-1)*3+m)+dump*(this%frame%charge(l)/rel_dist(l,j))*&
             this%frame%der_at_desc_dip(i)%desc((j-1)*3*(this%num_bisp-1)+(m-1)*(this%num_bisp-1)+k)&
             *this%beta_forces((this%frame%kind(j)-1)*(this%num_bisp-1)+k)
                
            end do
           end do
          
          end if       
         
         end do
        end do
       end do
       
       grad=-grad*F_conv

       !!!!!!deallocation part of the code 
       do i=1,this%frame%nats
        deallocate(this%frame%der_at_desc_dip(i)%desc)
       end do
       deallocate(this%frame%der_at_desc_dip)
       !!!!!!

       end subroutine get_Coul_force

       end module Coul_FF_class
