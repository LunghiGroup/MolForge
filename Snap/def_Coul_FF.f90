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
        real(kind=dbl), allocatable                  :: beta(:)
        integer                                      :: tot_kinds
        integer                                      :: num_bisp
        real(kind=dbl)                               :: r_screen

        contains

        procedure                        :: import_coeff => import_Coul_coeff
        procedure                        :: get_fval  => get_Coul_energy
        procedure                        :: get_fgrad => get_Coul_force
        procedure                        :: set_charges

        end type COUL_FF

        contains
        
        subroutine import_Coul_coeff(this)
        implicit none
        class(COUL_FF)                                :: this
        integer                                       :: i

        allocate(this%beta(this%num_bisp*this%tot_kinds))
        open(222,file='snapcoeff_dipoles',action='read')
         do i=1,this%num_bisp*this%tot_kinds
          read(222,*) this%beta(i)
         end do
        close(222)

        end subroutine import_Coul_coeff        
        
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

          this%frame%charge(j)=this%frame%charge(j)+this%frame%at_desc_dip(j)%desc(m)*this%beta&
                  ((this%frame%kind(j)-1)*this%num_bisp+m)

         end do
        end do

        !deallocation of the descriptors will not happen at this stage, but after coulomb energy is calculated in get_Coul_energy
        !subroutine

        !!!!!!!!!!!!!!!!!!!!!!Debugging line
        write(*,*) 'Total charge:',sum(this%frame%charge)

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

        do j=1,this%frame%nats
         do k=1,this%frame%nats
          if (k>j) then
        
           rel_dist(j,k)=A_to_B*this%frame%dist(j,1,k,1)
           write(*,*) rel_dist(j,k)
           if (rel_dist(j,k) <= this%r_screen) then
            screen=0.5*(1-dcos(PI*rel_dist(j,k)/this%r_screen))
           else
            screen=1.0
           end if

          val=val+screen*(this%frame%charge(j)*this%frame%charge(k)/rel_dist(j,k))
         end if
        end do
       end do
       deallocate(rel_dist)


       val=val*Har_to_Kc

       end subroutine get_Coul_energy

       subroutine get_Coul_force(this,vec,val,grad)
       implicit none
       class(COUL_FF)                                 :: this
       real(kind=dbl)                                  :: val
       real(kind=dbl), allocatable                     :: vec(:),grad(:)
        
       

       end subroutine get_Coul_force

       end module Coul_FF_class
