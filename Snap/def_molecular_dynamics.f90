        module md_class
        use SNAP_FF_class
        use lammps_class
        use potential_class
        use SNAP_FF_class
        use VdW_class
        use target_functions_class
        implicit none

        type                                    :: MD
         
        type(potential_list),allocatable        :: pot(:)
        integer                                 :: max_steps
        integer                                 :: step_size
        integer                                 :: num_pot
        integer                                 :: iter
        character(len=3)                        :: ensemble="nvt"
        real(kind=dbl)                          :: T_bath
        contains
        
        procedure                               :: initialize_vel        
        procedure                               :: link_potentials
        procedure                               :: propagate
        procedure                               :: velocity_verlet
        procedure                               :: bussi
        procedure                               :: print_info

        end type MD

        contains
        
        subroutine initialize_vel(this,frame,iseed,T_in)
        implicit none
        class(md),intent(inout)                    :: this
        real(kind=dbl), intent(in)                 :: T_in
        type(lammps_obj)                           :: frame
        integer                                    :: idist=3
        integer                                    :: N=3
        integer                                    :: i
        integer,dimension(4)                       :: iseed
        real(kind=dbl),dimension(3)                :: X
        real(kind=dbl),dimension(3)                :: lin_mom_tot
        real(kind=dbl)                             :: E_kin
        real(kind=dbl)                             :: T

        allocate(frame%v(frame%nats,3))
        
        do i=1,frame%nats
         call dlarnv(idist,iseed,N,X)
         frame%v(i,:)=dsqrt(boltz*T_in/(frame%mass(frame%kind(i))))*X
        end do

        E_kin=0.0d0
        lin_mom_tot=0.0d0

        do i=1,frame%nats
         lin_mom_tot=lin_mom_tot + frame%v(i,:)*frame%mass(frame%kind(i))
        end do

        do i=1,frame%nats
         frame%v(i,:)=frame%v(i,:)-(lin_mom_tot/(dble(frame%nats)*frame%mass(frame%kind(i))))
        end do

        do i=1,frame%nats
         E_kin=E_kin+0.5d0*sum((frame%v(i,:)**2))*(frame%mass(frame%kind(i)))
        end do

        T=(2.0d0*E_kin)/((3.0d0*frame%nats-3.0d0)*boltz)
        frame%v(:,:)=dsqrt(T_in/T)*frame%v(:,:)

        end subroutine initialize_vel

        subroutine link_potentials(this,keyword,SNAP,VdW)
        implicit none
        class(md)                       :: this
        character(len=*)                :: keyword
        type(SNAP_FF),target            :: SNAP
        type(VdW_FF),optional,target    :: VdW
        integer                         :: i

        allocate(this%pot(this%num_pot))
        
        if (keyword=="SNAP") then
         this%pot(1)%item=>SNAP
         do i=2,this%num_pot
          this%pot(i)%item=>NULL()
         end do      
        else if ((keyword=="SNAP_VdW").and.(present(VdW))) then 
         this%pot(1)%item=>SNAP
         this%pot(2)%item=>VdW
         do i=3,this%num_pot
          this%pot(i)%item=>NULL()
         end do
        
        end if
        end subroutine link_potentials

        subroutine propagate(this,frame)
        implicit none
        class(md), intent(inout) :: this
        type(lammps_obj)         :: frame
        integer                  :: i 
        
        this%iter=0
        call this%print_info(frame)
        do while (this%iter < this%max_steps)
         this%iter=this%iter+1
         call this%velocity_verlet(frame)
         if ((this%ensemble)=='nvt') then
         call this%bussi(frame)
         end if
         call this%print_info(frame)        
        end do

        end subroutine propagate

        subroutine velocity_verlet(this,frame)
        implicit none
        class(md), intent(inout)        :: this
        type(lammps_obj)                :: frame
        integer                         :: i,j
        real(kind=dbl),allocatable      :: sum_grad(:)
        real(kind=dbl),allocatable      :: vec(:)
        real(kind=dbl),allocatable      :: grad(:)
        real(kind=dbl),allocatable      :: acc(:,:)
        real(kind=dbl)                  :: val
        
        allocate(sum_grad(frame%nats*3))
        sum_grad=0.0d0

        if (this%iter==1) then 
         allocate(frame%acc(frame%nats,3))
         do i=1,this%num_pot
          call this%pot(i)%item%get_fgrad(vec,val,grad)
          sum_grad=sum_grad+grad
          deallocate(grad)
         end do
         sum_grad=sum_grad/F_conv

         do i=1,frame%nats
          do j=1,3

           frame%acc(i,j)=-sum_grad((i-1)*3+j)/frame%mass(frame%kind(i))

          end do
         end do

        end if
                
        frame%v=frame%v+0.5d0*frame%acc*this%step_size
        frame%x=(A_to_B)*frame%x+frame%v*this%step_size
        frame%x=frame%x/(A_to_B)
        sum_grad=0.0d0
         
        do i=1,this%num_pot
         this%pot(i)%item%frame=frame
         call this%pot(i)%item%get_fgrad(vec,val,grad)
         sum_grad=sum_grad+grad
         deallocate(grad)
        
        end do
        sum_grad=sum_grad/F_conv

        do i=1,frame%nats
         do j=1,3

          frame%acc(i,j)=-sum_grad((i-1)*3+j)/frame%mass(frame%kind(i))

         end do
        end do

        frame%v=frame%v+0.5d0*frame%acc*this%step_size
        deallocate(sum_grad)
        end subroutine velocity_verlet

        subroutine bussi(this,frame)
        implicit none
        class(md)                :: this
        type(lammps_obj)         :: frame
        real(kind=dbl)           :: E_kin
        real(kind=dbl)           :: E_kin_new
        real(kind=dbl)           :: resamplekin
        integer                  :: i,j

        E_kin=0.0

        do i=1,frame%nats
         do j=1,3

          E_kin=E_kin+0.5d0*(frame%v(i,j)**2)*frame%mass(frame%kind(i))

         end do
        end do

        E_kin_new= resamplekin(E_kin,((3.0d0*frame%nats-3.d0)/2.0d0)*this%T_bath &
        *boltz,(3*frame%nats)-3,100.0d0)

        frame%v=dsqrt(E_kin_new/E_kin)*frame%v

        end subroutine bussi

        subroutine print_info(this,frame)
        implicit none
        class(md)                       :: this
        type(lammps_obj)                :: frame
        real(kind=dbl)                  :: E_kin
        real(kind=dbl)                  :: val
        real(kind=dbl)                  :: sum_val
        real(kind=dbl)                  :: temp
        real(kind=dbl),allocatable      :: vec(:)
        integer                         :: i,j
        
        sum_val=0.0d0
        E_kin=0.0d0
        val=0.0d0

        do i=1,this%num_pot
        
         call this%pot(i)%item%get_fval(vec,val)
         sum_val=sum_val+val

        end do

        do i=1,frame%nats
         do j=1,3

          E_kin=E_kin+0.5d0*(frame%v(i,j)**2)*frame%mass(frame%kind(i))

         end do
        end do
        
        temp=(2.0d0*E_kin)/((3.0d0*frame%nats-3.0d0)*boltz)

        open(111, file="traj_MD_molforge.xyz", action="write",position='append')

        write(111,*)frame%nats
        write(111,*)'XXX'

        do i=1,frame%nats
         write(111,*) frame%label(frame%kind(i)),frame%x(i,:)
        end do
        
        close(111)

        open(111, file="etotal_kin_pot_temp_molforge.txt", action="write",position="append")
         write(111,*) E_kin*Har_to_kc+sum_val, E_kin*Har_to_kc,sum_val,temp
        close(111)
        
        end subroutine print_info

        end module md_class





! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the kinetic energy,
! which can be used to rescale the velocities.
! The procedure can be applied to all atoms or to smaller groups.
! If it is applied to intersecting groups in sequence, the kinetic energy
! that is given as an input (kk) has to be up-to-date with respect to the previous rescalings.
!
! When applied to the entire system, and when performing standard molecular dynamics (fixed c.o.m. (center of mass))
! the degrees of freedom of the c.o.m. have to be discarded in the calculation of ndeg,
! and the c.o.m. momentum HAS TO BE SET TO ZERO.
! When applied to subgroups, one can chose to:
! (a) calculate the subgroup kinetic energy in the usual reference frame, and count the c.o.m. in ndeg
! (b) calculate the subgroup kinetic energy with respect to its c.o.m. motion, discard the c.o.m. in ndeg
!     and apply the rescale factor with respect to the subgroup c.o.m. velocity.
! They should be almost equivalent.
! If the subgroups are expected to move one respect to the other, the choice (b) should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces to an istantaneous
! randomization of the kinetic energy, as described in paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator against discretization errors.
! The easiest recipe is:
! htilde = h + conint
! where h is the total energy (kinetic + potential)
! and conint is a quantity accumulated along the trajectory as minus the sum of all the increments of kinetic
! energy due to the thermostat.
!
!module rescale
!implicit none
!contains

function resamplekin(kk,sigma,ndeg,taut)
  implicit none
  double precision              :: resamplekin
  double precision,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  double precision,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  double precision,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  double precision :: factor,rr
  double precision, external :: gasdev

  if(taut>0.1) then
    factor=exp(-1.0/taut)
  else
    factor=0.0
  end if
  rr = gasdev()
  resamplekin = kk + (1.0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk) &
               + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor)

contains 

double precision function sumnoises(nn)
  implicit none
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
double precision, external :: gamdev,gasdev
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0*gamdev(nn/2)
  else
    sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

end function resamplekin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING ROUTINES ARE TRANSCRIBED FROM NUMERICAL RECIPES

double precision function gamdev(ia)
! gamma-distributed random number, implemented as described in numerical recipes

implicit none
integer, intent(in) :: ia
integer j
double precision am,e,s,v1,v2,x,y
double precision, external :: ran1
if(ia.lt.1)stop 'bad argument in gamdev'
if(ia.lt.6)then
  x=1.
  do 11 j=1,ia
    x=x*ran1()
11  continue
  x=-log(x)
else
1 v1=2.*ran1()-1.
    v2=2.*ran1()-1.
  if(v1**2+v2**2.gt.1.)goto 1
    y=v2/v1
    am=ia-1
    s=sqrt(2.*am+1.)
    x=s*y+am
  if(x.le.0.)goto 1
    e=(1.+y**2)*exp(am*log(x/am)-s*y)
  if(ran1().gt.e)goto 1
endif
gamdev=x
end function gamdev

double precision function gasdev()
! gaussian-distributed random number, implemented as described in numerical recipes

implicit none
integer, save :: iset = 0
double precision, save :: gset
double precision, external :: ran1
double precision fac,rsq,v1,v2
if(iset==0) then
1       v1=2.*ran1()-1.0d0
  v2=2.*ran1()-1.0d0
  rsq=v1**2+v2**2
  if(rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
end if
end function gasdev

FUNCTION ran1()
! random number generator
INTEGER IA,IM,IQ,IR,NTAB,NDIV
double precision ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
INTEGER, SAVE :: idum=0 !! ATTENTION: THE SEED IS HARDCODED
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
  do 11 j=NTAB+8,1,-1
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11  continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
      END function ran1
