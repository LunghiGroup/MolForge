        module control_variables
        use lists_class
        use general_types_class
        use pulses_class
        implicit none


         type(general_pulse),allocatable    :: pulse(:)
         type(general_pulse),allocatable    :: pulse_pi(:)
         type(general_pulse),allocatable    :: pulse_pi2(:)

         character(len=100)            :: operation

         double precision              :: spin_temp
         character(len=100)            :: type_rho0
         double precision, allocatable :: alpha0(:),beta0(:),gamma0(:)

         logical                       :: fulldiag=.true.
         logical                       :: tinv=.true.
         integer                       :: nex_max,ncorr_max,nexclude
         double precision, allocatable :: max_dist(:,:)
         double precision              :: dist_max
         type(sub_space), allocatable  :: active_space(:)
         integer, allocatable          :: ex_list(:,:)
         type(list)                    :: dist_kind

         logical                       :: dump_s=.true.
         character(len=20)             :: rho_restart_file
         integer                       :: dump_freq,s2print=0
         integer, allocatable          :: print_si(:)

         integer                       :: nsteps,step_nmult,start_step=0
         double precision              :: step,step_min,time=0.0d0

         integer                       :: lmax
         double precision              :: lambda
         logical                       :: compress=.false.

         integer                       :: echo_nsteps
         integer                       :: echo_spin
         double precision              :: echo_step

         double precision              :: max_phonon_ener=4000
         double precision              :: euler(3)=0.0d0

        contains

        subroutine bcast_method(spindy)
        use mpi
        use mpi_utils
        use hilbert_dist_class
        implicit none
        class(spins_hilbert)   :: spindy
        integer                :: l,nkinds,ncas,i,j,nelem
        double precision       :: val
        logical                :: do_cas=.false.
         
         call mpi_bcast(spindy%make_Rmat,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(spindy%make_R2mat,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(spindy%make_SA,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(spindy%make_PT2,1,mpi_logical,0,mpi_comm_world,err)
         call mpi_bcast(fulldiag,1,mpi_logical,0,mpi_comm_world,err) 

         if(fulldiag)then
          spindy%make_Heig=.true.
         else
          spindy%make_Heig=.false.
         endif

         call mpi_bcast(tinv,1,mpi_logical,0,mpi_comm_world,err) 
         call mpi_bcast(dump_s,1,mpi_logical,0,mpi_comm_world,err) 
         call mpi_bcast(nex_max,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(ncorr_max,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(dump_freq,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(s2print,1,mpi_integer,0,mpi_comm_world,err)
         call mpi_bcast(dist_max,1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(max_phonon_ener,1,mpi_double_precision,0,mpi_comm_world,err)

         call mpi_bcast(euler,3,mpi_double_precision,0,mpi_comm_world,err)

         call mpi_bcast(rho_restart_file,20,mpi_character,0,mpi_comm_world,err)

         if(.not.allocated(max_dist)) allocate(max_dist(spindy%nkinds,spindy%nkinds))
         max_dist=dist_max
         
         if(mpi_id.eq.0) nelem=dist_kind%nelem
         call mpi_bcast(nelem,1,mpi_integer,0,mpi_comm_world,err)
         if(mpi_id.eq.0) call dist_kind%reboot()
         do l=1,nelem 
          if(mpi_id.eq.0) call dist_kind%rd_val(i)
          if(mpi_id.eq.0) call dist_kind%skip()
          if(mpi_id.eq.0) call dist_kind%rd_val(j)
          if(mpi_id.eq.0) call dist_kind%skip()
          if(mpi_id.eq.0) call dist_kind%rd_val(val)
          if(mpi_id.eq.0) call dist_kind%skip()
          call mpi_bcast(i,1,mpi_integer,0,mpi_comm_world,err)
          call mpi_bcast(j,1,mpi_integer,0,mpi_comm_world,err)
          call mpi_bcast(val,1,mpi_double_precision,0,mpi_comm_world,err)
          max_dist(j,i)=val
          max_dist(i,j)=val
         enddo
         if(mpi_id.eq.0) call dist_kind%delete()

         call mpi_bcast(nexclude,1,mpi_integer,0,mpi_comm_world,err)
         if(nexclude.ne.0)then
          if(.not.allocated(ex_list)) allocate(ex_list(nexclude,2))
          do l=1,nexclude
           call mpi_bcast(ex_list(l,:),2,mpi_integer,0,mpi_comm_world,err)
          enddo
         endif
         
         if(.not.allocated(print_si)) allocate(print_si(s2print))
         call mpi_bcast(print_si,s2print,mpi_integer,0,mpi_comm_world,err)

         spindy%s2print=s2print
         allocate(spindy%print_si(spindy%s2print))
         spindy%print_si=print_si

         if(mpi_id.eq.0 .and. allocated(active_space)) do_cas=.true.
         call mpi_bcast(do_cas,1,mpi_logical,0,mpi_comm_world,err)
         if(do_cas)then
          if(allocated(active_space)) ncas=size(active_space)
          call mpi_bcast(ncas,1,mpi_integer,0,mpi_comm_world,err)
          allocate(spindy%active_space(ncas))
          do l=1,ncas
           if(allocated(active_space(l)%kind)) nkinds=size(active_space(l)%kind)
           call mpi_bcast(nkinds,1,mpi_integer,0,mpi_comm_world,err)
           allocate(spindy%active_space(l)%kind(nkinds))
           if(allocated(active_space)) then
            spindy%active_space(l)%kind=active_space(l)%kind
            spindy%active_space(l)%max_ex=active_space(l)%max_ex
           endif
           call mpi_bcast(spindy%active_space(l)%kind,nkinds,mpi_integer,0,mpi_comm_world,err)
           call mpi_bcast(spindy%active_space(l)%max_ex,1,mpi_integer,0,mpi_comm_world,err)
          enddo        
         endif

        return
        end subroutine bcast_method


        subroutine bcast_pulse
        use mpi
        use mpi_utils
        use pulses_class
        implicit none
        integer                :: i,npulses

         if(allocated(pulse)) npulses=size(pulse)
         call mpi_bcast(npulses,1,mpi_integer,0,mpi_comm_world,err)
         if(.not. allocated(pulse)) allocate(pulse(npulses))
         do i=1,size(pulse)
          call pulse(i)%bcast()
         enddo

        return
        end subroutine bcast_pulse


        subroutine bcast_echo
        use mpi
        use mpi_utils
        implicit none
        integer                :: nspins

         call mpi_bcast(echo_nsteps,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(echo_spin,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(echo_step,1,mpi_double_precision,0,mpi_comm_world,err)

        return
        end subroutine bcast_echo


        subroutine bcast_rho0
        use mpi
        use mpi_utils
        implicit none

         call mpi_bcast(type_rho0,100,mpi_character,0,mpi_comm_world,err) 
         call mpi_bcast(spin_temp,1,mpi_double_precision,0,mpi_comm_world,err)

        return
        end subroutine bcast_rho0


        subroutine bcast_propagator
        use mpi
        use mpi_utils
        implicit none

         call mpi_bcast(step_nmult,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(step,1,mpi_double_precision,0,mpi_comm_world,err)
         call mpi_bcast(step_min,1,mpi_double_precision,0,mpi_comm_world,err)

        return
        end subroutine bcast_propagator

        subroutine bcast_sh_mapping
        use mpi
        use mpi_utils
        implicit none
        integer                :: i

         call mpi_bcast(lmax,1,mpi_integer,0,mpi_comm_world,err) 
         call mpi_bcast(compress,1,mpi_logical,0,mpi_comm_world,err) 
         call mpi_bcast(lambda,1,mpi_double_precision,0,mpi_comm_world,err) 

        return
        end subroutine bcast_sh_mapping

        subroutine bcast_propagate
        use mpi
        use mpi_utils
        implicit none
        integer                :: i

         call mpi_bcast(nsteps,1,mpi_integer,0,mpi_comm_world,err) 

        return
        end subroutine bcast_propagate


        subroutine do_action(spindy,phondy,lattice,gsh)
        use spins_dist_rs_class
        use hilbert_dist_class
        use atoms_class
        use phonons_class
        use pulses_class
        use spinham_map_class
        use mpi_utils
        implicit none
        class(spins_hilbert)              :: spindy,gsh
        class(atoms_group)                :: lattice
        class(brillouin)                  :: phondy
        integer                           :: i,t1,t2,rate
        character(len=10)                 :: type_pulse

         select case (operation)

          case ('SET_SYSTEM')
           if(mpi_id.eq.0)    write(6,*) '  Setting Spin System'
           call spindy%spin_bcast()
           call spindy%dist_ij()

          case ('SET_SPINHAM')
           if(mpi_id.eq.0)    write(6,*) '  Setting Spin Hamiltonian'
           call spindy%SH%spinham_bcast()

          case ('SET_GSH_SYSTEM')
           if(mpi_id.eq.0)    write(6,*) '  Setting Spin System'
           call gsh%spin_bcast()
           call gsh%dist_ij()

          case ('SET_GSH_SPINHAM')
           if(mpi_id.eq.0)    write(6,*) '  Setting Spin Hamiltonian'
           call gsh%SH%spinham_bcast()

          case ('MAKE_SH_MAPP')
           if(mpi_id.eq.0)    write(6,*) '  MSH - GSH Mapping'
           call bcast_method(gsh)
           call gsh%make_basis(nex_max,ncorr_max,max_dist,tinv)
           call gsh%make_Hmat_nodes()
           call gsh%make_SH_rep(gsh%SH,-1,-1)
           call gsh%make_kbasis()
           call gsh%make_Hmat_2()
           call bcast_sh_mapping()
           call msh2gsh(spindy,gsh,lmax,compress,lambda) 

          case ('SET_SPH')
           if(mpi_id.eq.0)    write(6,*) '  Setting Spin-Phonon ', &
                                           'Hamiltonian'
           call lattice%atoms_bcast()
           call gen_vars_bcast() 
           call phondy%brillouin_bcast()
           call spindy%SPH%spinphonon_bcast()
           call spindy%SPH2%spinphonon_bcast()

          case ('MAKE_RHO0') 
           if(mpi_id.eq.0)    write(6,*) '  Building Density Matrix'
           call bcast_rho0()
           spindy%beta0=0.0d0
           spindy%beta0(2)=acos(-1.0d0)/2.0d0
           call spindy%make_rho0(type_rho0,spin_temp,rho_restart_file)

          case ('MAKE_HILBERT')
           if(mpi_id.eq.0)    write(6,*) '  Building Hilbert Space'
           call bcast_method(spindy)
           if(spindy%SH%make_dipolar) call spindy%set_dipolar(spindy%SH,ex_list)
           call spindy%make_basis(nex_max,ncorr_max,max_dist,tinv)
           call spindy%make_Hmat_nodes()
           call spindy%SH%rot(euler)
           call spindy%make_SH_rep(spindy%SH,-1,-1)
           call spindy%make_kbasis()
!            call spindy%make_Hmat()
           call spindy%make_Hmat_2()
           if(fulldiag) call spindy%diag_Hmat()
           call spindy%make_S()           
           if(dump_s) call spindy%dump_S()
           if(spindy%make_Rmat .or. spindy%make_R2mat) call spindy%make_Lbasis()
           if(spindy%make_Rmat .and. spindy%SPH%make_dipolar) then
            call spindy%set_sph_dipolar(spindy%SH,spindy%SPH,ex_list)
           endif
           if(spindy%make_R2mat .and. spindy%SPH2%make_dipolar) then
             call spindy%set_sph2_dipolar(spindy%SH,spindy%SPH2,ex_list)
           endif
           if(spindy%make_Rmat .or. spindy%make_R2mat)then 
            call phondy%calc_bands(lattice)
            if(read_fc3) call spindy%get_ph_lt(lattice,phondy,max_phonon_ener)
           endif

          case ('MAKE_U')
           call bcast_propagator()
           call spindy%make_U(step_min,step_nmult)
           if(spindy%make_Rmat .or. spindy%make_R2mat) then
            if(spindy%make_SA) then
             if(spindy%make_PT2)then
               call spindy%make_R2(lattice,phondy,step_min,step_nmult,max_phonon_ener,euler)
              else
               call spindy%make_R(lattice,phondy,step_min,step_nmult,max_phonon_ener,euler)
             endif
            else
             call spindy%make_RL(lattice,phondy,step_min,step_nmult,max_phonon_ener,euler)
            endif
           endif

          case ('PROPAGATE')
           call bcast_propagate()
           call spindy%propagate(start_step,time,nsteps,step,dump_freq)

          case ('MAKE_PULSE')
           if(mpi_id.eq.0)    write(6,*) '  Applying Pulse'
            call bcast_pulse()
            call spindy%make_rot(pulse)
            call spindy%rot_rho(pulse)

          case ('MAKE_ECHO')
           if(mpi_id.eq.0)    write(6,*) '  Calculating Echo'
            call bcast_echo()           
            type_pulse='PI2'
            call set_pulse(pulse_pi2,type_pulse,echo_spin)
            call spindy%make_rot(pulse_pi2)
            type_pulse='PI'
            call set_pulse(pulse_pi,type_pulse,echo_spin)
            call spindy%make_rot(pulse_pi)
            start_step=1
            if(fulldiag)then
             dump_freq=2
             nsteps=1
            else
            endif
            do i=1,echo_nsteps
!             if(i.gt.1) call spindy%make_rho0(type_rho0,spin_temp,rho_restart_file)
             if(i.gt.1) spindy%rho=spindy%rho0
             call spindy%rot_rho(pulse_pi2)
             call spindy%propagate(start_step,time,nsteps,step,dump_freq)
             call spindy%rot_rho(pulse_pi)
             call spindy%propagate(start_step,time,nsteps,step,dump_freq)
             if(fulldiag) then
              step_nmult=step_nmult+echo_step  !/step_min
              call spindy%make_U(step_min,step_nmult)
             else
             endif
            enddo

         end select

        return
        end subroutine do_action

        end module control_variables
