        module trajectory_class
        use kind_class
        use SNAP_FF_class
        use lammps_class
        use potential_class
        use SNAP_fit_class
        use VdW_class
        use target_functions_class
        implicit none
         
        type                                    :: trajectory
        type(potential_list),allocatable        :: pot(:)
        type(SNAP_fit)                          :: linear_fit
        integer                                 :: max_steps
        integer                                 :: num_pot
        integer                                 :: iter
          
        contains
        
        procedure                               :: init_potentials
        procedure                               :: import_linear_fits
        procedure                               :: link_potentials
        end type trajectory
         
        contains
        
        subroutine init_potentials(this,FF_SNAP,FF_VdW,frame,frame_file)
        implicit none
        class(trajectory)                                 :: this
        type(SNAP_FF),optional                    :: FF_SNAP
        type(VdW_FF),optional                     :: FF_VdW
        type(lammps_obj)                          :: frame
        character(len=*)                          :: frame_file

        if (present(FF_SNAP)) then
         call number_bispec(frame%twojmax,FF_SNAP%num_bisp)
         FF_SNAP%tot_kinds=frame%nkinds
         call FF_SNAP%import_coeff
         FF_SNAP%frame=frame
        end if

        if (present(FF_VdW)) then
         call FF_VdW%import_coeff
         call FF_VdW%import_geo(frame_file)
        end if

        end subroutine init_potentials

        subroutine import_linear_fits(this,active_learning,SNAP)
        implicit none
        class(trajectory)                       :: this
        logical                         :: active_learning
        type(SNAP_fit),optional         :: SNAP

        if (active_learning) then

         if (present(SNAP)) then
          this%linear_fit=SNAP
         end if

        end if
        end subroutine import_linear_fits

        subroutine link_potentials(this,keyword,SNAP,VdW)
        implicit none
        class(trajectory)               :: this
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

        end module trajectory_class
