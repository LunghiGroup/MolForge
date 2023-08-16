        module VdW_class        
        use atoms_class
        use parameters_class
        use potential_class
        implicit none
        
        type,extends(potential)                :: VdW_FF
        
        contains

        procedure                              :: import_coeff => import_coeff_VdW
        procedure                              :: import_geo => import_geo_VdW
        procedure                              :: get_fgrad => get_VdW_force
        procedure                              :: get_fval => get_VdW_en
        procedure                              :: grimme_d3

        end type VdW_FF

        contains

        subroutine get_VdW_en(this,vec,val)
        use dftd3_api
        implicit none
        class(VdW_FF)                              :: this
        real(kind=dbl),allocatable              :: vec(:)
        real(kind=dbl)                          :: val       
        
        call this%grimme_d3(flag_force=.false.)
        val=this%energy

        end subroutine get_VdW_en

        subroutine get_VdW_force(this,vec,val,grad)
        use dftd3_api
        implicit none
        class(VdW_FF)                :: this
        real(kind=dbl),allocatable      :: vec(:)
        real(kind=dbl)                  :: val
        real(kind=dbl),allocatable      :: grad(:)
        
        
        if (.not.allocated(this%grad)) allocate(this%grad(3*this%frame%nats))
        call this%grimme_d3(flag_force=.true.)

        end subroutine get_VdW_force

        subroutine import_coeff_VdW(this)
        implicit none
        class(VdW_FF)               :: this

        end subroutine import_coeff_VdW

        subroutine import_geo_VdW(this,filename)
        implicit none
        class(VdW_FF)               :: this
        character(len=50)         :: filename

        call this%frame%read_extended_xyz(3,trim(filename))

        end subroutine import_geo_VdW

        subroutine grimme_d3(this,flag_force)
        use dftd3_api
        implicit none
        
        class(VdW_FF)                                     :: this
        double precision, allocatable                  :: coords(:,:)
        logical                                        :: flag_force
        ! integer, parameter :: species(nAtoms) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        !  & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, &
        !  & 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

        ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
        ! They must be converted to Bohr before passed to dftd3
        !real(wp), parameter :: latVecs(3, 3) = reshape([&
        !    &  8.0000000000E+00,   0.0000000000E+00,   0.0000000000E+00, &
        !    &  0.0000000000E+00,   8.0000000000E+00,   0.0000000000E+00, &
        !    &  0.0000000000E+00,   0.0000000000E+00,   1.5000000000E+01  &
        !    & ] * AA__Bohr, [3, 3])

        !integer, parameter :: nSpecies = 4
        !character(2),  :: speciesNames(nSpecies) = [ 'N ', 'C ', 'O ', 'H ']
        real(kind=dbl),allocatable    :: grads(:,:)
        integer,allocatable           :: atnum(:)
        type(dftd3_input)             :: input
        type(dftd3_calc)              :: dftd3
        double precision              :: edisp
        double precision              :: stress(3, 3)
        integer                       :: i,j
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Initialize input
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! You can set input variables if you like, or just leave them on their
        ! defaults, which are the same as the dftd3 program uses.

        !! Threebody interactions (default: .false.)
        !input%threebody = .true.
        !
        !! Numerical gradients (default: .false.)
        !input%numgrad = .false.
        !
        !! Cutoffs (below you find the defaults)
        !input%cutoff = sqrt(9000.0_wp)
        !input%cutoff_cn = sqrt(1600.0_wp)

        ! Initialize dftd3
        call dftd3_init(dftd3, input)

        ! Choose functional. Alternatively you could set the parameters manually
        ! by the dftd3_set_params() function.
        call dftd3_set_functional(dftd3, func='pbe', version=4, tz=.false.)

        allocate(atnum(this%frame%nats))
        ! Convert species name to atomic number for each atom
        atnum(:) = get_atomic_number(this%frame%label(this%frame%kind))
        
        if (flag_force) allocate(grads(3,this%frame%nats))
        allocate(coords(3,this%frame%nats))
        coords=transpose(this%frame%x)*A_to_B
        
        !Calculate dispersion and gradients for non-periodic case

        call dftd3_dispersion(dftd3, coords,atnum,this%energy,grads)
        
        if (flag_force) then
         
         do i=1,this%frame%nats
          do j=1,3
           this%grad((i-1)*3+j)=grads(j,i)
          end do
         end do
        
         deallocate(grads)
        
        end if

        !this%en_VdW= edisp
        !this%grads_VdW=grads
        !write(*, "(A)") "*** Dispersion for non-periodic case"
        !open(111, file='VdW_ener_MolForge_20220225.txt',action='write',position='append')
        !write(111,*)  this%edisp*Har_to_kc
        !close(111)
        !write(*, "(A)") "Gradients [au]:"
        !write(*, "(3ES20.12)") this%grads_VdW
        !write(*, *) this%grads_VdW(1,:)
        !write(*, *)

        ! Calculate dispersion and gradients for periodic case
        !call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, edisp, grads, stress)
        !write(*, "(A)") "*** Dispersion for periodic case"
        !write(*, "(A,ES20.12)") "Energy [au]:", edisp
        !write(*, "(A)") "Gradients [au]:"
        !write(*, "(3ES20.12)") grads
        !write(*, "(A)") "Stress [au]:"
        !write(*, "(3ES20.12)") stress
        
        deallocate(coords,atnum)

        end subroutine grimme_d3

        end module VdW_class
