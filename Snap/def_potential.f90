        module potential_class
        use target_functions_class
        use lammps_class
        use atoms_class
        use kind_class
        implicit none
         
         type,extends(target_function)           :: potential
          type(lammps_obj)                       :: frame
          real(kind=dbl),allocatable             :: grad(:)
          real(kind=dbl)                         :: energy

          contains
          
          procedure                              :: import_coeff => import_generic
          procedure                              :: get_fgrad => get_forces
          procedure                              :: get_fval => get_energy
         end type potential
            
        contains

        subroutine get_energy(this,vec,val)
        implicit none
        class(potential)                :: this
        real(kind=dbl),allocatable      :: vec(:)
        real(kind=dbl)                  :: val


        end subroutine get_energy
        
        subroutine get_forces(this,vec,val,grad)
        implicit none
        class(potential)                :: this
        real(kind=dbl),allocatable      :: vec(:)
        real(kind=dbl)                  :: val
        real(kind=dbl),allocatable      :: grad(:)
        end subroutine get_forces

        subroutine import_generic(this)
        implicit none
        class(potential)               :: this
        
        end subroutine import_generic

         
        end module potential_class
