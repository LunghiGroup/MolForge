        module max_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use lammps_class
        use kind_class
        implicit none

        contains

        subroutine get_tot_kinds(set,tot_kinds)
        type(lammps_obj),allocatable,intent(in)   :: set(:)
        integer :: i
        integer,intent (out) :: tot_kinds

        tot_kinds=0

        do i=1,size(set)

         if (set(i)%nkinds > tot_kinds) then
          tot_kinds=set(i)%nkinds
         end if

        end do

        end subroutine get_tot_kinds

        subroutine get_ave_atoms(set,ave_atom,tot_atom)
        type(lammps_obj),allocatable,intent(in)                :: set(:)
        integer                                     :: part_sum
        integer                                     :: i
        integer,intent(out)                         :: tot_atom
        real(kind=dbl),intent(out)                  :: ave_atom

        part_sum=0

        do i=1,size(set)
          part_sum=part_sum+set(i)%nats
        end do
        
        tot_atom=part_sum
        ave_atom=dble(part_sum)/dble(size(set))

        end subroutine get_ave_atoms

        end module max_class
