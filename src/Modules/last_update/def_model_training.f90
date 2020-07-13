        module model_trainer_class
        use target_functions_class
        use nets_class
        use properties_class
        implicit none
        
        type, extends(target_function) :: chi2_prop
         type(property), pointer             :: prop
         type(property_ref_val), allocatable :: prop_tr(:)
         type(property_ref_val), allocatable :: prop_te(:)
         integer                             :: nparams 
         integer                             :: npoints_tr
         integer                             :: npoints_te
         double precision                    :: L2=0.001
         logical                             :: regular_L2=.true.
         logical                             :: norm_out=.false.
         contains
         procedure                           :: get_fval  => get_chi2_prop
         procedure                           :: get_fgrad => get_chi2_prop_grad
!         procedure                           :: set_nets
!         procedure                           :: map_nets
!         procedure                           :: out_results
        end type chi2_prop

        type, extends(chi2_prop)  :: chi2_mol_prop
         type(atoms_group), allocatable  :: tr(:)
         type(atoms_group), allocatable  :: te(:)
         contains
!         procedure                           :: read_sets
        end type chi2_mol_prop

        contains

         subroutine get_chi2_prop(this,vec,val)
         implicit none
         class(chi2_prop)               :: this
         double precision               :: val
         double precision, allocatable  :: vec(:)
         integer                        :: i,j

          val=0.0d0

          call this%prop%set_params(vec)
          do i=1,this%npoints_tr
           call this%prop%get_property(this%prop_tr(i)%desc,this%prop_tr(i)%type_desc)
           do j=1,this%prop%nout
            val=val+(this%prop_tr(i)%val(j)-this%prop%val(j))**2
           enddo          
          enddo

          val=val/this%npoints_tr

          ! Add L2 regularization
          if(this%regular_L2) val=val+this%L2*norm2(vec)**2

          val=sqrt(val)

         return
         end subroutine get_chi2_prop

         subroutine get_chi2_prop_grad(this,vec,val,grad)
         implicit none
         class(chi2_prop)               :: this
         double precision               :: val
         double precision, allocatable  :: vec(:),grad(:)
         integer                        :: i,j,v

          val=0.0d0
          grad=0.0d0

          call this%prop%set_params(vec)

          do i=1,this%npoints_tr
           call this%prop%get_property(this%prop_tr(i)%desc,this%prop_tr(i)%type_desc)

           do v=1,this%prop%nout

            val=val+(this%prop_tr(i)%val(v)-this%prop%val(v))**2

            do j=1,size(grad)
             grad(j)=grad(j)-2.0d0*(this%prop_tr(i)%val(v)-this%prop%val(v))*&
                this%prop%grad(j,v)
            enddo

           enddo
          enddo

          grad=grad/this%npoints_tr
          val=val/this%npoints_tr

          if(this%regular_L2)then
           do j=1,size(grad)
            grad(j)=grad(j)+2.0d0*this%L2*vec(j)
           enddo
          endif

          if(this%regular_L2) val=val+this%L2*norm2(vec)**2 
          val=sqrt(val)

         return
         end subroutine get_chi2_prop_grad


        end module model_trainer_class





