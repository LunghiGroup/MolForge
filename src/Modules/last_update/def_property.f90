        module properties_class
        use atoms_class
        use descriptors_class
        use nets_class
        implicit none

         type property_ref_val
          integer                         :: nout
          double precision, allocatable   :: val(:)
          type(descriptor), pointer       :: desc(:)
          integer, pointer                :: type_desc(:)
          contains
         end type property_ref_val

         type property
          integer                         ::  nnets
          integer                         ::  nout
          integer                         ::  nparams
          double precision, allocatable   ::  val(:)
          double precision, allocatable   ::  grad(:,:)
          type(net), pointer              ::  NN(:)
          integer, pointer                ::  types(:)
          double precision, pointer       ::  mean_out(:)
          double precision, pointer       ::  sigma_out(:)
          logical                         ::  norm_out=.false.
         contains
          procedure                       ::  get_property => get_property_plain
          procedure                       ::  get_grad => get_grad_plain
          procedure                       ::  set_params
         end type property

         type, extends(property) :: dipole
         contains
          procedure                       ::  get_property => get_dipole
          procedure                       ::  get_grad => get_grad_dipole
         end type dipole

        contains

         subroutine set_params(this,vec)
         implicit none
         class(property)                  :: this
         double precision, allocatable    :: vec(:),vec_loc(:)
         integer                          :: offset,i,npar

          offset=0
          do i=1,this%nnets
           npar=this%NN(i)%nparams
           if(allocated(vec_loc)) deallocate(vec_loc)
           vec_loc=vec(offset+1:offset+npar)
           call this%NN(i)%set_parameters(vec_loc)
           offset=offset+npar
          enddo
          if(allocated(vec_loc)) deallocate(vec_loc)

         return
         end subroutine set_params

         subroutine get_dipole(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(dipole)                    :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),inps(:),geo(:,:)
         integer                          :: at,de,i

          if(.not.allocated(this%val)) allocate(this%val(this%nout))
          this%val=0.0d0

          do at=1,size(at_desc)
           do i=1,size(this%types)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN(i)%get_output(inps,loc_prop)
             this%val=this%val+loc_prop
            endif
           enddo
          enddo       

         return
         end subroutine get_dipole

         subroutine get_grad_dipole(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(dipole)                    :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
          write(*,*) 'Dipole Gradient not available'
          stop
         return
         end subroutine get_grad_dipole

         subroutine get_property_plain(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(property)                  :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),inps(:)
         integer                          :: at,de,i

          if(.not.allocated(this%val)) allocate(this%val(this%nout))
          this%val=0.0d0

          do at=1,size(at_desc)
           do i=1,size(this%types)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN(i)%get_output(inps,loc_prop)
             this%val=this%val+loc_prop
            endif
           enddo
          enddo       

         return
         end subroutine get_property_plain

         subroutine get_grad_plain(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(property)                  :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: inps(:)
         double precision, allocatable    :: grad_loc(:,:),loc_prop(:)
         integer                          :: at,de,i,v,offset,npar

          if(.not.allocated(this%grad)) allocate(this%grad(this%nparams,this%nout))
          if(.not.allocated(this%val)) allocate(this%val(this%nout))
          
          this%grad=0.0d0
          this%val=0.0d0

          offset=0
          do i=1,this%nnets
           
           npar=this%NN(i)%nparams

           do at=1,size(at_desc)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then

             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc

             call this%NN(i)%get_output(inps,loc_prop)
             this%val=this%val+loc_prop

             call this%NN(i)%backprop(grad_loc)             
             do v=1,this%NN(i)%noutput
              this%grad(offset+1:offset+npar,v)=&
                         this%grad(offset+1:offset+npar,v)+grad_loc(:,v)
             enddo

            endif
           enddo
           offset=offset+npar
          enddo       

         return
         end subroutine get_grad_plain


        end module properties_class
