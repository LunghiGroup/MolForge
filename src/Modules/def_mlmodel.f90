        module mlmodels_class
        use nets_class
        implicit none

         type mlmodel
          integer                         ::  ndims
          integer                         ::  nparams
          double precision, allocatable   ::  output(:)
          double precision, allocatable   ::  grad(:,:)
          type(net), pointer              ::  NN
          double precision, allocatable   ::  mean_out(:)
          double precision, allocatable   ::  sigma_out(:)
          logical                         ::  norm_out=.false.
         contains
          procedure                       ::  set_params => set_model_params
          procedure                       ::  get_output => get_model_output
          procedure                       ::  get_grad => get_model_grad
         end type mlmodel

        contains

         subroutine set_model_params(this,vec)
         implicit none
         class(mlmodel)                   :: this
         double precision, allocatable    :: vec(:)
                 
           call this%NN%set_parameters(vec)

         return
         end subroutine set_model_params

         subroutine get_model_output(this,desc)         
         use descriptors_class
         implicit none
         class(mlmodel)                   :: this
         double precision, allocatable    :: desc(:)

          if(.not.allocated(this%output)) allocate(this%output(this%ndims))

          this%output=0.0d0

          call this%NN%get_output(desc,this%output)

          if(this%norm_out) this%output=this%sigma_out*this%output+this%mean_out

         return
         end subroutine get_model_output

         subroutine get_model_grad(this,desc)   
         implicit none
         class(mlmodel)                   :: this
         double precision, allocatable    :: desc(:)

          if(.not.allocated(this%output)) allocate(this%output(this%ndims))
          if(.not.allocated(this%grad)) allocate(this%grad(this%nparams,this%ndims))

          this%grad=0.0d0
          this%output=0.0d0

          call this%NN%get_output(desc,this%output)
          call this%NN%backprop(this%grad)

          if(this%norm_out) this%output=this%sigma_out*this%output+this%mean_out

         return
         end subroutine get_model_grad

        end module mlmodels_class
