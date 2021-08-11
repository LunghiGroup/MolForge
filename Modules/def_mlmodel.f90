        module mlmodels_class
        use nets_class
        implicit none

         type mlmodel
          integer                         ::  ndims
          integer                         ::  nparams
          double precision, allocatable   ::  output(:)
          double precision, allocatable   ::  grad(:,:)
          type(net)                       ::  NN
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
         double precision, allocatable    :: loc_prop(:),inps(:)


          if(.not.allocated(this%output)) allocate(this%output(this%ndims))

          this%output=0.0d0

          if(.not.allocated(loc_prop)) allocate(loc_prop(this%ndims))
          if(allocated(inps)) deallocate(inps)
          inps=desc

          call this%NN%get_output(inps,loc_prop)

          this%output=loc_prop

         return
         end subroutine get_model_output

         subroutine get_model_grad(this,desc)   
         implicit none
         class(mlmodel)                   :: this
         double precision, allocatable    :: desc(:)
         double precision, allocatable    :: loc_prop(:),inps(:),grad(:,:)
         integer                          :: v,l

          if(.not.allocated(this%output)) allocate(this%output(this%ndims))
          if(.not.allocated(this%grad)) allocate(this%grad(this%nparams,this%ndims))

          this%grad=0.0d0
          this%output=0.0d0

          if(.not.allocated(loc_prop)) allocate(loc_prop(this%ndims))
          if(allocated(inps)) deallocate(inps)
          if(allocated(grad)) deallocate(grad)
          inps=desc

          call this%NN%get_output(inps,loc_prop)
          call this%NN%backprop(grad)
          this%output=loc_prop
          this%grad=grad          

         return
         end subroutine get_model_grad

        end module mlmodels_class
