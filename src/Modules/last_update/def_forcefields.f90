        module ffs_class
        use nets_class
        implicit none

         type force_field
          double precision                ::  ener=0.0d0
          double precision                ::  disp_ener=0.0d0
          double precision                ::  coul_ener=0.0d0
          double precision                ::  local_ener=0.0d0
          double precision, allocatable   ::  charge(:)
          double precision, allocatable   ::  C6(:)
          double precision, allocatable   ::  forces(:)
          double precision, allocatable   ::  dipole(:)
          double precision, allocatable   ::  ffgrad(:)
          double precision, allocatable   ::  local_ffgrad(:)
          double precision, allocatable   ::  coul_ffgrad(:)
          double precision, allocatable   ::  disp_ffgrad(:)
          integer                         ::  tot_charge=0
          integer                         ::  nparams
          integer                         ::  local_nparams
          integer                         ::  coul_nparams
          integer                         ::  disp_nparams
          integer                         ::  nnets
          type(net), pointer              ::  NN(:)
          type(net), pointer              ::  NN_C6(:)
          type(net), pointer              ::  NN_charge(:)
          integer, pointer                ::  types(:)         
          logical                         ::  do_local_ener=.true.
          logical                         ::  do_coul_ener=.false.
          logical                         ::  do_disp_ener=.false.
          logical                         ::  fit_local_ener=.true.
          logical                         ::  fit_coul_ener=.false.
          logical                         ::  fit_disp_ener=.false.
          double precision                ::  mean_out
          double precision                ::  sigma_out
          logical                         ::  norm_out=.false.
         contains
          procedure                       ::  set_ff_params
          procedure                       ::  update_ff_params
          procedure                       ::  get_charges
          procedure                       ::  get_C6
          procedure                       ::  get_ener 
          procedure                       ::  get_local_ener
          procedure                       ::  get_disp_ener
          procedure                       ::  get_coul_ener
          procedure                       ::  get_ffgrad
          procedure                       ::  get_local_ffgrad
          procedure                       ::  get_coul_ffgrad
!          procedure                       ::  get_disp_ffgrad
         end type force_field

        contains

         subroutine update_ff_params(this,vec)
         implicit none
         class(force_field)               :: this
         double precision, allocatable    :: vec(:),vec_loc(:)
         integer                          :: offset,i,npar
                 
          offset=0

          if(this%fit_local_ener)then 
           do i=1,this%nnets
            npar=this%NN(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

          if(this%fit_coul_ener)then 
           do i=1,this%nnets
            npar=this%NN_charge(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN_charge(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

          if(this%fit_disp_ener)then 
           do i=1,this%nnets
            npar=this%NN_C6(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN_C6(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

         return
         end subroutine update_ff_params

         subroutine set_ff_params(this,vec)
         implicit none
         class(force_field)               :: this
         double precision, allocatable    :: vec(:),vec_loc(:)
         integer                          :: offset,i,npar
                 
          offset=0

          if(this%do_local_ener)then 
           do i=1,this%nnets
            npar=this%NN(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

          if(this%do_coul_ener)then 
           do i=1,this%nnets
            npar=this%NN_charge(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN_charge(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

          if(this%do_disp_ener)then 
           do i=1,this%nnets
            npar=this%NN_C6(i)%nparams
            if(allocated(vec_loc)) deallocate(vec_loc)
            vec_loc=vec(offset+1:offset+npar)
            call this%NN_C6(i)%set_parameters(vec_loc)
            offset=offset+npar
           enddo
           if(allocated(vec_loc)) deallocate(vec_loc)
          endif

         return
         end subroutine set_ff_params

         subroutine get_ffgrad(this,at_desc,at_type,rij)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         integer                          :: offset,npar
         double precision, allocatable    :: rij(:,:)

          if(.not.allocated(this%ffgrad)) allocate(this%ffgrad(this%nparams))

          this%ffgrad=0.0d0
          this%ener=0.0d0
          
          offset=0

          if(this%do_local_ener)then
           call this%get_local_ffgrad(at_desc,at_type)
           this%ener=this%ener+this%local_ener
           npar=this%local_nparams
           this%ffgrad(offset+1:offset+npar)=&
                         this%ffgrad(offset+1:offset+npar)+this%local_ffgrad
           offset=offset+npar
          endif

          if(this%do_coul_ener)then
           call this%get_charges(at_desc,at_type)         
           call this%get_coul_ffgrad(at_desc,at_type,rij)
           this%ener=this%ener+this%coul_ener
           npar=this%coul_nparams
           this%ffgrad(offset+1:offset+npar)=&
                         this%ffgrad(offset+1:offset+npar)+this%coul_ffgrad
           offset=offset+npar
          endif

          if(this%do_disp_ener)then
           call this%get_C6(at_desc,at_type)         
!           call this%get_disp_ffgrad(rij)
           this%ener=this%ener+this%disp_ener
           npar=this%disp_nparams
           this%ffgrad(offset+1:offset+npar)=&
                         this%ffgrad(offset+1:offset+npar)+this%disp_ffgrad
           offset=offset+npar
          endif

          if(this%norm_out) this%ener=this%sigma_out*this%ener+this%mean_out
          if(this%norm_out) this%ffgrad=this%sigma_out*this%ffgrad

         return
         end subroutine get_ffgrad

         subroutine get_ener(this,at_desc,at_type,rij)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: rij(:,:)
         
          this%ener=0.0d0
           
          if(this%do_local_ener)then
           call this%get_local_ener(at_desc,at_type)
           this%ener=this%ener+this%local_ener
          endif

          if(this%do_coul_ener)then
           call this%get_charges(at_desc,at_type)         
           call this%get_coul_ener(rij)
           this%ener=this%ener+this%coul_ener
          endif

          if(this%do_disp_ener)then
           call this%get_C6(at_desc,at_type)         
           call this%get_disp_ener(rij)
           this%ener=this%ener+this%disp_ener
          endif

          if(this%norm_out) this%ener=this%sigma_out*this%ener+this%mean_out

         return
         end subroutine get_ener 

         subroutine get_local_ener(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),vec(:),inps(:)
         integer                          :: at,de,i

          this%local_ener=0.0d0

          do at=1,size(at_desc)
           do i=1,size(this%types)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(.not.allocated(loc_prop)) allocate(loc_prop(1))
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN(i)%get_output(inps,loc_prop)
             this%local_ener=this%local_ener+loc_prop(1)       
            endif
           enddo
          enddo       

         return
         end subroutine get_local_ener

         subroutine get_local_ffgrad(this,at_desc,at_type)   
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),vec(:),inps(:)
         double precision, allocatable    :: grad_loc(:,:)
         integer                          :: at,de,i,offset,npar

          if(.not.allocated(this%local_ffgrad)) &
                  allocate(this%local_ffgrad(this%local_nparams))

          this%local_ffgrad=0.0d0
          this%local_ener=0.0d0

          offset=0

          do i=1,this%nnets

           npar=this%NN(i)%nparams
           if(allocated(grad_loc)) deallocate(grad_loc)
           allocate(grad_loc(npar,1))

           do at=1,size(at_desc)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(.not.allocated(loc_prop)) allocate(loc_prop(1))
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN(i)%get_output(inps,loc_prop)
             this%local_ener=this%local_ener+loc_prop(1)       
             call this%NN(i)%backprop(grad_loc)
             this%local_ffgrad(offset+1:offset+npar)=&
                        this%local_ffgrad(offset+1:offset+npar)+grad_loc(:,1)

            endif
           enddo
           offset=offset+npar
          enddo       

         return
         end subroutine get_local_ffgrad

         subroutine get_charges(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),vec(:),inps(:)
         integer                          :: at,de,i

          if(allocated(this%charge)) deallocate(this%charge)
          allocate(this%charge(size(at_desc)))
          this%charge=0.0d0

          do at=1,size(at_desc)
           do i=1,size(this%types)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(.not.allocated(loc_prop)) allocate(loc_prop(1))
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN_charge(i)%get_output(inps,loc_prop)
             this%charge(at)=this%charge(at)+loc_prop(1)  
            endif
           enddo
          enddo       

          this%charge=this%charge+(this%tot_charge-sum(this%charge))/size(this%charge)

         return
         end subroutine get_charges

         subroutine get_coul_ener(this,rij)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         double precision, allocatable    :: rij(:,:)
         integer                          :: i,j
         
          this%coul_ener=0.0d0

          do i=1,size(rij,1)
           do j=i+1,size(rij,2)
            this%coul_ener=this%coul_ener+0.5d0*(erf(3.1d0)+1)*this%charge(i)*this%charge(j)/rij(i,j)
           enddo
          enddo       

         return
         end subroutine get_coul_ener

         subroutine get_coul_ffgrad(this,at_desc,at_type,rij)   
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),vec(:),inps(:)
         double precision, allocatable    :: grad_loc(:,:),rij(:,:)
         integer                          :: at,de,i,offset,npar,j

          if(.not.allocated(this%coul_ffgrad)) &
                  allocate(this%coul_ffgrad(this%coul_nparams))

          this%coul_ffgrad=0.0d0
          this%coul_ener=0.0d0

          offset=0

          do i=1,this%nnets

           npar=this%NN_charge(i)%nparams
           if(allocated(grad_loc)) deallocate(grad_loc)
           allocate(grad_loc(npar,1))

           do at=1,size(at_desc)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(.not.allocated(loc_prop)) allocate(loc_prop(1))
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN_charge(i)%get_output(inps,loc_prop)
             call this%NN_charge(i)%backprop(grad_loc)
             do j=1,size(at_desc)
              if(j.ne.at) this%coul_ffgrad(offset+1:offset+npar)=&
                          this%coul_ffgrad(offset+1:offset+npar)+grad_loc(:,1)*&
                          this%charge(j)/rij(at,j)*0.5d0*(erf(3.1d0)+1)
             enddo
            endif
           enddo
           offset=offset+npar
          enddo       

          do i=1,size(rij,1)
           do j=i+1,size(rij,2)
            this%coul_ener=this%coul_ener+0.5d0*(erf(3.1d0)+1)*this%charge(i)*this%charge(j)/rij(i,j)
           enddo
          enddo       

         return
         end subroutine get_coul_ffgrad

         subroutine get_C6(this,at_desc,at_type)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         type(descriptor), pointer        :: at_desc(:)
         integer, pointer                 :: at_type(:)
         double precision, allocatable    :: loc_prop(:),vec(:),inps(:)
         integer                          :: at,de,i

          if(allocated(this%C6)) deallocate(this%C6)
          allocate(this%C6(size(at_desc)))
          this%C6=0.0d0

          do at=1,size(at_desc)
           do i=1,size(this%types)
            if(this%types(i).eq.at_type(at) .or. this%types(i).eq.0 )then
             if(.not.allocated(loc_prop)) allocate(loc_prop(1))
             if(allocated(inps)) deallocate(inps)
             inps=at_desc(at)%desc
             call this%NN_C6(i)%get_output(inps,loc_prop)
             this%C6(at)=this%C6(at)+loc_prop(1)  
            endif
           enddo
          enddo       

         return
         end subroutine get_C6

         subroutine get_disp_ener(this,rij)         
         use descriptors_class
         implicit none
         class(force_field)               :: this
         double precision, allocatable    :: rij(:,:)
         integer                          :: i,j
         
          this%disp_ener=0.0d0

          do i=1,size(rij,1)
           do j=i+1,size(rij,2)
            this%disp_ener=this%disp_ener+this%C6(i)*this%C6(j)/rij(i,j)**6
           enddo
          enddo

         return
         end subroutine get_disp_ener


        end module FFs_class
