        module mlmodel_trainer_class
        use target_functions_class
        use nets_class
        use atoms_class
        implicit none
        
        type, extends(target_function)   :: mlmodel_trainer
         type(mlmodel), pointer          :: ML
         integer                         :: npoints_tr,npoints_te
         type(atoms_group), allocatable  :: tr(:),te(:)                 
         double precision, allocatable   :: tr_val(:,:),te_val(:,:)  
         integer                         :: ndim=1
         double precision, allocatable   :: L2val(:)
         integer, allocatable            :: L2id(:)
         logical                         :: L2=.false.
         contains
         procedure                       :: get_fval  => get_chi2
         procedure                       :: get_fgrad => get_chi2_grad
         procedure                       :: set_nets
         procedure                       :: set_L2
         procedure                       :: map_nets
         procedure                       :: read_sets
         procedure                       :: std_sets
         procedure                       :: std_sets_inp
         procedure                       :: std_sets_out
         procedure                       :: out_results
        end type mlmodel_trainer

        contains

         subroutine out_results(this)                 
         implicit none
         class(mlmodel_trainer)         :: this
         integer                        :: i

          open(13,file='tr_rmse.dat')

          do i=1,size(this%tr)
           call this%tr(i)%ML%get_output(this%tr(i)%at_desc(1)%desc)           
           if(this%ML%norm_out)then
            this%tr_val(i,:)=this%tr_val(i,:)*this%ML%sigma_out+this%ML%mean_out
            this%tr(i)%ML%output=this%tr(i)%ML%output*this%ML%sigma_out+this%ML%mean_out           
           endif
           write(13,*) this%tr_val(i,:),this%tr(i)%ML%output
          enddo

          close(13)

          open(13,file='te_rmse.dat')

          do i=1,size(this%te)
           call this%te(i)%ML%get_output(this%te(i)%at_desc(1)%desc)
           if(this%ML%norm_out)then
            this%te_val(i,:)=this%te_val(i,:)*this%ML%sigma_out+this%ML%mean_out
            this%te(i)%ML%output=this%te(i)%ML%output*this%ML%sigma_out+this%ML%mean_out           
           endif
           write(13,*) this%te_val(i,:),this%te(i)%ML%output
          enddo

          close(13)

         return
         end subroutine out_results

         subroutine std_sets(this)
         implicit none
         class(mlmodel_trainer)                    :: this
           call this%std_sets_inp()
           call this%std_sets_out()
         return
         end subroutine  std_sets

         subroutine std_sets_inp(this)
         implicit none
         class(mlmodel_trainer)                    :: this
         integer                                   :: i,j,l,v

          this%ML%NN%norm_inp=.true.
          allocate(this%ML%NN%mean_inp(this%ML%NN%ninput))
          allocate(this%ML%NN%sigma_inp(this%ML%NN%ninput))
          this%ML%NN%mean_inp=0.0d0
          this%ML%NN%sigma_inp=0.0d0

          do i=1,this%npoints_tr
           do v=1,size(this%tr(i)%at_desc)
            do j=1,size(this%tr(i)%at_desc(v)%desc)
             this%ML%NN%mean_inp(j)=this%ML%NN%mean_inp(j)+this%tr(i)%at_desc(v)%desc(j)
            enddo
           enddo
          enddo

          this%ML%NN%mean_inp=this%ML%NN%mean_inp/this%npoints_tr

          do i=1,this%npoints_tr
           do v=1,size(this%tr(i)%at_desc)
            do j=1,size(this%tr(i)%at_desc(v)%desc)
              this%ML%NN%sigma_inp(j)=this%ML%NN%sigma_inp(j)+&
                     (this%tr(i)%at_desc(v)%desc(j)-this%ML%NN%mean_inp(j))**2/this%npoints_tr
            enddo
           enddo
          enddo

          this%ML%NN%sigma_inp=sqrt(this%ML%NN%sigma_inp)

         return
         end subroutine std_sets_inp

         subroutine std_sets_out(this)
         implicit none
         class(mlmodel_trainer)                 :: this
         integer                                :: i

         !! Normalize Training Set Output

          this%ML%norm_out=.true.
          allocate(this%ML%mean_out(this%ndim))
          allocate(this%ML%sigma_out(this%ndim))
          this%ML%mean_out=0.0d0
          this%ML%sigma_out=0.0d0

          do i=1,this%npoints_tr
           this%ML%mean_out=this%ML%mean_out+this%tr_val(i,:)/this%npoints_tr
          enddo

          do i=1,this%npoints_tr
           this%ML%sigma_out=this%ML%sigma_out+(this%tr_val(i,:)-this%ML%mean_out)**2/this%npoints_tr
          enddo

          this%ML%sigma_out=sqrt(this%ML%sigma_out)
           
          do i=1,this%npoints_tr
           this%tr_val(i,:)=(this%tr_val(i,:)-this%ML%mean_out)/this%ML%sigma_out
          enddo

          do i=1,this%npoints_te
           this%te_val(i,:)=(this%te_val(i,:)-this%ML%mean_out)/this%ML%sigma_out
          enddo

         return
         end subroutine std_sets_out


         subroutine read_sets(this,input_file)
         implicit none
         class(mlmodel_trainer)         :: this
         integer                        :: i,l,j
         character(len=100)             :: input_file
         character(len=100)             :: tr_file,te_file
         character(len=100)             :: tr_val_file,te_val_file
         character(len=10)              :: kind_desc='CART'

          open(13,file=input_file)
          read(13,*) this%npoints_tr,this%ndim,tr_file,tr_val_file
          read(13,*) this%npoints_te,this%ndim,te_file,te_val_file
          close(13)

          allocate(this%tr(this%npoints_tr))
          allocate(this%tr_val(this%npoints_tr,this%ndim))
          allocate(this%te(this%npoints_te))
          allocate(this%te_val(this%npoints_te,this%ndim))

          open(13,file=tr_file)          
          open(14,file=tr_val_file)

          do i=1,this%npoints_tr
!           call this%tr(i)%read_extended_xyz(IOid=13)
           call this%tr(i)%read_xyz(IOid=13)
           read(14,*) this%tr_val(i,:)
           call this%tr(i)%build_descriptors(kind_desc)
          enddo

          close(13)
          close(14)

          open(13,file=te_file)
          open(14,file=te_val_file)

          do i=1,this%npoints_te
!           call this%te(i)%read_extended_xyz(IOid=13)
           call this%te(i)%read_xyz(IOid=13)
           read(14,*) this%te_val(i,:)
           call this%te(i)%build_descriptors(kind_desc)
          enddo

          close(13)          
          close(14)
          
         return
         end subroutine read_sets

         subroutine set_L2(this,L2val,L2id,L2)
         use general_types_class
         use random_numbers_class 
         implicit none
         class(mlmodel_trainer)                              :: this
         integer                                  :: i,j
         double precision, allocatable, optional  :: L2val(:)
         integer, allocatable, optional           :: L2id(:)
         double precision, optional               :: L2
         
          if(present(L2)) this%L2=.true.
          if(present(L2id).and.present(L2val)) this%L2=.true.

          if(present(L2))then
           allocate(this%L2val(this%ML%nparams))
           allocate(this%L2id(this%ML%nparams))
           this%L2val=L2
           do i=1,size(this%L2id)
            this%L2id(i)=i
           enddo
          else if (this%L2) then                  
           this%L2val=L2val
           this%L2id=L2id
          endif

         return
         end subroutine set_L2

         subroutine set_nets(this,ninp,topo)
         use general_types_class
         use random_numbers_class 
         implicit none
         class(mlmodel_trainer)         :: this
         integer                        :: i,j
         integer                        :: ninp
         type(vector_int)               :: topo
         double precision, allocatable  :: vec(:)
        
          call init_random_seed()

          this%ML%nparams=0

          call this%ML%NN%set_topology(ninp,size(topo%v),topo%v)

          this%ML%nparams=this%ML%NN%nparams

          if(allocated(vec)) deallocate(vec)
          allocate(vec(this%ML%NN%nparams))
          do j=1,this%ML%NN%nparams
           call random_number(vec(j))
          enddo
          call this%ML%NN%set_parameters(vec)
          if(allocated(vec)) deallocate(vec)

          call this%map_nets()

         return
         end subroutine set_nets
        
         subroutine map_nets(this)
         implicit none
         class(mlmodel_trainer)                :: this
         integer                               :: i

          do i=1,size(this%tr)
           this%tr(i)%ML => this%ML
          enddo

          do i=1,size(this%te)
           this%te(i)%ML => this%ML
          enddo

         return
         end subroutine map_nets

         subroutine get_chi2(this,vec,val)
         implicit none
         class(mlmodel_trainer)         :: this
         double precision               :: val
         double precision, allocatable  :: vec(:),outs(:),inps(:),vec_loc(:)
         integer                        :: i,j,npar

          val=0.0d0

          call this%ML%set_params(vec)

          do i=1,size(this%tr)
           call this%tr(i)%ML%get_output(this%tr(i)%at_desc(1)%desc)
           do j=1,this%ndim
            val=val+(this%tr(i)%ML%output(j)-this%tr_val(i,j))**2
           enddo
          enddo
          
          val=val/size(this%tr)

          if( this%L2 )then
           do i=1,size(this%L2id)
            val=val+this%L2val(i)*vec(this%L2id(i))**2
           enddo
          endif

          val=sqrt(val)

         return
         end subroutine get_chi2

         subroutine get_chi2_grad(this,vec,val,grad)
         implicit none
         class(mlmodel_trainer)         :: this
         integer                        :: i,j,offset,npar,l
         double precision               :: val
         double precision, allocatable  :: vec(:),vec_loc(:),rij(:,:)
         double precision, allocatable  :: grad(:),grad_loc(:,:),grad2(:)


          call this%ML%set_params(vec)
          if(allocated(grad)) deallocate(grad)
          allocate(grad(this%ML%nparams))
          val=0.0d0
          grad=0.0d0

          do i=1,size(this%tr)

           call this%tr(i)%ML%get_grad(this%tr(i)%at_desc(1)%desc)

           do j=1,size(grad)
            do l=1,this%ndim
             grad(j)=grad(j)-2.0d0*(this%tr_val(i,l)-this%tr(i)%ML%output(l))*&
                 this%tr(i)%ML%grad(j,l)
            enddo
           enddo

           do j=1,this%ndim
            val=val+(this%tr(i)%ML%output(j)-this%tr_val(i,j))**2
           enddo

          enddo

          grad=grad/size(this%tr_val)
          val=val/size(this%tr_val)

          if(this%L2)then

           do j=1,size(this%L2id)
            val=val+this%L2val(j)*vec(this%L2id(j))**2
            grad(this%L2id(j))=grad(this%L2id(j))&
                +2.0e0*this%L2val(j)*vec(this%L2id(j))
           enddo

          endif

          val=sqrt(val)

         return
         end subroutine get_chi2_grad

        end module mlmodel_trainer_class





