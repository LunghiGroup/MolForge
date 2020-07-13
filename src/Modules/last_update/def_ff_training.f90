        module ff_trainer_class
        use target_functions_class
        use nets_class
        use atoms_class
        implicit none
        
        type, extends(target_function) :: chi2
         type(force_field), pointer      :: FF
         integer                         :: npoints_tr,npoints_te
         type(atoms_group), allocatable  :: tr(:),te(:)                 
         double precision, allocatable   :: tr_val(:),te_val(:)  
         double precision, allocatable   :: tr_fval(:,:),te_fval(:,:)  
         double precision, allocatable   :: tr_dval(:,:),te_dval(:,:)  
         integer                         :: ndim=1
         double precision                :: L2=0.001
         logical                         :: regular_L2=.true.
         logical                         :: fit_ener=.true.
         logical                         :: fit_forces=.false.
         logical                         :: fit_dipols=.false.
         contains
         procedure                       :: get_fval  => get_chi2
         procedure                       :: get_fgrad => get_chi2_grad
         procedure                       :: ridge
         procedure                       :: set_nets
         procedure                       :: map_nets
         procedure                       :: read_sets
         procedure                       :: std_sets
         procedure                       :: std_sets_inp
         procedure                       :: std_sets_out
         procedure                       :: out_results
        end type chi2

        contains

         subroutine ridge(this)
         implicit none
         class(chi2)                    :: this
         integer                        :: i,j,v,l
         integer                        :: dimA,dimB,lwork,inf
         integer                        :: npar,offset
         double precision, allocatable  :: A(:,:),B(:)
         double precision, allocatable  :: work(:),vec(:)
         double precision               :: ener

          write(*,*) 'Generating SNAP'

          dimA=this%FF%nparams
          dimB=size(this%tr_val)
          if(this%regular_L2) dimB=dimB+this%FF%nparams-1
          allocate(A(dimB,dimA))
          allocate(B(dimB))
          A=0.0d0
          B=0.0d0
          
          do i=1,size(this%tr_val)

           B(i)=this%tr_val(i)

           do v=1,this%tr(i)%nats
            offset=0
            do l=1,this%FF%nnets
             npar=this%FF%NN(l)%nparams
             if(this%tr(i)%kind(v).eq.this%FF%types(l)) exit
             offset=offset+this%FF%NN(l)%nparams
            enddo
            do j=1,npar-1
             A(i,offset+j)=A(i,offset+j)+this%tr(i)%at_desc(v)%desc(j)
            enddo
            A(i,offset+npar)=A(i,offset+npar)+1.0d0
           enddo

          enddo

          if(this%regular_L2)then
           offset=size(this%tr_val)
           do l=1,this%FF%nparams-1
            A(l+offset,l)=this%L2
            B(l+offset)=0.0d0
           enddo
          endif

          lwork=dimB+64*dimB+1000
          allocate(work(lwork))
          call dgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
          deallocate(work)
          if(inf.ne.0)then
           write(*,*) 'dgels crashed'
           stop
          endif

          offset=0

          ! Print Coefficients

          open(14,file='snapcoeff')
          write(14,*) this%FF%nnets,this%FF%NN(1)%nparams
          offset=0
          do i=1,this%FF%nnets 
           write(14,*) i,' 0.500',' 1.000'
           write(14,*) B(offset+this%FF%NN(i)%nparams)
           do l=1,this%FF%NN(i)%nparams-1
            write(14,*) B(offset+l)
           enddo
           offset=offset+this%FF%NN(i)%nparams
          enddo
          close(14)

          vec=B(1:this%FF%nparams)
          call this%FF%set_ff_params(vec)

          write(*,*) 'SNAP Generated Successfully'

         return
         end subroutine ridge

         subroutine out_results(this)                 
         implicit none
         class(chi2)                    :: this
         integer                        :: i
         double precision, allocatable  :: rij(:,:)

          open(13,file='tr_rmse.dat')
          open(14,file='tr_ener.dat')

          do i=1,size(this%tr)
           rij=this%tr(i)%dist(:,1,:,1)
           call this%tr(i)%FF%get_ener(this%tr(i)%at_desc,this%tr(i)%kind,rij)           
           write(13,*) this%tr_val(i),this%tr(i)%FF%ener
           write(14,*) this%FF%sigma_out*this%tr(i)%FF%local_ener,&
                this%FF%sigma_out*this%tr(i)%FF%coul_ener,this%FF%sigma_out*this%tr(i)%FF%disp_ener
          enddo

          close(13)
          close(14)

          open(13,file='te_rmse.dat')
          open(14,file='te_ener.dat')

          do i=1,size(this%te)
           rij=this%te(i)%dist(:,1,:,1)
           call this%te(i)%FF%get_ener(this%te(i)%at_desc,this%te(i)%kind,rij)
           write(13,*) this%te_val(i),this%te(i)%FF%ener
           write(14,*) this%FF%sigma_out*this%te(i)%FF%local_ener,&
                this%FF%sigma_out*this%te(i)%FF%coul_ener,this%FF%sigma_out*this%te(i)%FF%disp_ener
          enddo

          close(13)
          close(14)

         return
         end subroutine out_results

         subroutine std_sets(this)
         implicit none
         class(chi2)                    :: this
           call this%std_sets_inp()
           call this%std_sets_out()
         return
         end subroutine  std_sets

         subroutine std_sets_inp(this)
         implicit none
         class(chi2)                    :: this
         integer                        :: i,j,l,v
         integer, allocatable           :: nval(:)

          do l=1,this%FF%nnets
           if(this%FF%do_local_ener)then
            this%FF%NN(l)%norm_inp=.true.
            allocate(this%FF%NN(l)%mean_inp(this%FF%NN(l)%ninput))
            allocate(this%FF%NN(l)%sigma_inp(this%FF%NN(l)%ninput))
            this%FF%NN(l)%mean_inp=0.0d0
            this%FF%NN(l)%sigma_inp=0.0d0
           endif
           if(this%FF%do_coul_ener)then
            this%FF%NN_charge(l)%norm_inp=.true.
            allocate(this%FF%NN_charge(l)%mean_inp(this%FF%NN_charge(l)%ninput))
            allocate(this%FF%NN_charge(l)%sigma_inp(this%FF%NN_charge(l)%ninput))
            this%FF%NN_charge(l)%mean_inp=0.0d0
            this%FF%NN_charge(l)%sigma_inp=0.0d0
           endif
           if(this%FF%do_disp_ener)then
            this%FF%NN_C6(l)%norm_inp=.true.
            allocate(this%FF%NN_C6(l)%mean_inp(this%FF%NN_C6(l)%ninput))
            allocate(this%FF%NN_C6(l)%sigma_inp(this%FF%NN_C6(l)%ninput))
            this%FF%NN_C6(l)%mean_inp=0.0d0
            this%FF%NN_C6(l)%sigma_inp=0.0d0
           endif
          enddo

          allocate(nval(this%FF%nnets))
          nval=0

          do i=1,this%npoints_tr
           do v=1,size(this%tr(i)%at_desc)

            do l=1,this%FF%nnets
             if(this%tr(i)%kind(v).eq.this%FF%types(l) .or. this%FF%types(l).eq.0  )then
              do j=1,size(this%tr(i)%at_desc(v)%desc)
               if(this%FF%do_local_ener)then
                this%FF%NN(l)%mean_inp(j)=this%FF%NN(l)%mean_inp(j)+this%tr(i)%at_desc(v)%desc(j)
               endif
               if(this%FF%do_coul_ener)then
                this%FF%NN_charge(l)%mean_inp(j)=this%FF%NN_charge(l)%mean_inp(j)+this%tr(i)%at_desc(v)%desc(j)
               endif
               if(this%FF%do_disp_ener)then
                this%FF%NN_C6(l)%mean_inp(j)=this%FF%NN_C6(l)%mean_inp(j)+this%tr(i)%at_desc(v)%desc(j)
               endif
              enddo
              nval(l)=nval(l)+1
             endif
            enddo

           enddo
          enddo

          do l=1,this%FF%nnets
           if(this%FF%do_local_ener)then
            this%FF%NN(l)%mean_inp=this%FF%NN(l)%mean_inp/nval(l)
           endif
           if(this%FF%do_coul_ener)then
            this%FF%NN_charge(l)%mean_inp=this%FF%NN_charge(l)%mean_inp/nval(l)
           endif
           if(this%FF%do_disp_ener)then
            this%FF%NN_C6(l)%mean_inp=this%FF%NN_C6(l)%mean_inp/nval(l)
           endif
          enddo

          do i=1,this%npoints_tr
           do v=1,size(this%tr(i)%at_desc)

            do l=1,this%FF%nnets
             if(this%tr(i)%kind(v).eq.this%FF%types(l) .or. this%FF%types(l).eq.0  )then
              do j=1,size(this%tr(i)%at_desc(v)%desc)
               if(this%FF%do_local_ener)then
                this%FF%NN(l)%sigma_inp(j)=this%FF%NN(l)%sigma_inp(j)+&
                     (this%tr(i)%at_desc(v)%desc(j)-this%FF%NN(l)%mean_inp(j))**2/nval(l)
               endif
               if(this%FF%do_coul_ener)then
                this%FF%NN_charge(l)%sigma_inp(j)=this%FF%NN_charge(l)%sigma_inp(j)+&
                     (this%tr(i)%at_desc(v)%desc(j)-this%FF%NN_charge(l)%mean_inp(j))**2/nval(l)
               endif
               if(this%FF%do_disp_ener)then
                this%FF%NN_C6(l)%sigma_inp(j)=this%FF%NN_C6(l)%sigma_inp(j)+&
                     (this%tr(i)%at_desc(v)%desc(j)-this%FF%NN_C6(l)%mean_inp(j))**2/nval(l)
               endif
              enddo
             endif
            enddo

           enddo
          enddo

          do l=1,this%FF%nnets
           if(this%FF%do_local_ener)then
            this%FF%NN(l)%sigma_inp=sqrt(this%FF%NN(l)%sigma_inp)
           endif
           if(this%FF%do_coul_ener)then
            this%FF%NN_charge(l)%sigma_inp=sqrt(this%FF%NN_charge(l)%sigma_inp)
           endif
           if(this%FF%do_disp_ener)then
            this%FF%NN_C6(l)%sigma_inp=sqrt(this%FF%NN_C6(l)%sigma_inp)
           endif
          enddo

         return
         end subroutine std_sets_inp

         subroutine std_sets_out(this)
         implicit none
         class(chi2)                    :: this
         integer                        :: i,j,l,v
         integer, allocatable           :: nval(:)

         !! Normalize Training Set Output

          this%FF%norm_out=.true.

          this%FF%mean_out=0.0d0
          this%FF%sigma_out=0.0d0

          do i=1,this%npoints_tr
           this%FF%mean_out=this%FF%mean_out+this%tr_val(i)/this%npoints_tr
          enddo

          do i=1,this%npoints_tr
           this%FF%sigma_out=this%FF%sigma_out+(this%tr_val(i)-this%FF%mean_out)**2/this%npoints_tr
          enddo

          this%FF%sigma_out=sqrt(this%FF%sigma_out)

         return
         end subroutine std_sets_out


         subroutine read_sets(this,input_file)
         implicit none
         class(chi2)                    :: this
         integer                        :: i,l,j
         character(len=100)             :: input_file
         character(len=100)             :: tr_file,te_file
         character(len=100)             :: tr_val_file,te_val_file
         character(len=10)              :: kind_desc

          open(13,file=input_file)
          read(13,*) this%npoints_tr,tr_file,tr_val_file,kind_desc
          read(13,*) this%npoints_te,te_file,te_val_file,kind_desc
          close(13)

          allocate(this%tr(this%npoints_tr))
          allocate(this%tr_val(this%npoints_tr))
          allocate(this%te(this%npoints_te))
          allocate(this%te_val(this%npoints_te))

          open(13,file=tr_file)          
          open(14,file=tr_val_file)

          do i=1,this%npoints_tr
           call this%tr(i)%read_extended_xyz(IOid=13)
           read(14,*) this%tr_val(i)
           call this%tr(i)%build_descriptors(kind_desc,Jmax=8,r0=4.1d0)
          enddo

          close(13)          
          close(14)

          open(13,file=te_file)
          open(14,file=te_val_file)

          do i=1,this%npoints_te
           call this%te(i)%read_extended_xyz(IOid=13)
           read(14,*) this%te_val(i)
           call this%te(i)%build_descriptors(kind_desc,Jmax=8,r0=4.1d0)
          enddo

          close(13)          
          close(14)
          
         return
         end subroutine read_sets

         subroutine set_nets(this,nnets,types,ninp,topo)
         use general_types_class
         use random_numbers_class 
         implicit none
         class(chi2)                    :: this
         integer                        :: i,j,nnets
         integer, allocatable           :: types(:),ninp(:)
         type(vector_int), allocatable  :: topo(:)
         double precision, allocatable  :: vec(:)
        
          call init_random_seed()
  
          this%FF%nnets=nnets
          allocate(this%FF%types(this%FF%nnets))         
          this%FF%types=types
          this%FF%nparams=0

          if(this%FF%do_local_ener)then

           allocate(this%FF%NN(this%FF%nnets))
           this%FF%local_nparams=0

           do i=1,this%FF%nnets
            call this%FF%NN(i)%set_topology(ninp(i),size(topo(i)%v),topo(i)%v)
            if(allocated(vec)) deallocate(vec)
            allocate(vec(this%FF%NN(i)%nparams))
            do j=1,this%FF%NN(i)%nparams
             call random_number(vec(j))
            enddo
            this%FF%nparams=this%FF%nparams+size(vec)
            this%FF%local_nparams=this%FF%local_nparams+size(vec)
            call this%FF%NN(i)%set_parameters(vec)
           enddo

          endif

          if(this%FF%do_coul_ener)then

           allocate(this%FF%NN_charge(this%FF%nnets))
           this%FF%coul_nparams=0

           do i=1,this%FF%nnets
            call this%FF%NN_charge(i)%set_topology(ninp(i),size(topo(i)%v),topo(i)%v)
            if(allocated(vec)) deallocate(vec)
            allocate(vec(this%FF%NN_charge(i)%nparams))
            do j=1,this%FF%NN_charge(i)%nparams
             call random_number(vec(j))
            enddo
            this%FF%nparams=this%FF%nparams+size(vec)
            this%FF%coul_nparams=this%FF%coul_nparams+size(vec)
            call this%FF%NN_charge(i)%set_parameters(vec)
           enddo

          endif

          if(this%FF%do_disp_ener)then

           allocate(this%FF%NN_C6(this%FF%nnets))
           this%FF%disp_nparams=0

           do i=1,this%FF%nnets
            call this%FF%NN_C6(i)%set_topology(ninp(i),size(topo(i)%v),topo(i)%v)
            if(allocated(vec)) deallocate(vec)
            allocate(vec(this%FF%NN_C6(i)%nparams))
            do j=1,this%FF%NN_C6(i)%nparams
             call random_number(vec(j))
            enddo
            this%FF%nparams=this%FF%nparams+size(vec)
            this%FF%disp_nparams=this%FF%disp_nparams+size(vec)
            call this%FF%NN_C6(i)%set_parameters(vec)
           enddo

          endif

          call this%map_nets()

         return
         end subroutine set_nets
        
         subroutine map_nets(this)
         implicit none
         class(chi2)                    :: this
         integer                        :: i

          do i=1,size(this%tr)
           this%tr(i)%FF => this%FF
          enddo

          do i=1,size(this%te)
           this%te(i)%FF => this%FF
          enddo

         return
         end subroutine map_nets

         subroutine get_chi2(this,vec,val)
         implicit none
         class(chi2)                    :: this
         double precision               :: val
         double precision, allocatable  :: vec(:),outs(:),inps(:),vec_loc(:),rij(:,:)
         integer                        :: i,offset,npar

          val=0.0d0

          call this%FF%set_ff_params(vec)

          do i=1,size(this%tr)
           rij=this%tr(i)%dist(:,1,:,1)
           call this%tr(i)%FF%get_ener(this%tr(i)%at_desc,this%tr(i)%kind,rij)
           val=val+(this%tr(i)%FF%ener-this%tr_val(i))**2
          enddo

          val=val/size(this%tr)

          if(this%regular_L2) val=val+this%L2*norm2(vec)**2

          val=sqrt(val)

         return
         end subroutine get_chi2

         subroutine get_chi2_grad(this,vec,val,grad)
         implicit none
         class(chi2)                    :: this
         integer                        :: i,j,offset,npar
         double precision               :: val
         double precision, allocatable  :: vec(:),vec_loc(:)
         double precision, allocatable  :: grad(:),grad_loc(:,:),grad2(:)

          val=0.0d0
          grad=0.0d0
          
          call this%FF%set_ff_params(vec)

          do i=1,size(this%tr)

           call this%tr(i)%FF%get_ffgrad(this%tr(i)%at_desc,this%tr(i)%kind)

           do j=1,size(grad)
            grad(j)=grad(j)-2.0d0*(this%tr_val(i)-this%tr(i)%FF%ener)*&
                this%tr(i)%FF%ffgrad(j)
           enddo

           val=val+(this%tr(i)%FF%ener-this%tr_val(i))**2

          enddo

          grad=grad/size(this%tr_val)
          val=val/size(this%tr_val)

          if(this%regular_L2)then
           do j=1,size(grad)
            grad(j)=grad(j)+2.0d0*this%L2*vec(j)
           enddo
          endif

          if(this%regular_L2) val=val+this%L2*norm2(vec)**2
 
          val=sqrt(val)

         return
         end subroutine get_chi2_grad


        end module ff_trainer_class





