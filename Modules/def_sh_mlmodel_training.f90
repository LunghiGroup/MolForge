        module sh_mlmodel_trainer_class
        use mlmodel_trainer_class
        use proj_disp_class
        implicit none
        
        type, extends(mlmodel_trainer)   :: sh_mlmodel_trainer
         logical                         :: covariance=.false.
         type(molecule)                  :: ref_mol
         type(atoms_group), allocatable  :: tr0(:),te0(:)
         double precision, allocatable   :: tr_val0(:,:),te_val0(:,:)
         contains
         procedure                       :: out_results => out_results_sh
         procedure                       :: set_reference
         procedure                       :: set_rot
         procedure                       :: unset_rot
         procedure                       :: get_sph_1
         procedure                       :: get_sph_2
        end type sh_mlmodel_trainer

        contains

        subroutine out_results_sh(this)                 
        implicit none
        class(sh_mlmodel_trainer)         :: this
        integer                           :: i
        double precision                  :: gtens(3,3)

         open(13,file='tr_rmse.dat')

         do i=1,size(this%tr)
          call this%tr(i)%ML%get_output(this%tr(i)%at_desc(1)%desc)           
          if(this%ML%norm_out)then
           this%tr_val(i,:)=this%tr_val(i,:)*this%ML%sigma_out+this%ML%mean_out
           this%tr(i)%ML%output=this%tr(i)%ML%output*this%ML%sigma_out+this%ML%mean_out    
          endif

          gtens(1,:)=this%tr_val(i,1:3)
          gtens(2,:)=this%tr_val(i,4:6)
          gtens(3,:)=this%tr_val(i,7:9)

          gtens=matmul(gtens,transpose(this%tr(i)%at_desc(1)%rot))
          gtens=matmul(this%tr(i)%at_desc(1)%rot,gtens)

          this%tr_val(i,1:3)=gtens(1,:)
          this%tr_val(i,4:6)=gtens(2,:)
          this%tr_val(i,7:9)=gtens(3,:)

          gtens(1,:)=this%tr(i)%ML%output(1:3)
          gtens(2,:)=this%tr(i)%ML%output(4:6)
          gtens(3,:)=this%tr(i)%ML%output(7:9)

          gtens=matmul(gtens,transpose(this%tr(i)%at_desc(1)%rot))
          gtens=matmul(this%tr(i)%at_desc(1)%rot,gtens)

          this%tr(i)%ML%output(1:3)=gtens(1,:)
          this%tr(i)%ML%output(4:6)=gtens(2,:)
          this%tr(i)%ML%output(7:9)=gtens(3,:)

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

          gtens(1,:)=this%te_val(i,1:3)
          gtens(2,:)=this%te_val(i,4:6)
          gtens(3,:)=this%te_val(i,7:9)

          gtens=matmul(gtens,transpose(this%te(i)%at_desc(1)%rot))
          gtens=matmul(this%te(i)%at_desc(1)%rot,gtens)

          this%te_val(i,1:3)=gtens(1,:)
          this%te_val(i,4:6)=gtens(2,:)
          this%te_val(i,7:9)=gtens(3,:)

          gtens(1,:)=this%te(i)%ML%output(1:3)
          gtens(2,:)=this%te(i)%ML%output(4:6)
          gtens(3,:)=this%te(i)%ML%output(7:9)

          gtens=matmul(gtens,transpose(this%te(i)%at_desc(1)%rot))
          gtens=matmul(this%te(i)%at_desc(1)%rot,gtens)

          this%te(i)%ML%output(1:3)=gtens(1,:)
          this%te(i)%ML%output(4:6)=gtens(2,:)
          this%te(i)%ML%output(7:9)=gtens(3,:)

          write(13,*) this%te_val(i,:),this%te(i)%ML%output
         enddo

         close(13)

        return
        end subroutine out_results_sh

        subroutine set_reference(this,filename)
        use atoms_class
        implicit none
        class(sh_mlmodel_trainer)        :: this
        integer                          :: i,s,nats
        double precision,allocatable     :: mass(:),geo0(:,:)
        character(len=100)               :: filename
        character(len=2), allocatable    :: label(:)

         open(13,file=trim(adjustl(filename)))
         read(13,*) nats
         read(13,*)
         allocate(geo0(nats,3))
         allocate(mass(nats))
         allocate(label(nats))
         do i=1,nats
          read(13,*) label(i),geo0(i,1),geo0(i,2),geo0(i,3)
          call get_mass(label(i),mass(i))
         enddo
         close(13)

         call this%ref_mol%def_mol(geo0,mass)

        return
        end subroutine set_reference

        subroutine set_rot(this)
        implicit none
        class(sh_mlmodel_trainer)        :: this
        double precision, allocatable    :: geo(:,:),mass(:)
        double precision                 :: gtens(3,3)
        integer                          :: ii,v,i,s

         allocate(geo(this%ref_mol%nat,3))

         do ii=1,this%npoints_tr

          v=1
          do i=1,this%ref_mol%nat
           do s=1,3
            geo(i,s)=this%tr(ii)%at_desc(1)%desc(v)
            v=v+1
           enddo
          enddo

          call this%ref_mol%def_mol_dist(geo)
          call this%ref_mol%proj_disp()

          v=1
          do i=1,this%ref_mol%nat
           do s=1,3
            this%tr(ii)%at_desc(1)%desc(v)=this%ref_mol%cart_int(i,s)+this%ref_mol%cart_eq(i,s)
            v=v+1
           enddo
          enddo

          this%tr(ii)%at_desc(1)%rot=this%ref_mol%rot

          gtens(1,:)=this%tr_val(ii,1:3)
          gtens(2,:)=this%tr_val(ii,4:6)
          gtens(3,:)=this%tr_val(ii,7:9)

          gtens=matmul(gtens,this%ref_mol%rot)
          gtens=matmul(transpose(this%ref_mol%rot),gtens)

          this%tr_val(ii,1:3)=gtens(1,:)
          this%tr_val(ii,4:6)=gtens(2,:)
          this%tr_val(ii,7:9)=gtens(3,:)

         enddo

         do ii=1,this%npoints_te

          v=1
          do i=1,this%ref_mol%nat
           do s=1,3
            geo(i,s)=this%te(ii)%at_desc(1)%desc(v)
            v=v+1
           enddo
          enddo

          call this%ref_mol%def_mol_dist(geo)
          call this%ref_mol%proj_disp()

          v=1
          do i=1,this%ref_mol%nat
           do s=1,3
            this%te(ii)%at_desc(1)%desc(v)=this%ref_mol%cart_int(i,s)+this%ref_mol%cart_eq(i,s)
            v=v+1
           enddo
          enddo

          this%te(ii)%at_desc(1)%rot=this%ref_mol%rot

          gtens(1,:)=this%te_val(ii,1:3)
          gtens(2,:)=this%te_val(ii,4:6)
          gtens(3,:)=this%te_val(ii,7:9)

          gtens=matmul(gtens,this%ref_mol%rot)
          gtens=matmul(transpose(this%ref_mol%rot),gtens)

          this%te_val(ii,1:3)=gtens(1,:)
          this%te_val(ii,4:6)=gtens(2,:)
          this%te_val(ii,7:9)=gtens(3,:)

         enddo

        return
        end subroutine set_rot

        subroutine unset_rot(this)
        implicit none
        class(sh_mlmodel_trainer)        :: this
        integer                          :: i,s,ii,v,nats
        double precision                 :: gtens(3,3)
        double precision,allocatable     :: geo(:,:)

         do ii=1,this%npoints_tr

          gtens(1,:)=this%tr_val(ii,1:3)
          gtens(2,:)=this%tr_val(ii,4:6)
          gtens(3,:)=this%tr_val(ii,7:9)

          gtens=matmul(gtens,transpose(this%tr(ii)%at_desc(1)%rot))
          gtens=matmul(this%tr(ii)%at_desc(1)%rot,gtens)

          this%tr_val(ii,1:3)=gtens(1,:)
          this%tr_val(ii,4:6)=gtens(2,:)
          this%tr_val(ii,7:9)=gtens(3,:)

         enddo

         do ii=1,this%npoints_te

          gtens(1,:)=this%te_val(ii,1:3)
          gtens(2,:)=this%te_val(ii,4:6)
          gtens(3,:)=this%te_val(ii,7:9)

          gtens=matmul(gtens,transpose(this%te(ii)%at_desc(1)%rot))
          gtens=matmul(this%tr(ii)%at_desc(1)%rot,gtens)

          this%te_val(ii,1:3)=gtens(1,:)
          this%te_val(ii,4:6)=gtens(2,:)
          this%te_val(ii,7:9)=gtens(3,:)

         enddo

        return
        end subroutine unset_rot

        subroutine get_sph_1(this,vec0,step)
        implicit none
        class(sh_mlmodel_trainer)        :: this
        double precision, allocatable    :: vec0(:),vec(:),out0(:)
        double precision, allocatable    :: outp(:),outm(:)
        double precision, allocatable    :: sph(:)
        integer                          :: i,j,k,v
        double precision                 :: step
        character(len=1000)              :: num

         write(*,*) '     ','Computing First-Order Spin-Phonon Coupling'

         vec=vec0
         call this%ML%get_output(vec)
         out0=this%ML%output*this%ML%sigma_out+this%ML%mean_out
                 
         open(1314,file='SPh_1.dat')
         write(1314,*) this%ML%NN%ninput

         do i=1,this%ML%NN%ninput

         write(num,"(I1000)") i
         open(1313,file='SH_'//trim(adjustl(num))//'.dat')

          do k=10,-10,-1
           vec=vec0
           vec(i)=vec(i)+k*step
           call this%ML%get_output(vec)
           this%ML%output=this%ML%output*this%ML%sigma_out+this%ML%mean_out           
           if(k.eq.1) outp=this%ML%output
           if(k.eq.-1) outm=this%ML%output
           write(1313,*) k*step,(this%ML%output(j)-out0(j),j=1,this%ML%NN%noutput)
          enddo

          close(1313)

          sph=(outp-outm)/(2*step)         

          write(1314,*) i,'1'
!          write(1314,*) (sph(v),v=1,this%ML%NN%noutput)
          write(1314,*) (sph(v),v=1,3)
          write(1314,*) (sph(v),v=4,6)
          write(1314,*) (sph(v),v=7,9)

         enddo

         close(1314)

        return
        end subroutine get_sph_1

        subroutine get_sph_2(this,vec0,step)
        implicit none
        class(sh_mlmodel_trainer)        :: this                
        double precision, allocatable    :: vec0(:),vec(:),out0(:)
        double precision, allocatable    :: gradp(:,:),gradm(:,:)
        double precision, allocatable    :: sph(:),outpp(:),outmm(:),outmp(:),outpm(:)
        integer                          :: i1,i2,k1,k2,j,v
        double precision                 :: step
        character(len=10)                :: num1,num2

         write(*,*) '     ','Computing Second-Order Spin-Phonon Coupling'

         vec=vec0
         call this%ML%get_output(vec)
         out0=this%ML%output*this%ML%sigma_out+this%ML%mean_out

         open(1314,file='SPh_2.dat')
         write(1314,*) this%ML%NN%ninput*(this%ML%NN%ninput+1)/2

         do i1=1,this%ML%NN%ninput
          do i2=i1,this%ML%NN%ninput

           write(num1,"(I10)") i1
           write(num2,"(I10)") i2
           open(1313,file='SH_'//trim(adjustl(num1))//'_'//trim(adjustl(num2))//'.dat')

           do k1=5,-5,-1
            do k2=5,-5,-1
             vec=vec0
             vec(i1)=vec(i1)+k1*step
             vec(i2)=vec(i2)+k2*step
             call this%ML%get_output(vec)
             this%ML%output=this%ML%output*this%ML%sigma_out+this%ML%mean_out           
             if(k1.eq.1 .and. k2.eq.-1) outpm=this%ML%output
             if(k1.eq.-1 .and. k2.eq.1) outmp=this%ML%output
             if(k1.eq.-1 .and. k2.eq.-1) outmm=this%ML%output
             if(k1.eq.1 .and. k2.eq.1) outpp=this%ML%output
             write(1313,*) k1*step,k2*step,(this%ML%output(j)-out0(j),j=1,this%ML%NN%noutput)
            enddo
           enddo

           close(1313)

           sph=(outpp+outmm-outmp-outpm)/(4*step*step)         

           write(1314,*) i1,'1',i2,'1'
!           write(1314,*) (sph(v),v=1,this%ML%NN%noutput)
           write(1314,*) (sph(v),v=1,3)
           write(1314,*) (sph(v),v=4,6)
           write(1314,*) (sph(v),v=7,9)

          enddo
         enddo

         close(1314)

        return
        end subroutine get_sph_2

        end module sh_mlmodel_trainer_class

