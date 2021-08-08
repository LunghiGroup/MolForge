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
!         procedure                       :: read_sets => read_sets_sh
!         procedure                       :: out_results_sh
!         procedure                       :: rot_sets
!         procedure                       :: rotinv_sets
         procedure                       :: get_sph_1
         procedure                       :: get_sph_2
        end type sh_mlmodel_trainer

        contains

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

         do i=1,this%ML%NN%ninput

         write(num,"(I1000)") i
         open(1313,file='SH_'//trim(adjustl(num))//'.dat')

          do k=5,-5,-1
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
          write(1314,*) (sph(v),v=1,this%ML%NN%noutput)

         enddo

         close(1314)

        return
        end subroutine get_sph_1

        subroutine get_sph_2(this,vec0,step)
        implicit none
        class(sh_mlmodel_trainer)        :: this                
        double precision, allocatable    :: vec0(:),vec(:),out0(:)
        double precision, allocatable    :: gradp(:,:),gradm(:,:)
        double precision, allocatable    :: sph2(:,:)
        integer                          :: i1,i2,k1,k2,j,v
        double precision                 :: step
        character(len=1000)              :: num1,num2

         write(*,*) '     ','Computing Second-Order Spin-Phonon Coupling'

         vec=vec0
         call this%ML%get_output(vec)
         out0=this%ML%output*this%ML%sigma_out+this%ML%mean_out

!         open(1313,file='SPh_2.dat')
!         write(1313,*) this%ML%NN%ninput*(this%ML%NN%ninput+1)/2

!         do i=1,this%ML%NN%ninput

!          do k=1,-1,-2
!           vec=vec0
!           vec(i)=vec(i)+k*step
!           call this%ML%get_output(vec)
!           call this%ML%NN%get_gradinp()
!           if(k.eq.1) gradp=this%ML%NN%gradinp
!           if(k.eq.-1) gradm=this%ML%NN%gradinp
!          enddo
          
!          sph2=(gradp-gradm)/(2*step)         

!          do j=i,this%ML%NN%ninput
!           write(1313,*) i,'1',j,'1'
!           write(1313,*) (sph2(j,v),v=1,this%ML%NN%noutput)
!          enddo

!         enddo

!         close(1313)

         do i1=1,this%ML%NN%ninput
          do i2=i1,this%ML%NN%ninput

           write(num1,"(I1000)") i1
           write(num2,"(I1000)") i2
           open(1313,file='SH_'//trim(adjustl(num1))//'_'//trim(adjustl(num2))//'.dat')

           do k1=5,-5,-1
            do k2=5,-5,-1
             vec=vec0
             vec(i1)=vec(i1)+k1*step
             vec(i2)=vec(i2)+k2*step
             call this%ML%get_output(vec)
             this%ML%output=this%ML%output*this%ML%sigma_out+this%ML%mean_out           
             write(1313,*) k1*step,k2*step,(this%ML%output(j)-out0(j),j=1,this%ML%NN%noutput)
            enddo
           enddo

           close(1313)

          enddo
         enddo

        return
        end subroutine get_sph_2

        end module sh_mlmodel_trainer_class

