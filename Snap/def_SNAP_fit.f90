        module SNAP_fit_class

        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
        use LAMMPS
        use kind_class
        use max_class
        use atoms_class
        use lammps_class
        use parameters_class
        use lapack_inverse
        use VdW_class
        implicit none

        type                                         :: SNAP_fit
        type(lammps_obj), allocatable                :: set(:)
        real(kind=dbl),allocatable                   :: matrix(:,:)
        real(kind=dbl),allocatable                   :: target(:)
        real(kind=dbl), allocatable                  :: energies(:)
        real(kind=dbl), allocatable                  :: forces(:)
        real(kind=dbl)                               :: weight
        real(kind=dbl)                               :: lambda
        real(kind=dbl), allocatable                  :: beta(:)
        real(kind=dbl)                               :: s_z
        real(kind=dbl)                               :: cutoff
        integer                                      :: twojmax
        integer                                      :: num_bisp
        integer                                      :: nconfig
        logical, allocatable                         :: coeff_mask(:)
        logical                                      :: flag_energy
        logical                                      :: flag_forces
        character(len=120)                           :: energy_file,forces_file
        character(len=5)                             :: set_type
        
        contains
        
        procedure                                    :: import_set
        procedure                                    :: import_labels
        procedure                                    :: import_coeff
        procedure                                    :: add_sub_VdW
        procedure                                    :: LLS
        procedure                                    :: build_matrix
        procedure                                    :: build_target
        procedure                                    :: get_uncertainty

        end type SNAP_fit

        contains
        
       subroutine import_set(this,file_input,len_file_inp)
       implicit none
       class(SNAP_fit)                           :: this
       integer                                   :: nconfig, i,j,nats,ntypes
       character(len=100),allocatable            :: tmp(:,:)
       character(len=100),dimension(10)          :: tmp_cell_nkinds
       character(len=100)                        :: filename
       character(len=*)                          :: file_input
       integer,intent(in)                        :: len_file_inp
       character(len=150)                        :: file_inp
        
       call import_lammps_obj(nconfig=this%nconfig,file_input=file_input,len_file_inp=len_file_inp,set_array=this%set)

       do i=1,this%nconfig  

         this%set(i)%twojmax = this%twojmax
         this%set(i)%cutoff  = this%cutoff

        end do

        end subroutine import_set

        subroutine import_labels(this)
        implicit none
        
        class(SNAP_fit)                      :: this                    
        integer                              :: tot_atom
        real(kind=dbl)                       :: ave_atom
        real(kind=dbl),allocatable           :: gradients(:)
        integer                              :: i,j

        if (this%flag_energy) then

        allocate(this%energies(this%nconfig))
        open(1,file=trim(this%energy_file))

         do i=1,this%nconfig
           read(1,*) this%energies(i)
         end do
        
        close(1)
        end if

        if (this%flag_forces) then
        
        call get_ave_atoms(this%set,ave_atom,tot_atom)
        allocate(gradients(3*tot_atom))

        open(1,file=trim(this%forces_file))

        do j=1,3*tot_atom
         read(1,*) gradients(j)
        end do

        close(1)

        this%forces=-gradients
        deallocate(gradients)

        end if

        end subroutine import_labels
        
        subroutine add_sub_VdW(this,addsub)
        implicit none
        
        class(SNAP_fit)                                                       :: this
        integer                                                               :: i              
        character(len=3),intent(in)                                           :: addsub
        type(VdW_FF)                                                          :: FF_VdW
        real(kind=dbl),allocatable                                            :: vec(:)
        real(kind=dbl)                                                        :: val
        real(kind=dbl),allocatable                                            :: grad(:)

        do i=1,this%nconfig         

          FF_VdW%frame=this%set(i)
          call FF_VdW%get_fgrad(vec,val,grad)
          
          if ((this%flag_energy).and.(addsub=="sub")) then
            this%energies(i)=this%energies(i)-val
           else if ((this%flag_energy).and.(addsub=="add")) then
           this%energies(i)=this%energies(i)+val
          end if

          if ((this%flag_forces).and.(addsub=="sub")) then
             this%forces((i-1)*3*this%set(i)%nats+1:i*3*this%set(i)%nats)=&
             this%forces((i-1)*3*this%set(i)%nats+1:i*3*this%set(i)%nats)+grad
          else if ((this%flag_forces).and.(addsub=="add")) then
             this%forces((i-1)*3*this%set(i)%nats+1:i*3*this%set(i)%nats)=&
             this%forces((i-1)*3*this%set(i)%nats+1:i*3*this%set(i)%nats)-grad
          end if 
        
        end do

        end subroutine add_sub_VdW
        
        subroutine build_matrix(this)
        implicit none

        class(SNAP_fit)                                                       :: this
        integer                                                               :: counter,count_kinds_ezero,ez_cons_rows
        integer                                                               :: i,j,k,l,m,comp
        integer                                                               :: size_ref
        integer                                                               :: start_cycle,end_cycle,help_counter
        integer                                                               :: start_snap_force
        integer                                                               :: tot_atom
        integer                                                               :: tot_kinds
        real(kind=dbl)                                                        :: ave_atom
        
        do i=1,this%nconfig
         call this%set(i)%initialize()
         call this%set(i)%setup(this%set(i)%nkinds)
         if (this%flag_energy) then
          call this%set(i)%get_desc()
         end if
         if (this%flag_forces) then
          call this%set(i)%get_der_desc()
         end if
         call this%set(i)%finalize()

        end do

        call number_bispec(this%twojmax,this%num_bisp)
        call get_tot_kinds(this%set,tot_kinds)
        call get_ave_atoms(this%set,ave_atom,tot_atom)

        allocate(this%coeff_mask(tot_kinds))
        this%coeff_mask=.true.
        
        if (this%flag_energy) then
         this%coeff_mask(1)=.false.
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Calculation dimension of the matrix

        ez_cons_rows = count(this%coeff_mask)

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         size_ref=size(this%set) + ez_cons_rows

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         size_ref=3*tot_atom + ez_cons_rows

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         size_ref=size(this%set)+3*tot_atom+ez_cons_rows

        end if

        if (this%lambda.ne.0.0) then

         size_ref=size_ref+(this%num_bisp-1)*tot_kinds

        end if

        !!!!!!!!!!!!!!!!!!!!!!!! Allocation of the matrix

        allocate(this%matrix(size_ref,this%num_bisp*tot_kinds))
        this%matrix=0.0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! definizione variabili per fare i cicli

        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         start_cycle=size(this%set)+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         start_cycle=3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=3*tot_atom

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then

         start_cycle=size(this%set)+3*tot_atom+1
         end_cycle=size_ref-ez_cons_rows
         help_counter=size(this%set)+3*tot_atom

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! qui scriviamo proprio la matrice di SNAP        

        if (this%lambda.ne.0.0) then

         counter=0
         do j=start_cycle,end_cycle

          if (mod(j-help_counter,this%num_bisp-1)==1) then
           counter=counter+1
          end if
        
        this%matrix(j,j-help_counter+counter) = sqrt(this%lambda)

         end do

        end if

        if ((this%flag_energy).or.((this%flag_energy).and.(this%flag_forces))) then

         do j=1,size(this%set)
          do i=1,this%set(j)%nats
           do k=1,this%num_bisp

            this%matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) = &
                    this%matrix(j,(this%set(j)%kind(i)-1)*this%num_bisp+k) + this%set(j)%at_desc(i)%desc(k)

           end do
          end do
         end do

        this%matrix(1:size(this%set),:)=this%matrix(1:size(this%set),:)*this%weight

        end if

        counter=0

        do j=1,tot_kinds
         if (this%coeff_mask(j)) then

        counter=counter+1
        this%matrix(size_ref-ez_cons_rows+counter,(j-1)*this%num_bisp+1)=1

         end if
        end do

        if ((this%flag_forces).and.(.not.this%flag_energy)) then
         start_snap_force=0
        else if ((this%flag_forces).and.(this%flag_energy)) then
         start_snap_force=size(this%set)
        end if

        if (this%flag_forces) then

        do j=1,size(this%set)
          do i=1,this%set(j)%nkinds
           do l=1,this%set(j)%nats
            do comp=1,3
             do k=1,this%num_bisp-1

              this%matrix(start_snap_force+(j-1)*this%set(j)%nats*3+(l-1)*3+comp,this%num_bisp*(i-1)+1+k)=&
                this%set(j)%der_at_desc(l)%desc((i-1)*((this%num_bisp-1)*3)+(this%num_bisp-1)*(comp-1)+k)

             end do
            end do
           end do
          end do
         end do

        end if
        
       if (this%flag_energy) then
        do i=1,this%nconfig
         deallocate(this%set(i)%at_desc)
        end do
       end if

       if (this%flag_forces) then
        do i=1,this%nconfig       
         deallocate(this%set(i)%der_at_desc)
        end do
       end if
        
        end subroutine build_matrix
        
        subroutine build_target(this)
        implicit none
        class(SNAP_fit)                                                       :: this
        integer                                                               :: tot_atom
        real(kind=dbl)                                                        :: ave_atom

        call get_ave_atoms(this%set,ave_atom,tot_atom)

        allocate(this%target(size(this%matrix,1)))
        this%target=0.0
        
        if ((this%flag_energy).and.(.not.this%flag_forces)) then

         this%target(1:size(this%set))=this%energies*this%weight

        end if

        if ((this%flag_forces).and.(.not.this%flag_energy)) then

         this%target(1:3*tot_atom)=this%forces

        end if

        if ((this%flag_forces).and.(this%flag_energy)) then
         this%target(1:size(this%set))=this%weight*this%energies
         this%target(size(this%set)+1:size(this%set)+3*tot_atom)=this%forces

        end if

        end subroutine build_target

        subroutine LLS(this)
        implicit none
        
        class(SNAP_fit)                                                       :: this
        real(kind=dbl), dimension(:,:), allocatable                           :: C,F
        real(kind=dbl),dimension(:,:), allocatable                            :: temp_matrix
        real(kind=dbl), dimension(:), allocatable                             :: temp_target
        real(kind=dbl), dimension(:), allocatable                             :: ML_forces
        real(kind=dbl),dimension(:),allocatable                               :: ML_energies
        character(len=1)                                                      :: TRANS
        integer                                                               :: i
        integer                                                               :: M
        integer                                                               :: MF,N,NRHS,LDA,LDB,LWORK,INFO
        real(kind=dbl)                                                        :: ave_atom
        integer                                                               :: tot_atom
        integer                                                               :: start_snap_force
        real(kind=dbl),dimension(:),allocatable                               :: WORK
        
        if (this%set_type=="TRAIN") then
         TRANS='N'
         M=size(this%matrix,1)
         LDA=size(this%matrix,1)
         N=size(this%matrix,2)
         NRHS=1
         LDB=max(M,N)
         LWORK=2*min(M,N)
         allocate(WORK(LWORK))
        
         allocate(temp_matrix(size(this%matrix,1),size(this%matrix,2)))
         allocate(temp_target(size(this%matrix,1)))
         temp_matrix=this%matrix
         temp_target=this%target
         call dgels(TRANS,M,N,NRHS,temp_matrix,LDA,temp_target,LDB,WORK,LWORK,INFO)
         
         deallocate(WORK)
         deallocate(temp_matrix)
         if (info.ne.0) then
          write(*,*) 'Convergence issues: could not solve the linear least square problem'
         end if

         open(11,file='snapcoeff_energy',action='write')
          do i=1,size(this%matrix,2)
          write(11,*) temp_target(i)
          end do
         close(11)


         allocate(this%beta(size(this%matrix,2)))
         this%beta=temp_target(1:size(this%matrix,2))

         else

         call this%import_coeff

        end if

        if (this%flag_energy) then
         allocate(C(size(this%set),size(this%matrix,2)))
         
         C=this%matrix(1:size(this%set),:)/this%weight
         allocate(ML_energies(size(this%set)))
          ML_energies=matmul(C,this%beta)
         deallocate(C)

        end if

        if (this%flag_forces) then

         call get_ave_atoms(this%set,ave_atom,tot_atom) 
         
         if (this%flag_energy) then
           start_snap_force=size(this%set)
         else 
           start_snap_force=0
         end if
        
        allocate(F(3*tot_atom,size(this%matrix,2)))
        F=this%matrix(start_snap_force+1:start_snap_force+3*tot_atom,:)
        
         allocate(ML_forces(3*tot_atom))
         ML_forces=matmul(F,this%beta)         
         deallocate(F)
        
        end if
        
        call get_ave_atoms(this%set,ave_atom,tot_atom)

        if (this%flag_energy) then

        open(11,file='energy_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML energy','             ','DFT energy','           ','Error'

        do i=1,size(this%set)
         write(11,*) ' ', i, ML_energies(i), this%energies(i),ML_energies(i)- this%energies(i)
        end do

         write(11,*)'# RMS=',(1/ave_atom)*(sqrt(sum((ML_energies-this%energies)**2)/&
                 size(this%set))), 'kcal/mol/atom'
        close(11)

        end if

        if (this%flag_forces) then

        open(11,file='forces_rms_fpp.dat')

        write(11,*) '#','      ','Num config','       ','ML force','             ','DFT force','           ','Error'

        do i=1,3*tot_atom
         write(11,*) ' ', i, ML_forces(i), this%forces(i),ML_forces(i)-this%forces(i)
        end do

        write(11,*)'# RMS=',(sqrt(sum((ML_forces(:)-this%forces(:))**2)&
        /(dble(3*tot_atom)))), 'kcal/mol/angstrom'

        close(11)

        end if

        end subroutine LLS

        subroutine import_coeff(this)
        implicit none
        class(SNAP_fit)                                         :: this
        integer                                                 :: i

        allocate(this%beta(size(this%matrix,2)))

        open(11,file='snapcoeff_energy',action='read')
        do i=1,size(this%matrix,2)
         read(11,*) this%beta(i)
        end do
        close(11)

        end subroutine

        subroutine get_uncertainty(this,frame,calc_sz_flag,error)
        implicit none
        class(SNAP_fit)                                         :: this
        type(lammps_obj)                                        :: frame
        logical, intent(in)                                     :: calc_sz_flag
        real(kind=dbl),allocatable                              :: tmp(:)
        real(kind=dbl),allocatable                              :: y_pred(:)
        real(kind=dbl),allocatable                              :: error_array(:)
        real(kind=dbl),allocatable                              :: X(:,:),X_t(:,:)
        real(kind=dbl),allocatable                              :: U_matrix(:,:)
        real(kind=dbl),allocatable                              :: x_new(:,:)
        real(kind=dbl),allocatable                              :: K_mat(:,:)
        real(kind=dbl)                                          :: NUMERATOR,DENOMINATOR
        real(kind=dbl),intent(out)                              :: error
        integer                                                 :: i,j,k,l,comp
        integer                                                 :: N,size_ref
        integer                                                 :: start

        size_ref=size(this%matrix,1)
        N=size(this%matrix,2)

        allocate(X(size_ref,N))
        allocate(X_t(N,size_ref))
        allocate(U_matrix(N,N))
        allocate(y_pred(size_ref))
        
        if ((this%flag_energy).and.(.not.this%flag_forces)) then
         allocate(x_new(N,1))
         allocate(K_mat(N,1))
         allocate(tmp(1))
         allocate(error_array(1))
        else if ((this%flag_forces).and.(.not.this%flag_energy)) then
         allocate(x_new(N,frame%nats*3))
         allocate (K_mat(N,frame%nats*3))
         allocate (tmp(frame%nats*3))
         allocate(error_array(3*frame%nats))
        else if ((this%flag_forces).and.(this%flag_energy)) then
         allocate(x_new(N,1+frame%nats*3))
         allocate (K_mat(N,frame%nats*3+1))
         allocate (tmp(frame%nats*3+1))
         allocate(error_array(3*frame%nats+1))
        end if
        
        x_new=0.0d0
        start=0
        K_mat=0.0d0
        tmp=0.0d0
        NUMERATOR=0.0d0
        
        X=this%matrix
        X_t=transpose(X)
        U_matrix=matmul(X_t,X)
        call mat_inv(U_matrix,N)
        

        y_pred=matmul(this%matrix,this%beta)
        
        do j=1,size_ref
         NUMERATOR=NUMERATOR+(y_pred(j)-this%target(j))**2
        end do

        DENOMINATOR=dble(size_ref-N-1)
        this%s_z=dsqrt(NUMERATOR/DENOMINATOR)
        
        !!!!!!!!!!!!!!!!!!!!!!!fino a qui va fatto solo al primo ciclo di MD
        
        call frame%initialize()
        call frame%setup(frame%nkinds)

        if (this%flag_energy) then
         call frame%get_desc()
         
         start=1

         do i=1,frame%nats
          do k=1,this%num_bisp

           x_new((frame%kind(i)-1)*this%num_bisp+k,1) = x_new((frame%kind(i)-1)*this%num_bisp+k,1)&
           + frame%at_desc(i)%desc(k)

          end do
         end do

        end if

        if (this%flag_forces) then
         call frame%get_der_desc()
         
         do i=1,frame%nkinds
          do l=1,frame%nats
           do comp=1,3
            do k=1,this%num_bisp-1

             x_new((i-1)*(this%num_bisp)+k+1,(l-1)*3+comp+start) = frame&
             %der_at_desc(l)%desc((i-1)*((this%num_bisp-1)*3)+(this%num_bisp-1)*(comp-1)+k)

            end do
           end do
          end do
         end do
        end if
        
        call frame%finalize() 

        K_mat=matmul(U_matrix,x_new)

        if (this%flag_energy) then

         do i=1,N
          tmp(1)=tmp(1)+K_mat(i,1)*x_new(i,1)
         end do
        
        end if

        if (this%flag_forces) then

         do j=1,frame%nats*3+start
          do i=1,N
           tmp(j)=tmp(j)+K_mat(i,j)*x_new(i,j)
          end do
         end do

        end if

        error_array=this%s_z*dsqrt(1+tmp)
        error=maxval(error_array)
        
        deallocate(K_mat,x_new)
        deallocate(X)
        deallocate(X_t)
        deallocate(U_matrix)
        deallocate(y_pred)
        deallocate(tmp)
        deallocate(error_array)
        
        do i=1,frame%nats
         deallocate(frame%at_desc(i)%desc)
         deallocate(frame%der_at_desc(i)%desc)
        end do

        deallocate(frame%at_desc)
        deallocate(frame%der_at_desc)

        end subroutine get_uncertainty

        end module SNAP_fit_class
