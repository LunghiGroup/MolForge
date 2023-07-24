        module liuville_class
        use mpi
        use blacs_utils
        use general_types_class
        implicit none
                
        type :: liuville_space
         integer                            :: Hdim=0
         integer                            :: Ldim=0
         integer, allocatable               :: Lbasis(:,:)
         double precision, allocatable      :: Ener(:)
         type(dist_cmplx_mat)               :: H0
         type(dist_cmplx_mat)               :: H
         type(dist_cmplx_mat)               :: U
         type(dist_cmplx_mat)               :: R
         type(dist_cmplx_mat)               :: rho
         type(dist_cmplx_mat), allocatable  :: QMOP(:)
         contains
         procedure   ::  diag_Hmat
         procedure   ::  to_eigenbasis
         procedure   ::  make_unitary_propagator
         procedure   ::  make_basis_L
         procedure   ::  make_rho
         procedure   ::  diag_lindbladian
         procedure   ::  make_R21
         procedure   ::  make_R22
         procedure   ::  make_R41
         procedure   ::  propagate
         procedure   ::  make_expval
         procedure   ::  dump_expvals
         procedure   ::  rot_rho
        end type liuville_space

        contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   DIAGONALIZE THE HAMILTONIAN MATRIX
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine diag_Hmat(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        implicit none
        class(liuville_space)            :: this
        integer                         :: t1,t2,rate,i

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Diagonalizing the Spin Hamiltonian Matrix'
          flush(6)
         endif

         allocate(this%Ener(this%Hdim)) 

         call this%H%set(this%Hdim,this%Hdim,NB,MB)       

         call pzdiag(this%Hdim,this%H0%mat,this%Ener,this%H%mat, &
                      this%H0%desc,this%H%desc)

         this%Ener=this%Ener-this%Ener(1)

         if(mpi_id.eq.0)then
          open(11,file='eigenval.dat')
          do i=1,this%Hdim
           write(11,*) this%Ener(i)
          enddo
          close(11)
         endif

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine diag_Hmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   TRANSFOR A MATRIX FROM THE COMPUTATIONAL BASIS TO THE ONE OF THE H EIGENVECTORS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine to_eigenbasis(this,mat)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(liuville_space)      :: this
        type(dist_cmplx_mat)      :: AA,mat

         call AA%set(this%Hdim,this%Hdim,NB,MB)
         AA%mat=(0.0d0,0.0d0)

         call pzgemm('C','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),this%H%mat,1,1,this%H%desc,mat%mat,&
                     1,1,mat%desc,(0.0d0,0.0d0),AA%mat,1,1,AA%desc)


         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),AA%mat,1,1,AA%desc,this%H%mat,1,1,this%H%desc,&
                     (0.0d0,0.0d0),mat%mat,1,1,mat%desc) 

         call AA%dealloc()

        return
        end subroutine to_eigenbasis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE PROPAGATOR
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_unitary_propagator(this,step)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms
        implicit none
        class(liuville_space)    :: this
        double precision        :: step
        integer                 :: t1,t2,rate,v

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building unitary propagator operator'
          flush(6)
         endif
 
         if(allocated(this%U%mat)) call this%U%dealloc()
         allocate(this%U%mat(1,this%Hdim))
         this%U%mat=(0.0d0,0.0d0)

         do v=1,this%Hdim
          this%U%mat(1,v)=-step*cmplx(0.0d0,1.0d0,8)*this%Ener(v)*2*acos(-1.0d0)/hplank
          this%U%mat(1,v)=exp(this%U%mat(1,v))
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_unitary_propagator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   MAKE LIUVILLIAN BASIS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_basis_L(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(liuville_space)            :: this
        integer                         :: k,i,j

         this%Ldim=this%Hdim**2
         allocate(this%Lbasis(this%Ldim,2))

         k=1
         do i=1,this%Hdim
          this%Lbasis(k,1)=i
          this%Lbasis(k,2)=i
          k=k+1
         enddo
         do i=1,this%Hdim
          do j=1,this%Hdim
           if(i.eq.j)cycle
           this%Lbasis(k,1)=i
           this%Lbasis(k,2)=j
           k=k+1
          enddo
         enddo

         call this%R%set(this%Ldim,this%Ldim,NB,MB)

        return
        end subroutine make_basis_L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   MAKE DENSITY MATRIX
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_rho(this,type_rho0,temp,rho_restart_file)
        use mpi
        use mpi_utils
        use blacs_utils
        use sparse_class
        use units_parms 
        implicit none
        class(liuville_space)     :: this
        type(dist_cmplx_mat)      :: rho
        type(csr_mat_cmplx)       :: sprho
        double precision          :: temp,part_funct,valre,valim
        complex(8)                :: sum,val
        complex(8),allocatable    :: rho_i(:,:)
        integer                   :: i,ii,jj,j,k,v,l
        integer                   :: indxl2g,kpt,size_block
        integer                   :: t1,t2,rate
        character(len=100)        :: type_rho0
        character(len=20)         :: rho_restart_file

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building starting density matrix'
          flush(6)
         endif

         if(allocated(this%rho%mat)) call this%rho%dealloc()
         call this%rho%set(this%Hdim,this%Hdim,NB,MB)
         this%rho%mat=(0.0d0,0.0d0)

         select case (type_rho0)

          case ('READ_FILE')
        
           if(mpi_id.eq.0) open(121,file=rho_restart_file)
           do i=1,this%Hdim
            do j=1,this%Hdim            
             if(mpi_id.eq.0) read(121,*) ii,jj,valre,valim
             if(mpi_id.eq.0) val=cmplx(valre,valim,8)
             call mpi_bcast(val,1,mpi_double_complex,0,mpi_blacs_world,err)
             call pzelset(this%rho%mat,i,j,this%rho%desc,val)
            enddo
           enddo
           if(mpi_id.eq.0) close(121)

          case ('THERMAL_POPULATION') 
                      
           do i=1,this%Hdim
            call pzelset(this%rho%mat,i,i,this%rho%desc, &
                      cmplx(exp(-this%Ener(i)/(temp*kboltz)),0.0d0,8))
           enddo

           part_funct=0.0d0

           do ii=1,size(this%rho%mat,1)
            do jj=1,size(this%rho%mat,2)
             i=indxl2g(ii,NB,myrow,0,nprow)
             j=indxl2g(jj,MB,mycol,0,npcol)
             if(i.eq.j) part_funct=part_funct+this%rho%mat(ii,jj)
            enddo
           enddo

           call mpi_allreduce(part_funct,part_funct,1,mpi_double_precision,mpi_sum,mpi_blacs_world,err)

           this%rho%mat=this%rho%mat/part_funct

          case ('FULLY_POLARIZED')

           do ii=1,size(this%rho%mat,1)
            do jj=1,size(this%rho%mat,2)
             i=indxl2g(ii,NB,myrow,0,nprow)
             j=indxl2g(jj,MB,mycol,0,npcol)
             if(i.eq.1 .and. j.eq.1) this%rho%mat(ii,jj)=(1.0d0,0.0d0)
            enddo
           enddo

           call this%to_eigenbasis(this%rho)

         end select

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD THE LIMBLADIAN PROPAGATOR
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine diag_lindbladian(this,step)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_diag_asimm_cmplx
        use scalapack_inv
        use units_parms 
        implicit none
        class(liuville_space)          :: this
        type(dist_cmplx_mat)          :: AA,BB
        double complex, allocatable   :: Rval(:)
        double precision, allocatable :: rates(:)
        double precision              :: step,max_ener,min_ener
        integer                       :: t1,t2,rate,l,ii,jj,indxl2g
        character(len=30)             :: filename="R.dat"

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the lindbladian propagator by full diagonalization'
          flush(6)
         endif

         filename="R.dat"
         call this%R%print_mat(filename)

         call AA%set(this%Ldim,this%Ldim,NB,MB)          
         call BB%set(this%Ldim,this%Ldim,NB,MB)

         allocate(Rval(this%Ldim))

         call pzdiag2(this%Ldim,this%R,Rval,AA,BB)
         BB%mat=AA%mat
         call pzgeinv(this%Ldim,BB)
         filename="Reig.dat"
         call BB%print_mat(filename)

         if(mpi_id.eq.0) then

          write(*,*) '     Limbladian Matrix Eigenvalues (1/ms):'

          allocate(rates(this%Ldim))
          rates=abs(dble(Rval))
          call  order_array(rates)

          do l=1,this%Ldim
           write(*,*) '          ',l,rates(l)*1.0d9
          enddo

          deallocate(rates)

         endif

      ! Build Pop Propagator R= R exp(Rval) R^{\cross}

         Rval=exp(Rval*step)

         do ii=1,size(AA%mat,1)
          do jj=1,size(AA%mat,2)
           l=indxl2g(jj,MB,mycol,0,npcol)           
           AA%mat(ii,jj)=AA%mat(ii,jj)*Rval(l)
          enddo
         enddo

         call this%R%set(this%Ldim,this%Ldim,NB,MB)
         this%R%mat=(0.0d0,0.0d0)

         call pzgemm('N','N',this%Ldim,this%Ldim,this%Ldim,&
                     (1.0d0,0.0d0),AA%mat,1,1,AA%desc,BB%mat,&
                     1,1,BB%desc,(0.0d0,0.0d0),this%R%mat,1,1,this%R%desc)

         deallocate(Rval)
         call AA%dealloc()
         call BB%dealloc()
         
         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine diag_lindbladian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD SECOND-ORDER LIMBLADIAN OPERATOR WITH LINEAR SYSTEM-BATH COUPLING
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R21(this,Vmat,temp,freq,lw,type_smear)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms 
        implicit none
        class(liuville_space)         :: this
        double precision             :: Gf,DEner,freq,lw,temp,prefc,secular
        complex(8), allocatable      :: Vmat(:,:)       
        integer                      :: t1,t2,rate,indxl2g,type_smear
        integer                      :: ii,jj,kk,l,l2,la,lb,lc,ld

         if (.not.allocated(this%R%mat)) call this%R%set(this%Ldim,this%Ldim,NB,MB)

         prefc=pi*pi/hplank 

         do ii=1,size(this%R%mat,1)
          do jj=1,size(this%R%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)
 
           la=this%Lbasis(l,1)
           lb=this%Lbasis(l,2)

           lc=this%Lbasis(l2,1) 
           ld=this%Lbasis(l2,2) 

           Secular=this%Ener(la)-this%Ener(lc)+this%Ener(ld)-this%Ener(lb)
     
           if( abs(Secular).lt.1.0d-6 )then              

            Gf=0.0d0
            if( (this%Ener(lb)-this%Ener(ld)).gt.0.0d0 ) then
             DEner=this%Ener(lb)-this%Ener(ld)-freq
             Gf=Gf+bose(temp,freq)*delta(type_smear,DEner,lw)
            else
             DEner=this%Ener(lb)-this%Ener(ld)+freq
             Gf=Gf+(bose(temp,freq)+1.0d0)*delta(type_smear,DEner,lw)
            endif
            this%R%mat(ii,jj)=this%R%mat(ii,jj)+Vmat(la,lc)*conjg(Vmat(lb,ld))*Gf*prefc
 
            Gf=0.0d0
            if( (this%Ener(la)-this%Ener(lc)).gt.0.0d0) then
             DEner=this%Ener(la)-this%Ener(lc)-freq
             Gf=Gf+bose(temp,freq)*delta(type_smear,DEner,lw)
            else
             DEner=this%Ener(la)-this%Ener(lc)+freq
             Gf=Gf+(bose(temp,freq)+1.0d0)*delta(type_smear,DEner,lw)
            endif
            this%R%mat(ii,jj)=this%R%mat(ii,jj)+Vmat(la,lc)*conjg(Vmat(lb,ld))*Gf*prefc

            if(ld.eq.lb)then
             do kk=1,this%Hdim
              Gf=0.0d0
              if( (this%Ener(kk)-this%Ener(lc)).gt.0.0d0) then
               DEner=this%Ener(kk)-this%Ener(lc)-freq
               Gf=Gf+bose(temp,freq)*delta(type_smear,DEner,lw)
              else
               DEner=this%Ener(kk)-this%Ener(lc)+freq
               Gf=Gf+(bose(temp,freq)+1.0d0)*delta(type_smear,DEner,lw)
              endif
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Vmat(kk,la))*Vmat(kk,lc)*Gf*prefc 
             enddo
            endif

            if(lc.eq.la)then
             do kk=1,this%Hdim
              Gf=0.0d0
              if( (this%Ener(kk)-this%Ener(ld)).gt.0.0d0) then
               DEner=this%Ener(kk)-this%Ener(ld)-freq
               Gf=Gf+bose(temp,freq)*delta(type_smear,DEner,lw)
              else
               DEner=this%Ener(kk)-this%Ener(ld)+freq
               Gf=Gf+(bose(temp,freq)+1.0d0)*delta(type_smear,DEner,lw)
              endif
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-Vmat(kk,lb)*conjg(Vmat(kk,ld))*Gf*prefc
             enddo
            endif

           endif

          enddo  !! jj
         enddo !! ii

        return
        end subroutine make_R21

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD SECOND-ORDER LIMBLADIAN OPERATOR WITH QUADRATIC SYSTEM-BATH COUPLING
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R22(this,Vmat,temp,freq1,freq2,lw1,lw2,type_smear)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms 
        implicit none
        class(liuville_space)        :: this
        double precision             :: Gf,DEner,freq1,freq2,lw,lw1,lw2,temp,prefc,secular
        complex(8), allocatable      :: Vmat(:,:)       
        integer                      :: t1,t2,rate,indxl2g,type_smear
        integer                      :: ii,jj,kk,l,l2,la,lb,lc,ld
     
         if (.not.allocated(this%R%mat)) call this%R%set(this%Ldim,this%Ldim,NB,MB)

         prefc=pi*pi/hplank/2.0d0
         lw=lw1+lw2

         do ii=1,size(this%R%mat,1)
          do jj=1,size(this%R%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           la=this%Lbasis(l,1)
           lb=this%Lbasis(l,2)

           lc=this%Lbasis(l2,1)
           ld=this%Lbasis(l2,2)

           Secular=this%Ener(la)-this%Ener(lc)+this%Ener(ld)-this%Ener(lb)

           if( abs(Secular).lt.1.0e-6 )then              

            Gf=0.0d0
            DEner=this%Ener(l)-this%Ener(ld)-freq1-freq2
            Gf=bose(temp,freq1)*bose(temp,freq2)*delta(type_smear,DEner,lw)

            ! double conting something here?

            DEner=this%Ener(lb)-this%Ener(ld)+freq1-freq2
            Gf=Gf+(bose(temp,freq1)+1.0d0)*bose(temp,freq2)*delta(type_smear,DEner,lw)
 
            DEner=this%Ener(lb)-this%Ener(ld)-freq1+freq2
            Gf=Gf+bose(temp,freq1)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)
 
            !

            DEner=this%Ener(lb)-this%Ener(ld)+freq1+freq2  
            Gf=Gf+(bose(temp,freq1)+1.0d0)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)
 
            this%R%mat(ii,jj)=this%R%mat(ii,jj)+Vmat(la,lc)*conjg(Vmat(lb,ld))*Gf*prefc
 
            Gf=0.0d0
 
            DEner=this%Ener(la)-this%Ener(lc)-freq1-freq2  
            Gf=bose(temp,freq1)*bose(temp,freq2)*delta(type_smear,DEner,lw)
 
            ! double conting something here?

            DEner=this%Ener(la)-this%Ener(lc)+freq1-freq2
            Gf=Gf+(bose(temp,freq1)+1.0d0)*bose(temp,freq2)*delta(type_smear,DEner,lw)
  
            DEner=this%Ener(la)-this%Ener(lc)-freq1+freq2
            Gf=Gf+bose(temp,freq1)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)
 
            !

            DEner=this%Ener(la)-this%Ener(lc)+freq1+freq2  
            Gf=Gf+(bose(temp,freq1)+1.0d0)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)
 
            this%R%mat(ii,jj)=this%R%mat(ii,jj)+Vmat(la,lc)*conjg(Vmat(lb,ld))*Gf*prefc

            if(ld.eq.lb)then
             do kk=1,this%Hdim
              Gf=0.0d0
              DEner=this%Ener(kk)-this%Ener(lc)-freq1-freq2  
              Gf=bose(temp,freq1)*bose(temp,freq2)*delta(type_smear,DEner,lw)

              !

              DEner=this%Ener(kk)-this%Ener(lc)+freq1-freq2
              Gf=Gf+(bose(temp,freq1)+1.0d0)*bose(temp,freq2)*delta(type_smear,DEner,lw)

              DEner=this%Ener(kk)-this%Ener(lc)-freq1+freq2
              Gf=Gf+bose(temp,freq1)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)

              !

              DEner=this%Ener(jj)-this%Ener(lc)+freq1+freq2
              Gf=Gf+(bose(temp,freq1)+1.0d0)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-Vmat(kk,lc)*conjg(Vmat(kk,la))*Gf*prefc

             enddo
            endif

            if(lc.eq.la)then
             do kk=1,this%Hdim
              Gf=0.0d0
              DEner=this%Ener(kk)-this%Ener(ld)-freq1-freq2
              Gf=bose(temp,freq1)*bose(temp,freq2)*delta(type_smear,DEner,lw)

              !

              DEner=this%Ener(kk)-this%Ener(ld)+freq1-freq2
              Gf=Gf+(bose(temp,freq1)+1.0d0)*bose(temp,freq2)*delta(type_smear,DEner,lw)

              DEner=this%Ener(kk)-this%Ener(ld)-freq1+freq2
              Gf=Gf+bose(temp,freq1)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)

              !

              DEner=this%Ener(kk)-this%Ener(ld)+freq1+freq2
              Gf=Gf+(bose(temp,freq1)+1.0d0)*(bose(temp,freq2)+1.0d0)*delta(type_smear,DEner,lw)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-Vmat(kk,lb)*conjg(Vmat(kk,ld))*Gf*prefc

             enddo
            endif

           endif
             
          enddo  !! jj
         enddo !! ii

        return
        end subroutine make_R22

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   BUILD FOURTH-ORDER LIMBLADIAN OPERATOR WITH LINEAR SYSTEM-BATH COUPLING
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_R41(this,V1mat,V2mat,temp,freq1,freq2,lw1,lw2,type_smear)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms 
        implicit none
        class(liuville_space)        :: this
        double precision             :: Gf,DEner,freq1,freq2,lw,lw1,lw2,temp,prefc,secular
        double complex, allocatable  :: V1mat(:,:),V2mat(:,:)
        double complex, allocatable  :: Rabp(:,:),Rabm(:,:),Rbap(:,:),Rbam(:,:)
        integer                      :: t1,t2,rate,indxl2g,type_smear
        integer                      :: ii,jj,kk,l1,l2,la,lb,lc,ld
        double complex               :: val
              
         if (.not.allocated(this%R%mat)) call this%R%set(this%Ldim,this%Ldim,NB,MB)
         
         allocate(Rabp(this%Hdim,this%Hdim))
         allocate(Rabm(this%Hdim,this%Hdim))
         allocate(Rbap(this%Hdim,this%Hdim))
         allocate(Rbam(this%Hdim,this%Hdim))

         Rabp=(0.0d0,0.0d0)
         Rabm=(0.0d0,0.0d0)
         Rbap=(0.0d0,0.0d0)
         Rbam=(0.0d0,0.0d0)

         prefc=pi*pi/hplank
         lw=lw1+lw2

         do ii=1,this%Hdim
          do jj=1,this%Hdim          

           do kk=1,this%Hdim

             Rabp(ii,jj)=Rabp(ii,jj)+V1mat(ii,kk)*V2mat(kk,jj)&
                  /(this%Ener(kk)-this%Ener(jj)+freq2-cmplx(0.0d0,1.0d0,8)*lw2)

             Rabm(ii,jj)=Rabm(ii,jj)+V1mat(ii,kk)*V2mat(kk,jj)&
                  /(this%Ener(kk)-this%Ener(jj)-freq2-cmplx(0.0d0,1.0d0,8)*lw2)

             Rbap(ii,jj)=Rbap(ii,jj)+V2mat(ii,kk)*V1mat(kk,jj)&
                  /(this%Ener(kk)-this%Ener(jj)+freq1-cmplx(0.0d0,1.0d0,8)*lw1)

             Rbam(ii,jj)=Rbam(ii,jj)+V2mat(ii,kk)*V1mat(kk,jj)&
                  /(this%Ener(kk)-this%Ener(jj)-freq1-cmplx(0.0d0,1.0d0,8)*lw1)

           enddo ! kk
                
          enddo
         enddo

         do ii=1,size(this%R%mat,1)
          do jj=1,size(this%R%mat,2)         

           l1=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           la=this%Lbasis(l1,1)
           lb=this%Lbasis(l1,2)
           lc=this%Lbasis(l2,1)
           ld=this%Lbasis(l2,2)                 

           Secular=this%Ener(la)-this%Ener(lc)+this%Ener(ld)-this%Ener(lb)


           if( abs(Secular).lt.1.0e-6 )then              

            ! +- and -+ processes
           
            DEner=this%Ener(la)-this%Ener(lc)-freq1+freq2
            Gf=bose(temp,freq1)*(bose(temp,freq2)+1)*delta(type_smear,DEner,lw1)

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabp(lb,ld))*Rabp(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabp(lb,ld))*Rbam(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbam(lb,ld))*Rabp(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbam(lb,ld))*Rbam(la,lc)*Gf*prefc

            DEner=this%Ener(la)-this%Ener(lc)+freq1-freq2
            Gf=bose(temp,freq2)*(bose(temp,freq1)+1)*delta(type_smear,DEner,lw1)

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabm(lb,ld))*Rabm(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabm(lb,ld))*Rbap(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbap(lb,ld))*Rabm(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbap(lb,ld))*Rbap(la,lc)*Gf*prefc

            ! ++ processes

            DEner=this%Ener(la)-this%Ener(lc)-freq1-freq2
            Gf=bose(temp,freq1)*bose(temp,freq2)*delta(type_smear,DEner,lw1)

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabm(lb,ld))*Rabm(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbam(lb,ld))*Rbam(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbam(lb,ld))*Rabm(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabm(lb,ld))*Rbam(la,lc)*Gf*prefc

            ! -- processes

            DEner=this%Ener(la)-this%Ener(lc)+freq1+freq2
            Gf=(bose(temp,freq1)+1)*(bose(temp,freq2)+1)*delta(type_smear,DEner,lw1)

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabp(lb,ld))*Rabp(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbap(lb,ld))*Rbap(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rbap(lb,ld))*Rabp(la,lc)*Gf*prefc

            this%R%mat(ii,jj)=this%R%mat(ii,jj)+conjg(Rabp(lb,ld))*Rbap(la,lc)*Gf*prefc

            if(la.eq.lc)then

             do kk=1,this%Hdim

              DEner=this%Ener(kk)-this%Ener(ld)-freq1+freq2
              Gf=bose(temp,freq1)*(bose(temp,freq2)+1)*delta(type_smear,DEner,lw1)
  
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,ld))*Rabp(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,ld))*Rbam(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,ld))*Rabp(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,ld))*Rbam(kk,lb)*Gf*prefc*0.5d0

              DEner=this%Ener(kk)-this%Ener(ld)+freq1-freq2
              Gf=bose(temp,freq2)*(bose(temp,freq1)+1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,ld))*Rabm(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,ld))*Rbap(kk,lb)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,ld))*Rabm(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,ld))*Rbap(kk,lb)*Gf*prefc*0.5d0
       
              DEner=this%Ener(kk)-this%Ener(ld)-freq1-freq2
              Gf=bose(temp,freq2)*bose(temp,freq1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,ld))*Rabm(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,ld))*Rbam(kk,lb)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,ld))*Rabm(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,ld))*Rbam(kk,lb)*Gf*prefc*0.5d0

              DEner=this%Ener(kk)-this%Ener(ld)+freq1+freq2
              Gf=(bose(temp,freq2)+1)*(bose(temp,freq1)+1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,ld))*Rabp(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,ld))*Rbap(kk,lb)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,ld))*Rabp(kk,lb)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,ld))*Rbap(kk,lb)*Gf*prefc*0.5d0

             enddo

            endif

            if(lb.eq.ld)then

             do kk=1,this%Hdim

              DEner=this%Ener(kk)-this%Ener(lc)-freq1+freq2
              Gf=bose(temp,freq1)*(bose(temp,freq2)+1)*delta(type_smear,DEner,lw1)
  
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,la))*Rabp(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,la))*Rbam(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,la))*Rabp(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,la))*Rbam(kk,lc)*Gf*prefc*0.5d0

              DEner=this%Ener(kk)-this%Ener(lc)+freq1-freq2
              Gf=bose(temp,freq2)*(bose(temp,freq1)+1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,la))*Rabm(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,la))*Rbap(kk,lc)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,la))*Rabm(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,la))*Rbap(kk,lc)*Gf*prefc*0.5d0
         
              DEner=this%Ener(kk)-this%Ener(lc)-freq1-freq2
              Gf=bose(temp,freq2)*bose(temp,freq1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,la))*Rabm(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabm(kk,la))*Rbam(kk,lc)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,la))*Rabm(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbam(kk,la))*Rbam(kk,lc)*Gf*prefc*0.5d0

              DEner=this%Ener(kk)-this%Ener(lc)+freq1+freq2
              Gf=(bose(temp,freq2)+1)*(bose(temp,freq1)+1)*delta(type_smear,DEner,lw1)

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,la))*Rabp(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rabp(kk,la))*Rbap(kk,lc)*Gf*prefc*0.5d0
 
              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,la))*Rabp(kk,lc)*Gf*prefc*0.5d0

              this%R%mat(ii,jj)=this%R%mat(ii,jj)-conjg(Rbap(kk,la))*Rbap(kk,lc)*Gf*prefc*0.5d0            

             enddo

            endif

           endif

          enddo
         enddo

         deallocate(Rabp)
         deallocate(Rabm)
         deallocate(Rbap)
         deallocate(Rbam)

        return
        end subroutine make_R41


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   PROPAGATE THE DYNAMICS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine propagate(this,start_step,time,nsteps,step,dump_freq)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms
        implicit none
        class(liuville_space)         :: this
        type(dist_cmplx_vec)         :: pop,pop_new
        double precision             :: step,time
        double complex               :: val
        integer                      :: nsteps,dump_freq,start_step
        integer                      :: t1,t2,rate,i,ii,jj,l,v,indxl2g

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Propagation timestep  : ',step,' ps'
          write(*,*) '     Total propagation time: ',nsteps*step,' ps'
          flush(6)
         endif

         if(start_step.eq.0)then
          call dump_expvals(this,start_step,time)
          start_step=1
         endif

         if( allocated(this%R%mat) ) then
          call pop%set(this%Ldim,NB)
          pop%vec=(0.0d0,0.0d0)
          call pop_new%set(this%Ldim,NB)
          pop_new%vec=(0.0d0,0.0d0)
         endif

         do i=start_step,nsteps+start_step-1

          time=time+step

          do ii=1,size(this%rho%mat,1)
           do jj=1,size(this%rho%mat,2)

            l=indxl2g(ii,NB,myrow,0,nprow)
            v=indxl2g(jj,MB,mycol,0,npcol)
                    
            this%rho%mat(ii,jj)=this%rho%mat(ii,jj) &
                        *this%U%mat(1,l)*conjg(this%U%mat(1,v))
                          
           enddo
          enddo

          if( allocated(this%R%mat) ) then

           do l=1,this%Ldim
            ii=this%Lbasis(l,1)
            jj=this%Lbasis(l,2)
            call pzelget('A',' ',val,this%rho%mat,ii,jj,this%rho%desc)
            call pzelset(pop%vec,l,1,pop%desc,val)
           enddo

           call pzgemv('N',this%Ldim,this%Ldim,&
                        (1.0d0,0.0d0),this%R%mat,1,1,this%R%desc,pop%vec,1,1,pop%desc,&
                        1,(0.0d0,0d0),pop_new%vec,1,1,pop_new%desc,1)

           do l=1,this%Ldim
            ii=this%Lbasis(l,1)
            jj=this%Lbasis(l,2)
            call pzelget('A',' ',val,pop_new%vec,l,1,pop_new%desc)
            call pzelset(this%rho%mat,ii,jj,this%rho%desc,val)
           enddo   

          endif
         
          if ( mod(i,dump_freq).eq.0 ) then        
           call dump_expvals(this,i,time)
          endif

         enddo
         
         start_step=start_step+nsteps

         if( allocated(this%R%mat) ) then
          deallocate(pop%vec)
          deallocate(pop_new%vec)
         endif

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine propagate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   COMPUTE EXPECTATION VALUE OF AN HERMITIAN OPERATOR
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine make_expval(this,A,expval) 
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(liuville_space)    :: this
        type(dist_cmplx_mat)    :: A
        double precision        :: expval
        integer                 :: i,j

         expval=(0.0d0,0.0d0)

         do i=1,size(this%rho%mat,1)
          do j=1,size(this%rho%mat,2)
           expval=expval+this%rho%mat(i,j)*conjg(A%mat(i,j))
          enddo
         enddo
                      
         call mpi_allreduce(expval,expval,1,mpi_double_complex,mpi_sum,mpi_blacs_world,err)

        return
        end subroutine make_expval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!   PRINT OUT EXPECTATION VALUE OF A SET OF PRE-DEFINED HERMITIAN OPERATORS
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine dump_expvals(this,step,time)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none        
        class(liuville_space)                :: this
        integer                             :: i,k,step,unit_no
        double precision, allocatable       :: expval(:)
        double precision                    :: time
        double complex                      :: norm,val
      
         if(.not. allocated(this%QMOP))return       

         norm=(0.0d0,0.0d0)
         do k=1,this%Hdim
          val=(0.0d0,0.0d0)
          call pzelget(' ',' ',val,this%rho%mat,k,k,this%rho%desc)
          norm=norm+val
         enddo

         call mpi_allreduce(norm,norm,1,mpi_double_complex,mpi_sum,mpi_blacs_world,err)

         allocate(expval(size(this%QMOP)))
         do i=1,size(this%QMOP)
          call this%make_expval(this%QMOP(i),expval(i))
         enddo

         if(mpi_id.eq.0)then

          inquire(file='Dynamics.dat',number=unit_no)
          if(unit_no.eq.-1)  open(11,file='Dynamics.dat')

          write(11,*) step,time,expval(:),dble(norm)
          flush(11)

         endif

         deallocate(expval)

        return
        end subroutine dump_expvals

        subroutine rot_rho(this,rot)
        use mpi_utils
        use blacs_utils
        implicit none
        class(liuville_space)            :: this
        type(dist_cmplx_mat)             :: BB,rot

         call BB%set(this%Hdim,this%Hdim,NB,MB)
         BB%mat=(0.0d0,0.0d0)

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                       (1.0d0,0.0d0),rot%mat,1,1,rot%desc,this%rho%mat,&
                       1,1,this%rho%desc,(0.0d0,0.0d0),BB%mat,1,1,BB%desc)  
           
         call pzgemm('N','C',this%Hdim,this%Hdim,this%Hdim,&
                       (1.0d0,0.0d0),BB%mat,1,1,BB%desc,rot%mat,1,1,rot%desc,&
                       (0.0d0,0.0d0),this%rho%mat,1,1,this%rho%desc)

         call BB%dealloc()

        return
        end subroutine rot_rho

        end module liuville_class

