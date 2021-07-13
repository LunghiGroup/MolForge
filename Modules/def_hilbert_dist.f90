        module hilbert_dist_class
        use spins_dist_rs_class
        use spinham_class
        use lattice_class
        use sparse_class
        use blacs_utils
        use general_types_class
        use dist_class
        implicit none
                
        type, extends(spins_group)  :: spins_hilbert
         integer                            :: Hdim=0
         integer                            :: Ldim=0
         type(sub_space), allocatable       :: active_space(:)
         double precision, allocatable      :: basis(:,:)
         integer, allocatable               :: Lbasis(:,:)
         integer, allocatable               :: Hnodes(:,:)
         type(dist_cmplx_mat)               :: kbasis
         type(csr_mat_cmplx)                :: sp_kbasis
         double precision, allocatable      :: klist(:,:)
         integer, allocatable               :: kblc(:)
         type(vector_dbl), allocatable      :: Ener(:)
         type(dist1d)                       :: dos1p
         type(dist1d)                       :: dos2pm
         double precision                   :: smear=1.0d0
         type(dist_cmplx_mat)               :: rho
         type(csr_mat_cmplx)                :: sp_rho
         type(dist_cmplx_mat)               :: rho0
         type(csr_mat_cmplx)                :: sp_rho0
         type(dist_cmplx_mat), allocatable  :: H0(:)
         type(dist_cmplx_mat), allocatable  :: H(:)
         type(dist_cmplx_mat), allocatable  :: U(:)
         type(dist_cmplx_mat)               :: R
         double precision, allocatable      :: T1(:)
         complex(8), allocatable            :: Rval(:)
         type(dist_cmplx_mat)               :: T2
         type(dist_cmplx_mat), allocatable  :: Sx(:)
         type(dist_cmplx_mat), allocatable  :: Sy(:)
         type(dist_cmplx_mat), allocatable  :: Sz(:)
         type(csr_mat_cmplx), allocatable   :: sp_Sx(:)
         type(csr_mat_cmplx), allocatable   :: sp_Sy(:)
         type(csr_mat_cmplx), allocatable   :: sp_Sz(:)
         type(mat_cmplx), allocatable       :: SHrep(:,:)
         logical                            :: make_Heig=.false.
         logical                            :: make_Rmat=.false.
         logical                            :: make_R2mat=.false.
         logical                            :: make_SA=.false.
         logical                            :: make_PT2=.false.
         logical                            :: sparse=.false.
         logical                            :: printHeig=.false.
         logical                            :: printH0=.false.
         logical                            :: printRmat=.false.
         logical                            :: printRho=.false.
         integer, allocatable               :: print_si(:)
         integer                            :: s2print=0         
         contains
         procedure   ::  make_basis => make_basis_H
         procedure   ::  make_Lbasis => make_basis_L
         procedure   ::  make_kbasis
         procedure   ::  make_Hmat_nodes
         procedure   ::  make_SH_rep         
         procedure   ::  make_Hmat
         procedure   ::  make_Hmat_2
         procedure   ::  make_S => make_S_H
         procedure   ::  make_M => make_M_H
         procedure   ::  make_rho0 => make_rho0_H
         procedure   ::  make_U => make_propagator_H
         procedure   ::  make_R => make_R_H
         procedure   ::  make_R2 => make_R2_H
         procedure   ::  make_RL => make_RL_H
         procedure   ::  make_R01L => make_R01L_H
         procedure   ::  make_R02L => make_R02L_H
         procedure   ::  make_rot
         procedure   ::  get_ph_lt
         procedure   ::  rot_rho
         procedure   ::  set_dipolar
         procedure   ::  set_sph_dipolar
         procedure   ::  set_sph2_dipolar
         procedure   ::  diag_Hmat
         procedure   ::  get_Hij
         procedure   ::  get_Hij_2
         procedure   ::  get_Tij
         procedure   ::  to_eigenbasis
         procedure   ::  to_kbasis
         procedure   ::  propagate => propagate_H
         procedure   ::  dump_S => dump_S_H
         procedure   ::  dump_M => dump_M_H
         procedure   ::  dump_rho => dump_rho_H
         procedure   ::  dump_H0 => dump_H0_H
         procedure   ::  dump_Heig => dump_Heig_H
         procedure   ::  dump_kbasis => dump_kbasis_H
        end type spins_hilbert

        contains

        subroutine get_sph_dists(this,sys,phondy,max_ener)
        use mpi
        use mpi_utils
        use blacs_utils
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)         :: this
        class(brillouin)             :: phondy
        class(atoms_group)           :: sys
        double precision             :: sigma,step,max_ener
        integer                      :: ph,bn,nsteps,nloc,nstart,phx
        integer, allocatable         :: proc_grid(:)
        character(len=10)            :: file_name

         
         file_name='sph_dist'
         step=0.1d0
         nsteps=nint(3500.0/0.1)
         sigma=1.0d0

         call this%SPH%alloc_dists(sigma,step,nsteps)

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
          do bn=1,size(phondy%list(ph)%freq)

           if (ph.eq.1 .and. bn.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)

      ! Build SPH per phonon

           call this%SPH%cart2brill(this%rcell,sys,phondy,ph,bn)
           call this%SPH%get_dists(phondy%list(ph)%Freq(bn))

          enddo
          ph=ph+1
         enddo
        
         call this%SPH%merge_dists()
         if(mpi_id.eq.0) call this%SPH%dump_dists(trim(file_name))
         call this%SPH%remove_dists()

         deallocate(proc_grid)

        return
        end subroutine get_sph_dists

        subroutine get_sph2_dists(this,sys,phondy,max_ener)
        use mpi
        use mpi_utils
        use blacs_utils
        use atoms_class
        use phonons_class
        use units_parms 
        implicit none
        class(spins_hilbert)         :: this
        class(brillouin)             :: phondy
        class(atoms_group)           :: sys
        double precision             :: sigma,step,max_ener,diff
        integer                      :: ph,bn,nsteps,nloc,nstart,phx
        integer                      :: ph2,bn2,k1,k2
        integer, allocatable         :: proc_grid(:)
        character(len=10)            :: file_name

         file_name='sph2_dist'
         step=0.1d0
         nsteps=nint(3500.0/0.1)
         sigma=1.0d0

         call this%SPH2%alloc_dists(sigma,step,nsteps)

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
          do bn=1,size(phondy%list(ph)%freq)

           if (ph.eq.1 .and. bn.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)

           do ph2=1,phondy%ntot
            do bn2=1,size(phondy%list(ph2)%freq)

             k1=(ph-1)*size(phondy%list(1)%freq)+bn
             k2=(ph2-1)*size(phondy%list(1)%freq)+bn2
             if(k1.lt.k2)cycle

             if (ph2.eq.1 .and. bn2.le.3) cycle
             if (phondy%list(ph2)%freq(bn2).lt.0.0d0) cycle
             if (phondy%list(ph2)%freq(bn2).gt.max_ener) cycle

             if(.not.allocated(phondy%list(ph2)%hess) ) call phondy%list(ph2)%diagD(sys)

      ! Build SPH2 per phonon

             diff=abs(phondy%list(ph)%freq(bn)-phondy%list(ph2)%freq(bn2))*&
                     ((bose(temp(1),phondy%list(ph)%freq(bn))+1.0d0)*&
                     bose(temp(1),phondy%list(ph2)%freq(bn2))+&
                     (bose(temp(1),phondy%list(ph2)%freq(bn2))+1.0d0)*&
                     bose(temp(1),phondy%list(ph)%freq(bn)))
             call this%SPH2%cart2brill2(this%rcell,sys,phondy,ph,bn,ph2,bn2)
             call this%SPH2%get_dists(diff)

            enddo
           enddo
!           if(ph2.ne.ph .and. allocated(phondy%list(ph2)%hess)) deallocate(phondy%list(ph2)%hess)
          enddo
!          if(allocated(phondy%list(ph)%hess)) deallocate(phondy%list(ph2)%hess)
          ph=ph+1
         enddo
        
         call this%SPH2%merge_dists()
         if(mpi_id.eq.0) call this%SPH2%dump_dists(trim(file_name))
         call this%SPH2%remove_dists()

         deallocate(proc_grid)

        return
        end subroutine get_sph2_dists

        subroutine get_ph_lt(this,sys,phondy,max_ener)
        use mpi
        use mpi_utils
        use blacs_utils
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)         :: this
        class(brillouin)             :: phondy
        class(atoms_group)           :: sys
        double precision             :: val,norm,max_ener,avg_sph
        integer                      :: ph,bn

         do ph=1,phondy%ntot         
          do bn=1,size(phondy%list(ph)%freq)

      ! check spectrum overlap

           if (ph.eq.1 .and. bn.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle
           if (this%dos2pm%get_val(phondy%list(ph)%freq(bn)).lt.1.0d-6) cycle

!           call phondy%calc_linewidth_sp(sys,ph,bn)
!           write(*,*) ph,bn,phondy%list(ph)%width(bn,1)

          enddo
         enddo

        return
        end subroutine get_ph_lt

        subroutine make_R02L_H(this,sys,phondy,step_min,mult_fact,max_ener,R0,euler)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_diag_asimm_cmplx
        use scalapack_inv
        use lapack_diag_asimm
        use lapack_inverse
        use units_parms 
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)         :: this
        class(brillouin)             :: phondy
        class(atoms_group)           :: sys
        type(dist_cmplx_mat)         :: AA,BB,CC
        type(dist_cmplx_mat)         :: R0,R0inv
        double precision             :: Gf,DEner,step_min,coeff(3),DEner2
        double precision             :: val,norm,max_ener,avg_sph,val2
        double precision             :: euler(3)
        complex(8)                   :: valc,nodiag_sum,diag_sum,kcons,valc2
        complex(8), allocatable      :: Vmat(:,:)
        integer                      :: ph,bn,ii,jj,size_block,i1,ii_1,i,jj_1
        integer                      :: l,l2,la,lb,l2a,l2b,spin_id,spin_id2
        integer                      :: k,ia,ib,ka,kb,ic,kc,id,kd,vv
        integer                      :: nphonons,v,nze,bn2,ph2
        integer                      :: t1,t2,rate,indxl2g,px,py,pz
        integer                      :: mult_fact,s,t,j
        integer                      :: size_block_1,size_block_2,i2,l1
        integer                      :: nloc,nstart,phx
        integer, allocatable         :: proc_grid(:)
        logical                      :: check_SA

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Redfield matrix: 1st-order PT + 2nd-order of coupling strength'
          flush(6)
         endif

         call R0%set(this%Ldim,this%Ldim,NB,MB)
         R0%mat=(0.0d0,0.0d0)

         call AA%set(this%Hdim,this%Hdim,NB,MB)          

         allocate(Vmat(this%Hdim,this%Hdim))
         Vmat=(0.0d0,0.0d0)

         nphonons=0
         avg_sph=0.0d0
        
!         call get_sph2_dists(this,sys,phondy,max_ener)
!         stop
         ! make lists of phonons to be included and distribute them on
         ! processes        

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph)%width))then
            allocate(phondy%list(ph)%width(size(phondy%list(ph)%freq),1))
            if(phondy%effective_lt)then
              phondy%list(ph)%width=temp(1)*0.5d0
!             do i=1,size(phondy%list(ph)%freq)
!              phondy%list(ph)%width(i,1)=smear*&
!                (1+exp(phondy%list(ph)%freq(i)/(kboltz*temp(1)*2.0d0)))/&
!                (exp(phondy%list(ph)%freq(i)/(kboltz*temp(1)*2.0d0))-1)
!             enddo
            else
             phondy%list(ph)%width=smear
            endif
           endif
          endif

         do ph2=1,phondy%ntot
         
      ! Check for Ph lifetime

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph2)%width))then
            allocate(phondy%list(ph2)%width(size(phondy%list(ph2)%freq),1))
            phondy%list(ph2)%width=smear
           endif
          endif

          do bn=1,size(phondy%list(ph)%freq)
          do bn2=1,size(phondy%list(ph2)%freq)

           if( ((ph2-1)*size(phondy%list(ph2)%freq)+bn2) .lt. &
               ((ph-1)*size(phondy%list(ph)%freq)+bn) ) cycle

      ! check spectrum overlap

           if (ph.eq.1 .and. bn.le.3) cycle
           if (ph2.eq.1 .and. bn2.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph2)%freq(bn2).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle
           if (phondy%list(ph2)%freq(bn2).gt.max_ener) cycle
           DEner=abs(phondy%list(ph)%freq(bn)-phondy%list(ph2)%freq(bn2))
           DEner2=abs(phondy%list(ph)%freq(bn)+phondy%list(ph2)%freq(bn2))
           if (this%dos2pm%get_val(abs(DEner)).lt.1.0d-6 .and. &
               this%dos2pm%get_val(abs(DEner2)).lt.1.0d-6 ) cycle

           nphonons=nphonons+1

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)
           if(.not.allocated(phondy%list(ph2)%hess) ) call phondy%list(ph2)%diagD(sys)

      ! Build SPH per phonon

           call this%SPH2%cart2brill2(this%rcell,sys,phondy,ph,bn,ph2,bn2)

           do spin_id=1,this%nspins_pr
           do spin_id2=spin_id,this%nspins

      ! Make Vij per spin    

            AA%mat=(0.0d0,0.0d0)

            call make_SH_rep(this,this%SPH2,spin_id,spin_id2)

            do l=1,size(this%Hnodes,1)

             ii=this%Hnodes(l,1)
             jj=this%Hnodes(l,2)
            
             AA%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

            enddo ! l

      ! Rotate Vij

            if(this%ntot.gt.1)then
             call this%to_kbasis(AA)
            endif 

            if(this%make_Heig)then
             call this%to_eigenbasis(AA)
            endif

      ! Make Rij

            Vmat=(0.0d0,0.0d0)
            do l2=1,this%Hdim
             do l=1,this%Hdim
              call pzelget('A',' ',valc,AA%mat,l2,l,AA%desc)
              Vmat(l2,l)=valc
             enddo
            enddo                   

            do ii=1,size(R0%mat,1)
             do jj=1,size(R0%mat,2)

              l=indxl2g(ii,NB,myrow,0,nprow)
              l2=indxl2g(jj,MB,mycol,0,npcol)

              la=this%Lbasis(l,1)
              lb=this%Lbasis(l,2)

              l2a=this%Lbasis(l2,1)
              l2b=this%Lbasis(l2,2)

              do ii_1=1,this%ntot
               if(la.le.this%kblc(ii_1+1))then
                ka=ii_1
                ia=la-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(lb.le.this%kblc(ii_1+1))then
                kb=ii_1
                ib=lb-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(l2a.le.this%kblc(ii_1+1))then
                kc=ii_1
                ic=l2a-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(l2b.le.this%kblc(ii_1+1))then
                kd=ii_1
                id=l2b-this%kblc(ii_1)
                exit
               endif
              enddo

              if(this%ntot.gt.1)then
               stop
              else
               kcons=1.0d0
              endif

              Gf=0.0d0

              DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)  &
                         -phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
              Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*&
                 bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                 delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)  &
                         +phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                     bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)  &
                         -phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
              Gf=Gf+bose(temp(1),phondy%list(ph)%Freq(bn))*&
                    (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)  &
                         +phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                    (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*conjg(Vmat(lb,l2b))*Gf*kcons
 
              Gf=0.0d0
 
              DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)  &
                         -phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
              Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*&
                 bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                 delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)  &
                         +phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                     bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)  &
                         -phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
              Gf=Gf+bose(temp(1),phondy%list(ph)%Freq(bn))*&
                    (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)  &
                         +phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                    (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))
 
              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*conjg(Vmat(lb,l2b))*Gf*kcons

              if(l2b.eq.lb)then
!!!!!
              if(this%ntot.gt.1)then
               stop
              else
               kcons=1.0d0
              endif

              vv=0
              do ii_1=1,this%ntot
               do jj_1=1,(this%kblc(ii_1+1)-this%kblc(ii_1))
                vv=vv+1 

                Gf=0.0d0

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)  &
                           -phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
                Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*&
                   bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                   delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)  &
                           +phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
                Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                       bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                       delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)  &
                           -phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
                Gf=Gf+bose(temp(1),phondy%list(ph)%Freq(bn))*&
                      (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                      delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)  &
                           +phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
                Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                      (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                       delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(la,vv)*conjg(Vmat(l2a,vv))*Gf*kcons

               enddo
              enddo
             endif

             if(l2a.eq.la)then
!!!!!                     
              if(this%ntot.gt.1)then
               stop
              else
               kcons=1.0d0
              endif

              vv=0
              do ii_1=1,this%ntot
               do jj_1=1,(this%kblc(ii_1+1)-this%kblc(ii_1))
                vv=vv+1 

                Gf=0.0d0

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)  &
                           -phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
                Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*&
                   bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                   delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)  &
                           +phondy%list(ph)%Freq(bn)-phondy%list(ph2)%Freq(bn2)
                Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                       bose(temp(1),phondy%list(ph2)%Freq(bn2))*&
                       delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)  &
                           -phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
                Gf=Gf+bose(temp(1),phondy%list(ph)%Freq(bn))*&
                      (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                      delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)  &
                           +phondy%list(ph)%Freq(bn)+phondy%list(ph2)%Freq(bn2)
                Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*&
                      (bose(temp(1),phondy%list(ph2)%Freq(bn2))+1.0d0)*&
                       delta(type_smear,DEner,phondy%list(ph)%width(bn,1)+phondy%list(ph2)%width(bn2,1))

                R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(l2b,vv)*conjg(Vmat(lb,vv))*Gf*kcons

               enddo
              enddo
             endif
             
            enddo  !! jj
           enddo !! ii

          enddo ! spin_id2
          enddo ! spin_id

          enddo ! bn2
          enddo ! bn

          if(ph2.ne.ph)then
           if(allocated(phondy%list(ph2)%hess)) deallocate(phondy%list(ph2)%hess)
           if(allocated(phondy%list(ph2)%width)) deallocate(phondy%list(ph2)%width)
          endif

         enddo ! ph2

          if(allocated(phondy%list(ph)%hess)) deallocate(phondy%list(ph)%hess)
          if(allocated(phondy%list(ph)%width)) deallocate(phondy%list(ph)%width)

          ph=ph+1
         enddo ! phx

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(R0%mat(ii,jj),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           R0%mat(ii,jj)=valc
          enddo
         enddo
         
         call mpi_allreduce(nphonons,nze,1,&
              mpi_integer,mpi_sum,mpi_phonons_world,err)
         call mpi_allreduce(avg_sph,val,1,&
              mpi_double_precision,mpi_sum,mpi_phonons_world,err)

         nphonons=nze
         avg_sph=val
         
         if(mpi_id.eq.0) write(*,*) '     Total number of phonons included: ',nphonons

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           la=this%Lbasis(l,1)
           lb=this%Lbasis(l,2)

           l2a=this%Lbasis(l2,1)
           l2b=this%Lbasis(l2,2)

           do ii_1=1,this%ntot
            if(la.le.this%kblc(ii_1+1))then
             ka=ii_1
             ia=la-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(lb.le.this%kblc(ii_1+1))then
             kb=ii_1
             ib=lb-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(l2a.le.this%kblc(ii_1+1))then
             kc=ii_1
             ic=l2a-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(l2b.le.this%kblc(ii_1+1))then
             kd=ii_1
             id=l2b-this%kblc(ii_1)
             exit
            endif
           enddo

!!!!!!!!!!!!!!!!!!!!!!!

!           check_SA=.true.
!           if(la.eq.l2a .and. lb.eq.l2b) check_SA=.false.
!           if(la.eq.lb .and. l2a.eq.l2b) check_SA=.false.
!           if(check_SA) R0%mat(ii,jj)=(0.0d0,0.0d0)           

!!!!!!!!!!!!!!!!!!!!!!!

           DEner=this%Ener(ka)%v(ia)-this%Ener(kc)%v(ic)+  &
                 this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)

           if( abs(DEner).lt.1.0e-6 )then              

            R0%mat(ii,jj)=pi*pi*R0%mat(ii,jj)*step_min*mult_fact/hplank/2.0d0

           else

            R0%mat(ii,jj)=R0%mat(ii,jj)*hplank*2.0d0*acos(-1.0d0)*cmplx(0.0d0,-1.0d0,8)*  &
                          (exp(2.0d0*acos(-1.0d0)*cmplx(0.0d0,1.0d0,8)*DEner*step_min*mult_fact/hplank)-1) &
                          /DEner

            R0%mat(ii,jj)=pi*pi*R0%mat(ii,jj)/hplank/2.0d0

           endif
        
          enddo
         enddo

         do l2=1,this%Ldim
          do l=l2+1,this%Ldim
           call pzelget('A',' ',valc,R0%mat,l2,l,R0%desc)
           call pzelget('A',' ',valc2,R0%mat,l,l2,R0%desc)
           call pzelset(R0%mat,l2,l,R0%desc,conjg(valc2))
           call pzelset(R0%mat,l,l2,R0%desc,conjg(valc))
          enddo
         enddo                   

         call AA%dealloc()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_R02L_H

        subroutine make_R01L_H(this,sys,phondy,step_min,mult_fact,max_ener,R0,euler)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_diag_asimm_cmplx
        use scalapack_inv
        use lapack_diag_asimm
        use lapack_inverse
        use units_parms 
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)         :: this
        class(brillouin)             :: phondy
        class(atoms_group)           :: sys
        type(dist_cmplx_mat)         :: AA,BB,CC
        type(dist_cmplx_mat)         :: R0,R0inv
        double precision             :: Gf,DEner,step_min,coeff(3)
        double precision             :: val,norm,max_ener,avg_sph,val2,alpha,beta,gamma
        complex(8)                   :: valc,nodiag_sum,diag_sum,kcons,valc2,kcons2
        complex(8), allocatable      :: Vmat(:,:)
        double precision             :: euler(3)
        integer                      :: ph,bn,ii,jj,size_block,i1,ii_1,i,jj_1
        integer                      :: l,l2,la,lb,l2a,l2b,spin_id,spin_id2
        integer                      :: k,ia,ib,ka,kb,ic,kc,id,kd,vv,term
        integer                      :: nphonons,v,nze
        integer                      :: t1,t2,rate,indxl2g,px,py,pz,phx,nloc,nstart
        integer                      :: mult_fact,s,t,j
        integer                      :: size_block_1,size_block_2,i2,l1
        integer, allocatable         :: proc_grid(:)

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Redfield matrix: 1st-order PT + 1st-order of coupling strength'
          flush(6)
         endif

         call AA%set(this%Hdim,this%Hdim,NB,MB)          
         allocate(Vmat(this%Hdim,this%Hdim))
         Vmat=(0.0d0,0.0d0)

         nphonons=0
         avg_sph=0.0d0

         call get_sph_dists(this,sys,phondy,max_ener)

         ! make lists of phonons to be included and distribute them on
         ! processes        

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
         
      ! Check for Ph lifetime

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph)%width))then
            allocate(phondy%list(ph)%width(size(phondy%list(ph)%freq),1))            
            phondy%list(ph)%width=smear
            if(phondy%effective_lt)then
             do bn=1,size(phondy%list(ph)%freq)
              phondy%list(ph)%width(bn,1)=phondy%list(ph)%width(bn,1)+&
                     phondy%list(ph)%freq(bn)*&
                     exp(phondy%list(ph)%freq(bn)/kboltz/temp(1)/2)/&
                    (exp(phondy%list(ph)%freq(bn)/kboltz/temp(1))-1)
              if(isnan(phondy%list(ph)%width(bn,1))) phondy%list(ph)%width(bn,1)=0.0d0
             enddo
            endif
           endif
          endif

          do bn=1,size(phondy%list(ph)%freq)


      ! check spectrum overlap

           if (ph.eq.1 .and. bn.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle
!           if (this%dos2pm%get_val(phondy%list(ph)%freq(bn)).lt.1.0d-6) cycle

           nphonons=nphonons+1

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)

      ! Build SPH per phonon

           do spin_id=1,this%nspins_pr
           do spin_id2=spin_id,this%nspins
        
           call this%SPH%cart2brill(this%rcell,sys,phondy,ph,bn)

         ! DyCp NEV_FIC_RIJK_mllb
!           gamma=0.33007747868274917        
!           beta=-1.5453443590433216        
!           alpha=-1.8207758543913428

         ! Dyacac inverse euler of gmolcas
!            gamma=-2.2139203072213745
!            beta=-1.3102504949441061
!            alpha=2.0059623380955780
           
            if (any( abs(euler).gt.1.0d-4 )) call this%SPH%rot(euler)

      ! Make Vij per spin    

            AA%mat=(0.0d0,0.0d0)

            call make_SH_rep(this,this%SPH,spin_id,spin_id2)

            do l=1,size(this%Hnodes,1)

             ii=this%Hnodes(l,1)
             jj=this%Hnodes(l,2)
            
             AA%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

            enddo ! l

      ! Rotate Vij

            if(this%ntot.gt.1)then
             call this%to_kbasis(AA)
            endif 
 
            if(this%make_Heig)then
             call this%to_eigenbasis(AA)
            endif

      ! Make Rij

            Vmat=(0.0d0,0.0d0)
            do l2=1,this%Hdim
             do l=1,this%Hdim
              call pzelget('A',' ',valc,AA%mat,l2,l,AA%desc)
              Vmat(l2,l)=valc
             enddo
            enddo                   

            do ii=1,size(R0%mat,1)
             do jj=1,size(R0%mat,2)

              l=indxl2g(ii,NB,myrow,0,nprow)
              l2=indxl2g(jj,MB,mycol,0,npcol)

              la=this%Lbasis(l,1)
              lb=this%Lbasis(l,2)

              l2a=this%Lbasis(l2,1)
              l2b=this%Lbasis(l2,2)

              do ii_1=1,this%ntot
               if(la.le.this%kblc(ii_1+1))then
                ka=ii_1
                ia=la-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(lb.le.this%kblc(ii_1+1))then
                kb=ii_1
                ib=lb-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(l2a.le.this%kblc(ii_1+1))then
                kc=ii_1
                ic=l2a-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(l2b.le.this%kblc(ii_1+1))then
                kd=ii_1
                id=l2b-this%kblc(ii_1)
                exit
               endif
              enddo

 
              if(this%ntot.gt.1)then

               coeff(1)=phondy%list(ph)%k(1)-this%klist(ka,1)+this%klist(kc,1)   
               coeff(2)=phondy%list(ph)%k(2)-this%klist(ka,2)+this%klist(kc,2)   
               coeff(3)=phondy%list(ph)%k(3)-this%klist(ka,3)+this%klist(kc,3)   
               coeff(1)=coeff(1)-phondy%list(ph)%k(1)-this%klist(kd,1)+this%klist(kb,1)
               coeff(2)=coeff(2)-phondy%list(ph)%k(2)-this%klist(kd,2)+this%klist(kb,2)
               coeff(3)=coeff(3)-phondy%list(ph)%k(3)-this%klist(kd,3)+this%klist(kb,3)
               kcons=(0.0d0,0.0d0)
               do v=1,this%ntot
                kcons=kcons+exp(cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*&
                               (coeff(1)*this%rcell(v,1)+&
                                coeff(2)*this%rcell(v,2)+&
                                coeff(3)*this%rcell(v,3)))
               enddo

              else

               kcons=(1.0d0,0.0d0)

              endif

              Gf=0.0d0
              if(this%Ener(kd)%v(id).gt.this%Ener(kb)%v(ib))then
               DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)-phondy%list(ph)%Freq(bn)
               Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              else
               DEner=this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)+phondy%list(ph)%Freq(bn)
               Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              endif
!              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*Vmat(l2b,lb)*Gf*kcons*conjg(kcons)
              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*conjg(Vmat(lb,l2b))*Gf*kcons
 
              Gf=0.0d0
              if(this%Ener(kc)%v(ic).gt.this%Ener(ka)%v(ia))then
               DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)-phondy%list(ph)%Freq(bn)
               Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              else
               DEner=this%Ener(kc)%v(ic)-this%Ener(ka)%v(ia)+phondy%list(ph)%Freq(bn)
               Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              endif
!              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*Vmat(l2b,lb)*Gf*kcons*conjg(kcons)
              R0%mat(ii,jj)=R0%mat(ii,jj)+Vmat(la,l2a)*conjg(Vmat(lb,l2b))*Gf*kcons

              if(l2b.eq.lb)then
               vv=0
               do ii_1=1,this%ntot
                if(this%ntot.gt.1)then

                 coeff(1)=phondy%list(ph)%k(1)+this%klist(ka,1)-this%klist(ii_1,1)   
                 coeff(2)=phondy%list(ph)%k(2)+this%klist(ka,2)-this%klist(ii_1,2)   
                 coeff(3)=phondy%list(ph)%k(3)+this%klist(ka,3)-this%klist(ii_1,3)   
                 coeff(1)=coeff(1)-phondy%list(ph)%k(1)+this%klist(ii_1,1)-this%klist(kc,1)
                 coeff(2)=coeff(2)-phondy%list(ph)%k(2)+this%klist(ii_1,2)-this%klist(kc,2)
                 coeff(3)=coeff(3)-phondy%list(ph)%k(3)+this%klist(ii_1,3)-this%klist(kc,3)
                 kcons=(0.0d0,0.0d0)
                 do v=1,this%ntot
                  kcons=kcons+exp(cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*&
                              (coeff(1)*this%rcell(v,1)+&
                               coeff(2)*this%rcell(v,2)+&
                               coeff(3)*this%rcell(v,3)))
                 enddo
 
                else

                 kcons=(1.0d0,0.0d0)
 
                endif

                do jj_1=1,(this%kblc(ii_1+1)-this%kblc(ii_1))
                 vv=vv+1 
                 Gf=0.0d0
                 if(this%Ener(ii_1)%v(jj_1).gt.this%Ener(kc)%v(ic))then
                  DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)-phondy%list(ph)%Freq(bn)
                  Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                 else
                  DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kc)%v(ic)+phondy%list(ph)%Freq(bn)
                  Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                 endif
!                 R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(la,vv)*Vmat(vv,l2a)*Gf*kcons*conjg(kcons)
                 R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(la,vv)*conjg(Vmat(l2a,vv))*Gf*kcons
                enddo
               enddo

              endif

              if(l2a.eq.la)then

               vv=0
               do ii_1=1,this%ntot
                if(this%ntot.gt.1)then
                 coeff(1)=phondy%list(ph)%k(1)+this%klist(kd,1)-this%klist(ii_1,1)   
                 coeff(2)=phondy%list(ph)%k(2)+this%klist(kd,2)-this%klist(ii_1,2)   
                 coeff(3)=phondy%list(ph)%k(3)+this%klist(kd,3)-this%klist(ii_1,3)   
                 coeff(1)=coeff(1)-phondy%list(ph)%k(1)+this%klist(ii_1,1)-this%klist(kb,1)
                 coeff(2)=coeff(2)-phondy%list(ph)%k(2)+this%klist(ii_1,2)-this%klist(kb,2)
                 coeff(3)=coeff(3)-phondy%list(ph)%k(3)+this%klist(ii_1,3)-this%klist(kb,3)
                 kcons=(0.0d0,0.0d0)
                 do v=1,this%ntot
                  kcons=kcons+exp(cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*&
                               (coeff(1)*this%rcell(v,1)+&
                                coeff(2)*this%rcell(v,2)+&
                                coeff(3)*this%rcell(v,3)))
                 enddo
 
                else

                 kcons=(1.0d0,0.0d0)
 
                endif
                do jj_1=1,(this%kblc(ii_1+1)-this%kblc(ii_1))
                 vv=vv+1 
                 Gf=0.0d0
                 if(this%Ener(ii_1)%v(jj_1).gt.this%Ener(kd)%v(id))then
                  DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)-phondy%list(ph)%Freq(bn)
                  Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                 else
                  DEner=this%Ener(ii_1)%v(jj_1)-this%Ener(kd)%v(id)+phondy%list(ph)%Freq(bn)
                  Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                 endif
!                 R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(l2b,vv)*Vmat(vv,lb)*Gf*kcons*conjg(kcons)
                 R0%mat(ii,jj)=R0%mat(ii,jj)-Vmat(l2b,vv)*conjg(Vmat(lb,vv))*Gf*kcons
                enddo
               enddo
              endif

             enddo  !! jj
            enddo !! ii

           enddo ! spin_id2
           enddo ! spin_id

          enddo ! bn

          if(allocated(phondy%list(ph)%hess)) deallocate(phondy%list(ph)%hess)
          if(allocated(phondy%list(ph)%width)) deallocate(phondy%list(ph)%width)

          ph=ph+1
         enddo ! phx

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(R0%mat(ii,jj),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           R0%mat(ii,jj)=valc
          enddo
         enddo
         
         call mpi_allreduce(nphonons,nze,1,&
              mpi_integer,mpi_sum,mpi_phonons_world,err)

         nphonons=nze

         if(mpi_id.eq.0) write(*,*) '     Number of phonons included: ',nphonons

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           la=this%Lbasis(l,1)
           lb=this%Lbasis(l,2)

           l2a=this%Lbasis(l2,1)
           l2b=this%Lbasis(l2,2)

           do ii_1=1,this%ntot
            if(la.le.this%kblc(ii_1+1))then
             ka=ii_1
             ia=la-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(lb.le.this%kblc(ii_1+1))then
             kb=ii_1
             ib=lb-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(l2a.le.this%kblc(ii_1+1))then
             kc=ii_1
             ic=l2a-this%kblc(ii_1)
             exit
            endif
           enddo

           do ii_1=1,this%ntot
            if(l2b.le.this%kblc(ii_1+1))then
             kd=ii_1
             id=l2b-this%kblc(ii_1)
             exit
            endif
           enddo

           DEner=this%Ener(ka)%v(ia)-this%Ener(kc)%v(ic)+  &
                 this%Ener(kd)%v(id)-this%Ener(kb)%v(ib)

!          Enforce diagonal secular approximation

!           if ( ((ka.eq.kc .and. ia.eq.ic) .and. (kd.eq.kb .and. id.eq.ib)) &
!                .or. ((ka.eq.kb .and. ia.eq.ib) .and. (kd.eq.kc .and. id.eq.ic)) ) then
        
           if( abs(DEner).lt.1.0e-6 )then              

            R0%mat(ii,jj)=pi*pi*R0%mat(ii,jj)*step_min*mult_fact/hplank

!           if ( .not. ( ((ka.eq.kc .and. ia.eq.ic) .and. (kd.eq.kb .and. id.eq.ib)) &
!                .or. ((ka.eq.kb .and. ia.eq.ib) .and. (kd.eq.kc .and. id.eq.ic)) ) ) then
!             if( dble(R0%mat(ii,jj)*conjg(R0%mat(ii,jj))) .gt. 1.0e-12 ) then
!              write(*,*) ia,ib,ic,id,R0%mat(ii,jj)
!             endif
!            endif

           else

            R0%mat(ii,jj)=(0.0d0,0.0d0)

           endif
        
          enddo
         enddo

         ! Taking the Transpose conjugate of R0

         do l2=1,this%Ldim
          do l=l2+1,this%Ldim
           call pzelget('A',' ',valc,R0%mat,l2,l,R0%desc)
           call pzelget('A',' ',valc2,R0%mat,l,l2,R0%desc)
           call pzelset(R0%mat,l2,l,R0%desc,conjg(valc2))
           call pzelset(R0%mat,l,l2,R0%desc,conjg(valc))
          enddo
         enddo                   

         call AA%dealloc()

         if(this%printRmat)then
          if(mpi_id.eq.0) open(15,file='R.dat')
          do l=1,this%Ldim
           do l2=1,this%Ldim         
            call pzelget('A',' ',valc,R0%mat,l2,l,R0%desc)
            if(mpi_id.eq.0) write(15,*) l2,l,dble(valc),aimag(valc)
           enddo
           if(mpi_id.eq.0) write(15,*)
          enddo
          if(mpi_id.eq.0) close(15)
         endif

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_R01L_H

        subroutine make_RL_H(this,sys,phondy,step_min,mult_fact,max_ener,euler)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_diag_asimm_cmplx
        use scalapack_inv
        use lapack_diag_asimm
        use lapack_inverse
        use units_parms 
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)          :: this
        class(brillouin)              :: phondy
        class(atoms_group)            :: sys
        type(dist_cmplx_mat)          :: AA,BB,CC
        type(dist_cmplx_mat)          :: R0,R0inv,R02
        double precision              :: Gf,DEner,step_min,coeff(3)
        double precision              :: val,norm,max_ener,avg_sph,val2
        double precision, allocatable :: rates(:)
        double precision              :: euler(3)
        complex(8)                    :: valc,nodiag_sum,diag_sum,kcons,valc2
        complex(8), allocatable       :: Vmat(:,:)
        integer                       :: ph,bn,ii,jj,size_block,i1,ii_1,i,jj_1
        integer                       :: l,l2,la,lb,l2a,l2b,spin_id,spin_id2
        integer                       :: k,ia,ib,ka,kb,ic,kc,id,kd,vv
        integer                       :: nphonons,v,nze
        integer                       :: t1,t2,rate,indxl2g,px,py,pz
        integer                       :: mult_fact,s,t,j
        integer                       :: size_block_1,size_block_2,i2,l1

         this%make_Heig=.true.
         if(.not.this%make_Heig)then          
          write(*,*) 'Warning, open system symulation in the Sz basis', &
                     'is not implemented yet'
          return
         endif

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Initializing the calculation of spin-phonon relaxation with the Redfield theory'
          flush(6)
         endif

         call R0%set(this%Ldim,this%Ldim,NB,MB)
         R0%mat=(0.0d0,0.0d0)

         if(this%make_Rmat) call this%make_R01L(sys,phondy,step_min,mult_fact,max_ener,R0,euler)
         if(this%make_R2mat)then
          call this%make_R02L(sys,phondy,step_min,mult_fact,max_ener,R02,euler)
          R0%mat=R0%mat+R02%mat
          call R02%dealloc()
         endif       

         call AA%set(this%Ldim,this%Ldim,NB,MB)          
         AA%mat=(0.0d0,0.0d0)
         allocate(this%Rval(this%Ldim))

         call BB%set(this%Ldim,this%Ldim,NB,MB)
         BB%mat=(0.0d0,0.0d0)

         !!!!! Scalapack version

         call pzdiag2(this%Ldim,R0,this%Rval,AA,BB)
         BB%mat=AA%mat
         call pzgeinv(this%Ldim,BB)

         !!!!! Lapack version

!         AA%mat=R0%mat
!         call new_diag2(this%Ldim,AA%mat,this%Rval)
!         BB%mat=AA%mat
!         call mat_inv(BB%mat,this%Ldim) 

         if(mpi_id.eq.0) then
          write(*,*) '     Redfield Matrix Eigenvalues:'
          allocate(rates(this%Ldim))
          rates=abs(dble(this%Rval))
          call  order_array(rates)
          do l=1,this%Ldim
           write(*,*) '          ',l,rates(l)
          enddo
          deallocate(rates)
         endif

      ! Build Pop Propagator R= R exp(Rval) R^{\cross}

         this%Rval=exp(this%Rval)

         call this%R%set(this%Ldim,this%Ldim,NB,MB)
         this%R%mat=(0.0d0,0.0d0)

         do ii=1,size(AA%mat,1)
          do jj=1,size(AA%mat,2)
           l=indxl2g(jj,MB,mycol,0,npcol)
           AA%mat(ii,jj)=AA%mat(ii,jj)*this%Rval(l)
          enddo
         enddo

         call pzgemm('N','N',this%Ldim,this%Ldim,this%Ldim,&
                     (1.0d0,0.0d0),AA%mat,1,1,AA%desc,BB%mat,&
                     1,1,BB%desc,(0.0d0,0.0d0),this%R%mat,1,1,this%R%desc)


         call AA%dealloc()
         call BB%dealloc()
         
         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_RL_H

        subroutine make_R_H(this,sys,phondy,step_min,mult_fact,max_ener,euler)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_inv
        use lapack_diag_asimm
        use lapack_inverse
        use units_parms 
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)          :: this
        class(brillouin)              :: phondy
        class(atoms_group)            :: sys
        type(dist_cmplx_mat)          :: AA,BB,CC
        type(dist_dbl_mat)            :: R0,R0inv
        double precision              :: Gf,DEner,step_min,coeff(3)
        double precision              :: alpha,beta,gamma
        double precision              :: val,norm,max_ener,avg_sph,Gp,Gm
        complex(8),allocatable        :: Vii(:)
        double precision, allocatable :: rates(:)
        double precision              :: euler(3)
        complex(8)                    :: valc,nodiag_sum,diag_sum,kcons
        integer                       :: ph,bn,ii,jj,l,l2,size_block,i1,ii_1,i,vv
        integer                       :: k,ia,ib,ka,kb,nphonons,v,spin_id,spin_id2
        integer                       :: t1,t2,rate,indxl2g,px,py,pz,nloc,nstart
        integer                       :: mult_fact,s,t,j,term
        integer                       :: size_block_1,size_block_2,i2,l1,nze,phx
        integer, allocatable          :: proc_grid(:)

         this%make_Heig=.true.
         if(.not.this%make_Heig)then          
          write(*,*) 'Warning, open system symulation in the Sz basis', &
                     'is not implemented yet'
          return
         endif

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Redfield matrix: 1st-order PT + 1st-order of coupling strength'
          write(*,*) '     The diagonal approximation to the secular Redfield equations will be used'
          flush(6)
         endif

         call R0%set(this%Hdim,this%Hdim,NB,MB)
         R0%mat=0.0d0

         call this%T2%set(this%Hdim,this%Hdim,NB,MB)
         this%T2%mat=(0.0d0,0.0d0)

         allocate(this%T1(this%Hdim))
         allocate(Vii(this%Hdim))

         call AA%set(this%Hdim,this%Hdim,NB,MB)          

         nphonons=0
         avg_sph=0.0d0

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc
         
      ! Check for Ph lifetime

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph)%width))then
            allocate(phondy%list(ph)%width(size(phondy%list(ph)%freq),1))
            phondy%list(ph)%width=smear
!            read(12,*) val,(phondy%list(ph)%width(i,1),i=1,size(phondy%list(ph)%freq))
           endif
          endif

          do bn=1,size(phondy%list(ph)%freq)

      ! check spectrum overlap

           if (ph.eq.1 .and. bn.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle
!           if (this%dos2pm%get_val(phondy%list(ph)%freq(bn)).lt.1.0d-6) cycle

           nphonons=nphonons+1

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)

      ! Build SPH per phonon

           call this%SPH%cart2brill(this%rcell,sys,phondy,ph,bn)

         ! Eulers from gg of Dy_NEV_RIJX_X4_mlb.out ! This work 
!         alpha=-0.35182069521786968
!         beta=1.5438722669095704
!         gamma=1.3851648663175002

         ! Eulers from gg of Dy_NEV_RIJK_mlb.out ! This work 
!         alpha=-0.32993455234750912        
!         beta=1.5452772190818402       
!         gamma=0.88437557128663258

         ! Dyacac inverse euler of gmolcas
!           gamma=-2.2139203072213745
!           beta=-1.3102504949441061
!           alpha=2.0059623380955780

!           do v=1,size(this%SPH%O)
!            call this%SPH%O(v)%rot(alpha,beta,gamma)
!           enddo

          if (any( abs(euler).gt.1.0d-4 )) call this%SPH%rot(euler)

          do spin_id=1,this%nspins_pr
          do spin_id2=spin_id,this%nspins

      ! Make Vij

           AA%mat=(0.0d0,0.0d0)

           call make_SH_rep(this,this%SPH,spin_id,spin_id2)

           do l=1,size(this%Hnodes,1)

            ii=this%Hnodes(l,1)
            jj=this%Hnodes(l,2)
            
            AA%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

           enddo ! l

      ! Rotate Vij

           if(this%ntot.gt.1)then

           call this%to_kbasis(AA)
           
            do i=1,this%ntot
             do j=1,this%ntot
              coeff(1)=phondy%list(ph)%k(1)-this%klist(i,1)+this%klist(j,1)
              coeff(2)=phondy%list(ph)%k(2)-this%klist(i,2)+this%klist(j,2)
              coeff(3)=phondy%list(ph)%k(3)-this%klist(i,3)+this%klist(j,3)
              kcons=(0.0d0,0.0d0)
              do v=1,this%ntot
               kcons=kcons+exp(cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*&
                              (coeff(1)*this%rcell(v,1)+&
                               coeff(2)*this%rcell(v,2)+&
                               coeff(3)*this%rcell(v,3)))
              enddo
              coeff(1)=-phondy%list(ph)%k(1)-this%klist(i,1)+this%klist(j,1)
              coeff(2)=-phondy%list(ph)%k(2)-this%klist(i,2)+this%klist(j,2)
              coeff(3)=-phondy%list(ph)%k(3)-this%klist(i,3)+this%klist(j,3)
              do v=1,this%ntot
               kcons=kcons+exp(cmplx(0.0d0,1.0d0,8)*2.0d0*acos(-1.0d0)*&
                              (coeff(1)*this%rcell(v,1)+&
                               coeff(2)*this%rcell(v,2)+&
                               coeff(3)*this%rcell(v,3)))
              enddo
              kcons=kcons*0.5d0
              do ii=1,size(AA%mat,1)
               do jj=1,size(AA%mat,2)
                l=indxl2g(ii,NB,myrow,0,nprow)
                l2=indxl2g(jj,MB,mycol,0,npcol)
                if(l.le.this%kblc(i+1).and.l.gt.this%kblc(i) .and. &
                   l2.le.this%kblc(j+1).and.l2.gt.this%kblc(j)) then
                 AA%mat(ii,jj)=AA%mat(ii,jj)*kcons
                endif
               enddo
              enddo
             enddo ! j
            enddo ! i

           endif           

           if(this%make_Heig)then
            call this%to_eigenbasis(AA)
           endif
                     
           do l=1,this%Hdim
            call pzelget('A',' ',Vii(l),AA%mat,l,l,AA%desc)
           enddo

      ! Make Rij

           do ii=1,size(R0%mat,1)
            do jj=1,size(R0%mat,2)

             l=indxl2g(ii,NB,myrow,0,nprow)
             l2=indxl2g(jj,MB,mycol,0,npcol)

             do ii_1=1,this%ntot
              if(l.le.this%kblc(ii_1+1))then
               ka=ii_1
               ia=l-this%kblc(ii_1)
               exit
              endif
             enddo

             do ii_1=1,this%ntot
              if(l2.le.this%kblc(ii_1+1))then
               kb=ii_1
               ib=l2-this%kblc(ii_1)
               exit
              endif
             enddo

             Gf=0.0d0 
             Gm=0.0d0
             Gp=0.0d0
             if(this%Ener(ka)%v(ia).gt.this%Ener(kb)%v(ib))then
              DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)-phondy%list(ph)%Freq(bn)
              Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              Gm=Gf
             else
              DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)+phondy%list(ph)%Freq(bn)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1.0d0)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
              Gp=Gf
             endif

!             write(*,*) ii,jj,Gf,Gm,Gp,ph,bn,phondy%list(ph)%freq(bn)

             R0%mat(ii,jj)=R0%mat(ii,jj)+dble(AA%mat(ii,jj)*conjg(AA%mat(ii,jj)))*Gf

!             write(*,*) ph,bn,phondy%list(ph)%freq(bn),ii,jj,&
!                      dble(AA%mat(ii,jj)*conjg(AA%mat(ii,jj))),dble(AA%mat(ii,jj)),aimag(AA%mat(ii,jj))

            enddo
           enddo

        ! make T2 

           do ii=1,size(this%T2%mat,1)
            do jj=1,size(this%T2%mat,2)

              l=indxl2g(ii,NB,myrow,0,nprow)
              l2=indxl2g(jj,MB,mycol,0,npcol)

              DEner=phondy%list(ph)%Freq(bn)
              Gf=bose(temp(1),phondy%list(ph)%Freq(bn))*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))

              DEner=phondy%list(ph)%Freq(bn)
              Gf=Gf+(bose(temp(1),phondy%list(ph)%Freq(bn))+1)*delta(type_smear,DEner,phondy%list(ph)%width(bn,1))

              this%T2%mat(ii,jj)=this%T2%mat(ii,jj)+(2*Vii(l)*Vii(l2)-Vii(l)*Vii(l)-Vii(l2)*Vii(l2))*Gf

            enddo
           enddo

!          enddo ! term
          enddo ! spin_id2
          enddo ! spin_id
          enddo ! bn

          if(allocated(phondy%list(ph)%hess)) deallocate(phondy%list(ph)%hess)
          if(allocated(phondy%list(ph)%width)) deallocate(phondy%list(ph)%width)

         ph=ph+1
         enddo ! ph

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(R0%mat(ii,jj),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           R0%mat(ii,jj)=valc
          enddo
         enddo

         do ii=1,size(this%T2%mat,1)
          do jj=1,size(this%T2%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(this%T2%mat(ii,jj),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           this%T2%mat(ii,jj)=valc
          enddo
         enddo
         
         call mpi_allreduce(nphonons,nze,1,&
              mpi_integer,mpi_sum,mpi_phonons_world,err)
                  
         if(mpi_id.eq.0) write(*,*) '      Number of phonons included: ',nphonons

         R0%mat=2.0d0*pi*pi*R0%mat/hplank
         this%T2%mat=pi*pi*this%T2%mat/hplank   

         do l=1,this%Hdim
          norm=0.0d0
          do l2=1,this%Hdim         
           if(l.ne.l2)then  
            call pdelget('A',' ',val,R0%mat,l2,l,R0%desc)
            norm=norm-val
           endif
          enddo
          this%T1(l)=norm
          call pdelset(R0%mat,l,l,R0%desc,norm)
         enddo

         if(this%printRmat)then
          if(mpi_id.eq.0) open(15,file='R.dat')
          do l=1,this%Hdim
           do l2=1,this%Hdim         
            call pdelget('A',' ',val,R0%mat,l2,l,R0%desc)
            if(mpi_id.eq.0) write(15,*) l2,l,dble(val)
           enddo
           if(mpi_id.eq.0) write(15,*)
          enddo
          if(mpi_id.eq.0) close(15)
         endif
 
!         if(mpi_id.eq.0)then
!          write(*,*) '           T1 (ps)'
!          do i=1,this%Hdim
!           write(*,*) i,-1.0d0/this%T1(i)
!          enddo
!         endif

      ! Compute correlation propagator
       
         do ii=1,size(this%T2%mat,1)
          do jj=1,size(this%T2%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           this%T2%mat(ii,jj)=this%T2%mat(ii,jj)+ &
                                   0.50d0*(this%T1(l)+this%T1(l2))

           this%T2%mat(ii,jj)=exp(step_min*mult_fact*this%T2%mat(ii,jj))

          enddo
         enddo

      ! Diag Rij ! check routine

         AA%mat=(0.0d0,0.0d0)
         allocate(this%Rval(this%Hdim))

         call BB%set(this%Hdim,this%Hdim,NB,MB)
         BB%mat=(0.0d0,0.0d0)
          
         !!!!! Scalapack version

         call pddiag2(this%Hdim,R0,this%Rval,AA,BB)
         BB%mat=AA%mat
         call pzgeinv(this%Hdim,BB)

         if(mpi_id.eq.0) open(15,file='Reig.dat')
         do l=1,this%Hdim
          do l2=1,this%Hdim         
           call pzelget('A',' ',valc,AA%mat,l2,l,AA%desc)
           if(mpi_id.eq.0) write(15,*) l2,l,dble(valc),aimag(valc)
          enddo
          if(mpi_id.eq.0) write(15,*)
         enddo
         if(mpi_id.eq.0) close(15)


         !!!!! Lapack version

!         AA%mat=cmplx(R0%mat,0.0d0,8)
!         call new_diag2(this%Hdim,AA%mat,this%Rval)
!         BB%mat=AA%mat
!         call mat_inv(BB%mat,this%Hdim) 

         !!!!!

         call CC%set(this%Hdim,this%Hdim,NB,MB)
         CC%mat=(0.0d0,0.0d0)

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),BB%mat,1,1,BB%desc,AA%mat,&
                     1,1,AA%desc,(0.0d0,0.0d0),CC%mat,1,1,CC%desc)

         if(mpi_id.eq.0) write(*,*) '     Redfield Matrix Eigenvalues:'

         if(mpi_id.eq.0)then
          allocate(rates(this%Hdim))
          rates=abs(dble(this%Rval))
          call  order_array(rates)
          do l=1,this%Hdim
           write(*,*) '        ',l,rates(l) 
          enddo
          deallocate(rates)
         endif

         nodiag_sum=(0.0d0,0.0d0)
         diag_sum=(0.0d0,0.0d0)

         do jj=1,size(CC%mat,2)
          do ii=1,size(CC%mat,1)
           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           if(l.eq.l2)  diag_sum=diag_sum+CC%mat(ii,jj)
           if(l.ne.l2)  nodiag_sum=nodiag_sum+CC%mat(ii,jj)
!           if(l.eq.l2 .and. mpi_id.eq.0 )  write(*,*) l,this%Rval(l) 

          enddo
         enddo

         call mpi_allreduce(diag_sum,diag_sum,1,&
              mpi_double_precision,MPI_SUM,mpi_blacs_world,err)    
         call mpi_allreduce(nodiag_sum,nodiag_sum,1,&
              mpi_double_precision,MPI_SUM,mpi_blacs_world,err)    

         if(mpi_id.eq.0) then
          write(*,*) '     Diagonal R/L Overlap:',diag_sum/this%Hdim,' it should be close to one!'
          write(*,*) '     Off-diagonal R/L Overlap:',nodiag_sum,' it should be close to zero!'
         endif

      ! Build Pop Propagator R= R exp(Rval) R^{\cross}

         this%Rval=exp(step_min*mult_fact*this%Rval)

         call this%R%set(this%Hdim,this%Hdim,NB,MB)
         this%R%mat=(0.0d0,0.0d0)

         do ii=1,size(AA%mat,1)
          do jj=1,size(AA%mat,2)
           l=indxl2g(jj,MB,mycol,0,npcol)
           AA%mat(ii,jj)=AA%mat(ii,jj)*this%Rval(l)
          enddo
         enddo

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),AA%mat,1,1,AA%desc,BB%mat,&
                     1,1,BB%desc,(0.0d0,0.0d0),this%R%mat,1,1,this%R%desc)

         call AA%dealloc()
         call BB%dealloc()
         
         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_R_H

        subroutine make_R2_H(this,sys,phondy,step_min,mult_fact,max_ener,euler)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use scalapack_diag_asimm
        use scalapack_inv
        use lapack_diag_asimm
        use lapack_inverse
        use units_parms 
        use atoms_class
        use phonons_class
        implicit none
        class(spins_hilbert)          :: this
        class(brillouin)              :: phondy
        class(atoms_group)            :: sys
        type(dist_cmplx_mat)          :: AA,BB,CC,AA2
        type(dist_dbl_mat)            :: R0,R0inv
        double precision              :: Gf,DEner,DEner2,step_min,coeff(3)
        double precision              :: val,norm,max_ener,avg_sph,Gpp,Gmm,Gpm,Gmp
        double precision, allocatable :: rates(:)
        double precision              :: euler(3)
        complex(8),allocatable        :: Vii(:),R0mtmp(:),R0ptmp(:)
        complex(8)                    :: valc,nodiag_sum,diag_sum,kcons
        complex(8)                    :: R0p,R0m
        integer                       :: ph,bn,ii,jj,l,l2,size_block,i1,ii_1,i       
        integer                       :: k,ia,ib,ka,kb,nphonons,v,spin_id,spin_id2
        integer                       :: l3,l4,kk1,kk2,bn2,ph2,phx,ic,id,kc,kd
        integer                       :: t1,t2,rate,indxl2g,px,py,pz,nloc,nstart
        integer                       :: mult_fact,s,t,j
        integer                       :: size_block_1,size_block_2,i2,l1,nze
        integer, allocatable          :: proc_grid(:)

         if(.not.this%make_Heig)then          
          write(*,*) 'Warning, open system symulation in the Sz basis', &
                     'is not implemented yet'
          return
         endif

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)
          write(*,*) '     Building the Redfield matrix: 2nd-order PT + 1st-order of coupling strength'
          write(*,*) '     The diagonal approximation to the secular Redfield equations will be used'
          flush(6)
         endif

!         allocate(R0mtmp(this%Hdim))
!         allocate(R0ptmp(this%Hdim))

         call R0%set(this%Hdim,this%Hdim,NB,MB)
         R0%mat=0.0d0

         call this%T2%set(this%Hdim,this%Hdim,NB,MB)
         this%T2%mat=(0.0d0,0.0d0)

         allocate(this%T1(this%Hdim))
         allocate(Vii(this%Hdim))

         call AA%set(this%Hdim,this%Hdim,NB,MB)          
         call AA2%set(this%Hdim,this%Hdim,NB,MB)          

         nphonons=0

         if(allocated(proc_grid)) deallocate(proc_grid)
         call mpi_dist_nprocess(size(phondy%list,1),nloc,nstart,proc_grid,mpi_phonons_world)

         ph=nstart
         do phx=1,nloc

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph)%width))then
            allocate(phondy%list(ph)%width(size(phondy%list(ph)%freq),1))
             phondy%list(ph)%width=smear
           endif
          endif

         do ph2=1,phondy%ntot
         
      ! Check for Ph lifetime

          if(.not.read_fc3)then
           if(.not.allocated(phondy%list(ph2)%width))then
            allocate(phondy%list(ph2)%width(size(phondy%list(ph2)%freq),1))
            phondy%list(ph2)%width=smear
           endif
          endif

          do bn=1,size(phondy%list(ph)%freq)
          do bn2=1,size(phondy%list(ph2)%freq)

      ! check spectrum overlap

           if (ph.eq.1 .and. bn.le.3) cycle
           if (ph2.eq.1 .and. bn2.le.3) cycle
           if (phondy%list(ph)%freq(bn).lt.0.0d0) cycle
           if (phondy%list(ph2)%freq(bn2).lt.0.0d0) cycle
           if (phondy%list(ph)%freq(bn).gt.max_ener) cycle
           if (phondy%list(ph2)%freq(bn2).gt.max_ener) cycle

!           DEner=abs(phondy%list(ph)%freq(bn)-phondy%list(ph2)%freq(bn2))
!           DEner2=abs(phondy%list(ph)%freq(bn)+phondy%list(ph2)%freq(bn2))
!           if (this%dos2pm%get_val(abs(DEner)).lt.1.0d-6 .and. &
!               this%dos2pm%get_val(abs(DEner2)).lt.1.0d-6 ) cycle

           nphonons=nphonons+1

           if(.not.allocated(phondy%list(ph)%hess) ) call phondy%list(ph)%diagD(sys)
           if(.not.allocated(phondy%list(ph2)%hess) ) call phondy%list(ph2)%diagD(sys)

      ! Build SPH per phonon

           do spin_id=1,this%nspins_pr
           do spin_id2=spin_id,this%nspins

      ! Make Vij

            call this%SPH%cart2brill(this%rcell,sys,phondy,ph,bn)
            if (any( abs(euler).gt.1.0e-4 )) call this%SPH%rot(euler)

            AA%mat=(0.0d0,0.0d0)

            call make_SH_rep(this,this%SPH,spin_id,spin_id2)

            do l=1,size(this%Hnodes,1)

             ii=this%Hnodes(l,1)
             jj=this%Hnodes(l,2)
            
             AA%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

            enddo ! l

      ! Rotate Vij

            if(this%ntot.gt.1)then
             stop
            endif           

            if(this%make_Heig)then
            call this%to_eigenbasis(AA)
            endif
                      
            do l=1,this%Hdim
             call pzelget('A',' ',Vii(l),AA%mat,l,l,AA%desc)
            enddo

       ! make Vij for second phonon 

            call this%SPH%cart2brill(this%rcell,sys,phondy,ph2,bn2)
            if (any( abs(euler).gt.1.0e-4 )) call this%SPH%rot(euler)

            AA2%mat=(0.0d0,0.0d0)

            call make_SH_rep(this,this%SPH,spin_id,spin_id2)

            do l=1,size(this%Hnodes,1)

             ii=this%Hnodes(l,1)
             jj=this%Hnodes(l,2)
            
             AA2%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

            enddo ! l

      ! Rotate Vij for second phonon

            if(this%ntot.gt.1)then
             write(*,*) 'kpoints and Raman not implemented'
             stop
            endif           

            if(this%make_Heig)then
            call this%to_eigenbasis(AA2)
            endif
                      
            do l=1,this%Hdim
             call pzelget('A',' ',Vii(l),AA%mat,l,l,AA%desc)
            enddo
 
       ! Make Rij

            do ii=1,size(R0%mat,1)
             do jj=1,size(R0%mat,2)

              l=indxl2g(ii,NB,myrow,0,nprow)
              l2=indxl2g(jj,MB,mycol,0,npcol)

              do ii_1=1,this%ntot
               if(l.le.this%kblc(ii_1+1))then
                ka=ii_1
                ia=l-this%kblc(ii_1)
                exit
               endif
              enddo

              do ii_1=1,this%ntot
               if(l2.le.this%kblc(ii_1+1))then
                kb=ii_1
                ib=l2-this%kblc(ii_1)
                exit
               endif
              enddo

              R0p=(0.0d0,0.0d0)
              R0m=(0.0d0,0.0d0)
!              R0mtmp=(0.0d0,0.0d0)
!              R0ptmp=(0.0d0,0.0d0)
 
              do kk1=1,size(R0%mat,1)
               do kk2=1,size(R0%mat,2)

                l3=indxl2g(kk1,NB,myrow,0,nprow) 
                l4=indxl2g(kk2,MB,mycol,0,npcol)

                if(l3.eq.l4)then

                 do ii_1=1,this%ntot
                  if(l3.le.this%kblc(ii_1+1))then
                   kc=ii_1
                   ic=l3-this%kblc(ii_1)
                   exit
                  endif
                 enddo

                 do ii_1=1,this%ntot
                  if(l4.le.this%kblc(ii_1+1))then
                   kd=ii_1
                   id=l4-this%kblc(ii_1)
                   exit
                  endif
                 enddo
              

                 if(this%ener(kc)%v(ic).gt.this%Ener(kb)%v(ib))then  
                  R0m=R0m+AA%mat(ii,kk2)*AA2%mat(kk1,jj)&
                   /(this%Ener(kc)%v(ic)-this%Ener(kb)%v(ib)-phondy%list(ph2)%freq(bn2)&
                        +cmplx(0.0d0,1.0d0,8)*phondy%list(ph2)%width(bn2,1))

!                  R0mtmp(ic)=AA%mat(ii,kk2)*AA2%mat(kk1,jj)&
!                   /(this%Ener(kc)%v(ic)-this%Ener(kb)%v(ib)-phondy%list(ph2)%freq(bn2)&
!                        +cmplx(0.0d0,1.0d0,8)*phondy%list(ph2)%width(bn2,1))
                 else
                  R0p=R0p+AA%mat(ii,kk2)*AA2%mat(kk1,jj)&
                   /(this%Ener(kc)%v(ic)-this%Ener(kb)%v(ib)+phondy%list(ph2)%freq(bn2)&
                        +cmplx(0.0d0,1.0d0,8)*phondy%list(ph2)%width(bn2,1))

!                  R0ptmp(ic)=AA%mat(ii,kk2)*AA2%mat(kk1,jj)&
!                   /(this%Ener(kc)%v(ic)-this%Ener(kb)%v(ib)+phondy%list(ph2)%freq(bn2)&
!                        +cmplx(0.0d0,1.0d0,8)*phondy%list(ph2)%width(bn2,1))
                 endif

                endif              

               enddo ! kk1
              enddo ! kk2
              
              Gf=0.0d0 
              Gmp=0.0d0
              Gpp=0.0d0

              if(this%ener(ka)%v(ia).gt.this%Ener(kb)%v(ib) .and. &
                 phondy%list(ph2)%freq(bn2).gt.phondy%list(ph)%freq(bn)  )then
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)-phondy%list(ph2)%Freq(bn2)+phondy%list(ph)%Freq(bn)
                  Gf=bose(temp(1),phondy%list(ph2)%Freq(bn2))*(bose(temp(1),phondy%list(ph)%Freq(bn))+1)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
!                  Gmp=Gf
              endif

              if(this%ener(ka)%v(ia).lt.this%Ener(kb)%v(ib) .and. &
                 phondy%list(ph2)%freq(bn2).lt.phondy%list(ph)%freq(bn)  )then  
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)-phondy%list(ph2)%Freq(bn2)+phondy%list(ph)%Freq(bn)
                  Gf=Gf+bose(temp(1),phondy%list(ph2)%Freq(bn2))*(bose(temp(1),phondy%list(ph)%Freq(bn))+1)*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                  Gmp=Gf
              endif

              if(this%ener(ka)%v(ia).gt.this%Ener(kb)%v(ib))then
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)-phondy%list(ph2)%Freq(bn2)-phondy%list(ph)%Freq(bn)
                  Gf=Gf+bose(temp(1),phondy%list(ph2)%Freq(bn2))*bose(temp(1),phondy%list(ph)%Freq(bn))*&
                     delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                  Gpp=Gf-Gmp
              endif

              R0%mat(ii,jj)=R0%mat(ii,jj)+dble(R0m*conjg(R0m))*Gf
        
 !             if(Gf.gt.1.0d-8)then
!               if(ii.eq.1 .and. jj.eq.2)then
!               write(*,*) ii,jj,'-+,++',Gmp,Gpp,phondy%list(ph)%freq(bn),phondy%list(ph2)%freq(bn2)
!                do kk1=1,this%Hdim
!                 write(*,*) '-','+',kk1,dble(conjg(R0mtmp(kk1))*R0mtmp(kk1)),dble(conjg(R0ptmp(kk1))*R0ptmp(kk1))
!                enddo
!               endif
!               if(ii.eq.2 .and. jj.eq.1)then
!               write(*,*) ii,jj,'-+,++',Gmp,Gpp,phondy%list(ph)%freq(bn),phondy%list(ph2)%freq(bn2)
!               do kk1=1,this%Hdim
!                write(*,*) '-','+',kk1,dble(conjg(R0mtmp(kk1))*R0mtmp(kk1)),dble(conjg(R0ptmp(kk1))*R0ptmp(kk1))
!               enddo
!               endif
!              endif

              Gf=0.0d0
              Gpm=0.0d0
              Gmm=0.0d0
 
              if(this%ener(ka)%v(ia).lt.this%Ener(kb)%v(ib))then
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)+phondy%list(ph2)%Freq(bn2)+phondy%list(ph)%Freq(bn)
                  Gf=(bose(temp(1),phondy%list(ph2)%Freq(bn2))+1)*(bose(temp(1),phondy%list(ph)%Freq(bn))+1)*&
                      delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                  Gmm=Gf
              endif

              if(this%ener(ka)%v(ia).lt.this%Ener(kb)%v(ib) .and. &
                 phondy%list(ph2)%freq(bn2).gt.phondy%list(ph)%freq(bn)  )then  
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)+phondy%list(ph2)%Freq(bn2)-phondy%list(ph)%Freq(bn)
                  Gf=Gf+(bose(temp(1),phondy%list(ph2)%Freq(bn2))+1)*bose(temp(1),phondy%list(ph)%Freq(bn))*&
                      delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                  Gpm=Gf-Gmm
              endif

              if(this%ener(ka)%v(ia).gt.this%Ener(kb)%v(ib) .and. &
                 phondy%list(ph2)%freq(bn2).lt.phondy%list(ph)%freq(bn)  )then  
                  DEner=this%Ener(ka)%v(ia)-this%Ener(kb)%v(ib)+phondy%list(ph2)%Freq(bn2)-phondy%list(ph)%Freq(bn)
                  Gf=Gf+(bose(temp(1),phondy%list(ph2)%Freq(bn2))+1)*bose(temp(1),phondy%list(ph)%Freq(bn))*&
                      delta(type_smear,DEner,phondy%list(ph)%width(bn,1))
                  Gpm=Gf-Gmm
              endif
 
              R0%mat(ii,jj)=R0%mat(ii,jj)+dble(R0p*conjg(R0p))*Gf

!              if(Gf.gt.1.0d-8)then
!               if(ii.eq.1 .and. jj.eq.2)then
!                write(*,*) ii,jj,'--,+-',Gmm,Gpm,phondy%list(ph)%freq(bn),phondy%list(ph2)%freq(bn2)
!                do kk1=1,this%Hdim
!                 write(*,*) '-','+',kk1,dble(conjg(R0mtmp(kk1))*R0mtmp(kk1)),dble(conjg(R0ptmp(kk1))*R0ptmp(kk1))
!                enddo
!               endif
!               if(ii.eq.2 .and. jj.eq.1)then
!                write(*,*) ii,jj,'--,+-',Gmm,Gpm,phondy%list(ph)%freq(bn),phondy%list(ph2)%freq(bn2)
!                do kk1=1,this%Hdim
!                 write(*,*) '-','+',kk1,dble(conjg(R0mtmp(kk1))*R0mtmp(kk1)),dble(conjg(R0ptmp(kk1))*R0ptmp(kk1))
!                enddo
!               endif
!              endif

             enddo ! jj
            enddo ! ii

        ! make T2 To be implemented
                
          enddo ! spin_id2
          enddo ! spin_id

          enddo ! bn2
          enddo ! bn

          if(ph2.ne.ph)then
           if(allocated(phondy%list(ph2)%hess)) deallocate(phondy%list(ph2)%hess)
           if(allocated(phondy%list(ph2)%width)) deallocate(phondy%list(ph2)%width)
          endif

         enddo ! ph2

         if(allocated(phondy%list(ph)%hess)) deallocate(phondy%list(ph)%hess)
         if(allocated(phondy%list(ph)%width)) deallocate(phondy%list(ph)%width)

         ph=ph+1
         enddo ! ph

         do ii=1,size(R0%mat,1)
          do jj=1,size(R0%mat,2)
           valc=(0.0d0,0.0d0)
           call mpi_allreduce(R0%mat(ii,jj),valc,1,&
              mpi_double_complex,mpi_sum,mpi_phonons_world,err)
           R0%mat(ii,jj)=valc
          enddo
         enddo
         
         call mpi_allreduce(nphonons,nze,1,&
              mpi_integer,mpi_sum,mpi_phonons_world,err)

         nphonons=nze
         
         if(mpi_id.eq.0) write(*,*) '     Total number of phonons included: ',nphonons
         
         R0%mat=pi*pi*R0%mat/hplank
         this%T2%mat=pi*pi*this%T2%mat/hplank/2.0d0

         do l=1,this%Hdim
          norm=0.0d0
          do l2=1,this%Hdim         
           if(l.ne.l2)then  
            call pdelget('A',' ',val,R0%mat,l2,l,R0%desc)
            norm=norm-val
           endif
          enddo
          this%T1(l)=norm
          call pdelset(R0%mat,l,l,R0%desc,norm)
         enddo

         if(this%printRmat)then
          if(mpi_id.eq.0) open(15,file='R.dat')
          do l=1,this%Hdim
           do l2=1,this%Hdim         
            call pdelget('A',' ',val,R0%mat,l2,l,R0%desc)
            if(mpi_id.eq.0) write(15,*) l2,l,dble(val)
           enddo
           if(mpi_id.eq.0) write(15,*)
          enddo
          if(mpi_id.eq.0) close(15)
         endif

!         if(mpi_id.eq.0)then
!          write(*,*) '           T1 (ps)'
!          do i=1,this%Hdim
!           write(*,*) i,-1.0d0/this%T1(i)
!          enddo
!         endif

      ! Compute correlation propagator
       
         do ii=1,size(this%T2%mat,1)
          do jj=1,size(this%T2%mat,2)

           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           this%T2%mat(ii,jj)=this%T2%mat(ii,jj)+ &
                                   0.50d0*(this%T1(l)+this%T1(l2))

           this%T2%mat(ii,jj)=exp(step_min*mult_fact*this%T2%mat(ii,jj))

          enddo
         enddo

      ! Diag Rij ! check routine

         AA%mat=(0.0d0,0.0d0)
         allocate(this%Rval(this%Hdim))

         call BB%set(this%Hdim,this%Hdim,NB,MB)
         BB%mat=(0.0d0,0.0d0)
          
         !!!!! Scalapack version

         call pddiag2(this%Hdim,R0,this%Rval,AA,BB)
         BB%mat=AA%mat
         call pzgeinv(this%Hdim,BB)

         if(mpi_id.eq.0) open(15,file='Reig.dat')
         do l=1,this%Hdim
          do l2=1,this%Hdim         
           call pzelget('A',' ',valc,AA%mat,l2,l,AA%desc)
           if(mpi_id.eq.0) write(15,*) l2,l,dble(valc),aimag(valc)
          enddo
          if(mpi_id.eq.0) write(15,*)
         enddo
         if(mpi_id.eq.0) close(15)


         !!!!! Lapack version

!         AA%mat=cmplx(R0%mat,0.0d0,8)
!         call new_diag2(this%Hdim,AA%mat,this%Rval)
!         BB%mat=AA%mat
!         call mat_inv(BB%mat,this%Hdim) 

         !!!!!

         call CC%set(this%Hdim,this%Hdim,NB,MB)
         CC%mat=(0.0d0,0.0d0)

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),BB%mat,1,1,BB%desc,AA%mat,&
                     1,1,AA%desc,(0.0d0,0.0d0),CC%mat,1,1,CC%desc)

         if(mpi_id.eq.0) write(*,*) '     Redfield Matrix Eigenvalues:'

         nodiag_sum=(0.0d0,0.0d0)
         diag_sum=(0.0d0,0.0d0)

         if(mpi_id.eq.0)then
          allocate(rates(this%Hdim))
          rates=abs(dble(this%Rval))
          call  order_array(rates)
          do l=1,this%Hdim
           write(*,*) '          ',l,rates(l) 
          enddo
          deallocate(rates)
         endif

         do jj=1,size(CC%mat,2)
          do ii=1,size(CC%mat,1)
           l=indxl2g(ii,NB,myrow,0,nprow)
           l2=indxl2g(jj,MB,mycol,0,npcol)

           if(l.eq.l2)  diag_sum=diag_sum+CC%mat(ii,jj)
           if(l.ne.l2)  nodiag_sum=nodiag_sum+CC%mat(ii,jj)

          enddo
         enddo

         call mpi_allreduce(diag_sum,diag_sum,1,&
              mpi_double_precision,MPI_SUM,mpi_blacs_world,err)    
         call mpi_allreduce(nodiag_sum,nodiag_sum,1,&
              mpi_double_precision,MPI_SUM,mpi_blacs_world,err)    

         if(mpi_id.eq.0) then
          write(*,*) '     Diagonal Right-Left Overlap:',diag_sum/this%Hdim
          write(*,*) '     Out of Diagonal Right-Left Overlap:',nodiag_sum
         endif

      ! Build Pop Propagator R= R exp(Rval) R^{\cross}

         this%Rval=exp(step_min*mult_fact*this%Rval)

         call this%R%set(this%Hdim,this%Hdim,NB,MB)
         this%R%mat=(0.0d0,0.0d0)

         do ii=1,size(AA%mat,1)
          do jj=1,size(AA%mat,2)
           l=indxl2g(jj,MB,mycol,0,npcol)
           AA%mat(ii,jj)=AA%mat(ii,jj)*this%Rval(l)
          enddo
         enddo

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),AA%mat,1,1,AA%desc,BB%mat,&
                     1,1,BB%desc,(0.0d0,0.0d0),this%R%mat,1,1,this%R%desc)

         call AA%dealloc()
         call AA2%dealloc()
         call BB%dealloc()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '    Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_R2_H

        subroutine propagate_H(this,start_step,time,nsteps,step,dump_freq)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms
        implicit none
        class(spins_hilbert)         :: this
        type(dist_cmplx_mat)         :: AA
        type(dist_cmplx_vec)         :: pop,pop_new
        double precision             :: step,time
        complex(8)                   :: val
        integer                      :: nsteps,dump_freq,start_step
        integer                      :: l,v,k,i,ii,jj,indxl2g,j
        integer                      :: block1,block2,i1,i2,j1,j2
        integer                      :: kp1,kp2,kpt,ii_1,jj_1
        integer                      :: t1,t2,rate

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Propagation timestep  : ',step,' ps'
          write(*,*) '     Total propagation time: ',nsteps*step,' ps'
          flush(6)
         endif
         
         if(start_step.eq.0)then
          if(this%s2print.gt.0) call this%dump_M(0,time)
          if(this%printRho) call this%dump_rho()
          start_step=1
         endif

         if(this%make_Rmat .or. this%make_R2mat)then
          if(this%make_SA)then
           call pop%set(this%Hdim,NB)
           pop%vec=(0.0d0,0.0d0)
           call pop_new%set(this%Hdim,NB)
           pop_new%vec=(0.0d0,0.0d0)
          else
           call pop%set(this%Ldim,NB)
           pop%vec=(0.0d0,0.0d0)
           call pop_new%set(this%Ldim,NB)
           pop_new%vec=(0.0d0,0.0d0)
          endif
         endif

         do i=start_step,nsteps+start_step-1

          time=time+step

          if(this%make_Heig)then

           do ii=1,size(this%rho%mat,1)
            do jj=1,size(this%rho%mat,2)

             l=indxl2g(ii,NB,myrow,0,nprow)
             v=indxl2g(jj,MB,mycol,0,npcol)
        
             do ii_1=1,this%ntot
              if(l.le.this%kblc(ii_1+1))then
               kp1=ii_1
               j=l-this%kblc(ii_1)
               exit
              endif        
             enddo

             do ii_1=1,this%ntot             
              if(v.le.this%kblc(ii_1+1))then
               kp2=ii_1
               k=v-this%kblc(ii_1)
               exit
              endif        
             enddo
             
             this%rho%mat(ii,jj)=this%rho%mat(ii,jj) &
                        *this%U(kp1)%mat(1,j)*conjg(this%U(kp2)%mat(1,k))
                          
            enddo
           enddo

           if(this%make_Rmat .or. this%make_R2mat)then

            if(this%make_SA)then

            ! propagate populations

            do l=1,this%Hdim
             call pzelget('A',' ',val,this%rho%mat,l,l,this%rho%desc)
             call pzelset(pop%vec,l,1,pop%desc,val)
            enddo

            call pzgemv('N',this%Hdim,this%Hdim,&
                        (1.0d0,0.0d0),this%R%mat,1,1,this%R%desc,pop%vec,1,1,pop%desc,&
                        1,(0.0d0,0d0),pop_new%vec,1,1,pop_new%desc,1)

            do l=1,this%Hdim
             call pzelget('A',' ',val,pop_new%vec,l,1,pop_new%desc)
             call pzelset(this%rho%mat,l,l,this%rho%desc,val)
            enddo   
        
            ! propagate correlations

            do ii=1,size(this%rho%mat,1)
             do jj=1,size(this%rho%mat,2)

              l=indxl2g(ii,NB,myrow,0,nprow)
              v=indxl2g(jj,MB,mycol,0,npcol)

              if(l.ne.v)then
               this%rho%mat(ii,jj)=this%T2%mat(ii,jj)*this%rho%mat(ii,jj)
              endif

             enddo
            enddo

            else

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


            endif ! end SA if
           endif ! end Rmat if

          else

           if(this%make_Rmat)then
            write(*,*) 'Warning, open system symulation in the Szbasis',  &
                       'is not implemented yet'
           endif
          
           call AA%set(this%Hdim,this%Hdim,NB,MB)

           do kp1=1,this%ntot
            do kp2=1,this%ntot

             block1=this%kblc(kp1+1)-this%kblc(kp1)
             block2=this%kblc(kp2+1)-this%kblc(kp2)
             i1=this%kblc(kp1)+1
             j1=this%kblc(kp1)+1
             i2=this%kblc(kp2)+1
             j2=this%kblc(kp2)+1

             call pzgemm('N','N',block1,block2,block1,&
                         (1.0d0,0.0d0),this%U(kp1)%mat,1,1,this%U(kp1)%desc,this%rho%mat,&
                         i1,j2,this%rho%desc,(0.0d0,0.0d0),AA%mat,i1,j2,AA%desc)

            enddo
           enddo

           do kp1=1,this%ntot
            do kp2=1,this%ntot

             block1=this%kblc(kp1+1)-this%kblc(kp1)
             block2=this%kblc(kp2+1)-this%kblc(kp2)
             i1=this%kblc(kp1)+1
             j1=this%kblc(kp1)+1
             i2=this%kblc(kp2)+1
             j2=this%kblc(kp2)+1

             call pzgemm('N','C',block1,block2,block2,&
                         (1.0d0,0.0d0),AA%mat,i1,j2,AA%desc,this%U(kp2)%mat,1,1,this%U(kp2)%desc,&
                         (0.0d0,0.0d0),this%rho%mat,i1,j2,this%rho%desc) 

            enddo
           enddo

           call AA%dealloc()

          endif

          if ( mod(i,dump_freq).eq.0 ) then
           if(this%s2print .gt. 0) call this%dump_M(i,time)
           if(this%printRho)  call this%dump_rho()
          endif

         enddo

         if ( mod(i,dump_freq).eq.0 ) then
          if(this%printRho) call this%dump_rho()
         endif
         
         start_step=start_step+nsteps
         if(this%make_Rmat) call pop%dealloc()
         if(this%make_Rmat) call pop_new%dealloc()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine propagate_H


        subroutine dump_M_H(this,i,time)
        use mpi
        use mpi_utils
        implicit none
        class(spins_hilbert)   :: this
        integer                :: i,k,unit_no,nze
        double precision       :: time
        complex(8)             :: val,norm
        double precision       :: M(3),Mi(3,this%s2print)
        
         call this%make_M(M,Mi)

         norm=(0.0d0,0.0d0)
         do k=1,this%Hdim
          val=(0.0d0,0.0d0)
          call pzelget(' ',' ',val,this%rho%mat,k,k,this%rho%desc)
          norm=norm+val
         enddo

         call mpi_allreduce(norm,norm,1,mpi_double_complex,mpi_sum,mpi_blacs_world,err)

         nze=this%rho%get_nze(1.0d-4)

         if(mpi_id.eq.0)then

          if(this%s2print.gt.0)then

           inquire(file='Mx_dynamics.dat',number=unit_no)
           if(unit_no.eq.-1)  open(11,file='Mx_dynamics.dat')
           inquire(file='My_dynamics.dat',number=unit_no)
           if(unit_no.eq.-1)  open(12,file='My_dynamics.dat')
           inquire(file='Mz_dynamics.dat',number=unit_no)
           if(unit_no.eq.-1)  open(13,file='Mz_dynamics.dat')

           write(11,*) i,time,Mi(1,:),dble(norm),nze
           flush(11)
           write(12,*) i,time,Mi(2,:),dble(norm),nze
           flush(12)
           write(13,*) i,time,Mi(3,:),dble(norm),nze
           flush(13)

          endif

         endif

        return
        end subroutine dump_M_H

        subroutine dump_kbasis_H(this)
        use mpi
        use mpi_utils
        implicit none
        class(spins_hilbert)   :: this
        integer                :: l,j,k,size_block
        complex(8)             :: val

         if(mpi_id.eq.0) open(15,file='Kbasis.dat')

         do l=1,this%Hdim
          do j=1,this%Hdim
           call pzelget('A',' ',val,this%kbasis%mat,l,j,this%kbasis%desc)
           if(mpi_id.eq.0) write(15,*) k,l,j,dble(val),aimag(val)
          enddo
         enddo

         if(mpi_id.eq.0) close(15)

        return
        end subroutine dump_kbasis_H

        subroutine dump_Heig_H(this)
        use mpi
        use mpi_utils
        implicit none
        class(spins_hilbert)   :: this
        integer                :: l,j,k,size_block
        complex(8)             :: val

         if(mpi_id.eq.0) open(15,file='Heig.dat')

         do k=1,this%ntot
          size_block=this%kblc(k+1)-this%kblc(k)
          do j=1,size_block
           do l=1,size_block
            call pzelget('A',' ',val,this%H(k)%mat,l,j,this%H(k)%desc)
            if(mpi_id.eq.0) write(15,*) k,l,j,dble(val),aimag(val),dble(val*conjg(val))
           enddo
            if(mpi_id.eq.0) write(15,*) 
          enddo
         enddo

         if(mpi_id.eq.0) close(15)

        return
        end subroutine dump_Heig_H

        subroutine dump_H0_H(this)
        use mpi
        use mpi_utils
        implicit none
        class(spins_hilbert)   :: this
        integer                :: l,j,k,size_block
        complex(8)             :: val

         if(mpi_id.eq.0) open(15,file='H0.dat')

         do k=1,this%ntot
          size_block=this%kblc(k+1)-this%kblc(k)
          do l=1,size_block
           do j=1,size_block
            call pzelget('A',' ',val,this%H0(k)%mat,l,j,this%H0(k)%desc)
            if(mpi_id.eq.0) write(15,*) k,l,j,dble(val),aimag(val)
           enddo
           if(mpi_id.eq.0) write(15,*) 
          enddo
         enddo

         if(mpi_id.eq.0) flush(15)
         if(mpi_id.eq.0) close(15)

        return
        end subroutine dump_H0_H

        subroutine dump_rho_H(this)
        use mpi
        use mpi_utils
        implicit none
        class(spins_hilbert)   :: this
        integer                :: l,j
        complex(8)             :: val

         if(mpi_id.eq.0) open(15,file='rho.dat')

         do l=1,this%Hdim
          do j=1,this%Hdim         
           call pzelget('A',' ',val,this%rho%mat,l,j,this%rho%desc)
           if(mpi_id.eq.0) write(15,*) l,j,dble(val),aimag(val)
          enddo
          if(mpi_id.eq.0) write(15,*) 
         enddo

         if(mpi_id.eq.0) close(15)

        return
        end subroutine dump_rho_H

        subroutine make_propagator_H(this,step_min,mult_fact)
        use mpi
        use mpi_utils
        use blacs_utils
        use units_parms
        implicit none
        class(spins_hilbert)    :: this
        integer                 :: mult_fact,indxl2g
        integer                 :: i,k,l,ll,v,j
        integer                 :: ii,jj,kpt,size_block
        double precision        :: step_min
        type(dist_cmplx_mat)    :: AA,BB
        complex(8)              :: val,coeff,val_tot
        integer                 :: t1,t2,rate
        logical, allocatable    :: check(:)


         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building unitary propagator operator'
          flush(6)
         endif

         if(allocated(this%U)) deallocate(this%U)
         allocate(this%U(this%ntot))
 
         if(this%make_Heig)then

          do kpt=1,this%ntot

           size_block=this%kblc(kpt+1)-this%kblc(kpt)
           if(allocated(this%U(kpt)%mat)) call this%U(kpt)%dealloc()
           allocate(this%U(kpt)%mat(1,size_block))
           this%U(kpt)%mat=(0.0d0,0.0d0)

           do v=1,size_block
            this%U(kpt)%mat(1,v)=-step_min*mult_fact*cmplx(0.0d0,1.0d0,8)*this%Ener(kpt)%v(v)*2*acos(-1.0d0)/hplank
            this%U(kpt)%mat(1,v)=exp(this%U(kpt)%mat(1,v))
           enddo

          enddo

         else

          allocate(check(mpi_nproc))

          do kpt=1,this%ntot

           size_block=this%kblc(kpt+1)-this%kblc(kpt)
           if(allocated(this%U(kpt)%mat)) call this%U(kpt)%dealloc()
           call this%U(kpt)%set(size_block,size_block,NB,MB)
           

           call AA%set(size_block,size_block,NB,MB)
           call BB%set(size_block,size_block,NB,MB)

           do ii=1,size(this%U(kpt)%mat,1)
            do jj=1,size(this%U(kpt)%mat,2)
             i=indxl2g(ii,NB,myrow,0,nprow)
             j=indxl2g(jj,MB,mycol,0,npcol)
             if(i.ne.j)then
              this%U(kpt)%mat(ii,jj)=(0.0d0,0.0d0)
              AA%mat(ii,jj)=(0.0d0,0.0d0)
             else
              this%U(kpt)%mat(ii,jj)=(1.0d0,0.0d0)
              AA%mat(ii,jj)=(1.0d0,0.0d0)
             endif
            enddo
           enddo
                

        !  Taylor expansion for the starting Propagator               

           i=0
           check=.false.
           do while(.not. all(check) )
           
            i=i+1
            coeff=(1.0d0/dble(i))*(-1.0d0*step_min*cmplx(0.0d0,1.0d0,8)*2.0d0*pi/hplank)

            call pzgemm('N','N',size_block,size_block,size_block,&
                         coeff,AA%mat,1,1,AA%desc,this%H0(kpt)%mat,1,1,this%H0(kpt)%desc,&
                         (0.0d0,0.0d0),BB%mat,1,1,BB%desc) 

            do ii=1,size(this%U(kpt)%mat,1)
             do jj=1,size(this%U(kpt)%mat,2)
              this%U(kpt)%mat(ii,jj)=this%U(kpt)%mat(ii,jj)+BB%mat(ii,jj)
             enddo
            enddo

            i=i+1
            coeff=(1.0d0/dble(i))*(-1.0d0*step_min*cmplx(0.0d0,1.0d0,8)*2.0d0*pi/hplank)

            call pzgemm('N','N',size_block,size_block,size_block,&
                         coeff,BB%mat,1,1,BB%desc,this%H0(kpt)%mat,1,1,this%H0(kpt)%desc,&
                         (0.0d0,0.0d0),AA%mat,1,1,AA%desc) 
        

            do ii=1,size(this%U(kpt)%mat,1)
             do jj=1,size(this%U(kpt)%mat,2)
              this%U(kpt)%mat(ii,jj)=this%U(kpt)%mat(ii,jj)+AA%mat(ii,jj)
             enddo
            enddo

            if( maxval(abs(dble(AA%mat))).lt.1.0d-18 .and. &
                maxval(abs(aimag(AA%mat))).lt.1.0d-18 ) then
             check(mpi_id+1)=.true.
            else
             check(mpi_id+1)=.false.
            endif

            do l=0,mpi_nproc-1
             call mpi_bcast(check(l+1),1,mpi_logical,l,mpi_blacs_world,err)
            enddo

           enddo

           if(mpi_id.eq.0) write(*,*) '     Taylor expansion converged in ',i,'steps'

        !  Quadrature of the propagator

           BB%mat=this%U(kpt)%mat

           do i=1,mult_fact-1

            call pzgemm('N','N',size_block,size_block,size_block,&
                         (1.0d0,0.0d0),BB%mat,1,1,BB%desc,this%U(kpt)%mat,1,1,this%U(kpt)%desc,&
                         (0.0d0,0.0d0),AA%mat,1,1,AA%desc) 

            BB%mat=AA%mat
            
           enddo

           this%U(kpt)%mat=BB%mat

           call AA%dealloc()
           call BB%dealloc()

           ii=this%U(kpt)%get_nze(1.0d-9)
           if(mpi_id.eq.0)  &
           write(*,*) '     Propagator Sparsity: ',10*(1-ii/dble(size_block**2)),'%'

          enddo ! end on kpt

         endif
         

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_propagator_H


        subroutine make_rho0_i(si,rho_i,alpha0,beta0,gamma0)
        use rotations_class
        implicit none
        complex(8), allocatable   :: rho_i(:,:),rot(:,:)
        double precision          :: alpha0,beta0,gamma0,si
        integer                   :: ni,i

         ni=nint(2*si+1)
         if(allocated(rho_i)) deallocate(rho_i)
         allocate(rho_i(ni,ni))

         rho_i=(0.0d0,0.0d0)
         rho_i(1,1)=(1.0d0,0.0d0)                

         call rot_wig(si,alpha0,beta0,gamma0,rot)

         rho_i=matmul(rot,rho_i)
         rho_i=matmul(rho_i,conjg(transpose(rot)))

        return
        end subroutine make_rho0_i

        subroutine make_rho0_H(this,type_rho0,temp,rho_restart_file)
        use mpi
        use mpi_utils
        use blacs_utils
        use sparse_class
        use units_parms 
        implicit none
        class(spins_hilbert)      :: this
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

         if(allocated(this%rho0%mat)) call this%rho0%dealloc()
         call this%rho0%set(this%Hdim,this%Hdim,NB,MB)
         this%rho0%mat=(0.0d0,0.0d0)

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
               
           if(this%make_Heig)then
        
            do kpt=1,this%ntot
             size_block=this%kblc(kpt+1)-this%kblc(kpt)
             do l=1,size_block
              i=this%kblc(kpt)+l 
              call pzelset(this%rho%mat,i,i,this%rho%desc, &
                       cmplx(exp(-this%Ener(kpt)%v(l)/(temp*kboltz)),0.0d0,8))
             enddo
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

           else

         !! implement for non eigenbasis

           endif

          case ('FULLY_POLARIZED')

           do ii=1,size(this%rho%mat,1)
            do jj=1,size(this%rho%mat,2)
             i=indxl2g(ii,NB,myrow,0,nprow)
             j=indxl2g(jj,MB,mycol,0,npcol)
             if(i.eq.1 .and. j.eq.1) this%rho%mat(ii,jj)=(1.0d0,0.0d0)
            enddo
           enddo

           if(this%sparse)then
            call this%sp_rho%tosparse(this%Hdim,this%rho)
            if(this%ntot.gt.1) then
             call mult2_sparse_cmplx(.false.,this%sp_rho,this%sp_kbasis,sprho)
             call this%sp_rho%delete()
             call mult2_sparse_cmplx(.true.,this%sp_kbasis,sprho,this%sp_rho)
             call sprho%delete()
            endif
             call this%rho%dealloc()
             call this%sp_rho%todense(this%rho) 
            else
             if(this%ntot.gt.1) then
              call this%to_kbasis(this%rho)
             endif
            endif
            if(this%make_Heig) call this%to_eigenbasis(this%rho)

          case ('GENERAL_ENTANGLE')

          case ('SPIN_BATHS')

          case ('EIGENSTATES_MIX')

           this%rho%mat=(0.0d0,0.0d0)
           this%rho%mat(1,1)=cmplx(0.5d0,0.0d0,8)
           this%rho%mat(2,1)=cmplx(0.5d0,0.0d0,8)
           this%rho%mat(1,2)=cmplx(0.5d0,0.0d0,8)
           this%rho%mat(2,2)=cmplx(0.5d0,0.0d0,8)

          case ('NO_ENTANGLE')

           do v=1,this%nspins                       
            call make_rho0_i(this%spin(this%kind(v)),rho_i,this%alpha0(v),this%beta0(v),this%gamma0(v))
             write(*,*) v,rho_i(1,1),rho_i(1,2)
             write(*,*) v,rho_i(2,1),rho_i(2,2)
             do ii=1,size(this%rho%mat,1)
              do jj=1,size(this%rho%mat,2)
               i=indxl2g(ii,NB,myrow,0,nprow)
               j=indxl2g(jj,MB,mycol,0,npcol)
                k=nint(  this%spin(this%kind(v))+1+this%basis(i,v) )
                l=nint(  this%spin(this%kind(v))+1+this%basis(j,v) )
                if(v.eq.1)then          
                 this%rho%mat(ii,jj)=rho_i(k,l)
                else
                 this%rho%mat(ii,jj)=this%rho%mat(ii,jj)*rho_i(k,l)
                endif
               enddo
              enddo
            enddo
         
           if(this%sparse)then
            call this%sp_rho%tosparse(this%Hdim,this%rho)
            if(this%ntot.gt.1) then
             call mult2_sparse_cmplx(.false.,this%sp_rho,this%sp_kbasis,sprho)
             call this%sp_rho%delete()
             call mult2_sparse_cmplx(.true.,this%sp_kbasis,sprho,this%sp_rho)
             call sprho%delete()
            endif
             call this%rho%dealloc()
             call this%sp_rho%todense(this%rho) 
           else
            if(this%ntot.gt.1) then
             call this%to_kbasis(this%rho)
            endif
           endif
           if(this%make_Heig) call this%to_eigenbasis(this%rho)

           deallocate(rho_i)

         end select

         this%rho0%mat=this%rho%mat

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_rho0_H

        subroutine to_eigenbasis(this,mat)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)      :: this
        type(dist_cmplx_mat)      :: AA,mat
        integer                   :: kp1,kp2,block1,block2
        integer                   :: i1,i2,j1,j2,k

         call AA%set(this%Hdim,this%Hdim,NB,MB)
         AA%mat=(0.0d0,0.0d0)

         do kp1=1,this%ntot
          do kp2=1,this%ntot

           block1=this%kblc(kp1+1)-this%kblc(kp1)
           block2=this%kblc(kp2+1)-this%kblc(kp2)
           i1=this%kblc(kp1)+1
           j1=this%kblc(kp1)+1
           i2=this%kblc(kp2)+1
           j2=this%kblc(kp2)+1

           call pzgemm('C','N',block1,block2,block1,&
                       (1.0d0,0.0d0),this%H(kp1)%mat,1,1,this%H(kp1)%desc,mat%mat,&
                       i1,j2,mat%desc,(0.0d0,0.0d0),AA%mat,i1,j2,AA%desc)

          enddo
         enddo


         do kp1=1,this%ntot
          do kp2=1,this%ntot

           block1=this%kblc(kp1+1)-this%kblc(kp1)
           block2=this%kblc(kp2+1)-this%kblc(kp2)
           i1=this%kblc(kp1)+1
           j1=this%kblc(kp1)+1
           i2=this%kblc(kp2)+1
           j2=this%kblc(kp2)+1


           call pzgemm('N','N',block1,block2,block2,&
                       (1.0d0,0.0d0),AA%mat,i1,j2,AA%desc,this%H(kp2)%mat,1,1,this%H(kp2)%desc,&
                       (0.0d0,0.0d0),mat%mat,i1,j2,mat%desc) 


          enddo
         enddo


         call AA%dealloc()

        return
        end subroutine to_eigenbasis


        subroutine to_kbasis(this,mat)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)      :: this
        type(dist_cmplx_mat)      :: BB,mat
        integer                   :: k
        
         call BB%set(this%Hdim,this%Hdim,NB,MB)

         call pzgemm('C','N',this%Hdim,this%Hdim,this%Hdim,&
                     (1.0d0,0.0d0),this%kbasis%mat,1,1,this%kbasis%desc,mat%mat,&
                     1,1,mat%desc,(0.0d0,0.0d0),BB%mat,1,1,BB%desc) 

         call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                      (1.0d0,0.0d0),BB%mat,1,1,BB%desc,this%kbasis%mat,1,1,this%kbasis%desc,&
                      (0.0d0,0.0d0),mat%mat,1,1,mat%desc)

         call BB%dealloc()

        return
        end subroutine to_kbasis

        subroutine make_M_H(this,Mr,Mir) 
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)    :: this
        complex(8)              :: M(3),Mi(3,this%s2print),val
        double precision        :: Mr(3),Mir(3,this%s2print)
        integer                 :: i,j,v,s,info,loc_id,kpt,blacs_pnum
        integer                 :: mpi_status(MPI_STATUS_SIZE)

         Mi=(0.0d0,0.0d0)
         M=(0.0d0,0.0d0)

         do i=1,size(this%rho%mat,1)
          do j=1,size(this%rho%mat,2)
           do v=1,this%s2print
            Mi(1,v)=Mi(1,v)+this%rho%mat(i,j)*conjg(this%Sx(v)%mat(i,j))
            Mi(2,v)=Mi(2,v)+this%rho%mat(i,j)*conjg(this%Sy(v)%mat(i,j))
            Mi(3,v)=Mi(3,v)+this%rho%mat(i,j)*conjg(this%Sz(v)%mat(i,j))
           enddo
          enddo
         enddo
                      
         do v=1,this%s2print
          do s=1,3

           if(mpi_id.eq.0)then
            Mir(s,v)=dble(Mi(s,v))
            do i=0,nprow-1
             do j=0,npcol-1
              loc_id=blacs_pnum(context,i,j)
              if(loc_id.ne.0)then
               info=2001
               call mpi_recv(val,1,mpi_double_complex,loc_id,info,mpi_blacs_world,mpi_status,err)
               Mir(s,v)=Mir(s,v)+dble(val)
              endif
             enddo
            enddo
           else
            info=2001
            call mpi_send(Mi(s,v),1,mpi_double_complex,0,info,mpi_blacs_world,err)
           endif

          enddo
         enddo

        return
        end subroutine make_M_H


        subroutine diag_Hmat(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        use lapack_diag_simm
        implicit none
        class(spins_hilbert)            :: this
        integer                         :: i,j,info,nze,k,k1,i1,ii,jj
        integer                         :: t1,t2,rate,size_block
        complex(8)                      :: s_tmp,AA(4,4)
        double precision                :: max_eigen,min_eigen,coeff


         if(mpi_id.eq.0)then
          call system_clock(t1,rate)       
          write(*,*) '' 
          write(*,*) '     Diagonalizing the Spin Hamiltonian Matrix'
          flush(6)
         endif

         this%make_Heig=.true.
         
         allocate(this%Ener(this%ntot)) 
         do k=1,this%ntot
          size_block=this%kblc(k+1)-this%kblc(k)
          allocate(this%Ener(k)%v(size_block)) 
         enddo
        
        ! allocating eigenvectors matrix on the 2d grid
         
         allocate(this%H(this%ntot))
         do k=1,this%ntot
          size_block=this%kblc(k+1)-this%kblc(k)
          call this%H(k)%set(size_block,size_block,NB,MB)
         enddo

         if( this%printH0) call this%dump_H0()                
         
         
        ! diagonalizing hamiltonian with scalapack

         do k=1,this%ntot
          size_block=this%kblc(k+1)-this%kblc(k)
          call pzdiag(size_block,this%H0(k)%mat,this%Ener(k)%v,this%H(k)%mat, &
                      this%H0(k)%desc,this%H(k)%desc)
         enddo


        ! diagonalizing hamiltonian with lapack

!         do k=1,this%ntot
!          size_block=this%kblc(k+1)-this%kblc(k)
!          this%H(k)%mat=this%H0(k)%mat
!          call new_diag(size_block,this%H(k)%mat,this%Ener(k)%v)
!         enddo

        if(this%printHeig) call this%dump_Heig()

        ! find the largest/minimum eigenvalue

         min_eigen=this%Ener(1)%v(1)
         do k=1,this%ntot
          if(min_eigen.gt.this%Ener(k)%v(1))then
           min_eigen=this%Ener(k)%v(1)
          endif
         enddo

         do k=1,this%ntot         
          this%Ener(k)%v=this%Ener(k)%v-min_eigen        
         enddo

         min_eigen=0.0d0
         max_eigen=0.0d0
         do k=1,this%ntot
          if(max_eigen.lt.this%Ener(k)%v(size(this%Ener(k)%v)))then
           max_eigen=this%Ener(k)%v(size(this%Ener(k)%v))
          endif
         enddo         

        ! allocate dos1p and dos2pm

         this%dos1p%sigma=1.0d0
         this%dos1p%step=0.1d0
         this%dos1p%nsteps=nint(max_eigen/this%dos1p%step)+10*nint(this%dos1p%sigma/this%dos1p%step)
         call this%dos1p%alloc_dist()

         do k=1,this%ntot
          do i=1,size(this%ener(k)%v)
           call this%dos1p%update_dist(this%ener(k)%v(i),1.0d0)
          enddo
         enddo

         if(this%make_Rmat .or. this%make_R2mat)then
          this%dos2pm%sigma=5.d0
          this%dos2pm%step=0.1d0
          this%dos2pm%nsteps=nint(max_eigen/this%dos2pm%step)+10*nint(this%dos2pm%sigma/this%dos2pm%step)
          call this%dos2pm%alloc_dist()

          ii=0
          do k=1,this%ntot
           do i=1,size(this%ener(k)%v)
            ii=ii+1
            jj=0
            do k1=1,this%ntot
             do i1=1,size(this%ener(k1)%v)
              jj=jj+1
              if(jj.le.ii)cycle
              call this%dos2pm%update_dist(abs(this%ener(k)%v(i)-this%ener(k1)%v(i1)),1.0d0)
             enddo
            enddo
           enddo
          enddo         
         endif

         if(mpi_id.eq.0)then
          open(11,file='eigenval.dat')
          do k=1,this%ntot
           do i=1,size(this%ener(k)%v)
            write(11,*) k,i,this%Ener(k)%v(i)
           enddo
          enddo
          close(11)
          open(11,file='spin_dos1p.dat')
          do k=1,this%dos1p%nsteps
           write(11,*) this%dos1p%step*(k-1),this%dos1p%dist(k)
          enddo
          close(11)
!          if(this%make_Rmat)then
!           open(11,file='spin_dos2pm.dat')
!           do k=1,this%dos2pm%nsteps
!            write(11,*) this%dos2pm%step*(k-1),this%dos2pm%dist(k)
!           enddo
!           close(11)
!          endif
         endif

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
         endif

        return
        end subroutine diag_Hmat


        subroutine dump_S_H(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use scalapack_diag_simm
        implicit none
        class(spins_hilbert)            :: this
        integer                         :: i,j,k,info
        integer                         :: t1,t2,rate
        complex(8), allocatable         :: s_tmp(:,:)

         allocate(s_tmp(3,size(this%print_si)))
         
         if(mpi_id.eq.0) open(11,file='Sx_diag.dat')  
         if(mpi_id.eq.0) open(12,file='Sy_diag.dat')  
         if(mpi_id.eq.0) open(13,file='Sz_diag.dat')  
         do i=1,this%Hdim          
          do j=1,this%s2print
           call pzelget('A',' ',s_tmp(3,j),this%Sz(j)%mat,i,i,this%Sz(j)%desc)
           call pzelget('A',' ',s_tmp(2,j),this%Sy(j)%mat,i,i,this%Sy(j)%desc)
           call pzelget('A',' ',s_tmp(1,j),this%Sx(j)%mat,i,i,this%Sx(j)%desc)
          enddo
          if(mpi_id.eq.0)  write(11,*) i,dble(s_tmp(1,:))
          if(mpi_id.eq.0)  write(12,*) i,dble(s_tmp(2,:))
          if(mpi_id.eq.0)  write(13,*) i,dble(s_tmp(3,:))
         enddo
         if(mpi_id.eq.0) close(11)
         if(mpi_id.eq.0) close(12)
         if(mpi_id.eq.0) close(13)

         deallocate(s_tmp)

        return
        end subroutine dump_S_H


        subroutine make_S_H(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)                :: this
        type(list)                          :: AJ,Aval
        type(csr_mat_cmplx), allocatable    :: Sz(:),Sy(:),Sx(:)
        logical,allocatable                 :: a(:),skip
        integer                             :: i,j,v,k,s,t,ii,jj,kk,l
        integer                             :: indxl2g,kpt,size_block
        integer                             :: t1,t2,rate,t1_2
        double precision                    :: aaa
        complex(8), allocatable             :: Mtmp(:)
        complex(8)                          :: sum,val
         
         if(allocated(this%Sx)) return         

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Calculation of the Spin operators matrix elements'
          flush(6)
         endif

         if(this%s2print.gt.0)then
       
          allocate(this%Sx(this%s2print))
          allocate(this%Sy(this%s2print))
          allocate(this%Sz(this%s2print))

          !!!!
          allocate(this%sp_Sx(this%s2print))
          allocate(this%sp_Sy(this%s2print))
          allocate(this%sp_Sz(this%s2print))
          allocate(Sz(this%s2print))
          allocate(Sy(this%s2print))
          allocate(Sx(this%s2print))
          !!!!

          do i=1,this%s2print
           call this%Sx(i)%set(this%Hdim,this%Hdim,NB,MB)
           call this%Sy(i)%set(this%Hdim,this%Hdim,NB,MB)
           call this%Sz(i)%set(this%Hdim,this%Hdim,NB,MB)
           this%Sx(i)%mat=(0.0d0,0.0d0)
           this%Sy(i)%mat=(0.0d0,0.0d0)
           this%Sz(i)%mat=(0.0d0,0.0d0)
          enddo

         endif
                 
         allocate(a(this%nspins))

         allocate(Mtmp(3))
         Mtmp=(0.0d0,0.0d0)

         do ii=1,size(this%Sz(1)%mat,1)
          do jj=1,size(this%Sz(1)%mat,2)
           s=indxl2g(ii,NB,myrow,0,nprow)
           t=indxl2g(jj,MB,mycol,0,npcol)
           if (s.eq.t) then                   
            do v=1,this%nspins

             do k=1,this%s2print
              if(this%print_si(k).eq.v)then
               this%Sz(k)%mat(ii,jj)=this%basis(s,v)
              endif
              if(this%print_si(k).eq.-1)then
               this%Sz(k)%mat(ii,jj)=this%Sz(k)%mat(ii,jj)+this%basis(s,v)
              endif
             enddo
 
            enddo
           endif
          enddo
         enddo
          
         do ii=1,size(this%Sx(1)%mat,1)
          do jj=1,size(this%Sx(1)%mat,2)

           s=indxl2g(ii,NB,myrow,0,nprow)
           t=indxl2g(jj,MB,mycol,0,npcol)

           do v=1,this%nspins

            skip=.true.
            do k=1,this%s2print
             if(this%print_si(k).eq.v) skip=.false.
             if(this%print_si(k).eq.-1) skip=.false.
            enddo

            if(skip) cycle

            a=.true.
            do k=1,this%nspins
             if(k.ne.v)then
              if(ABS(this%basis(s,k)-this%basis(t,k)).gt.1.0d-8)then
               a(k)=.false.
               exit
              endif
             endif
            enddo

            if(ALL(a))then

! S+
             if(abs(this%basis(s,v)-this%basis(t,v)-1).lt.1.0E-06)then

              aaa=dsqrt((this%spin(this%kind(v))-this%basis(t,v))* &
                        (this%spin(this%kind(v))+this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)
              Mtmp(2)=cmplx(0.0d0,-1.0d0*aaa,8)

              do k=1,this%s2print
               if(this%print_si(k).eq.v)then
                this%Sx(k)%mat(ii,jj)=Mtmp(1)
                this%Sy(k)%mat(ii,jj)=Mtmp(2)
               endif
               if(this%print_si(k).eq.-1)then
                this%Sx(k)%mat(ii,jj)=this%Sx(k)%mat(ii,jj)+Mtmp(1)
                this%Sy(k)%mat(ii,jj)=this%Sy(k)%mat(ii,jj)+Mtmp(2)
               endif
              enddo

             endif

! S-
             if(abs(this%basis(s,v)-this%basis(t,v)+1).lt.1.0E-06)then

              aaa=dsqrt((this%spin(this%kind(v))+this%basis(t,v))* &
                        (this%spin(this%kind(v))-this%basis(t,v)+1) )
              aaa=aaa/2.0d0

              Mtmp(1)=cmplx(aaa,0.0d0,8)
              Mtmp(2)=cmplx(0.0d0,aaa,8)

              do k=1,this%s2print
               if(this%print_si(k).eq.v)then
                this%Sx(k)%mat(ii,jj)=Mtmp(1)
                this%Sy(k)%mat(ii,jj)=Mtmp(2)
               endif
               if(this%print_si(k).eq.-1)then
                this%Sx(k)%mat(ii,jj)=this%Sx(k)%mat(ii,jj)+Mtmp(1)
                this%Sy(k)%mat(ii,jj)=this%Sy(k)%mat(ii,jj)+Mtmp(2)
               endif
              enddo

             endif

            endif

           enddo
          enddo
         enddo

         if(this%sparse)then

          do s=1,this%s2print

           call this%sp_Sz(s)%tosparse(this%Hdim,this%Sz(s))
           call this%sp_Sy(s)%tosparse(this%Hdim,this%Sy(s))
           call this%sp_Sx(s)%tosparse(this%Hdim,this%Sx(s))

           if(this%ntot.gt.1)then
            call mult2_sparse_cmplx(.false.,this%sp_Sz(s),this%sp_kbasis,Sz(s))          
            call this%sp_Sz(s)%delete()
            call mult2_sparse_cmplx(.true.,this%sp_kbasis,Sz(s),this%sp_Sz(s))
            call Sz(s)%delete()

            call mult2_sparse_cmplx(.false.,this%sp_Sy(s),this%sp_kbasis,Sy(s))
            call this%sp_Sy(s)%delete()
            call mult2_sparse_cmplx(.true.,this%sp_kbasis,Sy(s),this%sp_Sy(s))
            call Sy(s)%delete()

            call mult2_sparse_cmplx(.false.,this%sp_Sx(s),this%sp_kbasis,Sx(s))
            call this%sp_Sx(s)%delete()
            call mult2_sparse_cmplx(.true.,this%sp_kbasis,Sx(s),this%sp_Sx(s))
            call Sx(s)%delete()
           endif

          call this%Sz(s)%dealloc()
          call this%Sy(s)%dealloc()
          call this%Sx(s)%dealloc()

          call this%sp_Sz(s)%todense(this%Sz(s)) 
          call this%sp_Sy(s)%todense(this%Sy(s)) 
          call this%sp_Sx(s)%todense(this%Sx(s)) 

         enddo ! print
        
         else

          if(this%ntot.gt.1)then

           do s=1,this%s2print
            call this%to_kbasis(this%Sx(s))
            call this%to_kbasis(this%Sy(s))
            call this%to_kbasis(this%Sz(s))             
           enddo

          endif

         endif  ! sparse

         if(this%make_Heig)then

          do s=1,this%s2print
           call this%to_eigenbasis(this%Sx(s))
           call this%to_eigenbasis(this%Sy(s))
           call this%to_eigenbasis(this%Sz(s))
          enddo
        
         endif

         deallocate(Mtmp)
         deallocate(a)

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return 
        end subroutine make_S_H


        subroutine check_couple(nspins,mapper,mapper2,a,t,v)
        implicit none
        integer :: k,t,v
        double precision :: mapper(nspins)
        double precision :: mapper2(nspins)
        double precision :: ZERO=1.0E-8
        logical          :: a(nspins)
        integer          :: nspins

         a=.true.

         do k=1,nspins
          if(k.ne.t .and. k.ne.v)then
           if(ABS(mapper(k)-mapper2(k)).gt.ZERO)then
            a(k)=.false.
            exit
           endif
          endif
         enddo

        return
        end subroutine check_couple


        subroutine make_Hmat_nodes(this)
        use lists_class   
        use blacs_utils
        implicit none
        class(spins_hilbert)       :: this
        integer                    :: l,l2,ii,jj,s1,s2,indxl2g,i
        integer                    :: Nloc_row,Nloc_col,numroc
        integer                    :: spin(2),ndiff
        integer                    :: t1,t2,rate
        logical                    :: save_it
        type(list)                 :: nodes(4)


         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Calculation of non-zero elements of H0'
          flush(6)
         endif
         
         call nodes(1)%init()
         call nodes(2)%init()
         call nodes(3)%init()
         call nodes(4)%init()

!         do i=1,size(context)
!          call blacs_gridinfo(context(i),nprow,npcol,myrow,mycol)
!          if(myrow.ne.-1)then

           Nloc_row = NUMROC(this%Hdim,NB,myrow,0,nprow)
           Nloc_col = NUMROC(this%Hdim,MB,mycol,0,npcol)

           do ii=1,Nloc_row
            do jj=1,Nloc_col

             l=indxl2g(ii,NB,myrow,0,nprow)
             l2=indxl2g(jj,MB,mycol,0,npcol)
         
             ndiff=0
             save_it=.true.
             spin=0

             do s1=1,this%nspins
              if(ABS(this%basis(l,s1)-this%basis(l2,s1)).gt.1.0d-8)then
               ndiff=ndiff+1
               if(ndiff.gt.2)then
                save_it=.false.
                exit
               endif
               spin(ndiff)=s1
              endif
             enddo

             if(save_it)then               
              
              s1=min(spin(1),spin(2))
              s2=max(spin(1),spin(2))

              if(s1.gt.this%nspins_pr) goto 12

              call nodes(1)%add_node(ii)
              call nodes(2)%add_node(jj)
              call nodes(3)%add_node(s1)
              call nodes(4)%add_node(s2)

             endif

12           continue           

            enddo
           enddo

!          endif
!         enddo

         allocate(this%Hnodes(nodes(1)%nelem,4))

         call nodes(1)%reboot()
         call nodes(2)%reboot()
         call nodes(3)%reboot()
         call nodes(4)%reboot()

         do ii=1,nodes(1)%nelem
          call nodes(1)%rd_val(this%Hnodes(ii,1))
          call nodes(2)%rd_val(this%Hnodes(ii,2))
          call nodes(3)%rd_val(this%Hnodes(ii,3))
          call nodes(4)%rd_val(this%Hnodes(ii,4))
          call nodes(1)%skip()
          call nodes(2)%skip()
          call nodes(3)%skip()
          call nodes(4)%skip()
         enddo 

         call nodes(1)%delete()
         call nodes(2)%delete()
         call nodes(3)%delete()
         call nodes(4)%delete()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_Hmat_nodes



        function get_Hij_2(this,node) result(val)
        implicit none
        class(spins_hilbert)       :: this
        integer                    :: s1,s2,l,l2,i,j
        integer                    :: ii,jj,indxl2g
        complex(8)                 :: val
        integer                    :: node(4)
        
          
!         do cntx=1,size(context)
!         call blacs_gridinfo(context(cntx),nprow,npcol,myrow,mycol)
!          if(myrow.ne.-1)then
           l=indxl2g(node(1),NB,myrow,0,nprow)
           l2=indxl2g(node(2),MB,mycol,0,npcol)
!          endif
!         enddo

         val=(0.0d0,0.0d0)

         if(node(3).ne.0 .and. node(4).ne.0)then

          s1=node(3)
          s2=node(4)

          i=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
          j=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
          ii=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
          i=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
          j=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
          jj=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j

          val=val+this%SHrep(s1,s2)%mat(ii,jj)

         endif

         if(node(3).eq.0 .and. node(4).ne.0)then

          s2=node(4)
          if(s2.gt.this%nspins_pr)then
           do s1=1,this%nspins_pr
            if(s1.lt.s2)then
             i=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
            if(s1.gt.s2)then
             i=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
             j=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             ii=(i-1)*(nint(2*this%spin(this%kind(s1)))+1)+j
             i=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             j=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
             jj=(i-1)*(nint(2*this%spin(this%kind(s1)))+1)+j
             val=val+this%SHrep(s2,s1)%mat(ii,jj)
            endif
            if(s1.eq.s2)then
             ii=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
           enddo
          else
           do s1=1,this%nspins
            if(s1.lt.s2)then
             i=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
            if(s1.gt.s2)then
             i=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
             j=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             ii=(i-1)*(nint(2*this%spin(this%kind(s1)))+1)+j
             i=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             j=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
             jj=(i-1)*(nint(2*this%spin(this%kind(s1)))+1)+j
             val=val+this%SHrep(s2,s1)%mat(ii,jj)
            endif
            if(s1.eq.s2)then
             ii=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             val=val+this%SHrep(s1,s2)%mat(ii,jj)
            endif
           enddo
          endif

         endif
        
         
         if(node(3).eq.0 .and. node(4).eq.0)then

          do s1=1,this%nspins_pr
           do s2=s1,this%nspins
            if(s1.ne.s2)then
             i=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l,s2)+this%spin(this%kind(s2))+1)
             ii=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
             i=nint(this%basis(l2,s1)+this%spin(this%kind(s1))+1)
             j=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
             jj=(i-1)*(nint(2*this%spin(this%kind(s2)))+1)+j
            else
             ii=nint(this%basis(l,s1)+this%spin(this%kind(s1))+1)
             jj=nint(this%basis(l2,s2)+this%spin(this%kind(s2))+1)
            endif
            val=val+this%SHrep(s1,s2)%mat(ii,jj)
           enddo
          enddo

         endif

        return
        end function get_Hij_2


        subroutine make_SH_rep(this,SH,spin_id,spin_id2)
        use spinham_class
        implicit none
        class(spins_hilbert)       :: this
        class(SpinHamiltonian)     :: SH
        integer                    :: i,j,l2,l,s1,s2,l1
        integer                    :: ii,jj,ii1,ii2,jj1,jj2
        integer                    :: t1,t2,rate,a,b,c,d,spin_id,spin_id2
        double precision           :: spin2(2),psi2(2,2)
        double precision           :: spin1,psi1(2)

!         if(mpi_id.eq.0)then
!          call system_clock(t1,rate)        
!          write(*,*) '   Calculation of SH representaion'
!          flush(6)
!         endif         

         if(.not.allocated(this%SHrep))then
          allocate(this%SHrep(this%nspins_pr,this%nspins))
         endif

         do s1=1,this%nspins_pr  
          if(allocated(this%SHrep(s1,s1)%mat)) deallocate(this%SHrep(s1,s1)%mat)
          ii=nint(2*this%spin(this%kind(s1)))+1
          allocate(this%SHrep(s1,s1)%mat(ii,ii))
          this%SHrep(s1,s1)%mat=(0.0d0,0.0d0)
         enddo
         do s1=1,this%nspins_pr
          do s2=s1+1,this%nspins
           if(allocated(this%SHrep(s1,s2)%mat)) deallocate(this%SHrep(s1,s2)%mat)
           ii=nint(2*this%spin(this%kind(s1)))+1
           ii=ii*(nint(2*this%spin(this%kind(s2)))+1)
           allocate(this%SHrep(s1,s2)%mat(ii,ii))
           this%SHrep(s1,s2)%mat=(0.0d0,0.0d0)
          enddo
         enddo

         do s1=1,this%nspins_pr  
          if(spin_id.ne.spin_id2) cycle
          if(spin_id.ne.-1 .and. spin_id.ne.s1) cycle

          if(allocated(this%SHrep(s1,s1)%mat)) deallocate(this%SHrep(s1,s1)%mat)
          ii=nint(2*this%spin(this%kind(s1)))+1
          allocate(this%SHrep(s1,s1)%mat(ii,ii))
          this%SHrep(s1,s1)%mat=(0.0d0,0.0d0)

          do ii=1,nint(2*this%spin(this%kind(s1)))+1
           do jj=1,nint(2*this%spin(this%kind(s1)))+1

            psi1(1)=ii-this%spin(this%kind(s1))-1
            psi1(2)=jj-this%spin(this%kind(s1))-1
            spin1=this%spin(this%kind(s1))

        !  Stevens operators matrix elements
                   
            do i=1,SH%nO
             if(this%kind(s1).eq.SH%O(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%O(i)%mat_elem(psi1,spin1)
             endif
            enddo

        !  Zeeman Stevens operators matrix elements

            do i=1,SH%nG
             if(this%kind(s1).eq.SH%G(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%G(i)%mat_elem(psi1,spin1,this%Bfield,this%bohr_mag(this%kind(s1)))
             endif
            enddo

        !  Single Ion Anisotropy operators matrix elements

            do i=1,SH%nDSI
             if(this%kind(s1).eq.SH%DSI(i)%kind)then
              this%SHrep(s1,s1)%mat(ii,jj)=this%SHrep(s1,s1)%mat(ii,jj)&
                +SH%DSI(i)%mat_elem(psi1,spin1)
             endif
            enddo


           enddo ! 1 spin basis
          enddo ! 1 spin basis

         enddo ! spins

         do s1=1,this%nspins_pr
          if(spin_id.ne.-1 .and. spin_id.ne.s1) cycle
          do s2=s1+1,this%nspins
           if(spin_id2.ne.-1 .and. spin_id2.ne.s2) cycle

           if(allocated(this%SHrep(s1,s2)%mat)) deallocate(this%SHrep(s1,s2)%mat)
           ii=nint(2*this%spin(this%kind(s1)))+1
           ii=ii*(nint(2*this%spin(this%kind(s2)))+1)
           allocate(this%SHrep(s1,s2)%mat(ii,ii))
           this%SHrep(s1,s2)%mat=(0.0d0,0.0d0)

           ii=1
           do ii1=1,nint(2*this%spin(this%kind(s1)))+1
            do ii2=1,nint(2*this%spin(this%kind(s2)))+1

             jj=1
             do jj1=1,nint(2*this%spin(this%kind(s1)))+1
              do jj2=1,nint(2*this%spin(this%kind(s2)))+1

               psi2(1,1)=ii1-this%spin(this%kind(s1))-1
               psi2(1,2)=ii2-this%spin(this%kind(s2))-1
               psi2(2,1)=jj1-this%spin(this%kind(s1))-1
               psi2(2,2)=jj2-this%spin(this%kind(s2))-1
               spin2(1)=this%spin(this%kind(s1))
               spin2(2)=this%spin(this%kind(s2))

        !  Isotropic exchange operators matrix elements

               do i=1,SH%nJ
                if(this%kind(s1).eq.SH%J(i)%kind(1) .and. this%kind(s2).eq.SH%J(i)%kind(2) )then
                 if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)
                 endif
                endif
                if(this%kind(s1).eq.SH%J(i)%kind(2) .and. this%kind(s2).eq.SH%J(i)%kind(1) .and. & 
                   SH%J(i)%kind(1).ne.SH%J(i)%kind(2) ) then
                 if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%J(i)%mat_elem(psi2,spin2)
                 endif
                endif
               enddo

        !  Exchange Anisotry operators matrix elements

               do i=1,SH%nD2S

                b=MOD((s2-1),this%nspins_pr)+1
                d=INT((s2-1)/this%nspins_pr)+1
!                if(this%dist(s1,1,b,d).le.SH%D2S(i)%cutoff)then
                if(this%dist(s1,1,b,d).le.0.5d0)then


                if(this%kind(s1).eq.SH%D2S(i)%kind(1) .and. this%kind(s2).eq.SH%D2S(i)%kind(2)) then
                 if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)
                 endif
                endif
                if(this%kind(s1).eq.SH%D2S(i)%kind(2) .and. this%kind(s2).eq.SH%D2S(i)%kind(1) .and. & 
                   SH%D2S(i)%kind(1).ne.SH%D2S(i)%kind(2) ) then
                 if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%D2S(i)%mat_elem(psi2,spin2)
                 endif
                endif              

                endif
               enddo


        !  Dipolar Exchange operators matrix elements

               do i=1,SH%nDdip
                if( s1.eq.SH%Ddip(i)%kind(1) .and. s2.eq.SH%Ddip(i)%kind(2) ) then
                 if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%Ddip(i)%mat_elem(psi2,spin2)*0.5d0
                 else
                  this%SHrep(s1,s2)%mat(ii,jj)=this%SHrep(s1,s2)%mat(ii,jj)&
                        +SH%Ddip(i)%mat_elem(psi2,spin2)
                 endif
                endif
      
               enddo

               jj=jj+1
              enddo
             enddo

             ii=ii+1
            enddo
           enddo

          enddo ! s2           
         enddo ! s1

        return
        end subroutine make_SH_rep


        function get_Hij(this,SH,l,l2) result(val)
        use spinham_class
        implicit none
        class(spins_hilbert)       :: this
        class(SpinHamiltonian)     :: SH
        integer                    :: i,j,l2,l,s1,s2,cell1,cell2,ii,jj
        double precision           :: spin2(2),psi2(2,2),sum
        double precision           :: spin1,psi1(2)
        logical, allocatable       :: a(:)
        complex(8)                 :: val
         
         allocate(a(this%nspins))

         val=(0.0d0,0.0d0)

         do s1=1,this%nspins_pr
          
!          call check_couple(this%nspins,this%basis(l,:),this%basis(l2,:),a,s1,s1)
          a=.true.
          do ii=1,this%nspins
           if(ii.ne.s1 )then
            if(ABS(this%basis(l,ii)-this%basis(l2,ii)).gt.1.0d-8)then
             a(ii)=.false.
             exit
            endif
           endif
          enddo

           if(all(a))then

        !  Stevens operators matrix elements

            do i=1,SH%nO
             if(this%kind(s1).eq.SH%O(i)%kind)then
              spin1=this%spin(SH%O(i)%kind)
              psi1(1)=this%basis(l,s1)
              psi1(2)=this%basis(l2,s1)
              val=val+SH%O(i)%mat_elem(psi1,spin1)
             endif
            enddo

        !  Zeeman Stevens operators matrix elements

            do i=1,SH%nG
             if(this%kind(s1).eq.SH%G(i)%kind)then
              spin1=this%spin(SH%G(i)%kind)
              psi1(1)=this%basis(l,s1)
              psi1(2)=this%basis(l2,s1)
              val=val+SH%G(i)%mat_elem(psi1,spin1,this%Bfield,this%bohr_mag(this%kind(s1)))
             endif
            enddo

        !  Single Ion Anisotropy operators matrix elements

            do i=1,SH%nDSI
             if(this%kind(s1).eq.SH%DSI(i)%kind)then
              spin1=this%spin(SH%DSI(i)%kind)
              psi1(1)=this%basis(l,s1)
              psi1(2)=this%basis(l2,s1)
              val=val+SH%DSI(i)%mat_elem(psi1,spin1)
             endif
            enddo

           endif

          enddo

          do s1=1,this%nspins_pr
           do s2=s1+1,this%nspins

!            call check_couple(this%nspins,this%basis(l,:),this%basis(l2,:),a,s1,s2)

            a=.true.
            do ii=1,this%nspins
             if(ii.ne.s1 .and. ii.ne.s2)then
              if(ABS(this%basis(l,ii)-this%basis(l2,ii)).gt.1.0d-8)then
               a(ii)=.false.
               exit
              endif
             endif
            enddo

            if(all(a))then

        !  Isotropic exchange operators matrix elements

             do i=1,SH%nJ
              if(this%kind(s1).eq.SH%J(i)%kind(1) .and. this%kind(s2).eq.SH%J(i)%kind(2) )then
               spin2(1)=this%spin(SH%J(i)%kind(1))
               spin2(2)=this%spin(SH%J(i)%kind(2))
               psi2(1,1)=this%basis(l,s1)
               psi2(1,2)=this%basis(l,s2)
               psi2(2,1)=this%basis(l2,s1)
               psi2(2,2)=this%basis(l2,s2)
               if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                val=val+SH%J(i)%mat_elem(psi2,spin2)*0.5d0
               else
                val=val+SH%J(i)%mat_elem(psi2,spin2)
               endif
              endif
              if(this%kind(s1).eq.SH%J(i)%kind(2) .and. this%kind(s2).eq.SH%J(i)%kind(1) .and. & 
                 SH%J(i)%kind(1).ne.SH%J(i)%kind(2) ) then
               spin2(1)=this%spin(SH%J(i)%kind(2))
               spin2(2)=this%spin(SH%J(i)%kind(1))
               psi2(1,1)=this%basis(l,s1)
               psi2(1,2)=this%basis(l,s2)
               psi2(2,1)=this%basis(l2,s1)
               psi2(2,2)=this%basis(l2,s2)
               if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                val=val+SH%J(i)%mat_elem(psi2,spin2)*0.5d0
               else
                val=val+SH%J(i)%mat_elem(psi2,spin2)
               endif
              endif
             enddo

        !  Exchange Anisotry operators matrix elements

             do i=1,SH%nD2S
              if(this%kind(s1).eq.SH%D2S(i)%kind(1) .and. this%kind(s2).eq.SH%D2S(i)%kind(2)) then
               spin2(1)=this%spin(SH%D2S(i)%kind(1))
               spin2(2)=this%spin(SH%D2S(i)%kind(2))
               psi2(1,1)=this%basis(l,s1)
               psi2(1,2)=this%basis(l,s2)
               psi2(2,1)=this%basis(l2,s1)
               psi2(2,2)=this%basis(l2,s2)
               if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                val=val+SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
               else
                val=val+SH%D2S(i)%mat_elem(psi2,spin2)
               endif
              endif
              if(this%kind(s1).eq.SH%D2S(i)%kind(2) .and. this%kind(s2).eq.SH%D2S(i)%kind(1) .and. & 
                 SH%D2S(i)%kind(1).ne.SH%D2S(i)%kind(2) ) then
               spin2(1)=this%spin(SH%D2S(i)%kind(2))
               spin2(2)=this%spin(SH%D2S(i)%kind(1))
               psi2(1,1)=this%basis(l,s1)
               psi2(1,2)=this%basis(l,s2)
               psi2(2,1)=this%basis(l2,s1)
               psi2(2,2)=this%basis(l2,s2)
               if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                val=val+SH%D2S(i)%mat_elem(psi2,spin2)*0.5d0
               else
                val=val+SH%D2S(i)%mat_elem(psi2,spin2)
               endif
              endif              
             enddo


        !  Dipolar Exchange operators matrix elements

              do i=1,SH%nDdip
               
               if( s1.eq.SH%Ddip(i)%kind(1) .and. s2.eq.SH%Ddip(i)%kind(2) ) then
                spin2(1)=this%spin(this%kind(s1))
                spin2(2)=this%spin(this%kind(s2))
                psi2(1,1)=this%basis(l,s1)
                psi2(1,2)=this%basis(l,s2)
                psi2(2,1)=this%basis(l2,s1)
                psi2(2,2)=this%basis(l2,s2)     
                if(this%ntot.gt.1 .and. s2.gt.this%nspins_pr)then
                 val=val+SH%Ddip(i)%mat_elem(psi2,spin2)*0.5d0
                else
                 val=val+SH%Ddip(i)%mat_elem(psi2,spin2)
                endif
               endif
       
              enddo

            endif ! if on mapp
           enddo ! s2           
          enddo ! s1

         deallocate(a)

        return
        end function get_Hij

        subroutine set_sph2_dipolar(this,SH,SPH,ex_list)
        use spinham_class
        use mpi_utils
        implicit none
        class(spins_hilbert)             :: this
        class(SpinHamiltonian)           :: SH
        class(SpinPhononHamiltonian)     :: SPH
        integer                          :: i,j,s1,s2,ii,jj,v,l,m,ex,is
        integer                          :: ii_1,jj_1,vv
        integer, allocatable             :: ex_list(:,:)
        integer                          :: celli,cellj,Hdim1,Hdim2,Hdim12
        double precision                 :: dist0(3)
        double precision                 :: spin2(2),spin,psi(2),psi2(2,2)
        integer                          :: t1,t2,rate
        logical                          :: skip

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building 2nd-order Spin-Phonon Dipolar Network'
          flush(6)
         endif

         SPH%nDdip=0

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif

             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SPH%dipolar_thr)then

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then                 
               SPH%nDdip=SPH%nDdip+1
              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then
               SPH%nDdip=SPH%nDdip+1
              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         allocate(SPH%Ddip(SPH%nDdip))
         allocate(SPH%Ddip_t(SPH%nDdip))

         if(mpi_id.eq.0)  &
         write(*,*) '     Total Number of Spin-Phonon '&
                                 'Dipolar interactions: ',SH%nDdip

         v=1

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif

             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SPH%dipolar_thr)then

              dist0=this%dist_vec_pbc(this%x(s1,:),this%x(s2,:))  

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then

                call SPH%Ddip_t(v)%make_dD(SH%G(i)%G,SH%G(j)%G, &
                      this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
                SPH%Ddip_t(v)%nderiv=36
                allocate(SPH%Ddip_t(v)%map_s2a(36,4))
                vv=0
                do ii_1=1,6
                 do jj_1=1,6
                  vv=vv+1
                  if(ii_1.le.3)then
                   SPH%Ddip_t(v)%map_s2a(vv,2)=celli
                   SPH%Ddip_t(v)%map_s2a(vv,1)=(SPH%mapp(ii)-1)*3+mod(ii_1-1,3)+1
                  else
                   SPH%Ddip_t(v)%map_s2a(vv,2)=cellj
                   SPH%Ddip_t(v)%map_s2a(vv,1)=(SPH%mapp(jj)-1)*3+mod(ii_1-1,3)+1
                  endif
                  if(jj_1.le.3)then
                   SPH%Ddip_t(v)%map_s2a(vv,4)=celli
                   SPH%Ddip_t(v)%map_s2a(vv,3)=(SPH%mapp(ii)-1)*3+mod(jj_1-1,3)+1
                  else
                   SPH%Ddip_t(v)%map_s2a(vv,4)=cellj
                   SPH%Ddip_t(v)%map_s2a(vv,3)=(SPH%mapp(jj)-1)*3+mod(jj_1-1,3)+1
                  endif
                 enddo
                enddo
                SPH%Ddip_t(v)%kind(1)=s1
                SPH%Ddip_t(v)%kind(2)=s2
                 
                v=v+1

              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then

                call SPH%Ddip_t(v)%make_ddD(SH%G(j)%G,SH%G(i)%G,  &
                     this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
                SPH%Ddip_t(v)%nderiv=36
                allocate(SPH%Ddip_t(v)%map_s2a(36,4))
                vv=0
                do ii_1=1,6
                 do jj_1=1,6
                  vv=vv+1
                  if(ii_1.le.3)then
                   SPH%Ddip_t(v)%map_s2a(vv,2)=celli
                   SPH%Ddip_t(v)%map_s2a(vv,1)=(SPH%mapp(ii)-1)*3+mod(ii_1-1,3)+1
                  else
                   SPH%Ddip_t(v)%map_s2a(vv,2)=cellj
                   SPH%Ddip_t(v)%map_s2a(vv,1)=(SPH%mapp(jj)-1)*3+mod(ii_1-1,3)+1
                  endif
                  if(jj_1.le.3)then
                   SPH%Ddip_t(v)%map_s2a(vv,4)=celli
                   SPH%Ddip_t(v)%map_s2a(vv,3)=(SPH%mapp(ii)-1)*3+mod(jj_1-1,3)+1
                  else
                   SPH%Ddip_t(v)%map_s2a(vv,4)=cellj
                   SPH%Ddip_t(v)%map_s2a(vv,3)=(SPH%mapp(jj)-1)*3+mod(jj_1-1,3)+1
                  endif
                 enddo
                enddo
                SPH%Ddip_t(v)%kind(1)=s1
                SPH%Ddip_t(v)%kind(2)=s2

                v=v+1

              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine set_sph2_dipolar

        subroutine set_sph_dipolar(this,SH,SPH,ex_list)
        use spinham_class
        use mpi_utils
        implicit none
        class(spins_hilbert)             :: this
        class(SpinHamiltonian)           :: SH
        class(SpinPhononHamiltonian)     :: SPH
        integer                          :: i,j,s1,s2,ii,jj,v,l,m,ex,is
        integer, allocatable             :: ex_list(:,:)
        integer                          :: celli,cellj,Hdim1,Hdim2,Hdim12
        double precision                 :: dist0(3)
        double precision                 :: spin2(2),spin,psi(2),psi2(2,2)
        integer                          :: t1,t2,rate
        logical                          :: skip

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building Spin-Phonon Dipolar Network'
          flush(6)
         endif

         SPH%nDdip=0

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif

             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SPH%dipolar_thr)then

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then                 
               SPH%nDdip=SPH%nDdip+1
              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then
               SPH%nDdip=SPH%nDdip+1
              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         allocate(SPH%Ddip(SPH%nDdip))
         allocate(SPH%Ddip_t(SPH%nDdip))

         if(mpi_id.eq.0)  &
         write(*,*) '     Total Number of Spin-Phonon '&
                                 'Dipolar interactions: ',SH%nDdip

         v=1

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif


             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SPH%dipolar_thr)then

              dist0=this%dist_vec_pbc(this%x(s1,:),this%x(s2,:))  

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then

                call SPH%Ddip_t(v)%make_dD(SH%G(i)%G,SH%G(j)%G, &
                      this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
                SPH%Ddip_t(v)%nderiv=6
                allocate(SPH%Ddip_t(v)%map_s2a(6,2))
                SPH%Ddip_t(v)%map_s2a(1:3,2)=celli
                SPH%Ddip_t(v)%map_s2a(4:6,2)=cellj
                SPH%Ddip_t(v)%map_s2a(1,1)=(SPH%mapp(ii)-1)*3+1
                SPH%Ddip_t(v)%map_s2a(2,1)=(SPH%mapp(ii)-1)*3+2
                SPH%Ddip_t(v)%map_s2a(3,1)=(SPH%mapp(ii)-1)*3+3
                SPH%Ddip_t(v)%map_s2a(4,1)=(SPH%mapp(jj)-1)*3+1
                SPH%Ddip_t(v)%map_s2a(5,1)=(SPH%mapp(jj)-1)*3+2
                SPH%Ddip_t(v)%map_s2a(6,1)=(SPH%mapp(jj)-1)*3+3
                SPH%Ddip_t(v)%kind(1)=s1
                SPH%Ddip_t(v)%kind(2)=s2
                 
                v=v+1

              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then

                call SPH%Ddip_t(v)%make_dD(SH%G(j)%G,SH%G(i)%G,  &
                     this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
                SPH%Ddip_t(v)%nderiv=6
                allocate(SPH%Ddip_t(v)%map_s2a(6,2))
                SPH%Ddip_t(v)%map_s2a(1:3,2)=celli
                SPH%Ddip_t(v)%map_s2a(4:6,2)=cellj
                SPH%Ddip_t(v)%map_s2a(1,1)=(SPH%mapp(ii)-1)*3+1
                SPH%Ddip_t(v)%map_s2a(2,1)=(SPH%mapp(ii)-1)*3+2
                SPH%Ddip_t(v)%map_s2a(3,1)=(SPH%mapp(ii)-1)*3+3
                SPH%Ddip_t(v)%map_s2a(4,1)=(SPH%mapp(jj)-1)*3+1
                SPH%Ddip_t(v)%map_s2a(5,1)=(SPH%mapp(jj)-1)*3+2
                SPH%Ddip_t(v)%map_s2a(6,1)=(SPH%mapp(jj)-1)*3+3
                SPH%Ddip_t(v)%kind(1)=s1
                SPH%Ddip_t(v)%kind(2)=s2

                v=v+1

              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo


         ! compute matrix representation of SH%Ddips

!         do i=1,SH%nDdip
                
!          s1=SH%Ddip(i)%kind(1)
!         s2=SH%Ddip(i)%kind(2)
                
!          spin2(1)=this%spin(this%kind(s1))
!          spin2(2)=this%spin(this%kind(s2))

!          Hdim1=(2*spin2(1)+1)
!          Hdim2=(2*spin2(2)+1)
!          Hdim12=Hdim1*Hdim2
!          allocate(SH%Ddip(i)%mat(Hdim12,Hdim12))
!          SH%Ddip(i)%mat=(0.0d0,0.0d0)

!          do l=1,Hdim12
!           do m=1,Hdim12
!            psi2(1,1)=-this%spin(this%kind(s1))+int((l-1)/Hdim2)
!            psi2(1,2)=-this%spin(this%kind(s2))+mod((l-1),Hdim2)
!            psi2(2,1)=-this%spin(this%kind(s1))+int((m-1)/Hdim2)
!            psi2(2,2)=-this%spin(this%kind(s2))+mod((m-1),Hdim2)
!            SH%Ddip(i)%mat(l,m)=SH%Ddip(i)%mat_elem(psi2,spin2)
!           enddo
!          enddo

!         enddo


         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine set_sph_dipolar


        subroutine set_dipolar(this,SH,ex_list)
        use spinham_class
        use mpi_utils
        implicit none
        class(spins_hilbert)             :: this
        class(SpinHamiltonian)           :: SH
        integer                          :: i,j,s1,s2,ii,jj,v,l,m,ex,is
        integer, allocatable             :: ex_list(:,:)
        integer                          :: celli,cellj,Hdim1,Hdim2,Hdim12
        double precision                 :: dist0(3)
        double precision                 :: spin2(2),spin,psi(2),psi2(2,2)
        integer                          :: t1,t2,rate
        logical                          :: skip

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building Dipolar Network'
          flush(6)
         endif

         SH%nDdip=0

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif

             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SH%dipolar_thr)then

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then                 
               SH%nDdip=SH%nDdip+1
              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then
               SH%nDdip=SH%nDdip+1
              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo

         allocate(SH%Ddip(SH%nDdip))

         if(mpi_id.eq.0)  &
         write(*,*) '     Total Number of spin-spin dipolar interactions: ',SH%nDdip

         v=1

         do i=1,SH%nG
          do j=i,SH%nG
           do ii=1,this%nspins_pr 
           do celli=1,1!this%ntot
            do jj=1,this%nspins_pr 
            do cellj=1,this%ntot 

             s1=this%nspins_pr*(celli-1)+ii
             s2=this%nspins_pr*(cellj-1)+jj

             skip=.false.
             if(allocated(ex_list))then
              do ex=1,size(ex_list,1)
               if(s1.eq.ex_list(ex,1) .and. s2.eq.ex_list(ex,2) ) skip=.true.
               if(s1.eq.ex_list(ex,2) .and. s2.eq.ex_list(ex,1) ) skip=.true.
              enddo
             endif


             if(s2.le.s1 .or. skip) cycle

             if(this%dist(ii,celli,jj,cellj).le.SH%dipolar_thr)then

              dist0=this%dist_vec_pbc(this%x(s1,:),this%x(s2,:))  

              if(this%kind(s1).eq.SH%G(i)%kind .and. this%kind(s2).eq.SH%G(j)%kind )then

               call SH%Ddip(v)%make_D(SH%G(i)%G,SH%G(j)%G, &
                    this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
               SH%Ddip(v)%kind(1)=s1
               SH%Ddip(v)%kind(2)=s2

               v=v+1

              endif

              if(this%kind(s1).eq.SH%G(j)%kind .and. this%kind(s2).eq.SH%G(i)%kind .and. & 
                 SH%G(j)%kind .ne. SH%G(i)%kind  )then

               call SH%Ddip(v)%make_D(SH%G(j)%G,SH%G(i)%G,  &
                    this%bohr_mag(this%kind(s1)),this%bohr_mag(this%kind(s2)),dist0,this%dist(ii,celli,jj,cellj))
               SH%Ddip(v)%kind(1)=s1
               SH%Ddip(v)%kind(2)=s2

               v=v+1

              endif

             endif

            enddo
            enddo
           enddo
           enddo
          enddo
         enddo


         ! compute matrix representation of SH%Ddips

!         do i=1,SH%nDdip
                
!          s1=SH%Ddip(i)%kind(1)
!         s2=SH%Ddip(i)%kind(2)
                
!          spin2(1)=this%spin(this%kind(s1))
!          spin2(2)=this%spin(this%kind(s2))

!          Hdim1=(2*spin2(1)+1)
!          Hdim2=(2*spin2(2)+1)
!          Hdim12=Hdim1*Hdim2
!          allocate(SH%Ddip(i)%mat(Hdim12,Hdim12))
!          SH%Ddip(i)%mat=(0.0d0,0.0d0)

!          do l=1,Hdim12
!           do m=1,Hdim12
!            psi2(1,1)=-this%spin(this%kind(s1))+int((l-1)/Hdim2)
!            psi2(1,2)=-this%spin(this%kind(s2))+mod((l-1),Hdim2)
!            psi2(2,1)=-this%spin(this%kind(s1))+int((m-1)/Hdim2)
!            psi2(2,2)=-this%spin(this%kind(s2))+mod((m-1),Hdim2)
!            SH%Ddip(i)%mat(l,m)=SH%Ddip(i)%mat_elem(psi2,spin2)
!           enddo
!          enddo

!         enddo


         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine set_dipolar



        subroutine make_Hmat_2(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use spinham_class
        implicit none
        class(spins_hilbert)       :: this
        type(dist_cmplx_mat)       :: AA,BB
        integer                    :: l2,l,size_block
        integer                    :: ii,jj,i1
        integer                    :: t1,t2,rate,indxl2g
        complex(8)                 :: val

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '     Building Hamiltonian Matrix'
          flush(6)
         endif

        !  Run over Hamiltonian rows for each process and fill them with
        !  matrix elements
       
         allocate(this%H0(this%ntot))
         do l=1,this%ntot
          size_block=this%kblc(l+1)-this%kblc(l)
          call this%H0(l)%set(size_block,size_block,NB,MB)  
         enddo

         call AA%set(this%Hdim,this%Hdim,NB,MB)
         AA%mat=(0.0d0,0.0d0)

         do l=1,size(this%Hnodes,1)
             
          ii=this%Hnodes(l,1)
          jj=this%Hnodes(l,2)
           
          AA%mat(ii,jj)=this%get_Hij_2(this%Hnodes(l,1:4))

         enddo ! l

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

         if(this%ntot.gt.1)then

          if(mpi_id.eq.0)then
           call system_clock(t1,rate)        
           write(*,*) '     Rotating Hamiltonian in T invariant basis'
           flush(6)
          endif
         
         ! rotate H in Tinvariant basis e block-factorize it
                 
          do l=1,this%ntot

           size_block=this%kblc(l+1)-this%kblc(l)
           i1=this%kblc(l)+1

           call BB%set(size_block,this%Hdim,NB,MB)

           call pzgemm('C','N',size_block,this%Hdim,this%Hdim,&
                       (1.0d0,0.0d0),this%kbasis%mat,1,i1,this%kbasis%desc,AA%mat,&
                       1,1,AA%desc,(0.0d0,0.0d0),BB%mat,1,1,BB%desc) 
 
           call pzgemm('N','N',size_block,size_block,this%Hdim,&
                       (1.0d0,0.0d0),BB%mat,1,1,BB%desc,this%kbasis%mat,1,i1,this%kbasis%desc,&
                       (0.0d0,0.0d0),this%H0(l)%mat,1,1,this%H0(l)%desc)

           call BB%dealloc()

           this%H0(l)%mat=this%H0(l)%mat*this%ntot
           
          enddo

         else

          this%H0(1)%mat=AA%mat

         endif

         if(mpi_phonons_id.eq.0)then
          do l=1,this%ntot
           size_block=this%kblc(l+1)-this%kblc(l)
           jj=this%H0(l)%get_nze(1.0D-9)
           if(mpi_blacs_id.eq.0) writE(*,*) '     Spin Hamiltonian Matrix Sparsity:',100*(1-jj/dble(size_block**2)),'%'
          enddo
         endif

         call AA%dealloc()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_Hmat_2


        subroutine make_Hmat(this)
        use mpi
        use mpi_utils
        use blacs_utils
        use spinham_class
        implicit none
        class(spins_hilbert)       :: this
        type(dist_cmplx_mat)       :: AA,BB
        integer                    :: l2,l,size_block
        integer                    :: ii,jj,i1,i
        integer                    :: t1,t2,rate,indxl2g
        complex(8)                 :: val

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '   Building Hamiltonian Matrix'
          flush(6)
         endif

        !  Run over Hamiltonian rows for each process and fill them with
        !  matrix elements
       
         allocate(this%H0(this%ntot))
         do l=1,this%ntot
          size_block=this%kblc(l+1)-this%kblc(l)
          call this%H0(l)%set(size_block,size_block,NB,MB)  
         enddo

         call AA%set(this%Hdim,this%Hdim,NB,MB)
         AA%mat=(0.0d0,0.0d0)

         open(1212,file='h0.dat')

         do ii=1,size(AA%mat,1)
          do jj=1,size(AA%mat,2)
 
!           do icntx=1,size(context)
!            call blacs_gridinfo(context(icntx),nprow,npcol,myrow,mycol)
!            if(myrow.ne.-1)then
             l=indxl2g(ii,NB,myrow,0,nprow)
             l2=indxl2g(jj,MB,mycol,0,npcol)
!            endif
!           enddo
           
           AA%mat(ii,jj)=this%get_Hij(this%SH,l,l2)

          enddo ! l2
         enddo ! l

         close(1212)

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

         if(this%ntot.gt.1)then

          if(mpi_id.eq.0)then
           call system_clock(t1,rate)        
           write(*,*) '          Rotating Hamiltonian in T invariant basis'
           flush(6)
          endif
         
         ! rotate H in Tinvariant basis e block-factorize it
                 
          do l=1,this%ntot

           size_block=this%kblc(l+1)-this%kblc(l)
           i1=this%kblc(l)+1

           call BB%set(size_block,this%Hdim,NB,MB)

           call pzgemm('C','N',size_block,this%Hdim,this%Hdim,&
                       (1.0d0,0.0d0),this%kbasis%mat,1,i1,this%kbasis%desc,AA%mat,&
                       1,1,AA%desc,(0.0d0,0.0d0),BB%mat,1,1,BB%desc) 
 
           call pzgemm('N','N',size_block,size_block,this%Hdim,&
                       (1.0d0,0.0d0),BB%mat,1,1,BB%desc,this%kbasis%mat,1,i1,this%kbasis%desc,&
                       (0.0d0,0.0d0),this%H0(l)%mat,1,1,this%H0(l)%desc)

           call BB%dealloc()

           this%H0(l)%mat=this%H0(l)%mat*this%ntot
           
          enddo

         else

          this%H0(1)%mat=AA%mat

         endif

         do l=1,this%ntot
          size_block=this%kblc(l+1)-this%kblc(l)
          jj=this%H0(l)%get_nze(1.0D-9)
          if(mpi_id.eq.0) writE(*,*) '          H0 Matrix Sparsity:',100*(1-jj/dble(size_block**2)),'%'
         enddo

         call AA%dealloc()

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_Hmat
        
        subroutine make_basis_L(this)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)            :: this
        integer                         :: k,i,j


         this%Ldim=this%Hdim**2
         allocate(this%Lbasis(this%Ldim,2))

         k=1
         do i=1,this%hdim
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

        return
        end subroutine make_basis_L

        subroutine make_basis_H(this,max_ex,max_corr,max_dist,tinv)
        use mpi
        use mpi_utils
        use blacs_utils
        implicit none
        class(spins_hilbert)            :: this
        integer                         :: si,v,i,sj,ii,jj,indxl2g
        integer                         :: max_ex,max_corr
        double precision, allocatable   :: max_dist(:,:)
        integer                         :: t1,t2,rate
        integer, allocatable            :: basis_tmp(:)
        logical                         :: tinv


         if (max_ex.eq.-1 ) max_ex=this%nspins
         if (max_ex.ge.this%nspins/2) max_ex=this%nspins
                 
         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) '' 
          write(*,*) '     Building Basis Set'
          write(*,*) '     Maximum number of spin excitations: ',max_ex
          flush(6)
         endif

         allocate(basis_tmp(this%nspins))

         si=1
         v=1

         call get_Hdim(this,basis_tmp,max_ex,max_dist,v,si)

         this%Hdim=(v-1)
         if (max_ex.lt.this%nspins/2.and.tinv) this%Hdim=this%Hdim*2

         allocate(this%basis(this%Hdim,this%nspins))

         si=1
         v=1

         call build_basis(this,basis_tmp,max_ex,max_dist,v,si)

         if (max_ex.lt.this%nspins/2.and.tinv) then
          si=this%Hdim
          do v=1,this%Hdim/2
           this%basis(si,:)=-1.0d0*this%basis(v,:)
           si=si-1
          enddo
         endif

         deallocate(basis_tmp)

         if(mpi_id.eq.0)then
!          open(15,file='basis')
          do v=1,this%Hdim
!           write(*,*) this%basis(v,:)
          enddo
!          close(15)
          write(*,*) '     Total Hilbert space size: ',this%Hdim
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif

        return
        end subroutine make_basis_H

        recursive subroutine get_Hdim(this,basis_tmp,max_ex,max_dist,v,si)
        implicit none
        class(spins_hilbert)          :: this
        integer                       :: i,j,si,sj,v,max_ex,nex
        integer                       :: celli,cellj,k1,k2
        double precision,allocatable  :: max_dist(:,:)
        integer, allocatable          :: basis_tmp(:)
        logical                       :: skip_spin=.false.

        if ( v.eq.1 .and. si.eq.1 ) then
         do j=si+1,this%nspins
          basis_tmp(j)=-nint(2*this%spin(this%kind(j)))
         enddo
        endif

        skip_spin=.false.

        do j=1,si-1
         celli=int((si-1)/this%nspins_pr)+1
          cellj=int((j-1)/this%nspins_pr)+1
          k1=mod((si-1),this%nspins_pr)+1
          k2=mod((j-1),this%nspins_pr)+1
         if (basis_tmp(j).ne.-nint(2*this%spin(this%kind(j)))) then 
          if(abs(this%dist(k1,celli,k2,cellj)).gt.max_dist(this%kind(j),this%kind(si))) skip_spin=.true.        
         endif
        enddo

        if(skip_spin)then

         basis_tmp(si)=-nint(2*this%spin(this%kind(si)))
         if(si.lt.this%nspins)then
          sj=si+1  
          call get_Hdim(this,basis_tmp,max_ex,max_dist,v,sj)
         else
          v=v+1
         endif

        else

         if(si.lt.this%nspins)then

          do i=-nint(2*this%spin(this%kind(si))),nint(2*this%spin(this%kind(si))),2

           basis_tmp(si)=i

           do j=si+1,this%nspins
            basis_tmp(j)=-nint(2*this%spin(this%kind(j)))
           enddo

           call check_nex(this,basis_tmp,max_ex,skip_spin)

           if( .not. skip_spin )then
            sj=si+1  
            call get_Hdim(this,basis_tmp,max_ex,max_dist,v,sj)
           else
            v=v+1
           endif

          enddo

         else

          do i=-nint(2*this%spin(this%kind(si))),nint(2*this%spin(this%kind(si))),2
           basis_tmp(si)=i
           v=v+1
          enddo

         endif
        endif

        return
        end subroutine get_Hdim


        recursive subroutine build_basis(this,basis_tmp,max_ex,max_dist,v,si)
        use blacs_utils        
        implicit none
        class(spins_hilbert)          :: this
        integer                       :: i,j,si,sj,v,max_ex
        integer                       :: celli,cellj,k1,k2
        integer                       :: nex,l
        double precision,allocatable  :: max_dist(:,:)
        integer, allocatable          :: basis_tmp(:)
        logical                       :: skip_spin=.false.

        if ( v.eq.1 .and. si.eq.1 ) then
         do j=si+1,this%nspins
          basis_tmp(j)=-nint(2*this%spin(this%kind(j)))
         enddo
        endif

        skip_spin=.false.

        do j=1,si-1
         celli=int((si-1)/this%nspins_pr)+1
          cellj=int((j-1)/this%nspins_pr)+1
          k1=mod((si-1),this%nspins_pr)+1
          k2=mod((j-1),this%nspins_pr)+1
         if (basis_tmp(j).ne.-nint(2*this%spin(this%kind(j)))) then 
          if(abs(this%dist(k1,celli,k2,cellj)).gt.max_dist(this%kind(j),this%kind(si))) skip_spin=.true.        
         endif
        enddo

        if(skip_spin)then

          basis_tmp(si)=-nint(2*this%spin(this%kind(si)))
          if(si.lt.this%nspins)then
           sj=si+1  
           call build_basis(this,basis_tmp,max_ex,max_dist,v,sj)
          else
           this%basis(v,:)=basis_tmp(:)/2.0d0
           v=v+1
          endif

        else

        if(si.lt.this%nspins)then

         do i=-nint(2*this%spin(this%kind(si))),nint(2*this%spin(this%kind(si))),2
                   
          basis_tmp(si)=i

          do j=si+1,this%nspins
           basis_tmp(j)=-nint(2*this%spin(this%kind(j)))
          enddo

          call check_nex(this,basis_tmp,max_ex,skip_spin)

          if( .not. skip_spin )then
           sj=si+1  
           call build_basis(this,basis_tmp,max_ex,max_dist,v,sj)
          else
           this%basis(v,:)=basis_tmp(:)/2.0d0
           v=v+1
          endif

         enddo

        else

          do i=-nint(2*this%spin(this%kind(si))),nint(2*this%spin(this%kind(si))),2
           basis_tmp(si)=i
           this%basis(v,:)=basis_tmp(:)/2.0d0
           v=v+1
          enddo
       
        endif
        endif

        return

        end subroutine build_basis


        subroutine check_nex(this,basis_tmp,max_ex,skip_spin)
        implicit none
        class(spins_hilbert)  :: this
        integer               :: i,j,k,nex,max_ex
        integer,allocatable   :: basis_tmp(:),cas_nex(:)
        logical               :: skip_spin

         skip_spin=.false.
        
         if( allocated(this%active_space) )then            
          allocate(cas_nex(size(this%active_space)))
          cas_nex=0
         endif

         nex=0
         do i=1,this%nspins
          if((basis_tmp(i)+nint(2*this%spin(this%kind(i)))).ne.0) then

         ! total nex

           nex=nex+1

           if(nex.ge.max_ex)then
            skip_spin=.true.
            return
           endif            

         ! compute the sub_space nex if any

           if( allocated(this%active_space) )then            
            do j=1,size(this%active_space)
             do k=1,size(this%active_space(j)%kind)
            
              if(this%active_space(j)%kind(k).eq.this%kind(i))then
               cas_nex(j)=cas_nex(j)+1
              endif

             enddo

             if(cas_nex(j).ge.this%active_space(j)%max_ex)then
              skip_spin=.true.
              return
             endif            

            enddo                   
           endif


          endif
         enddo

         if( allocated(cas_nex) ) deallocate(cas_nex)

        return
        end subroutine check_nex
       
        function get_Tij(this,i,j,v) result (Tij)
        implicit none 
        class(spins_hilbert) :: this
        complex(8)           :: Tij
        integer              :: i,j,l,m,v

         Tij=(1.0d0,0.0d0)
               
         do l=1,this%nspins

          m=this%tr_map(l,v)

          if(abs(this%basis(i,l)-this%basis(j,m)).gt.1.0e-4) then
           Tij=(0.0d0,0.d0)
           goto 20
          endif

         enddo

20       continue                

        return
        end function get_Tij

        subroutine make_kbasis(this)
        use lists_class
        use sparse_class
        use mpi_utils
        use units_parms
        use general_types_class
        use lapack_diag_asimm
        implicit none
        class(spins_hilbert)    :: this
        type(sarc),allocatable  :: list_i(:,:),list_tmp(:)
        type(list)              :: AJ,Aval
        type(csr_mat_cmplx)     :: Tmat
        integer                 :: i,j,v,ii,jj,l,k
        integer                 :: vv,ll,Ti_old
        integer                 :: t1,t2,rate
        integer                 :: size_block,kk,ktot,Ti,nblc
        integer, allocatable    :: mapp(:),blc(:)
        complex(8), allocatable :: kval_tmp(:),Tmat_tmp(:,:)
        complex(8)              :: v_tmp,norm,mat_elem
        double precision        :: kx,ky,kz,ki(3),kj(3)
        logical, allocatable    :: check(:)
        integer                 :: nloc,nstart,istart,nelem_tot
        integer,allocatable     :: proc_grid(:)


         if(this%ntot.eq.1)then

          allocate(this%kblc(2))
          allocate(this%klist(1,3))
          this%klist(1,:)=0.0d0
          this%kblc(1)=0
          this%kblc(2)=this%Hdim
          return

         endif

         if(mpi_id.eq.0)then
          call system_clock(t1,rate)        
          write(*,*) ''
          write(*,*) '     Building translational invariant basis set'
          flush(6)
         endif
        
         allocate(list_i(this%Hdim,3))
         do i=1,this%Hdim
          do j=1,3
           list_i(i,j)%k(:)=0.0d0
          enddo
         enddo

         Ti=1

         if(Ti.eq.1 .and. this%nx.eq.1 ) Ti=Ti+1
         if(Ti.eq.2 .and. this%ny.eq.1 ) Ti=Ti+1
         if(Ti.eq.3 .and. this%nz.eq.1 ) stop

         Ti_old=Ti
          
         ! parallelizza sulle righe i

         call mpi_dist_nprocess(this%Hdim,nloc,nstart,proc_grid,mpi_blacs_world)

         call AJ%init()
         call Aval%init()

         allocate(Tmat%AI(this%Hdim+1))
         Tmat%AI(1)=0

         i=nstart
         do l=1,nloc
          Tmat%AI(i+1)=0   
          do j=1,this%Hdim             
           
           mat_elem=this%get_Tij(i,j,Ti) 

           if(abs(dble(mat_elem)).gt.1.0d-10.or.abs(aimag(mat_elem)).gt.1.0d-10)then
            Tmat%AI(i+1)=Tmat%AI(i+1)+1  
            call AJ%add_node(j)
            call Aval%add_node(mat_elem)
           endif

          enddo
          i=i+1
         enddo

         nelem_tot=AJ%nelem
         call mpi_allreduce(nelem_tot,nelem_tot,1,mpi_integer,mpi_sum,mpi_blacs_world,err)
         allocate(Tmat%AJ(nelem_tot))
         allocate(Tmat%A(nelem_tot))
         Tmat%A=(0.0d0,0.0d0)
         Tmat%AJ=0

         do i=1,this%Hdim
          call mpi_bcast(Tmat%AI(i+1),1,mpi_integer,proc_grid(i),mpi_blacs_world,err)
         enddo
         do i=2,this%Hdim+1
          Tmat%AI(i)=Tmat%AI(i-1)+Tmat%AI(i)
         enddo

         call AJ%reboot()
         call Aval%reboot()

         istart=Tmat%AI(nstart)
         do i=1,AJ%nelem
                     
          call AJ%rd_val(Tmat%AJ(i+istart))
          call Aval%rd_val(Tmat%A(i+istart))

          call AJ%skip()
          call Aval%skip()
        
         enddo

         do i=1,nelem_tot
          call mpi_allreduce(Tmat%AJ(i),Tmat%AJ(i),1,mpi_integer,mpi_sum,mpi_blacs_world,err)
          call mpi_allreduce(Tmat%A(i),Tmat%A(i),1,mpi_double_complex,mpi_sum,mpi_blacs_world,err)
         enddo

         call AJ%delete()
         call Aval%delete()                        

         call Tmat%do_coord()

         if(allocated(blc)) deallocate(blc)
         if(allocated(mapp)) deallocate(mapp)
         call Tmat%block(blc,mapp)

         ktot=1
        
         do j=1,size(blc)-1
         
          size_block=blc(j+1)-blc(j)
          allocate(Tmat_tmp(size_block,size_block))
          allocate(kval_tmp(size_block))
          Tmat_tmp=(0.0d0,0.0d0)

          do ii=1,size_block
           do jj=1,size_block             
            Tmat_tmp(ii,jj)=this%get_Tij(mapp(blc(j)+ii),mapp(blc(j)+jj),Ti)
           enddo
          enddo
        
          call new_diag2(size_block,Tmat_tmp,kval_tmp) 

          do ii=1,size_block
           kx=dble(log(kval_tmp(ii))/(2.0d0*pi*cmplx(0.0d0,1.0d0,8)))
           if((kx+0.50000).lt.1.0D-4) kx=0.50000000
           allocate(list_i(ktot,Ti)%v(size_block))
           allocate(list_i(ktot,Ti)%vi(size_block))
           list_i(ktot,Ti)%k(Ti)=kx
           do jj=1,size_block
            list_i(ktot,Ti)%v(jj)=Tmat_tmp(jj,ii)
            list_i(ktot,Ti)%vi(jj)=mapp(blc(j)+jj)
           enddo
           ktot=ktot+1
          enddo
          
          deallocate(Tmat_tmp)
          deallocate(kval_tmp)

         enddo

         if(allocated(Tmat%A)) deallocate(Tmat%A)
         if(allocated(Tmat%AI)) deallocate(Tmat%AI)
         if(allocated(Tmat%AJ)) deallocate(Tmat%AJ)
         if(allocated(Tmat%AC)) deallocate(Tmat%AC)

        ! Scorri su le dimensioni Ti rimanenti calcolando 
        ! gli elemenenti di matrice nella base di T_old

          Ti_old=Ti
          Ti=Ti_old+1

         do while (Ti.le.3)

          if(Ti.eq.2 .and. this%ny.eq.1 )then
           Ti=Ti+1
           cycle
          endif
          if(Ti.eq.3 .and. this%nz.eq.1 )then
           Ti=Ti+1
           cycle
          endif

          call AJ%init()
          call Aval%init()

          allocate(Tmat%AI(this%Hdim+1))
          Tmat%AI(1)=0

          i=nstart
          do k=1,nloc
           Tmat%AI(i+1)=0
           do j=1,this%Hdim             

            if(abs(list_i(i,Ti_old)%k(Ti_old)-list_i(j,Ti_old)%k(Ti_old)).gt.1.0e-3) cycle
           
            mat_elem=(0.0d0,0.0d0)
            do ii=1,size(list_i(i,Ti_old)%vi) 
             do jj=1,size(list_i(j,Ti_old)%vi) 
              l=list_i(i,Ti_old)%vi(ii)
              v=list_i(j,Ti_old)%vi(jj)
              mat_elem=mat_elem+this%get_Tij(l,v,Ti)&
                       *conjg(list_i(i,Ti_old)%v(ii))*list_i(j,Ti_old)%v(jj)

             enddo
            enddo

            if(abs(dble(mat_elem)).gt.1.0d-12.or.abs(aimag(mat_elem)).gt.1.0d-12)then
             Tmat%AI(i+1)=Tmat%AI(i+1)+1
             call AJ%add_node(j)
             call Aval%add_node(mat_elem)
            endif

           enddo
           i=i+1
          enddo

          nelem_tot=AJ%nelem
          call mpi_allreduce(nelem_tot,nelem_tot,1,mpi_integer,mpi_sum,mpi_blacs_world,err)
          allocate(Tmat%AJ(nelem_tot))
          allocate(Tmat%A(nelem_tot))
          Tmat%A=(0.0d0,0.0d0)
          Tmat%AJ=0

          do i=1,this%Hdim
           call mpi_bcast(Tmat%AI(i+1),1,mpi_integer,proc_grid(i),mpi_blacs_world,err)
          enddo
          do i=2,this%Hdim+1
           Tmat%AI(i)=Tmat%AI(i-1)+Tmat%AI(i)
          enddo

          call AJ%reboot()
          call Aval%reboot()

          istart=Tmat%AI(nstart)
          do i=1,AJ%nelem
                     
          call AJ%rd_val(Tmat%AJ(i+istart))
          call Aval%rd_val(Tmat%A(i+istart))

          call AJ%skip()
          call Aval%skip()
        
         enddo

         do i=1,nelem_tot
          call mpi_allreduce(Tmat%AJ(i),Tmat%AJ(i),1,mpi_integer,mpi_sum,mpi_blacs_world,err)          
          call mpi_allreduce(Tmat%A(i),Tmat%A(i),1,mpi_double_complex,mpi_sum,mpi_blacs_world,err)
         enddo
                  
          call Tmat%do_coord()

          if(allocated(blc)) deallocate(blc)
          if(allocated(mapp)) deallocate(mapp)
          call Tmat%block(blc,mapp)


          ktot=1

          do j=1,size(blc)-1
         
           size_block=blc(j+1)-blc(j)
           allocate(Tmat_tmp(size_block,size_block))
           allocate(kval_tmp(size_block))
           Tmat_tmp=(0.0d0,0.0d0)

           do ii=1,size_block
            do jj=1,size_block             
             do ll=1,size(list_i(mapp(blc(j)+ii),Ti_old)%vi) 
              do vv=1,size(list_i(mapp(blc(j)+jj),Ti_old)%vi) 

               l=list_i(mapp(blc(j)+ii),Ti_old)%vi(ll)
               v=list_i(mapp(blc(j)+jj),Ti_old)%vi(vv)

               Tmat_tmp(ii,jj)=Tmat_tmp(ii,jj)+this%get_Tij(l,v,Ti)&
                        *conjg(list_i(mapp(blc(j)+ii),Ti_old)%v(ll))  &
                        *list_i(mapp(blc(j)+jj),Ti_old)%v(vv)

              enddo
             enddo
            enddo
           enddo
       
           call new_diag2(size_block,Tmat_tmp,kval_tmp) 

           do ii=1,size_block
            kx=dble(log(kval_tmp(ii))/(2.0d0*pi*cmplx(0.0d0,1.0d0,8)))
            if((kx+0.50000).lt.1.0D-4) kx=0.50000000
            list_i(ktot,Ti)%k(Ti_old)=list_i(mapp(blc(j)+1),Ti_old)%k(Ti_old)
            list_i(ktot,Ti)%k=list_i(mapp(blc(j)+1),Ti_old)%k
            list_i(ktot,Ti)%k(Ti)=dble(kx) 
            allocate(list_i(ktot,Ti)%v(this%Hdim))
            list_i(ktot,Ti)%v(:)=(0.0d0,0.0d0)
            do jj=1,size_block
             do ll=1,size(list_i(mapp(blc(j)+jj),Ti_old)%vi)
              l=list_i(mapp(blc(j)+jj),Ti_old)%vi(ll)
              list_i(ktot,Ti)%v(l)=list_i(ktot,Ti)%v(l)+               &
                Tmat_tmp(jj,ii)*list_i(mapp(blc(j)+jj),Ti_old)%v(ll)
             enddo
            enddo

            ! compress the vector

            call AJ%init()
            call Aval%init()
            do jj=1,this%Hdim
             if( abs(dble(list_i(ktot,Ti)%v(jj))).gt.1.0d-10 .or.  &
                 abs(aimag(list_i(ktot,Ti)%v(jj))).gt.1.0d-10 )then
              call AJ%add_node(jj)
              call Aval%add_node(list_i(ktot,Ti)%v(jj))
             endif
            enddo
    
            call aj%reboot()
            call aval%reboot()
            deallocate(list_i(ktot,Ti)%v)
            allocate(list_i(ktot,Ti)%vi(AJ%nelem))
            allocate(list_i(ktot,Ti)%v(AJ%nelem))
            do jj=1,AJ%nelem                   
             call AJ%rd_val(list_i(ktot,Ti)%vi(jj)) 
             call Aval%rd_val(list_i(ktot,Ti)%v(jj)) 
             call AJ%skip()
             call Aval%skip()
            enddo
            call AJ%delete()
            call Aval%delete()

            ! end compression

            ktot=ktot+1
           enddo
          
           deallocate(Tmat_tmp)
           deallocate(kval_tmp)

          enddo

          if(allocated(Tmat%A)) deallocate(Tmat%A)
          if(allocated(Tmat%AI)) deallocate(Tmat%AI)
          if(allocated(Tmat%AJ)) deallocate(Tmat%AJ)
          if(allocated(Tmat%AC)) deallocate(Tmat%AC)

          ! deallocate Told

          Ti_old=Ti
          Ti=Ti+1

         enddo

         Ti=Ti_old

         ! allocate kbasis

         if(allocated(mapp)) deallocate(mapp)
         allocate(this%kblc(this%ntot+1))
         allocate(this%klist(this%ntot,3))
         allocate(check(this%Hdim))
         allocate(mapp(this%Hdim))

         l=1
         i=1
         check=.false.
         this%kblc(1)=0
         nblc=0

         do while(.not. all(check))

          if(check(i))then
           i=i+1
           cycle
          endif
          
          nblc=nblc+1
          this%kblc(nblc+1)=this%kblc(nblc)+1
          mapp(l)=i
          l=l+1
          check(i)=.true.
          
          ki(1)=list_i(i,Ti)%k(1)
          ki(2)=list_i(i,Ti)%k(2)
          ki(3)=list_i(i,Ti)%k(3)
          this%klist(nblc,1)=ki(1)
          this%klist(nblc,2)=ki(2)
          this%klist(nblc,3)=ki(3)

          do j=1,this%Hdim
           
           if(check(j))cycle

           kj(1)=list_i(j,Ti)%k(1)
           kj(2)=list_i(j,Ti)%k(2)
           kj(3)=list_i(j,Ti)%k(3)

           if(abs(kj(1)-ki(1)).lt.1.0d-4  .and. &
              abs(kj(2)-ki(2)).lt.1.0d-4  .and. &
              abs(kj(3)-ki(3)).lt.1.0d-4  )then

            this%kblc(nblc+1)=this%kblc(nblc+1)+1
            mapp(l)=j
            l=l+1
            check(j)=.true.
           endif

          enddo

         enddo ! end while

         ! set kbasis

         call this%kbasis%set(this%Hdim,this%Hdim,NB,MB)
         this%kbasis%mat=(0.0d0,0.0d0)

         do i=1,this%ntot
          size_block=this%kblc(i+1)-this%kblc(i)
          do j=1,size_block
           jj=this%kblc(i)+j

           do l=1,size(list_i(mapp(this%kblc(i)+j),Ti)%vi)
            ii=list_i(mapp(this%kblc(i)+j),Ti)%vi(l)
            call pzelset(this%kbasis%mat,ii,jj,this%kbasis%desc, &
                         list_i(mapp(this%kblc(i)+j),Ti)%v(l))

           enddo

          enddo
         enddo

         ! build kbasis as sparse

         if(this%sparse) call this%sp_kbasis%tosparse(this%Hdim,this%kbasis)

         !!!!!

         jj=this%kbasis%get_nze(1.0D-9)
         if(mpi_id.eq.0) writE(*,*) '     K-basis Sparsity:',100*(1-jj/dble(this%Hdim**2)),'%'
          
         ! free memory

         do Ti=1,3
          do i=1,this%Hdim
           if(allocated(list_i(i,Ti)%vi)) deallocate(list_i(i,Ti)%vi)
           if(allocated(list_i(i,Ti)%v)) deallocate(list_i(i,Ti)%v)
          enddo
         enddo

         deallocate(list_i)

         if(mpi_id.eq.0)then
          call system_clock(t2)
          write(*,*) '     Task completed in ',real(t2-t1)/real(rate),'s'
          flush(6)
         endif
         
        return
        end subroutine make_kbasis

        subroutine make_rot(this,pulse)
        use mpi_utils
        use blacs_utils
        use rotations_class
        use pulses_class
        use units_parms
        use exponential_class
        implicit none
        class(spins_hilbert)            :: this
        type(general_pulse),allocatable :: pulse(:)
        integer                         :: i,ii,k,indxl2g
        complex(8)                      :: angle
        type(dist_cmplx_mat)            :: Sn
        integer                         :: spin


         do i=1,size(pulse)
       
          call Sn%set(this%Hdim,this%Hdim,NB,MB)
          Sn%mat=(0.0d0,0.0d0)

          if(allocated(pulse(i)%rot%mat)) call pulse(i)%rot%dealloc()
          call pulse(i)%rot%set(this%Hdim,this%Hdim,NB,MB)
          pulse(i)%rot%mat=(0.0d0,0.0d0)

          do k=1,this%s2print
           if(this%print_si(k).eq.pulse(i)%spin)then
            spin=k
            exit
           endif
          enddo

          pulse(i)%n=pulse(i)%n/  &
                sqrt(pulse(i)%n(1)**2+pulse(i)%n(2)**2+pulse(i)%n(3)**2)

          Sn%mat=this%Sx(spin)%mat*pulse(i)%n(1)+     &
                 this%Sy(spin)%mat*pulse(i)%n(2)+     &
                 this%Sz(spin)%mat*pulse(i)%n(3)

          angle=cmplx(pulse(i)%beta,0.0d0,8)

          call exp_taylor(this%Hdim,Sn,pulse(i)%rot,angle,100) 
          call Sn%dealloc()

          ii=pulse(i)%rot%get_nze(1.0d-9)
          if(mpi_id.eq.0) write(*,*) '     Pulse Sparsity:',100*(1-ii/dble(this%Hdim**2)),'%'

         enddo

        return
        end subroutine make_rot

        subroutine rot_rho(this,pulse)
        use mpi_utils
        use blacs_utils
        use pulses_class
        implicit none
        class(spins_hilbert)             :: this
        type(general_pulse),allocatable  :: pulse(:)
        integer                          :: si,i,j,ii,jj,indxl2g
        type(dist_cmplx_mat)             :: AA,BB,newrho

         call newrho%set(this%Hdim,this%Hdim,NB,MB)          
         call AA%set(this%Hdim,this%Hdim,NB,MB) 
         newrho%mat=(0.0d0,0.0d0)
        
         do i=1,size(pulse)

          AA%mat=(0.0d0,0.0d0)

          if(pulse(i)%sx .and. (.not. pulse(i)%dx) )then

           call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                        (1.0d0,0.0d0),pulse(i)%rot%mat,1,1,pulse(i)%rot%desc,this%rho%mat,&
                         1,1,this%rho%desc,(0.0d0,0.0d0),AA%mat,1,1,AA%desc)  
          
          endif


          if( pulse(i)%dx .and. (.not. pulse(i)%sx) )then

           call pzgemm('N','C',this%Hdim,this%Hdim,this%Hdim,&
                        (1.0d0,0.0d0),this%rho%mat,1,1,this%rho%desc,pulse(i)%rot%mat,1,1,pulse(i)%rot%desc,&
                        (0.0d0,0.0d0),AA%mat,1,1,AA%desc) 

          endif


          if( pulse(i)%dx .and. pulse(i)%sx )then
 
           call BB%set(this%Hdim,this%Hdim,NB,MB)
           BB%mat=(0.0d0,0.0d0)

           call pzgemm('N','N',this%Hdim,this%Hdim,this%Hdim,&
                         (1.0d0,0.0d0),pulse(i)%rot%mat,1,1,pulse(i)%rot%desc,this%rho%mat,&
                         1,1,this%rho%desc,(0.0d0,0.0d0),BB%mat,1,1,BB%desc)  
           
           call pzgemm('N','C',this%Hdim,this%Hdim,this%Hdim,&
                         (1.0d0,0.0d0),BB%mat,1,1,BB%desc,pulse(i)%rot%mat,1,1,pulse(i)%rot%desc,&
                         (0.0d0,0.0d0),AA%mat,1,1,AA%desc)

           call BB%dealloc()

          endif

          newrho%mat=newrho%mat+pulse(i)%weight*AA%mat

         enddo
         
         this%rho%mat=newrho%mat

         call AA%dealloc()
         call newrho%dealloc()

        return
        end subroutine rot_rho

        end module hilbert_dist_class

