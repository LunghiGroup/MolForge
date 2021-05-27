        module spindy_parser

        contains

        subroutine parse_input(input,spindy,phondy,lattice,gsh)
        use spins_dist_rs_class
        use hilbert_dist_class
        use atoms_class
        use phonons_class
        use parser_class
        use pulses_class
        use control_variables
        implicit none
        class(spins_hilbert)             :: spindy,gsh
        class(brillouin)                 :: phondy
        class(atoms_group)               :: lattice
        character(len=20)                :: input
        character(len=:), allocatable    :: line,word
        integer            :: l
        logical            :: eof=.false.

         do

          call get_line(10,line,eof)
          if (eof) then
           operation='END'
           exit
          endif

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('&HILBERT_SPACE')

            call parse_global(spindy)
            operation='MAKE_HILBERT'
            exit

           case('&SPIN_H')

            call parse_spinham(spindy)
            operation='SET_SPINHAM'
            exit

           case('&GS_SPIN_H')

            call parse_spinham(gsh)
            operation='SET_GSH_SPINHAM'
            exit

           case('&SPH_H')

            call parse_sph(spindy,phondy,lattice)
            operation='SET_SPH'
            exit

           case('&SYSTEM')

            call parse_system(spindy)
            operation='SET_SYSTEM'
            exit

           case('&GS_SYSTEM')

            call parse_system(gsh)
            operation='SET_GSH_SYSTEM'
            exit

           case('&PULSE')

            call parse_pulse(spindy)
            operation='MAKE_PULSE'
            exit

           case('&SPIN_DYNAMICS')

            call parse_dynamics
            operation='SET_SPINDYN'
            exit

           case('PROPAGATE')

            call get_word(line,word,2)
            read(word,*) nsteps
            operation='PROPAGATE'
            exit

           case('BUILD_PROPAGATOR')

            call get_word(line,word,2)
            read(word,*) step_min
            call get_word(line,word,3)
            read(word,*) step_nmult

            step=step_min*step_nmult

            operation='MAKE_U'
            exit

           case('&DENSITY_MATRIX')

            call parse_rho0
            operation='MAKE_RHO0'
            exit

           case('&SPIN_ECHO')

            call parse_spin_echo
            operation='MAKE_ECHO'
            exit

           case('SH_MAPP')

            call get_word(line,word,2)
            read(word,*) lmax
            call get_word(line,word,3)
            call To_upper(word)
            if(trim(word).eq.'COMPRESS')then
             compress=.true.
             call get_word(line,word,4)
             read(word,*) lambda
            endif
            operation='MAKE_SH_MAPP'
            exit


           case default
           write(*,*) 'Word',trim(word),&
                      'do not recognized in subroutine parse_input'

          end select

         enddo

        return
        end subroutine parse_input

        subroutine parse_sph(spindy,phondy,lattice)
        use spins_dist_rs_class
        use hilbert_dist_class
        use atoms_class
        use phonons_class
        use spin_phonon_class
        use parser_class
        use control_variables
        implicit none
        class(spins_hilbert)             :: spindy
        class(brillouin)                 :: phondy
        class(atoms_group)               :: lattice
        type(Othermos)                   :: Otmp 
        type(Gthermos)                   :: Gtmp 
        type(DSIthermos)                 :: DSItmp 
        type(D2Sthermos)                 :: D2Stmp 
        type(Jthermos)                   :: Jtmp 
        type(Ot_list)                    :: Olist,O2list
        type(Jt_list)                    :: Jlist,J2list
        type(Gt_list)                    :: Glist,G2list
        type(DSIt_list)                  :: DSIlist,DSI2list
        type(D2St_list)                  :: D2Slist,D2S2list
        character(len=:), allocatable    :: line,word
        integer                          :: l,n,i
        logical                          :: eof=.false.


         call Glist%init()
         call G2list%init()

         call D2Slist%init()
         call D2S2list%init()

         call DSIlist%init()
         call DSI2list%init()

         call Jlist%init()

         call Olist%init()

         do

          call get_line(10,line,eof)
          if (eof) then
           operation='END'
           exit
          endif
          
          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('&PHONDY')
            call parse_phondy(phondy,lattice)

           case('&G_BATH')
            call get_word(line,word,2)
            read(word,*) Gtmp%kind
            call parse_Gbath(Gtmp)
            if (Gtmp%norder.eq.1)then
             spindy%make_Rmat=.true.
             spindy%SPH%nG=spindy%SPH%nG+1
             call Glist%add_node(Gtmp)
            endif
            if (Gtmp%norder.eq.2)then
             spindy%make_R2mat=.true.
             spindy%SPH2%nG=spindy%SPH2%nG+1
             call G2list%add_node(Gtmp)
            endif

           case('&O_BATH')
            spindy%make_Rmat=.true.
            spindy%SPH%nO=spindy%SPH%nO+1
            call get_word(line,word,2)
            read(word,*) Otmp%kind
            call get_word(line,word,3)
            read(word,*) Otmp%k
            call parse_Obath(Otmp)
            call Olist%add_node(Otmp)

           case('&J_BATH')
            spindy%make_Rmat=.true.
            spindy%SPH%nJ=spindy%SPH%nJ+1
            call get_word(line,word,2)
            read(word,*) Jtmp%kind(1)
            call get_word(line,word,3)
            read(word,*) Jtmp%kind(2)
            call parse_Jbath(Jtmp)
            call Jlist%add_node(Jtmp)

           case('&DSI_BATH')
            call get_word(line,word,2)
            read(word,*) DSItmp%kind
            call parse_DSIbath(DSItmp)
            if (DSItmp%norder.eq.1)then
             spindy%make_Rmat=.true.
             spindy%SPH%nDSI=spindy%SPH%nDSI+1
             call DSIlist%add_node(DSItmp)
            endif
            if (DSItmp%norder.eq.2)then
             spindy%make_R2mat=.true.
             spindy%SPH2%nDSI=spindy%SPH2%nDSI+1
             call DSI2list%add_node(DSItmp)
            endif

           case('&D2S_BATH')
            call get_word(line,word,2)
            read(word,*) D2Stmp%kind(1)
            call get_word(line,word,3)
            read(word,*) D2Stmp%kind(2)
            call parse_D2Sbath(D2Stmp)
            if (D2Stmp%norder.eq.1)then
             spindy%make_Rmat=.true.
             spindy%SPH%nD2S=spindy%SPH%nD2S+1
             call D2Slist%add_node(D2Stmp)
            endif
            if (D2Stmp%norder.eq.2)then
             spindy%make_R2mat=.true.
             spindy%SPH2%nD2S=spindy%SPH2%nD2S+1
             call D2S2list%add_node(D2Stmp)
            endif

           case('SPINSPIN')
            spindy%make_Rmat=.true.
            spindy%SPH%make_dipolar=.true.
            call get_word(line,word,2)
            read(word,*) spindy%SPH%dipolar_thr

           case('SPINSPIN2')
            spindy%make_R2mat=.true.
            spindy%SPH2%make_dipolar=.true.
            call get_word(line,word,2)
            read(word,*) spindy%SPH2%dipolar_thr

           case('MAPPING')
            call get_word(line,word,2)
            read(word,*) n
            allocate(spindy%SPH%mapp(n))
            allocate(spindy%SPH2%mapp(n))
            do i=1,n
             call get_word(line,word,i+2)
             read(word,*) spindy%SPH%mapp(i)  
             read(word,*) spindy%SPH2%mapp(i)  
            enddo

           case('SECULAR')
            spindy%make_SA=.true.

           case('PT2')
            spindy%make_PT2=.true.

           case('&END')

            if( spindy%SPH%nO.gt.0) then
             allocate(spindy%SPH%O_t(spindy%SPH%nO))
             call Olist%reboot()
             do l=1,spindy%SPH%nO
              call Olist%rd_node(spindy%SPH%O_t(l))
              call Olist%skip()
             enddo
            endif
            call Olist%delete()

            if( spindy%SPH%nG.gt.0) then
             allocate(spindy%SPH%G_t(spindy%SPH%nG))
             call Glist%reboot()
             do l=1,spindy%SPH%nG
              call Glist%rd_node(spindy%SPH%G_t(l))
              call Glist%skip()
             enddo
            endif
            call Glist%delete()

            if( spindy%SPH2%nG.gt.0) then
             allocate(spindy%SPH2%G_t(spindy%SPH2%nG))
             call G2list%reboot()
             do l=1,spindy%SPH2%nG
              call G2list%rd_node(spindy%SPH2%G_t(l))
              call G2list%skip()
             enddo
            endif
            call G2list%delete()

            if( spindy%SPH%nD2S.gt.0) then
             allocate(spindy%SPH%D2S_t(spindy%SPH%nD2S))
             call D2Slist%reboot()
             do l=1,spindy%SPH%nD2S
              call D2Slist%rd_node(spindy%SPH%D2S_t(l))
              call D2Slist%skip()
             enddo
            endif
            call D2Slist%delete()

            if( spindy%SPH2%nD2S.gt.0) then
             allocate(spindy%SPH2%D2S_t(spindy%SPH2%nD2S))
             call D2S2list%reboot()
             do l=1,spindy%SPH2%nD2S
              call D2S2list%rd_node(spindy%SPH2%D2S_t(l))
              call D2S2list%skip()
             enddo
            endif
            call D2S2list%delete()

            if( spindy%SPH%nDSI.gt.0) then
             allocate(spindy%SPH%DSI_t(spindy%SPH%nDSI))
             call DSIlist%reboot()
             do l=1,spindy%SPH%nDSI
              call DSIlist%rd_node(spindy%SPH%DSI_t(l))
              call DSIlist%skip()
             enddo
            endif
            call DSIlist%delete()

            if( spindy%SPH2%nDSI.gt.0) then
             allocate(spindy%SPH2%DSI_t(spindy%SPH2%nDSI))
             call DSI2list%reboot()
             do l=1,spindy%SPH2%nDSI
              call DSI2list%rd_node(spindy%SPH2%DSI_t(l))
              call DSI2list%skip()
             enddo
            endif
            call DSI2list%delete()
             
            return

          end select

         enddo

        stop
        end subroutine parse_sph


        subroutine parse_DSIbath(DSI)
        use spin_phonon_class
        use parser_class
        implicit none
        class(DSIthermos)                :: DSI
        character(len=:), allocatable    :: line,word
        double precision                 :: Dr(3) 
        integer                          :: l,s
        logical                          :: eof=.false.

         if(allocated(DSI%Dcart)) deallocate(DSI%Dcart)
         if(allocated(DSI%map_s2a)) deallocate(DSI%map_s2a)

         do

          call get_line(10,line,eof)

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('FILENAME')
            call get_word(line,word,2)
            open(12,file=trim(word))

           case('NORDER')
            call get_word(line,word,2)
            read(word,*) DSI%norder

           case('&END')
            read(12,*) DSI%nderiv
            allocate(DSI%Dcart(DSI%nderiv))
            allocate(DSI%map_s2a(DSI%nderiv,DSI%norder*2))
            do l=1,DSI%nderiv
             read(12,*) DSI%map_s2a(l,:)
             do s=1,3
              read(12,*) Dr(1),Dr(2),Dr(3)
              DSI%Dcart(l)%D(s,1:3)=Dr(1:3)
              DSI%Dcart(l)%kind=DSI%kind
             enddo
            enddo
            close(12)
            return

          end select

         enddo

        stop
        end subroutine parse_DSIbath


        subroutine parse_Obath(O)
        use spin_phonon_class
        use parser_class
        implicit none
        class(Othermos)                  :: O
        character(len=:), allocatable    :: line,word
        integer                          :: l,s
        logical                          :: eof=.false.
        double precision                 :: coeff

         if(allocated(O%Ocart)) deallocate(O%Ocart)
         if(allocated(O%map_s2a)) deallocate(O%map_s2a)

         do

          call get_line(10,line,eof)

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('FILENAME')
            call get_word(line,word,2)
            open(12,file=trim(word))

           case('NORDER')
            call get_word(line,word,2)
            read(word,*) O%norder

           case('&END')
            read(12,*) O%nderiv
            allocate(O%Ocart(O%nderiv))
            allocate(O%map_s2a(O%nderiv,O%norder*2))
            do l=1,O%nderiv
             O%Ocart(l)%k=O%k
             O%Ocart(l)%kind=O%kind
             read(12,*) O%map_s2a(l,:)
             allocate(O%Ocart(l)%B(2*O%Ocart(l)%k+1))
             allocate(O%Ocart(l)%q(2*O%Ocart(l)%k+1))
             do s=1,2*O%Ocart(l)%k+1
              read(12,*) O%Ocart(l)%q(s),coeff
              O%Ocart(l)%B(s)=cmplx(coeff,0.0d0,8)
             enddo
            enddo
            close(12)
            return

          end select

         enddo

        stop
        end subroutine parse_Obath


        subroutine parse_Gbath(G)
        use spin_phonon_class
        use parser_class
        implicit none
        class(Gthermos)                  :: G
        double precision                 :: Gr(3)
        character(len=:), allocatable    :: line,word
        integer                          :: l,s,t
        logical                          :: eof=.false.

         if(allocated(G%Gcart)) deallocate(G%Gcart)
         if(allocated(G%map_s2a)) deallocate(G%map_s2a)         

         do

          call get_line(10,line,eof)

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('FILENAME')
            call get_word(line,word,2)
            open(12,file=trim(word))

           case('NORDER')
            call get_word(line,word,2)
            read(word,*) G%norder

           case('&END')
            read(12,*) G%nderiv
            allocate(G%Gcart(G%nderiv))
            allocate(G%map_s2a(G%nderiv,G%norder*2))
            do l=1,G%nderiv
             read(12,*) G%map_s2a(l,:)
             do s=1,3
              read(12,*) Gr(1),Gr(2),Gr(3)
              G%Gcart(l)%G(s,1:3)=Gr(1:3)
              G%Gcart(l)%kind=G%kind
             enddo
            enddo
            close(12)
            return

          end select

         enddo

        stop
        end subroutine parse_Gbath


        subroutine parse_D2Sbath(D2S)
        use spin_phonon_class
        use parser_class
        implicit none
        class(D2Sthermos)                :: D2S
        double precision                 :: Dr(3)
        character(len=:), allocatable    :: line,word
        integer                          :: l,s
        logical                          :: eof=.false.

         if(allocated(D2S%Dcart)) deallocate(D2S%Dcart)
         if(allocated(D2S%map_s2a)) deallocate(D2S%map_s2a)

         do

          call get_line(10,line,eof)

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('FILENAME')
            call get_word(line,word,2)
            open(12,file=trim(word))

           case('NORDER')
            call get_word(line,word,2)
            read(word,*) D2S%norder

           case('&END')
            read(12,*) D2S%nderiv
            allocate(D2S%Dcart(D2S%nderiv))
            allocate(D2S%map_s2a(D2S%nderiv,D2S%norder*2))
            do l=1,D2S%nderiv
             read(12,*) D2S%map_s2a(l,:)
             do s=1,3
              read(12,*) Dr(1),Dr(2),Dr(3)
              D2S%Dcart(l)%D(s,1)=cmplx(Dr(1),0.0d0,8)
              D2S%Dcart(l)%D(s,2)=cmplx(Dr(2),0.0d0,8)
              D2S%Dcart(l)%D(s,3)=cmplx(Dr(3),0.0d0,8)
              D2S%Dcart(l)%kind=D2S%kind
             enddo
!             call D2S%Dcart(l)%traceless() 
            enddo
            close(12)
            return

          end select

         enddo

        stop
        end subroutine parse_D2Sbath


        subroutine parse_Jbath(J)
        use spin_phonon_class
        use parser_class
        implicit none
        class(Jthermos)                  :: J
        character(len=:), allocatable    :: line,word
        double precision                 :: Jr
        integer                          :: l,s
        logical                          :: eof=.false.

         if(allocated(J%Jcart)) deallocate(J%Jcart)
         if(allocated(J%map_s2a)) deallocate(J%map_s2a)

         do

          call get_line(10,line,eof)

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('FILENAME')
            call get_word(line,word,2)
            open(12,file=trim(word))

           case('NORDER')
            call get_word(line,word,2)
            read(word,*) J%norder

           case('&END')
            read(12,*) J%nderiv
            allocate(J%Jcart(J%nderiv))
            allocate(J%map_s2a(J%nderiv,J%norder*2))
            do l=1,J%nderiv
             read(12,*) J%map_s2a(l,:)
             read(12,*) Jr
             J%Jcart(l)%J=Jr
             J%Jcart(l)%kind=J%kind
            enddo
            close(12)
            return

          end select

         enddo

        stop
        end subroutine parse_Jbath

        subroutine parse_phondy(phondy,lattice)
        use atoms_class
        use phonons_class
        use parser_class
        use control_variables
        implicit none
        class(brillouin)                 :: phondy
        class(atoms_group)               :: lattice
        character(len=:), allocatable    :: line,word
        integer                          :: l
        logical                          :: eof=.false.

         do

          call get_line(10,line,eof)
          if (eof) then
           operation='END'
           exit
          endif

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('MAX_ENER')
            call get_word(line,word,2)
            read(word,*)  max_phonon_ener

           case('EFFECTIVE_LT')
            phondy%effective_lt=.true.
            write(*,*) 'Effective phonon linewidth on'

           case('TEMP')
            ntemps=1
            allocate(temp(1))
            call get_word(line,word,2)
            read(word,*) temp(1) 

           case('K_MESH')
            call get_word(line,word,2)
            read(word,*) phondy%nx
            call get_word(line,word,3)
            read(word,*) phondy%ny
            call get_word(line,word,4)
            read(word,*) phondy%nz
            phondy%ntot=phondy%nx*phondy%ny*phondy%nz

           case('SMEAR')
            call get_word(line,word,2)
            read(word,*) smear

           case('TYPE_SMEAR')
            call get_word(line,word,2)
            read(word,*) type_smear

           case('FC2')
            call get_word(line,word,2)
            read(word,*) fc2_file

           case('FC3')
            call get_word(line,word,2)
            read_fc3=.true.
            read(word,*) fc3_file
            
           case('&END')

           return

          end select


         enddo

        stop
        end subroutine parse_phondy


        subroutine parse_spin_echo
        use control_variables
        use parser_class
        implicit none
        character(len=:),allocatable   :: line,word
        integer            :: l
        logical            :: eof=.false.

         echo_nsteps=500
         echo_step=500
         echo_spin=1

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)
    
           case('S')
            call get_word(line,word,2)
            read(word,*) echo_spin

           case('NSTEPS')
            call get_word(line,word,2)
            read(word,*) echo_nsteps

           case('STEP')
            call get_word(line,word,2)
            read(word,*) echo_step

           case('&END')
            return

          end select

         enddo

        return
        end subroutine parse_spin_echo
    
        subroutine parse_rho0
        use control_variables
        use parser_class
        implicit none
        character(len=:),allocatable   :: line,word
        integer            :: l
        logical            :: eof=.false.

         type_rho0='FULLY_POLARIZED'
         spin_temp=1.0d0

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)
    
           case('TYPE')
            call get_word(line,word,2)
            call To_upper(word)
            read(word,*) type_rho0

           case('TEMP')
            call get_word(line,word,2)
            read(word,*) spin_temp

           case('RESTART')
            call get_word(line,word,2)
            read(word,*) rho_restart_file

           case('&END')
            return
           
          end select

         enddo

        return
        end subroutine parse_rho0


        subroutine parse_global(spindy)
        use control_variables
        use parser_class
        use hilbert_dist_class
        implicit none
        class(spins_hilbert)           :: spindy
        character(len=:),allocatable   :: line,word
        integer                        :: l,i,int_val
        double precision               :: dbl_val
        logical                        :: eof=.false.

         fulldiag=.true.
         dump_s=.false.
         nex_max=-1
         ncorr_max=-1
         nexclude=0
         dist_max=100.0d0
         dump_freq=1
         tinv=.false.
         allocate(print_si(1))
         s2print=1
         print_si(1)=-1

         call dist_kind%init()

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)
    
           case('FULLDIAG')
            fulldiag=.true.
      
           case('NODIAG')
            fulldiag=.false.

           case('DUMP_H0')
            spindy%printH0=.true.

           case('DUMP_HEIG')
            spindy%printHeig=.true.

           case('DUMP_RHO')
            spindy%printRho=.true.

           case('DUMP_RMAT')
            spindy%printRmat=.true.

           case('DUMP_S')
            dump_s=.true.

           case('MAX_DIST')
            call get_word(line,word,2)
            read(word,*) dist_max

           case('MAX_DIST_KIND')
            call get_word(line,word,2)
            read(word,*) int_val
            call dist_kind%add_node(int_val)
            call get_word(line,word,3)
            read(word,*) int_val
            call dist_kind%add_node(int_val)
            call get_word(line,word,4)
            read(word,*) dbl_val
            call dist_kind%add_node(dbl_val)

           case('EX_LIST')
            call get_word(line,word,2)
            read(word,*) nexclude
            allocate(ex_list(nexclude,2))
            i=1
            do l=1,nexclude
             call get_word(line,word,2+i)
             read(word,*) ex_list(l,1)
             i=i+1
             call get_word(line,word,2+i)
             read(word,*) ex_list(l,2)
             i=i+1
            enddo

           case('T_INV')
            tinv=.true.

           case('MAX_EX')
            call get_word(line,word,2)
            read(word,*) nex_max

           case('MAX_CORR')
            call get_word(line,word,2)
            read(word,*) ncorr_max

           case('DUMP_FREQ')
            call get_word(line,word,2)
            read(word,*) dump_freq

           case('DUMP_MI')
            call get_word(line,word,2)
            read(word,*) s2print
            deallocate(print_si)
            allocate(print_si(s2print))
            do l=3,s2print+2
             call get_word(line,word,l)
             read(word,*) print_si(l-2)
            enddo

           case('&ACTIVE_SPACE')
            call parse_active_space
    
           case('&END')
            return

           case default
            write(6,*) 'Err: Keyword ',trim(word),' &
                        in subroutine parse_global not recognized'
            stop

          end select

         enddo

         write(*,*) 'Err: Unexpected end of file &
                     in subroutine parse_global'
         stop

        return
        end subroutine parse_global


        subroutine parse_active_space
        use control_variables
        use spins_dist_rs_class
        use parser_class
        use general_types_class
        implicit none
        character(len=:),allocatable   :: line,word
        integer                        :: l,size_old=0,nkinds
        logical                        :: eof=.false.
        type(sub_space), allocatable   :: active_tmp(:)


         if(allocated(active_space))then
          size_old=size(active_space)
         endif

         allocate(active_tmp(size_old+1))

         if(allocated(active_space))then
          do l=1,size_old
           active_tmp(l)%max_ex=active_space(l)%max_ex
           active_tmp(l)%kind=active_space(l)%kind
          enddo
         endif

         do 

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case(word)

           case('S')
            call get_word(line,word,2)
            read(word,*) nkinds
            allocate(active_tmp(size_old+1)%kind(nkinds))
            do l=1,nkinds
             call get_word(line,word,l+2)
             read(word,*) active_tmp(size_old+1)%kind(l)
            enddo
            
           case('MAX_EX')
             call get_word(line,word,2)
             read(word,*) active_tmp(size_old+1)%max_ex

           case('&END')

           if(allocated(active_space))then
            deallocate(active_space)
           endif

           allocate(active_space(size_old+1))

           do l=1,size(active_tmp)
            active_space(l)%max_ex=active_tmp(l)%max_ex
            active_space(l)%kind=active_tmp(l)%kind
           enddo

           return

          end select

         enddo

         write(6,*) 'Err: Unexpected end of file ,'&
                    'in subroutine parse_active_space'
         stop

        return
        end subroutine parse_active_space

        subroutine parse_system(spindy)
        use control_variables
        use spins_dist_rs_class
        use parser_class
        implicit none
        class(spins_group)             :: spindy
        type(SpinBath), allocatable    :: spin_bath(:)
        character(len=:),allocatable   :: line,word
        integer                        :: l
        logical                        :: eof=.false.

         spindy%bfield=0.0d0

         do 

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)


          select case(word)

           case('&DEF_SPINS')
            call parse_def_spins(spindy)

           case('B')
            call get_word(line,word,2)
            read(word,*) spindy%bfield(1)
            call get_word(line,word,3)
            read(word,*) spindy%bfield(2)
            call get_word(line,word,4)
            read(word,*) spindy%bfield(3)

           case('&CELL')
            call parse_cell(spindy)

           case('&END')
            return

           case default
           write(6,*) 'Err: Keyword ',trim(word),&
                      ' in subroutine parse_system not recognized'
           stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of file ,'&
                    'in subroutine parse_system'
         stop

        return
        end subroutine parse_system


        subroutine parse_spinham(spindy)
        use control_variables
        use spins_dist_rs_class
        use spinham_class
        use parser_class
        use lists_class
        implicit none
        class(spins_group)   :: spindy
        character(len=:),allocatable   :: line,word
        integer              :: l,i
        logical              :: eof=.false.
        type(Jiso_list)      :: J 
        type(DSItensor_list) :: DSI
        type(OSItensor_list) :: O
        type(D2Stensor_list) :: D2S
        type(Gtensor_list)   :: G
        type(Jiso)           :: J_tmp
        type(DSItensor)      :: DSI_tmp
        type(D2Stensor)      :: D2S_tmp
        type(OSItensor)      :: O_tmp
        type(Gtensor)        :: G_tmp
        class(*),pointer     :: bho

         call G%init()
         call J%init()
         call O%init()
         call DSI%init()
         call D2S%init()

         do 

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)


          select case(word)

           case('&DEF_G')
            spindy%SH%nG=spindy%SH%nG+1
            call get_word(line,word,2)
            read(word,*) G_tmp%kind
            call read_3x3mat(G_tmp%G)
            call G%add_node(G_tmp)

           case('&DEF_J')
            spindy%SH%nJ=spindy%SH%nJ+1
            call get_word(line,word,2)
            read(word,*) J_tmp%kind(1)
            call get_word(line,word,3)
            read(word,*) J_tmp%kind(2)
            call read_scalar(J_tmp%J)
            call J%add_node(J_tmp)

           case('&DEF_O')
            spindy%SH%nO=spindy%SH%nO+1
            call get_word(line,word,2)
            read(word,*) O_tmp%kind
            call get_word(line,word,3)
            read(word,*) O_tmp%k
            call read_otens(O_tmp)
            call O%add_node(O_tmp)

           case('&DEF_DSI')
            spindy%SH%nDSI=spindy%SH%nDSI+1
            call get_word(line,word,2)
            read(word,*) DSI_tmp%kind
            call read_3x3mat(DSI_tmp%D)
            call DSI%add_node(DSI_tmp)

           case('&DEF_D2S')
            spindy%SH%nD2S=spindy%SH%nD2S+1
            call get_word(line,word,2)
            read(word,*) D2S_tmp%kind(1)
            call get_word(line,word,3)
            read(word,*) D2S_tmp%kind(2)
            call read_3x3mat(D2S_tmp%D)
            call D2S%add_node(D2S_tmp)
            call D2S%add_node(D2S_tmp)

           case('SPINSPIN')
            spindy%SH%make_dipolar=.true.
            call get_word(line,word,2)
            read(word,*) spindy%SH%dipolar_thr

           case('EULER')
            call get_word(line,word,2)
            read(word,*) euler(1)
            call get_word(line,word,3)
            read(word,*) euler(2)
            call get_word(line,word,4)
            read(word,*) euler(3)

!            euler=euler*acos(-1.0d0)/180.0d0

           case('&END')

            if( spindy%SH%nG.gt.0) then
             allocate(spindy%SH%G(spindy%SH%nG))
             call G%reboot()             
             do i=1,spindy%SH%nG
              call G%rd_node(spindy%SH%G(i))
              call G%skip()
             enddo
            endif

            if( spindy%SH%nJ.gt.0) then
             allocate(spindy%SH%J(spindy%SH%nJ))
             call J%reboot()
             do i=1,spindy%SH%nJ
              call J%rd_node(spindy%SH%J(i))
              call J%skip()
             enddo
            endif

            if( spindy%SH%nO.gt.0) then
             allocate(spindy%SH%O(spindy%SH%nO))
             call O%reboot()
             do i=1,spindy%SH%nO
              call O%rd_node(spindy%SH%O(i))
              call O%skip()
             enddo
            endif

            if( spindy%SH%nDSI.gt.0) then
             allocate(spindy%SH%DSI(spindy%SH%nDSI))
             call DSI%reboot()
             do i=1,spindy%SH%nDSI
              call DSI%rd_node(spindy%SH%DSI(i))
              call DSI%skip()
             enddo
            endif

            if( spindy%SH%nD2S.gt.0) then
             allocate(spindy%SH%D2S(spindy%SH%nD2S))
             call D2S%reboot()
             do i=1,spindy%SH%nD2S              
              call D2S%rd_node(spindy%SH%D2S(i))
              call D2S%skip()
             enddo
            endif

            call G%delete()
            call J%delete()
            call O%delete()
            call DSI%delete()
            call D2S%delete()

            return

           case default

            write(6,*) 'Err: Keyword ',trim(word),&
                       ' in subroutine parse_spinham not recognized'
            stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of ',&
                    'file in subroutine parse_spinham'
         stop

        return
        end subroutine parse_spinham


        subroutine read_otens(O)
        use parser_class
        use spinham_class
        implicit none
        type(OSItensor)    :: O
        character(len=:),allocatable   :: line,word
        integer            :: l,i
        logical            :: eof=.false.
        double precision   :: B

         if(allocated(O%q))deallocate(O%q)
         if(allocated(O%B))deallocate(O%B)
         allocate(O%q(2*O%k+1))
         allocate(O%B(2*O%k+1))
         O%B=(0.0d0,0.0d0)

         do i=1,2*O%k+1

          call get_line(10,line,eof)
          if (eof) exit
          
          call get_word(line,word,1)
          if(word.eq.'&END') return
          read(word,*) O%q(i)
          call get_word(line,word,2)
          read(word,*) B
          O%B(i)=cmplx(B,0.0d0,8)

         enddo

         call get_line(10,line,eof)
         if (eof) stop
             
         call get_word(line,word,1)
         if(word.eq.'&END') return

        stop
        end subroutine read_otens

        subroutine read_3x3mat(G)
        use parser_class
        implicit none
        complex(8)         :: G(3,3)
        double precision   :: Gr(3,3)
        character(len=:),allocatable   :: line,word
        integer            :: l,i
        logical            :: eof=.false.

         do i=1,3

          call get_line(10,line,eof)
          if (eof) exit
          
          call get_word(line,word,1)
          read(word,*) Gr(i,1)
          call get_word(line,word,2)
          read(word,*) Gr(i,2)
          call get_word(line,word,3)
          read(word,*) Gr(i,3)

         enddo

         G=CMPLX(Gr,0.0d0,8)

         call get_line(10,line,eof)
         if (eof) stop
             
         call get_word(line,word,1)

         if(word.eq.'&END') return

        stop
        end subroutine read_3x3mat


        subroutine read_scalar(J)
        use parser_class
        implicit none
        complex(8)         :: J
        double precision   :: Jr
        character(len=:),allocatable   :: line,word
        integer            :: l
        logical            :: eof=.false.

         call get_line(10,line,eof)
         if (eof) stop
          
         call get_word(line,word,1)
         read(word,*) Jr
         J=CMPLX(Jr,0.0d0,8)

         call get_line(10,line,eof)
         if (eof) stop
             
         call get_word(line,word,1)

         if(word.eq.'&END') return

        stop
        end subroutine read_scalar


    
        subroutine parse_dynamics
        use control_variables
        use parser_class
        implicit none
        character(len=:),allocatable   :: line,word
        integer            :: l
        logical            :: eof=.false.

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('NSTEPS')
            call get_word(line,word,2)
            read(word,*) nsteps

           case('STEP_MIN')
            call get_word(line,word,2)
            read(word,*) step_min

           case('STEP_MULT')
            call get_word(line,word,2)
            read(word,*) step_nmult

           case('&END')

            step=step_min*(step_nmult+1)
            return

           case default
            write(*,*) 'Err: Keyword ',trim(word),&
                       ' in subroutine parse_dynamics not recognized'
            stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of ',&
                    'file in subroutine parse_dynamics'
         stop

        return
        end subroutine parse_dynamics


        subroutine parse_cell(spindy)
        use spins_dist_rs_class
        use parser_class
        implicit none
        class(spins_group) :: spindy
        character(len=:),allocatable   :: line,word
        integer                        :: l,repx,repy,repz,v,k,s,i,m
        integer                        :: nspins,i2,l2
        logical                        :: eof=.false.,check(3)
        double precision               :: cell(3,3),cell_inv(3,3),diff(3)
        double precision, allocatable  :: x(:,:)
        integer, allocatable           :: kind(:),coord3d(:,:)        

         repx=1
         repy=1
         repz=1
         nspins=0         

         do 

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)
        
          select case (word)

           case('A')            
            call get_word(line,word,2)
            read(word,*) cell(1,1)
            call get_word(line,word,3)
            read(word,*) cell(1,2)
            call get_word(line,word,4)
            read(word,*) cell(1,3)

           case('B')
            call get_word(line,word,2)
            read(word,*) cell(2,1)
            call get_word(line,word,3)
            read(word,*) cell(2,2)
            call get_word(line,word,4)
            read(word,*) cell(2,3)

           case('C')
            call get_word(line,word,2)
            read(word,*) cell(3,1)
            call get_word(line,word,3)
            read(word,*) cell(3,2)
            call get_word(line,word,4)
            read(word,*) cell(3,3)

           case('NREP')
            call get_word(line,word,2)
            read(word,*) repx
            call get_word(line,word,3)
            read(word,*) repy
            call get_word(line,word,4)
            read(word,*) repz

           case('&COORD')
            call def_coord(nspins,x,kind)

           case('&END')
            
            call spindy%read_cell(cell)
            spindy%ntot=repx*repy*repz
            spindy%nx=repx
            spindy%ny=repy
            spindy%nz=repz
            call spindy%do_supercell()

            spindy%nspins_pr=nspins
            spindy%nspins=nspins*repx*repy*repz
            allocate(spindy%x(spindy%nspins,3))
            allocate(spindy%kind(spindy%nspins))
                     
            spindy%x(1:nspins,:)=x(1:nspins,:)

            call spindy%cart2frac()

            if(spindy%ntot.gt.1) allocate(coord3d(nspins*spindy%ntot,3))

            v=1
            do k=1,repz
             do s=1,repy
              do i=1,repx
               do l=1,nspins
                spindy%x(v,1)=spindy%x(l,1)+i-1
                spindy%x(v,2)=spindy%x(l,2)+s-1
                spindy%x(v,3)=spindy%x(l,3)+k-1
                spindy%kind(v)=kind(l)
                if(spindy%ntot.gt.1) coord3d(v,1)=i
                if(spindy%ntot.gt.1) coord3d(v,2)=s
                if(spindy%ntot.gt.1) coord3d(v,3)=k
                v=v+1
               enddo
              enddo
             enddo
            enddo

            if(spindy%ntot .gt. 1)then

             allocate(spindy%tr_map(spindy%nspins,3))

             do i=1,spindy%ntot
              do l=1,spindy%nspins_pr

               k=(i-1)*spindy%nspins_pr+l
               spindy%tr_map(k,:)=k

               do i2=1,spindy%ntot
               do l2=1,spindy%nspins_pr

                m=(i2-1)*spindy%nspins_pr+l2

                if(l2.ne.l) cycle

                check=.false.

                if (spindy%nx.eq.1 .or.  coord3d(k,1)+1.eq.coord3d(m,1)) check(1)=.true.
                if (spindy%nx.ne.1 .and. coord3d(m,1).eq.1 .and. coord3d(k,1).eq.spindy%nx) check(1)=.true.
                if (spindy%ny.eq.1 .or.  coord3d(k,2).eq.coord3d(m,2)) check(2)=.true.
                if (spindy%nz.eq.1 .or.  coord3d(k,3).eq.coord3d(m,3)) check(3)=.true.

                if (all(check)) spindy%tr_map(k,1)=m


                check=.false.

                if (spindy%ny.eq.1 .or.  coord3d(k,2)+1.eq.coord3d(m,2)) check(2)=.true.
                if (spindy%ny.ne.1 .and. coord3d(m,2).eq.2 .and. coord3d(k,2).eq.spindy%ny) check(2)=.true.
                if (spindy%nx.eq.1 .or.  coord3d(k,1).eq.coord3d(m,1)) check(1)=.true.
                if (spindy%nz.eq.1 .or.  coord3d(k,3).eq.coord3d(m,3)) check(3)=.true.

                if (all(check)) spindy%tr_map(k,2)=m


                check=.false.

                if (spindy%nz.eq.1 .or.  coord3d(k,3)+1.eq.coord3d(m,3)) check(3)=.true.
                if (spindy%nz.ne.1 .and. coord3d(m,3).eq.1 .and. coord3d(k,3).eq.spindy%nz) check(3)=.true.
                if (spindy%ny.eq.1 .or.  coord3d(k,2).eq.coord3d(m,2)) check(2)=.true.
                if (spindy%nx.eq.1 .or.  coord3d(k,1).eq.coord3d(m,1)) check(1)=.true.

                if (all(check)) spindy%tr_map(k,3)=m

               enddo
               enddo

              enddo              
             enddo

             deallocate(coord3d)

            endif

            
            call spindy%frac2cart()

            return

           case default
            write(*,*) 'Err: Keyword ',trim(word),&
                       ' in subroutine parse_cell not recognized'
            stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of ',&
                    'file in subroutine parse_cell'
         stop

        return
        end subroutine parse_cell


        subroutine parse_def_spins(spindy)
        use spins_dist_rs_class
        use lists_class
        use parser_class
        implicit none
        class(spins_group) :: spindy
        character(len=:),allocatable   :: line,word
        integer            :: l,val_int,i,j
        double precision   :: val_dbl
        logical            :: eof=.false.
        type(list)         :: bohr,spin,kind
        class(*), pointer  :: bho

         call bohr%init()
         call spin%init()
         call kind%init()

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('S')

            call get_word(line,word,2)
            read(word,*) val_int
            call kind%add_node(val_int)

            call get_word(line,word,3)
            read(word,*) val_dbl
            call spin%add_node(val_dbl)

            call get_word(line,word,4)
            read(word,*) val_dbl
            call bohr%add_node(val_dbl)

            spindy%nkinds=spindy%nkinds+1

           case('&END')

            allocate(spindy%bohr_mag(spindy%nkinds))
            allocate(spindy%spin(spindy%nkinds))
            call bohr%reboot()
            call spin%reboot()
            call kind%reboot()

            do i=1,spindy%nkinds

             call kind%rd_val(j)
             call bohr%rd_val(spindy%bohr_mag(j))
             call spin%rd_val(spindy%spin(j))

             call bohr%skip()
             call spin%skip()
             call kind%skip()

            enddo

            call bohr%delete()
            call spin%delete()
            call kind%delete()

            return

           case default
            write(*,*) 'Err: Keyword ',trim(word),&
                       ' in subroutine parse_def_spins not recognized'
            stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of ,'&
                    'file in subroutine parse_def_spins'
         stop

        return
        end subroutine parse_def_spins


        subroutine def_coord(nspins,x,kinds)
        use lists_class
        use parser_class
        implicit none
        character(len=:),allocatable   :: line,word
        integer            :: l,val_int,i,j,nspins
        double precision   :: val_dbl
        double precision, allocatable :: x(:,:)
        integer, allocatable :: kinds(:)
        logical            :: eof=.false.
        type(list)         :: coord(3),kind
        class(*), pointer  :: bho
        

         call coord(1)%init()
         call coord(2)%init()
         call coord(3)%init()
         call kind%init()

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('S')

            call get_word(line,word,2)
            read(word,*) val_int
            call kind%add_node(val_int)

            call get_word(line,word,3)
            read(word,*) val_dbl
            call coord(1)%add_node(val_dbl)

            call get_word(line,word,4)
            read(word,*) val_dbl
            call coord(2)%add_node(val_dbl)

            call get_word(line,word,5)
            read(word,*) val_dbl
            call coord(3)%add_node(val_dbl)

            nspins=nspins+1

           case('&END')

            allocate(x(nspins,3))
            allocate(kinds(nspins))

            call coord(1)%reboot()
            call coord(2)%reboot()
            call coord(3)%reboot()
            call kind%reboot()

            do i=1,nspins

             select type (bho=>kind%node%key)
              type is (integer)
               kinds(i)=bho
             end select

             select type (bho=>coord(1)%node%key)
              type is (double precision)
               x(i,1)=bho
             end select

             select type (bho=>coord(2)%node%key)
              type is (double precision)
               x(i,2)=bho
             end select

             select type (bho=>coord(3)%node%key)
              type is (double precision)
               x(i,3)=bho
             end select

             call coord(1)%skip()
             call coord(2)%skip()
             call coord(3)%skip()
             call kind%skip()

            enddo

            call coord(1)%delete()
            call coord(2)%delete()
            call coord(3)%delete()
            call kind%delete()

            return

           case default

            write(*,*) 'Err: Keyword ',trim(word),&
                       ' in subroutine def_coord not recognized'
            stop

          end select

         enddo

         write(6,*) 'Err: Unexpected end of ,'&
                    'file in subroutine def_coord'
         stop

        return
        end subroutine def_coord



        subroutine parse_pulse(spindy)
        use parser_class
        use spins_dist_rs_class
        use hilbert_dist_class
        use pulses_class
        use control_variables
        implicit none
        character(len=:),allocatable     :: line,word
        logical                          :: eof=.false.
        class(spins_hilbert)             :: spindy
        type(general_pulse)              :: rot
        type(general_pulse_list)         :: pulse_list
        integer                          :: k
        
         call pulse_list%init()

         if(allocated(pulse)) deallocate(pulse)

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)


          select case (word)

           case('&ROT') 

            call parse_rot(rot)
            call pulse_list%add_node(rot)

           case('&END')

            call pulse_list%reboot()
           
            allocate(pulse(pulse_list%nelem))
            do k=1,pulse_list%nelem

             call pulse_list%rd_node(pulse(k))
             call pulse_list%skip()

            enddo
            
            call pulse_list%delete()

            return
            
          end select         

         enddo

        return
        end subroutine parse_pulse

        subroutine parse_rot(rot)
        use parser_class
        use spins_dist_rs_class
        use hilbert_dist_class
        use pulses_class
        implicit none
        character(len=:),allocatable   :: line,word
        logical                        :: eof=.false.
        integer                        :: i
        type(general_pulse)            :: rot
        
         rot%sx=.true.
         rot%dx=.true.
         rot%weight=1.0d0
         rot%spin=1
         rot%beta=acos(-1.0d0)/2.0d0
         rot%n(1)=0.0d0
         rot%n(2)=1.0d0
         rot%n(3)=0.0d0

         do

          call get_line(10,line,eof)
          if (eof) exit

          call get_word(line,word,1)
          call To_upper(word)

          select case (word)

           case('DX')
            call get_word(line,word,2)
            call To_upper(word)
            read(word,*) rot%dx

           case('SX')
            call get_word(line,word,2)
            call To_upper(word)
            read(word,*) rot%sx

           case('WEIGHT')
            call get_word(line,word,2)
            read(word,*) rot%weight

           case('S')

             call get_word(line,word,2)
             read(word,*) rot%spin

             call get_word(line,word,3)
             read(word,*) rot%beta
             rot%beta=rot%beta*acos(-1.0d0)/180.d0

             call get_word(line,word,4)
             read(word,*) rot%n(1)

             call get_word(line,word,5)
             read(word,*) rot%n(2)

             call get_word(line,word,6)
             read(word,*) rot%n(3)

           case('&END')

            return

          end select

         enddo
         
        return
        end subroutine parse_rot

        end module spindy_parser
