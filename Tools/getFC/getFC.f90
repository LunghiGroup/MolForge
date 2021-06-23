        program generate_disps
        use mpi
        use atoms_class
        implicit none
        type(atoms_group)              :: sys
        integer                        :: nats,nx,ny,nz,err,l1,s1,icell1,i2,s2,k2
        integer                        :: i,s,k,j,l,v,icell,k1,v0,v1,v2,kx,ky,kz
        double precision               :: step,cell(3,3),sum
        double precision               :: mux,muy,muz,dDip(3)
        double precision, allocatable  :: x0(:,:),x(:,:),fx(:),fy(:),fz(:)
        character(len=100)             :: forces_file_2,dipole_file
        character(len=100)             :: arg,cell_file,inp_file
        character(len=10), allocatable :: label(:)
        logical                        :: gen_disps_2=.false.
        logical                        :: read_forces_2=.false.
        logical                        :: read_dipole=.false.,new_type,isol=.false.
        
        call MPI_INIT(err)


         if(iargc().eq.0)then
          write(*,*) 'Usage:'
          write(*,*) '-xyz           = Unit cell xyz'
          write(*,*) '-cell          = Unic cell 3x3 matrx '
          write(*,*) '-n             = Supercell nx ny nz'
          write(*,*) '-step          = Diff. step in Angstrom'
          write(*,*) '-gen_disps     = '
          write(*,*) '-read_forces   = '
          write(*,*) '-read_dipole   = '
          stop
         endif

         do i=1,iargc()

          call getarg(i,arg)

          if (arg.eq.'-xyz')then
           call getarg(i+1,arg)
           inp_file=arg
          endif

          if (arg.eq.'-cell')then
           call getarg(i+1,arg)
           cell_file=arg
          endif

          if (arg.eq.'-step')then
           call getarg(i+1,arg)
           read(arg,*) step
           step=step*1.889725989d0
          endif

          if (arg.eq.'-n')then
           call getarg(i+1,arg)
           read(arg,*) nx
           call getarg(i+2,arg)
           read(arg,*) ny
           call getarg(i+3,arg)
           read(arg,*) nz
          endif

          if (arg.eq.'-isol') isol=.true.

          if (arg.eq.'-gen_disps') gen_disps_2=.true.

          if (arg.eq.'-read_forces') then
           read_forces_2=.true.
           call getarg(i+1,arg)
           forces_file_2=arg
          endif

          if (arg.eq.'-read_dipoles') then
           read_dipole=.true.
           call getarg(i+1,arg)
           dipole_file=arg
          endif

         enddo

         open(11,file=cell_file)

         read(11,*) cell(1,:)
         read(11,*) cell(2,:)
         read(11,*) cell(3,:)
         close(11)

         call sys%read_cell(cell)
         sys%nx=nx
         sys%ny=ny
         sys%nz=nz
         sys%ntot=nz*ny*nx
         call sys%do_supercell()

         call sys%read_xyz(10,inp_file)

         if(gen_disps_2)then

          j=1
          do i=1,sys%nats
           do s=1,3
            do k=1,-1,-2

             if(j.lt.10) write(arg,"(I1)") j
             if(j.lt.100 .and. j.gt.9 ) write(arg,"(I2)") j
             if(j.lt.1000 .and. j.gt.99) write(arg,"(I3)") j
             if(j.lt.10000 .and. j.gt.999) write(arg,"(I4)") j
             open(10,file=trim(arg)//'_2nd_'//trim(inp_file))

             x=x0
             x(i,s)=x0(i,s)+k*step*0.529177249d0

             write(10,*) sys%nats*sys%ntot
             write(10,*)  
             do l=1,sys%nats*sys%ntot
              write(10,*) label(l),x(l,:)     
             enddo

             close(10)
         
             j=j+1
            enddo
           enddo
          enddo        

         endif

         if(read_forces_2) call get_fc2(forces_file_2,sys,step,isol)

         if(read_dipole)then

          open(14,file=dipole_file)
          open(15,file='dDipole.dat')

          do l=1,sys%nats
           do s=1,3
            dDip=0.0d0
            do k=-1,1,2

              read(14,*) mux,muy,muz
    
              dDip(1)=dDip(1)+k*mux/(2.0d0*step) 
              dDip(2)=dDip(2)+k*muy/(2.0d0*step) 
              dDip(3)=dDip(3)+k*muz/(2.0d0*step) 

              write(15,*) dDip(:)

            enddo ! k
           enddo ! s 
          enddo ! l

          close(14)
          close(15)

         endif

         call MPI_FINALIZE(err)

        return
        end program generate_disps


        subroutine get_fc2(forces_file_2,sys,step,isol)
        use lists_class
        use atoms_class
        use lapack_inverse
        use lapack_diag_simm
        implicit none
        type(atoms_group)                    :: sys
        integer                              :: l,s,k,icell,i,v,j,t
        integer                              :: kz,ky,kx
        double precision                     :: step,sum,sumx,sumy,sumz,inertia(3,3),coeff,mass1
        double precision                     :: cmass_eq(3),eig(3),vec(3),inertia_inv(3,3),M(3,3,3)
        double precision, allocatable        :: fz(:),fx(:),fy(:),x0(:,:),AA(:,:),H(:,:),Heig(:)
        double precision, allocatable        :: Dmat(:,:),Amat(:,:),Bmat(:,:),mass(:)
        character(len=10)                    :: label
        character(len=100)                   :: forces_file_2        
        logical                              :: isol
        
         write(*,*) 'Reading Forces from: ',trim(forces_file_2)

         allocate(x0(sys%ntot*sys%nats,3))

         call sys%cart2frac()

         k=1
         do kz=1,sys%nz
          do ky=1,sys%ny
           do kx=1,sys%nx
            do i=1,sys%nats
             j=(k-1)*sys%nats+i
             x0(j,1)=sys%x(i,1)+kx
             x0(j,2)=sys%x(i,2)+ky
             x0(j,3)=sys%x(i,3)+kz            
             x0(j,:)=matmul(sys%J,x0(j,:))
            enddo
            k=k+1
           enddo
          enddo
         enddo

         call sys%frac2cart()

         open(14,file=forces_file_2)
         allocate(fx(sys%ntot*sys%nats))                   
         allocate(fy(sys%ntot*sys%nats))                   
         allocate(fz(sys%ntot*sys%nats))                   
         call sys%init_fcs2()          

         do l=1,sys%nats
          do s=1,3
           do k=1,-1,-2

            read(14,*) 
            read(14,*) 
            do i=1,sys%ntot*sys%nats
             read(14,*) fx(i),fy(i),fz(i)
            enddo

            do icell=1,sys%ntot
             do i=1,sys%nats
              v=(icell-1)*sys%nats+i
              sys%fcs2(icell,l,s,i,1)=sys%fcs2(icell,l,s,i,1)-k*fx(v)/(2.0d0*step)
              sys%fcs2(icell,l,s,i,2)=sys%fcs2(icell,l,s,i,2)-k*fy(v)/(2.0d0*step)
              sys%fcs2(icell,l,s,i,3)=sys%fcs2(icell,l,s,i,3)-k*fz(v)/(2.0d0*step)
             enddo
            enddo

           enddo ! k
          enddo ! s 
         enddo ! l

        if(.not. isol)then

         ! Symmetrise by average of out-of-diag elements 

!         do s=1,3
!          do l=1,sys%nats
!           do k=1,3
!            do i=1,sys%nats
!             v=(l-1)*3+s
!             j=(i-1)*3+k
!             if(j.lt.v)cycle 
!             sum=sys%fcs2(1,l,s,i,k)+sys%fcs2(1,i,k,l,s)
!             sys%fcs2(1,l,s,i,k)=sum/2.0d0
!             sys%fcs2(1,i,k,l,s)=sum/2.0d0
!            enddo
!           enddo
!          enddo
!         enddo

         ! Acoustic Sum Rule for solids

         do s=1,3
          do k=1,3
           do i=1,sys%nats
            sum=0.0d0
            sys%fcs2(1,i,s,i,k)=0.0d0
            do j=1,sys%nats
             do icell=1,sys%ntot
              sum=sum+sys%fcs2(icell,i,s,j,k)
             enddo
            enddo

            sys%fcs2(1,i,s,i,k)=-sum

           enddo
          enddo
         enddo

         endif
         
         call sys%write_restart_file('FC2') 
         close(14)

        return
        end subroutine get_fc2

