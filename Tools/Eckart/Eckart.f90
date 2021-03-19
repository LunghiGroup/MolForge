        program eckart
        use proj_disp_class
        implicit none
        integer                         :: i,j,l,s,v,nats,nframes

        type(molecule)                  :: mol
        double precision,allocatable    :: geo(:,:),mass(:),geo0(:,:)
        character(len=5), allocatable   :: label(:)

        character(len=100)              :: traj_filename,filename,option
        character(len=100)              :: int_filename='traj_int.xyz'
        character(len=100)              :: rot_filename='traj_rot.xyz'
        character(len=100)              :: tr_filename='traj_tr.xyz'
        
        if(iargc().eq.0)then
         write(*,*) '-ref_geo     : xyz file containing the reference geometry'
         write(*,*) '-traj        : xyz trajectory file containing distorted frames of ref_geo'
         write(*,*) '-frames      : Number of distorted frames in traj'
         write(*,*) '-int_out     : xyz filename where to print internal dynamics (optional)'
         write(*,*) '-rot_out     : xyz filename where to print rotational dynamics (optional)'
         write(*,*) '-tr_out      : xyz filename where to print translational dynamics (optional)'
         stop
        endif

        do l=1,iargc()

         call  getarg(l,option)

         select case (option)

         case('-ref_geo')
          call getarg(l+1,option)
          read(option,*) filename

          open(unit=13,file=filename)

          read(13,*) nats
          read(13,*) 

          allocate(geo(nats,3))
          allocate(geo0(nats,3))
          allocate(mass(nats))
          allocate(label(nats))

          do i=1,nats
           read(13,*) label(i),geo0(i,:)          
          ! set masses
           if(trim(label(i)).eq.'C') mass(i)=12.0
           if(trim(label(i)).eq.'H') mass(i)=1.0
           if(trim(label(i)).eq.'Se') mass(i)=79.9
          enddo
          
          call mol%def_mol(geo0,mass)
          close(13)

         case('-traj')
          call getarg(l+1,option)
          read(option,*) traj_filename

         case('-frames')
          call getarg(l+1,option)
          read(option,*) nframes

         case('-int_out')
          call getarg(l+1,option)
          read(option,*) int_filename

         case('-rot_out')
          call getarg(l+1,option)
          read(option,*) rot_filename

         case('-tr_out')
          call getarg(l+1,option)
          read(option,*) tr_filename

         end select

        enddo
            
        open(13,file=traj_filename)
        open(14,file=int_filename)         
        open(15,file=rot_filename)         
        open(16,file=tr_filename)         

        do v=1,nframes

         read(13,*) nats
         read(13,*) 
         do i=1,nats
          read(13,*) label(i),geo(i,1),geo(i,2),geo(i,3)
         enddo

         call mol%def_mol_dist(geo)
         call mol%proj_disp()
         
         write(14,*) nats
         writE(14,*) 
         do i=1,nats
          write(14,*) label(i),(mol%cart_int(i,s)+mol%cart_eq(i,s),s=1,3)
         enddo   

         write(15,*) nats
         writE(15,*) 
         do i=1,nats
          write(15,*) label(i),(mol%cart_rot(i,s),s=1,3)
         enddo       

         write(16,*) nats
         writE(16,*) 
         do i=1,nats
          write(16,*) label(i),(mol%tr(s)+mol%cart_eq(i,s),s=1,3)
         enddo       

        enddo

        close(13)
        close(14)
        close(15)
        close(16)
                 
        return
        end program eckart
