        program dipole2charges
        implicit none
        type dataset
         double precision, allocatable :: x(:,:)
         double precision              :: cm(3)=0.0d0
         double precision              :: dipole(3)=0.0d0
        end type dataset
        integer                        :: tr_frames,te_frames,dimA,dimB,lwork,inf
        integer                        :: i,j,v,nats
        double precision, allocatable  :: B(:),A(:,:),work(:)
        double precision               :: L2coeff=1.0d0,tot_charge=0.0d0,dipole(3)
        type(dataset), allocatable     :: trdata(:),tedata(:)
        character(len=100)             :: word,tr_dipoles,tr_xyz,te_dipoles,te_xyz,inp_filename
        character(len=10)              :: label,origin_type='zero'
        logical                        :: skip_fit=.false.,regularization=.false.,normalization=.false.
        

         if( iargc().eq.0)then
          write(*,*) 'fitcharges Usage:'
          write(*,*) '-inp          : Total number of frames in the xyz file'
          write(*,*) '-origin       : "com" form center of mass or "zero" for r=(0,0,0)'
          write(*,*) '-tot_charge   : Constraint total charge'
          write(*,*) '-compress     : L2 normalization parameter'
          stop
         endif

         do i=1,iargc()
          call getarg(i,word)
          
          select case (trim(word))

             case ('-inp')
                 call getarg(i+1,word)
                 read(word,*) inp_filename

             case ('-origin')
                 call getarg(i+1,word)
                 read(word,*) origin_type

             case ('-tot_charge')
                 call getarg(i+1,word)
                 read(word,*) tot_charge
                 normalization=.true.

             case ('-compress')
                 call getarg(i+1,word)
                 read(word,*) L2coeff
                 regularization=.true.

          end select

         enddo

         ! Read training and test set

         open(11,file=inp_filename)
         read(11,*) tr_frames,tr_xyz,tr_dipoles
         allocate(trdata(tr_frames))
         read(11,*) te_frames,te_xyz,te_dipoles
         allocate(tedata(te_frames))
         close(11)

         open(13,file=tr_xyz)
         open(14,file=tr_dipoles)

         v=1
         do i=1,tr_frames
          
          read(14,*) trdata(i)%dipole

          read(13,*) nats
          read(13,*)
          allocate(trdata(i)%x(nats,3))

          do j=1,nats
           read(13,*) label,trdata(i)%x(j,:)
           trdata(i)%cm=trdata(i)%cm+trdata(i)%x(j,:)
          enddo

         enddo

         close(13)
         close(14)

         open(13,file=te_xyz)
         open(14,file=te_dipoles)

         v=1
         do i=1,te_frames
          
          read(14,*) tedata(i)%dipole

          read(13,*) nats
          read(13,*)

          allocate(tedata(i)%x(nats,3))

          do j=1,nats
           read(13,*) label,tedata(i)%x(j,:)
           tedata(i)%cm=tedata(i)%cm+tedata(i)%x(j,:)
          enddo

         enddo

         close(13)
         close(14)

         ! Build and solve system B=Ax , B=list of dipoles , A=coordinates, x=charges

         dimB=3*tr_frames
         dimA=nats
         if(regularization) dimB=dimB+nats
         if(normalization) dimB=dimB+1
         allocate(A(dimB,dimA))
         allocate(B(dimB))
         A=0.0d0
         B=0.0d0

         v=1
         do i=1,tr_frames
          B(v)=trdata(i)%dipole(1)
          B(v+1)=trdata(i)%dipole(2)
          B(v+2)=trdata(i)%dipole(3)
          if(trim(origin_type).eq.'com')then
           do j=1,nats
            A(v,j)=trdata(i)%x(j,1)-trdata(i)%cm(1)
            A(v+1,j)=trdata(i)%x(j,2)-trdata(i)%cm(2)
            A(v+2,j)=trdata(i)%x(j,3)-trdata(i)%cm(3)
           enddo
          endif
          if(trim(origin_type).eq.'zero')then
           do j=1,nats
            A(v,j)=trdata(i)%x(j,1)
            A(v+1,j)=trdata(i)%x(j,2)
            A(v+2,j)=trdata(i)%x(j,3)
           enddo
          endif
          v=v+3
         enddo
         if(regularization)then
          do j=1,nats
           A(v,j)=1.0d0
           B(v)=L2coeff
           v=v+1
          enddo
         endif
         if(normalization)then
          do j=1,nats
           A(v,j)=1.0d3
          enddo
          B(v)=tot_charge*1.0d3
          v=v+1
         endif
         
         open(15,file='charges.dat')

         lwork=dimB+64*dimB+1000
         allocate(work(lwork))
         call dgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
         deallocate(work)
         if(inf.ne.0)then
          write(*,*) 'zgels failed',inf
          stop
         else
          do i=1,nats
           write(15,*) B(i)
          enddo
         endif

         close(15)

         ! Predict dipoles for both 

         open(11,file='dipoles_tr.dat')
         do i=1,tr_frames
          if(trim(origin_type).eq.'com')then
           dipole=0.0d0
           do j=1,nats 
            dipole(1)=dipole(1)+(trdata(i)%x(j,1)-trdata(i)%cm(1))*B(j)
            dipole(2)=dipole(2)+(trdata(i)%x(j,2)-trdata(i)%cm(2))*B(j)
            dipole(3)=dipole(3)+(trdata(i)%x(j,3)-trdata(i)%cm(3))*B(j)
           enddo
           write(11,*) i,trdata(i)%dipole,dipole
          endif
          if(trim(origin_type).eq.'zero')then
           dipole=0.0d0
           do j=1,nats 
            dipole(1)=dipole(1)+trdata(i)%x(j,1)*B(j)
            dipole(2)=dipole(2)+trdata(i)%x(j,2)*B(j)
            dipole(3)=dipole(3)+trdata(i)%x(j,3)*B(j)
           enddo
           write(11,*) i,trdata(i)%dipole,dipole
          endif
         enddo
         close(11)

         open(11,file='dipoles_te.dat')
         do i=1,te_frames
          if(trim(origin_type).eq.'com')then
           dipole=0.0d0
           do j=1,nats 
            dipole(1)=dipole(1)+(tedata(i)%x(j,1)-tedata(i)%cm(1))*B(j)
            dipole(2)=dipole(2)+(tedata(i)%x(j,2)-tedata(i)%cm(2))*B(j)
            dipole(3)=dipole(3)+(tedata(i)%x(j,3)-tedata(i)%cm(3))*B(j)
           enddo
           write(11,*) i,tedata(i)%dipole,dipole
          endif
          if(trim(origin_type).eq.'zero')then
           dipole=0.0d0
           do j=1,nats 
            dipole(1)=dipole(1)+tedata(i)%x(j,1)*B(j)
            dipole(2)=dipole(2)+tedata(i)%x(j,2)*B(j)
            dipole(3)=dipole(3)+tedata(i)%x(j,3)*B(j)
           enddo
           write(11,*) i,tedata(i)%dipole,dipole
          endif
         enddo
         close(11)
                 

        return
        end program dipole2charges


