        program dipole2charges
        implicit none
        integer                        :: nframes,nats,dimA,dimB,lwork,inf
        integer                        :: i,j,v
        double precision               :: cm(3)
        double precision, allocatable  :: B(:),A(:,:),work(:),x(:,:)
        character(len=100)             :: word,dipoles_filename,xyz_filename
        character(len=10)              :: label,origin_type='zero'
        

         if( iargc().eq.0)then
          write(*,*) 'fitcharges Usage:'
          write(*,*) '-nframes      : Total number of frames in the xyz file'
          write(*,*) '-nats         : Number of atoms in each frame'
          write(*,*) '-xyz          : Name of the xyz file to read from'
          write(*,*) '-dipoles      : Name of the file with dipole vectors'
          write(*,*) '-origin       : "com" form center of mass or "zero" for r=(0,0,0)'
          stop
         endif

         do i=1,iargc()
          call getarg(i,word)
          
          select case (trim(word))

             case ('-xyz')
                 call getarg(i+1,word)
                 read(word,*) xyz_filename

             case ('-dipoles')
                 call getarg(i+1,word)
                 read(word,*) dipoles_filename

             case ('-frames')
                 call getarg(i+1,word)
                 read(word,*) nframes

             case ('-nats')
                 call getarg(i+1,word)
                 read(word,*) nats

             case ('-origin')
                 call getarg(i+1,word)
                 read(word,*) origin_type

          end select

         enddo

         dimB=3*nframes
         dimA=nats
         allocate(A(dimB,dimA))
         allocate(B(dimB))
         A=0.0d0
         B=0.0d0
         cm=0.0d0

         open(13,file=xyz_filename)
         open(14,file=dipoles_filename)

         allocate(x(nats,3))

         v=1
         do i=1,nframes
          read(14,*) B(v),B(v+1),B(v+2)
          read(13,*)
          read(13,*)
          cm=0.0d0
          do j=1,nats
           read(13,*) label,x(j,:)
           cm=cm+x(j,:)
          enddo
          if(trim(origin_type).eq.'com')then
           do j=1,nats
            A(v,j)=x(j,1)-cm(1)
            A(v+1,j)=x(j,2)-cm(2)
            A(v+2,j)=x(j,3)-cm(3)
           enddo
          endif
          if(trim(origin_type).eq.'zero')then
           do j=1,nats
            A(v,j)=x(j,1)
            A(v+1,j)=x(j,2)
            A(v+2,j)=x(j,3)
           enddo
          endif
          v=v+3
         enddo

         lwork=dimB+64*dimB+1000
         allocate(work(lwork))
         call dgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
         deallocate(work)
         if(inf.ne.0)then
          write(*,*) 'zgels failed',inf
          stop
         else
          do i=1,nats
           write(*,*) B(i)
          enddo
         endif

         close(13)
         close(14)

        return
        end program dipole2charges


