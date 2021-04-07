        program make_cells
        use atoms_class
        use variables
        implicit none
        DOUBLE PRECISION              :: ZERO=1.0E-12,mat_h(3,3),mat_hinv(3,3),aaa,bbb,pi,deg2rad,lat
        double precision, pointer     :: pos(:)
        double precision              :: mass
        integer                       :: s,k,i,ialloc,v,l,shift,nkind
        integer,allocatable           :: label_num(:)
        character*10, allocatable     :: label(:),kind(:),old_kind(:)
        character(len=2)              :: lab
        character*100                 :: word,unit_cell
        logical                       :: lammps_file_gen=.false. 

         if( iargc().eq.0)then
          write(*,*) 'SuperCellGen Usage:'               
          write(*,*) '-unit    : Name of the xyz file w unit cell'               
          write(*,*) '-box     : a,b,c,alpha,beta,gamma:'               
          write(*,*) '-n       : nx ny nz'               
          write(*,*) '-lammps  : Generate lammps data file'               
          stop
         endif

         do l=1,100
          call getarg(l,word)

          if(trim(word).eq.'-unit')then               
           call getarg(l+1,word)
           read(word,*) unit_cell
          endif

          if(trim(word).eq.'-n')then               
           call getarg(l+1,word)
           read(word,*) nx
           call getarg(l+2,word)
           read(word,*) ny
           call getarg(l+3,word)
           read(word,*) nz
          endif

          if(trim(word).eq.'-box')then
           
           call getarg(l+1,word)
           read(word,*) cell(1)
           call getarg(l+2,word)
           read(word,*) cell(2)
           call getarg(l+3,word)
           read(word,*) cell(3)
           call getarg(l+4,word)
           read(word,*) cell(4)
           call getarg(l+5,word)
           read(word,*) cell(5)
           call getarg(l+6,word)
           read(word,*) cell(6)
 
          endif

          if(trim(word).eq.'-lammps') lammps_file_gen=.true.

         enddo

         open(9,file=trim(unit_cell))
         open(12,file='coord.xyz')


         read(9,*) nmolxcell
         read(9,*) 

         allocate(pos0(nmolxcell,3),stat=ialloc)
         allocate(label(nmolxcell),stat=ialloc)
         allocate(label_num(nmolxcell),stat=ialloc)
         allocate(kind(1),stat=ialloc)
         nkind=1

         do l=1,nmolxcell
          read(9,*) label(l),pos0(l,1),pos0(l,2),pos0(l,3)
          if (l.gt.1)then
           if( .not.any(label(l).eq.label(1:l-1) ) )then
            nkind=nkind+1
            allocate(old_kind(nkind-1),stat=ialloc)
            old_kind=kind
            deallocate(kind)
            allocate(kind(nkind),stat=ialloc)
            kind(1:nkind-1)=old_kind(1:nkind-1)
            deallocate(old_kind)
            kind(nkind)=trim(label(l)) 
            endif
           else
            kind(1)=trim(label(l)) 
          endif
         enddo 

         do i=1,nmolxcell
          do k=1,nkind
          if(trim(label(i)).eq.trim(kind(k)))then
            label_num(i)=k
          endif
          enddo
         enddo

         nmol=(nx*ny*nz)*nmolxcell
         allocate(pos(3),stat=ialloc)
   
         pos=0.0d0

         call xyztoxray(mat_hinv)
         call xraytoxyz(mat_h)

         do i=1,nmolxcell
          pos0(i,1:3)=MATMUL(mat_hinv,pos0(i,1:3))
         enddo

         write(12,*) nmol
         write(12,*) 

         pi=acos(-1.0d0)
         deg2rad=pi/180.0d0
         cell(1)=cell(1)*nx
         cell(2)=cell(2)*ny
         cell(3)=cell(3)*nz
         cell(4)=cell(4)*deg2rad
         cell(5)=cell(5)*deg2rad
         cell(6)=cell(6)*deg2rad
 
         aaa=cell(3)*cos(cell(5))
         bbb=(cell(2)*cell(3)*cos(cell(4))-cell(2)*cos(cell(6))*cell(3)*cos(cell(5)))/(cell(2)*sin(cell(6)))

         v=1
         do k=1,nz
          do s=1,ny
           do i=1,nx
            do l=1,nmolxcell
             pos(1)=pos0(l,1)+i-1
             pos(2)=pos0(l,2)+s-1
             pos(3)=pos0(l,3)+k-1
             pos(1:3)=MATMUL(mat_h,pos(1:3))
             write(12,*) label(l),pos(1),pos(2),pos(3)
             v=v+1
            enddo
           enddo
          enddo
         enddo

         close(12)

         if(lammps_file_gen)then

         open(11,file='data')

          write(11,*) 'LAMMPS data file produced by map_phonon'
          write(11,*) 
          write(11,*) nmol,'atoms'
          write(11,*) nkind,'atom types'
          write(11,*)
          write(11,*) 0.0,cell(1),'xlo xhi'
          write(11,*) 0.0,cell(2)*sin(cell(6)),'ylo yhi'
          write(11,*) 0.0,sqrt(cell(3)**2-aaa**2-bbb**2),'zlo zhi'
          write(11,*) cell(2)*cos(cell(6)),cell(3)*cos(cell(5)), &
                    (cell(2)*cell(3)*cos(cell(4))-cell(2)*cos(cell(6))*cell(3)*cos(cell(5)))/(cell(2)*sin(cell(6))) &
                    ,'xy xz yz' 
          write(11,*) 
          write(11,*) 'Masses'
          write(11,*)
          do i=1,nkind
           lab=trim(kind(i))
           call get_mass(lab,mass)
           write(11,*)  i,mass!,trim(kind(i))
          enddo
          write(11,*)
          write(11,*) 'Atoms'
          write(11,*)

          v=1
          do k=1,nz
           do s=1,ny
            do i=1,nx
             do l=1,nmolxcell
              pos(1)=pos0(l,1)+i-1
              pos(2)=pos0(l,2)+s-1
              pos(3)=pos0(l,3)+k-1
              pos(1:3)=MATMUL(mat_h,pos(1:3))
              write(11,*) v,label_num(l),pos(1),pos(2),pos(3)
              v=v+1
             enddo
            enddo
           enddo
          enddo

         close(11)

        endif        

       return
       end program


