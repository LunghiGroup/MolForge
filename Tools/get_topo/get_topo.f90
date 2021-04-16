        program topo_ale
        use atoms_class, only : get_cov_radius,get_mass
        implicit none
        integer                          :: aa,ba,i,s,j,nat,ialloc,nbonds,nangles,max_options,nkind,nbond_kinds,nangle_kinds
        integer                          :: num_exclude=0,l
        integer,allocatable              :: label_num(:),bond_list(:,:),old_bond_list(:,:),angle_list(:,:),old_angle_list(:,:)
        integer,allocatable              :: bond_num(:),angle_num(:),ncoor(:)
        double precision,allocatable     :: coord(:,:),dist(:,:),charge(:)
        double precision                 :: bond_thr,a,b,c,alpha,beta,gamma,cell(6),mat_h(3,3),mat_hinv(3,3),&
                                            aaa,bbb,deg2rad,pi,thr,mass
        character*200                    :: geo_in,topo_out,reading,charge_in
        character*30, allocatable        :: label(:),kind(:),bond_kind(:),angle_kind(:),exclA(:),exclB(:)
        character*30, allocatable        :: old_kind(:),old_bond_kind(:),old_angle_kind(:)
        character(len=2)                 :: lab
        logical                          :: wrap_pbc=.false.,charges=.false.,list=.false.,nbnd=.false.,flag=.true.

         if( iargc().eq.0)then
          write(*,*) 'Get_Topo Usage:'               
          write(*,*) '-i    '               
          write(*,*) '-o    '               
          write(*,*) '-c    '               
          write(*,*) '-nbnd '               
          write(*,*) '-box  '               
          write(*,*) '-ex  '               
          stop
         endif

         do i=1,iargc()

          call getarg(i,reading)

          if (reading.eq.'-nbnd')then
           nbnd=.true.
          endif

          if (reading.eq.'-ex')then
           call getarg(i+1,reading)
           read(reading,*) num_exclude
           allocate(exclA(num_exclude))
           allocate(exclB(num_exclude))
           j=2
           do s=1,num_exclude
            call getarg(i+j,reading)
            read(reading,*) exclA(s) 
            call getarg(i+j+1,reading)
            read(reading,*) exclB(s) 
            j=j+2
           enddo

          endif

          if (reading.eq.'-i')then
           call getarg(i+1,reading)
           read(reading,*) geo_in
           open(10,file=trim(geo_in))
           write(*,*) 'reading coordinate file from ',trim(geo_in)
          endif

          if (reading.eq.'-c')then
           call getarg(i+1,reading)
           read(reading,*) charge_in
           open(9,file=trim(charge_in))
           write(*,*) 'charges read from file ',trim(charge_in)
           charges=.true.
           list=.true.
          endif

          if (reading.eq.'-o')then
           call getarg(i+1,reading)
           read(reading,*) topo_out
           open(11,file=trim(topo_out))
           write(*,*) 'final topology file will be written on ',trim(topo_out)
          endif


          if (reading.eq.'-box')then

           wrap_pbc=.true.
           call getarg(i+1,reading)
           read(reading,*) cell(1)
           call getarg(i+2,reading)
           read(reading,*) cell(2)
           call getarg(i+3,reading)
           read(reading,*) cell(3)
           call getarg(i+4,reading)
           read(reading,*) cell(4)
           call getarg(i+5,reading)
           read(reading,*) cell(5)
           call getarg(i+6,reading)
           read(reading,*) cell(6)
          endif

         enddo

         read(10,*) nat
         read(10,*)

         allocate(coord(nat,3),stat=ialloc)
         allocate(label(nat),stat=ialloc)
         allocate(label_num(nat),stat=ialloc)
         allocate(kind(1),stat=ialloc)
         nkind=1


         do i=1,nat
          read(10,*) label(i),coord(i,1),coord(i,2),coord(i,3)
         
          if (i.gt.1)then
            if( .not.any(label(i).eq.label(1:i-1) ) )then
             nkind=nkind+1
                 allocate(old_kind(nkind-1),stat=ialloc)
                 old_kind=kind
                 deallocate(kind)
                 allocate(kind(nkind),stat=ialloc)
                 kind(1:nkind-1)=old_kind(1:nkind-1)
                 deallocate(old_kind)
             kind(nkind)=trim(label(i)) 
            endif
           else
             kind(1)=trim(label(i)) 
          endif

         enddo
      

         do i=1,nat
          do j=1,nkind
          if(trim(label(i)).eq.trim(kind(j)))then
            label_num(i)=j
          endif
          enddo
         enddo


         if(.not. nbnd)then


         allocate(charge(nat),stat=ialloc)
         charge=0.0d0

         if(charges .and. list)then
          do i=1,nat
           read(9,*) charge(i)
          enddo
         endif

         if(wrap_pbc)then
 
           call xyztoxray(mat_hinv,cell)
           call xraytoxyz(mat_h,cell)

           allocate(dist(nat,nat),stat=ialloc)
           dist=0.0d0

           do i=1,nat
            do s=i+1,nat
             call dist_pbc(coord(i,1:3),coord(s,1:3),mat_h,mat_hinv,dist(i,s),1)
             dist(s,i)=dist(i,s)
            enddo
           enddo

         else

         allocate(dist(nat,nat),stat=ialloc)
         dist=0.0d0

         do i=1,nat
          do j=i+1,nat

           dist(i,j)=(coord(i,1)-coord(j,1))**2
           dist(i,j)=dist(i,j)+(coord(i,2)-coord(j,2))**2
           dist(i,j)=dist(i,j)+(coord(i,3)-coord(j,3))**2
           dist(i,j)=dsqrt(dist(i,j))
           dist(j,i)=dist(i,j)
 
          enddo
         enddo

        endif

        s=0

        allocate(bond_list(1,2),stat=ialloc)
        allocate(bond_kind(1),stat=ialloc)

        do i=1,nat
         do j=i+1,nat
          flag=.true.  
          do l=1,num_exclude          
           if ( (trim(label(i)).eq.exclA(l) .and. trim(label(j)).eq.exclB(l)) ) flag=.false.
           if ( (trim(label(i)).eq.exclB(l) .and. trim(label(j)).eq.exclA(l)) ) flag=.false.
          enddo

          if ( flag ) then

           lab=trim(label(i))
           call get_cov_radius(lab,thr)
           bond_thr=thr
           lab=trim(label(j))
           call get_cov_radius(lab,thr)
           bond_thr=bond_thr+thr

           if (dist(i,j).lt.bond_thr+0.3)then

           s=s+1
            if(s.gt.1)then        
     
               allocate(old_bond_list(s-1,2),stat=ialloc)
               old_bond_list=bond_list
               deallocate(bond_list)       
               allocate(bond_list(s,2),stat=ialloc)
               bond_list(1:s-1,1:2)=old_bond_list(1:s-1,1:2)
               deallocate(old_bond_list)       

            endif

            bond_list(s,1)=i
            bond_list(s,2)=j

                 if (s.gt.1)then


                   if( (.not. any( trim(label(i))//'-'//trim(label(j)).eq.bond_kind(1:nbond_kinds) )) .and. &
                       (.not. any( trim(label(j))//'-'//trim(label(i)).eq.bond_kind(1:nbond_kinds) ))   )  then

                        nbond_kinds=nbond_kinds+1
                        allocate(old_bond_kind(nbond_kinds-1),stat=ialloc)

                        old_bond_kind=bond_kind

                        deallocate(bond_kind)
                        allocate(bond_kind(nbond_kinds),stat=ialloc)
                        bond_kind(1:nbond_kinds-1)=old_bond_kind(1:nbond_kinds-1)
                        deallocate(old_bond_kind)
                        bond_kind(nbond_kinds)=trim(label(i))//'-'//trim(label(j))


                   endif

                  else

                        bond_kind(1)=trim(label(i))//'-'//trim(label(j))
                        nbond_kinds=1

                endif
        
           endif  
        
          endif
         enddo
        enddo



        nbonds=s
        write(*,*) nbonds,'bonds detected'

        allocate(bond_num(nbonds),stat=ialloc)
        allocate(ncoor(nat),stat=ialloc)
        ncoor=0

        do i=1,nbonds

         ncoor(bond_list(i,1))=ncoor(bond_list(i,1))+1
         ncoor(bond_list(i,2))=ncoor(bond_list(i,2))+1

                 do j=1,nbond_kinds
                   if( trim(label(bond_list(i,1)))//'-'//trim(label(bond_list(i,2))).eq. trim(bond_kind(j)) .or. &
                     trim(label(bond_list(i,2)))//'-'//trim(label(bond_list(i,1))).eq. trim(bond_kind(j))  )then
                     bond_num(i)=j
                   endif
                 enddo

        enddo

!        do i=1,nat   
!        if(label(i).eq.'Si' .and. ncoor(i).ne.4)then
!         write(*,*) label(i),i-1,ncoor(i)
!        endif
!
!        if(label(i).eq.'O' .and. ncoor(i).ne.2)then
!         write(*,*) label(i),i-1,ncoor(i)
!        endif
!        enddo

        s=0
        allocate(angle_list(1,3),stat=ialloc)
        allocate(angle_kind(1),stat=ialloc)

        do i=1,nbonds
         do j=i+1,nbonds

          if (bond_list(i,1).eq.bond_list(j,1) )then
          s=s+1
             if(s.gt.1)then
               allocate(old_angle_list(s-1,3),stat=ialloc)
               old_angle_list=angle_list
               deallocate(angle_list)       
               allocate(angle_list(s,3),stat=ialloc)
               angle_list(1:s-1,1:3)=old_angle_list(1:s-1,1:3)
               deallocate(old_angle_list)       
             endif

                angle_list(s,1)=bond_list(i,2)
                angle_list(s,2)=bond_list(i,1)
                angle_list(s,3)=bond_list(j,2)
                 if (s.gt.1)then


        if( (.not. any( trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3))) & 
              .eq. angle_kind(1:nangle_kinds) )) .and. &
            (.not. any( trim(label(angle_list(s,3)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,1))) & 
              .eq. angle_kind(1:nangle_kinds) ))    )  then

                        nangle_kinds=nangle_kinds+1
                        allocate(old_angle_kind(nangle_kinds-1),stat=ialloc)

                        old_angle_kind=angle_kind

                        deallocate(angle_kind)
                        allocate(angle_kind(nangle_kinds),stat=ialloc)
                        angle_kind(1:nangle_kinds-1)=old_angle_kind(1:nangle_kinds-1)
                        deallocate(old_angle_kind)
      angle_kind(nangle_kinds)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))


         endif

                  else

                 angle_kind(1)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))
                        nangle_kinds=1

                endif
          endif

          if (bond_list(i,1).eq.bond_list(j,2) )then

          s=s+1
             if(s.gt.1)then
               allocate(old_angle_list(s-1,3),stat=ialloc)
               old_angle_list=angle_list
               deallocate(angle_list)       
               allocate(angle_list(s,3),stat=ialloc)
               angle_list(1:s-1,1:3)=old_angle_list(1:s-1,1:3)
               deallocate(old_angle_list)
             endif

                angle_list(s,1)=bond_list(i,2)
                angle_list(s,2)=bond_list(i,1)
                angle_list(s,3)=bond_list(j,1)
                 if (s.gt.1)then


        if( (.not. any( trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3))) & 
              .eq. angle_kind(1:nangle_kinds) )) .and. &
            (.not. any( trim(label(angle_list(s,3)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,1))) & 
              .eq. angle_kind(1:nangle_kinds) ))    )  then

                        nangle_kinds=nangle_kinds+1
                        allocate(old_angle_kind(nangle_kinds-1),stat=ialloc)

                        old_angle_kind=angle_kind

                        deallocate(angle_kind)
                        allocate(angle_kind(nangle_kinds),stat=ialloc)
                        angle_kind(1:nangle_kinds-1)=old_angle_kind(1:nangle_kinds-1)
                        deallocate(old_angle_kind)
      angle_kind(nangle_kinds)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))


         endif

                  else

                 angle_kind(1)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))
                        nangle_kinds=1

                endif
          endif

          if (bond_list(i,2).eq.bond_list(j,1) )then

          s=s+1

             if(s.gt.1)then
               allocate(old_angle_list(s-1,3),stat=ialloc)
               old_angle_list=angle_list
               deallocate(angle_list)       
               allocate(angle_list(s,3),stat=ialloc)
               angle_list(1:s-1,1:3)=old_angle_list(1:s-1,1:3)
               deallocate(old_angle_list)
             endif

                angle_list(s,1)=bond_list(i,1)
                angle_list(s,2)=bond_list(i,2)
                angle_list(s,3)=bond_list(j,2)
                 if (s.gt.1)then


        if( (.not. any( trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3))) & 
              .eq. angle_kind(1:nangle_kinds) )) .and. &
            (.not. any( trim(label(angle_list(s,3)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,1))) & 
              .eq. angle_kind(1:nangle_kinds) ))    )  then

                        nangle_kinds=nangle_kinds+1
                        allocate(old_angle_kind(nangle_kinds-1),stat=ialloc)

                        old_angle_kind=angle_kind

                        deallocate(angle_kind)
                        allocate(angle_kind(nangle_kinds),stat=ialloc)
                        angle_kind(1:nangle_kinds-1)=old_angle_kind(1:nangle_kinds-1)
                        deallocate(old_angle_kind)
      angle_kind(nangle_kinds)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))


         endif

                  else

                 angle_kind(1)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))
                        nangle_kinds=1

                endif
          endif

          if (bond_list(i,2).eq.bond_list(j,2) )then

          s=s+1
            
             if(s.gt.1)then
               allocate(old_angle_list(s-1,3),stat=ialloc)
               old_angle_list=angle_list
               deallocate(angle_list)       
               allocate(angle_list(s,3),stat=ialloc)
               angle_list(1:s-1,1:3)=old_angle_list(1:s-1,1:3)
               deallocate(old_angle_list)
             endif

                angle_list(s,1)=bond_list(i,1)
                angle_list(s,2)=bond_list(i,2)
                angle_list(s,3)=bond_list(j,1)

                 if (s.gt.1)then


        if( (.not. any( trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3))) & 
              .eq. angle_kind(1:nangle_kinds) )) .and. &
            (.not. any( trim(label(angle_list(s,3)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,1))) & 
              .eq. angle_kind(1:nangle_kinds) ))    )  then

                        nangle_kinds=nangle_kinds+1
                        allocate(old_angle_kind(nangle_kinds-1),stat=ialloc)

                        old_angle_kind=angle_kind

                        deallocate(angle_kind)
                        allocate(angle_kind(nangle_kinds),stat=ialloc)
                        angle_kind(1:nangle_kinds-1)=old_angle_kind(1:nangle_kinds-1)
                        deallocate(old_angle_kind)
      angle_kind(nangle_kinds)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))


         endif

                  else

                 angle_kind(1)=trim(label(angle_list(s,1)))//'-'//trim(label(angle_list(s,2)))//'-'//trim(label(angle_list(s,3)))
                        nangle_kinds=1

                endif

          endif




         enddo
        enddo


        nangles=s
        write(*,*) nangles,'angles detected'
        

        allocate(angle_num(nangles),stat=ialloc)

        do i=1,nangles
         do j=1,nangle_kinds
           if( trim(label(angle_list(i,1)))//'-'//trim(label(angle_list(i,2)))//'-'//trim(label(angle_list(i,3))) & 
          .eq. trim(angle_kind(j)) .or. &
               trim(label(angle_list(i,3)))//'-'//trim(label(angle_list(i,2)))//'-'//trim(label(angle_list(i,1))) &
          .eq. trim(angle_kind(j))  )then
             angle_num(i)=j
           endif
         enddo
        enddo

        endif

        write(*,*) '##################################################'



        write(11,*) 'LAMMPS data file created by topoale'
        write(11,*) nat,' atoms'
        if(.not.nbnd)then
                write(11,*) nbonds,' bonds'
                write(11,*) nangles,' angles'       
                write(11,*) '0 dihedrals'
                write(11,*) '0 impropers'
        endif
        write(11,*) nkind,' atom types'
        if(.not.nbnd)then
                write(11,*) nbond_kinds,' bond types'
                write(11,*) nangle_kinds,' angle types'
                write(11,*) '0 dihedral types'
                write(11,*) '0 improper types'
        endif
! metti cella triclina per lammps

        pi=acos(-1.0d0)
        deg2rad=pi/180.0d0


        cell(4)=cell(4)*deg2rad
        cell(5)=cell(5)*deg2rad
        cell(6)=cell(6)*deg2rad

        aaa=cell(3)*cos(cell(5))
        bbb=(cell(2)*cell(3)*cos(cell(4))-cell(2)*cos(cell(6))*cell(3)*cos(cell(5)))/(cell(2)*sin(cell(6)))

        write(11,*) 
        write(11,*) '0.0',cell(1),' xlo xhi '
        write(11,*) '0.0',cell(2)*sin(cell(6)),' ylo yhi '
        write(11,*) '0.0',sqrt(cell(3)**2-aaa**2-bbb**2),' zlo zhi '
        write(11,*) cell(2)*cos(cell(6)),cell(3)*cos(cell(5)),&
                    (cell(2)*cell(3)*cos(cell(4))-cell(2)*cos(cell(6))*cell(3)*cos(cell(5)))/(cell(2)*sin(cell(6))),' xy xz yz'
        

        write(11,*) 
        write(11,*) 'Masses'
        write(11,*) 
        do i=1,nkind
         lab=trim(kind(i))
         call get_masS(lab,mass)
         write(11,*)  i,mass!'insert_mass',trim(kind(i))
        enddo
        
        write(11,*) 
        write(11,*) 'Atoms'
        write(11,*) 

        do i=1,nat

        if(.not.nbnd)then
         write(11,*) i,'1',label_num(i),charge(i),coord(i,1),coord(i,2),coord(i,3),'#',trim(label(i))
        else
         write(11,*) i,label_num(i),coord(i,1),coord(i,2),coord(i,3),'#',trim(label(i))
        endif

        enddo

        if(.not.nbnd)then
        write(11,*) 
        write(11,*) 'Bonds'
        write(11,*) 

        do i=1,nbonds
         write(11,*) i,bond_num(i),bond_list(i,1),bond_list(i,2),'#',trim(label(bond_list(i,1)))//'-'//trim(label(bond_list(i,2)))
        enddo

        write(11,*) 
        write(11,*) 'Angles'
        write(11,*) 

        do i=1,nangles
           write(11,*) i,angle_num(i),angle_list(i,1),angle_list(i,2),angle_list(i,3),'#', &
           trim(label(angle_list(i,1)))//'-'//trim(label(angle_list(i,2)))//'-'//trim(label(angle_list(i,3)))

        enddo
        endif

        return
        end




