        module spinham_map_class
        implicit none

        contains

         subroutine msh2gsh(msh,gsh,lmax,compress,lambda)
         use hilbert_dist_class
         use orthonorm_class
         use lapack_diag_simm
         type(spins_hilbert)            :: msh,gsh
         integer                        :: lmax,i,k1,k2
         double precision, allocatable  :: ener(:),Ms(:)
         double complex, allocatable    :: O(:),Sx(:,:),Sz(:,:),Sx2(:,:)
         double precision               :: lambda,phi
         logical                        :: compress
                       
!          call ProjMG(msh,gsh)
!          call orthonorm_set(gsh%Hdim,gsh%H0(1)%mat)

          allocate(Sz(gsh%Hdim,gsh%Hdim))
          allocate(Ms(gsh%Hdim))
          allocate(Sx(gsh%Hdim,gsh%Hdim))
          allocate(Sx2(gsh%Hdim,gsh%Hdim))
          Sz=msh%Sz(1)%mat(1:gsh%Hdim,1:gsh%Hdim)
          Sx=msh%Sx(1)%mat(1:gsh%Hdim,1:gsh%Hdim)
          Sx2=(0.0d0,0.0d0)

          call new_diag(gsh%Hdim,Sz,Ms)

          do i=1,gsh%Hdim-1
           do k1=1,gsh%Hdim
            do k2=1,gsh%Hdim
             Sx2(i,i+1)=Sx2(i,i+1)&
                            +Sx(k1,k2)*conjg(Sz(k1,i))*Sz(k2,i+1)
            enddo
           enddo
           phi=atan2(aimag(Sx2(i,i+1)),dble(Sx2(i,i+1)))
           Sz(:,i+1)=Sz(:,i+1)*exp(cmplx(0.0d0,-phi,8))
          enddo         

          Sz=transpose(conjg(Sz))

          allocate(ener(1:gsh%Hdim))
          ener=msh%Ener(1)%v(1:gsh%Hdim)
!          call projH(lmax,gsh%Hdim,gsh%spin(1),ener,gsh%H0(1)%mat,O,compress,lambda)          
          call projH(lmax,gsh%Hdim,gsh%spin(1),ener,Sz,O,compress,lambda)          
          call diagH(gsh%Hdim,gsh%spin(1),lmax,O)

         return
         end subroutine msh2gsh

         subroutine ProjMG(msh,gsh)
         use hilbert_dist_class
         implicit none
         integer                       :: l,s,i
         type(spins_hilbert)           :: msh,gsh

          gsh%H0(1)%mat=(0.0d0,0.0d0)

          do l=1,gsh%Hdim
           do s=1,gsh%Hdim
            do i=1,msh%Hdim
             gsh%H0(1)%mat(s,l)=gsh%H0(1)%mat(s,l)+msh%H(1)%mat(i,l)*&
                               newProjSzS(i,msh%spin(1),msh%basis(i,1),2,gsh%spin(1),gsh%basis(s,1),msh%spin,msh%basis)
            enddo
           enddo
          enddo


        return
        end subroutine projMG

        recursive function newProjSzS(i,si,szi,counter,S,Ms,spin,base) result(AP)
        use anglib
        implicit none
        double precision, allocatable   :: base(:,:),spin(:)
        double precision                :: AP,si,szi,S,MS,sj,szj,smin,smax
        integer                         :: i,counter,v,t

         sj=spin(counter)        
         szj=base(i,counter)
         smax=si+sj
         smin=ABS(si-sj)

         AP=0.0

         if(counter.lt.size(base,2))then

          do v=NINT(2*smin),NINT(2*smax),2
           do t=-v,v,2
            AP=AP+cleb(NINT(2*si),NINT(2*szi),NINT(2*sj),NINT(2*szj),v,t)*newProjSzS(i,v/2.0d0,t/2.0d0,counter+1,S,Ms,spin,base)
           enddo
          enddo
         else
            AP=cleb(NINT(2*si),NINT(2*szi),NINT(2*sj),NINT(2*szj),NINT(2*S),NINT(2*Ms))
         endif

        return
        end function newProjSzS

        subroutine diagH(Hdim,jj,lmax,O)
        use lapack_diag_simm
        use stevens_class
        integer                        :: lmax,Hdim,i,j,k,l,q,s
        double precision               :: j1,j2,jj
        double precision, allocatable  :: Ener(:)
        double complex, allocatable    :: Hmat(:,:),O(:)
        double complex                 :: val_tmp,val

         allocate(Hmat(Hdim,Hdim))         
         allocate(Ener(Hdim))         

         j1=-jj
         do i=1,Hdim
          j2=-jj
          do j=1,Hdim

           val=(0.0d0,0.0d0)
           s=1
           do l=2,lmax,2
            do q=-l,l
             val_tmp=(0.0d0,0.0d0)
             call stevens_mat_elem(l,q,jj,j1,jj,j2,val_tmp)
             val=val+val_tmp*O(s)
             s=s+1
            enddo
           enddo

           Hmat(i,j)=val

           j2=j2+1.0d0
          enddo
          j1=j1+1.0d0
         enddo

         call new_diag(Hdim,Hmat,Ener)

         write(*,*) (Ener-Ener(1))

        return
        end subroutine diagH

        subroutine projH(lmax,Nj,jj,Ener,Jz,O,compress,lambda)
        use stevens_class
        use particles_swarm_class
        implicit none
        integer                        :: lmax,Odim,Nj,dimB,dimA
        integer                        :: i,j,l,s,k,i1,i2,v,lwork,inf
        double precision, allocatable  :: Ener(:)
        double precision               :: jj,q1,q2,avg_ener,lambda
        double complex, allocatable    :: Jz(:,:),B(:),A(:,:),O(:),work(:)
        double complex                 :: mat_elem

        logical                        :: compress
        character(len=10)              :: lsmf_opt


         lsmf_opt='RIDGE'

         Odim=0
         do i=2,lmax,2
          Odim=Odim+(2*i+1)
         enddo

         avg_ener=0.0d0
         do i=1,Nj
          avg_ener=avg_ener+Ener(i)
         enddo
         avg_ener=avg_ener/Nj

         allocate(O(Odim))
         if(compress)then
          dimB=Nj*Nj+Odim
          dimA=Odim+1
          allocate(B(dimB))
          allocate(A(dimB,dimA))
         else
          dimB=Nj*Nj
          dimA=Odim+1
          allocate(B(dimB))
          allocate(A(dimB,dimA))
         endif

         A=(0.0d0,0.0d0)        
         B=(0.0d0,0.0d0)        
         O=(0.0d0,0.0d0)        

         k=1
         do i1=1,Nj
          do i2=1,Nj
          
           do v=1,Nj
            B(k)=B(k)+(Ener(v))*conjg(Jz(i2,v))*Jz(i1,v)
           enddo

           q1=i1-((Nj-1)/2.0d0)-1
           q2=i2-((Nj-1)/2.0d0)-1
           if(q1.eq.q2)then
            A(k,1)=1.0d0
           else
            A(k,1)=0.0d0
           endif
           s=2
           do l=2,lmax,2
            do j=-l,l
             call stevens_mat_elem(l,j,jj,q1,jj,q2,mat_elem)
             A(k,s)=mat_elem
             s=s+1
            enddo
           enddo

           k=k+1
          enddo
         enddo

         if(compress)then
          do s=2,Odim+1
           A(k,s)=lambda        
           k=k+1
          enddo 
         endif


         select case (lsmf_opt)

         case ("LASSO")
          write(*,*) 'LASSO not implemented yet'
          stop

         case ("RIDGE")

          lwork=dimB+64*dimB+1000
          allocate(work(lwork))

          call zgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
           if(inf.ne.0)then
           write(*,*) 'zgels failed',inf
           stop
          endif

         end select

         write(*,*) 'Crystal Field Parameters:'
         write(*,*) 0,0,dble(B(1)),aimag(B(1))
         s=2
         do l=2,lmax,2
          do j=-l,l
           write(*,*) l,j,dble(B(s)),aimag(B(s))
           O(s-1)=B(s)
           s=s+1
          enddo
         enddo
 
        return
        end subroutine projH



        end module spinham_map_class



