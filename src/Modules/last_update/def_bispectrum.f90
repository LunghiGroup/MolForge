        module bispectrum_class
        implicit none

         type hyper_spherical_harmonic
          double complex, allocatable   :: Y(:)
          integer                       :: order !! 2J
          contains
          procedure                     :: setup => setup_hsh
          procedure                     :: reset => reset_hsh
          procedure                     :: get_val => get_val
          procedure                     :: delete => delete_hsh
         end type hyper_spherical_harmonic 

         type bispectrum_components
          type(hyper_spherical_harmonic), allocatable  :: SH(:)
          double complex, allocatable                  :: B(:)
          integer                                      :: max_order !! 2J Max
          integer                                      :: tot_bi
          double precision, allocatable                :: Cmat(:,:,:,:,:,:)
          contains
          procedure                                    :: setup => setup_bi
          procedure                                    :: reset => reset_bi
          procedure                                    :: get_bi => get_bi
          procedure                                    :: delete => delete_bi
          procedure                                    :: get_CG
         end type bispectrum_components

        contains

         subroutine delete_bi(this)
         implicit none
         class(bispectrum_components) :: this
         integer                      :: i
          deallocate(this%B)
          deallocate(this%Cmat)
          do i=1,size(this%SH)
           call this%SH(i)%delete()
          enddo
          deallocate(this%SH)
         return
         end subroutine delete_bi

         subroutine delete_hsh(this)
         implicit none
         class(hyper_spherical_harmonic) :: this
         integer                         :: i
          deallocate(this%Y)
         return
         end subroutine delete_hsh

         subroutine cart2spher(vec1,vec2)
         implicit none
         double precision        :: vec1(3),vec2(3),a,b

         ! vec2(1)=|r| vec2(2)=phi vec2(3)=theta

          vec2=0.0d0
          vec2(1)=vec2(1)+vec1(1)**2
          vec2(1)=vec2(1)+vec1(2)**2
          vec2(1)=vec2(1)+vec1(3)**2
          vec2(1)=sqrt(vec2(1))

          vec2(2)=atan2(vec1(2),vec1(1))
          vec2(3)=acos(vec1(3)/vec2(1))

         return
         end subroutine cart2spher

         subroutine setup_bi(this)
         implicit none
         class(bispectrum_components) :: this
         integer                      :: i

          if(allocated(this%SH)) call this%delete()
          allocate(this%SH(this%max_order+1))

          do i=0,this%max_order 
           this%SH(i+1)%order=i
           call this%SH(i+1)%setup()
          enddo

          this%tot_bi=nint((this%max_order/2.0d0+1)*(this%max_order/2.0d0+2)*&
          (this%max_order/2.0d0+1.50d0)/3.0d0)

          if(allocated(this%B)) deallocate(this%B)
          allocate(this%B(this%tot_bi))
          this%B=(0.0d0,0.0d0)

          call this%get_CG()

         return
         end subroutine setup_bi

         subroutine setup_hsh(this)
         implicit none
         class(hyper_spherical_harmonic) :: this
         integer                         :: i,j,indx
         
          if(allocated(this%Y)) deallocate(this%Y)
          allocate(this%Y((this%order+1)**2))
          this%Y=(0.0d0,0.0d0)
          indx=1
          do i=-this%order,this%order,2
           do j=-this%order,this%order,2
            if(i.eq.j) this%Y(indx)=(1.0d0,0.0d0)
            if(i.ne.j) this%Y(indx)=(0.0d0,0.0d0)
            indx=indx+1
           enddo
          enddo

         return
         end subroutine setup_hsh

         subroutine reset_bi(this)
         implicit none
         class(bispectrum_components) :: this
         integer                      :: i

          this%B=(0.0d0,0.0d0)
          do i=0,this%max_order 
           call this%SH(i+1)%reset()
          enddo

         return
         end subroutine reset_bi

         subroutine reset_hsh(this)
         implicit none
         class(hyper_spherical_harmonic) :: this
         integer                         :: i,j,indx
         
          this%Y=(0.0d0,0.0d0)
          indx=1
          do i=-this%order,this%order,2
           do j=-this%order,this%order,2
            if(i.eq.j) this%Y(indx)=(1.0d0,0.0d0)
            if(i.ne.j) this%Y(indx)=(0.0d0,0.0d0)
            indx=indx+1
           enddo
          enddo

         return
         end subroutine reset_hsh

         subroutine get_CG(this)
         use anglib
         implicit none
         class(bispectrum_components)      :: this
         integer                           :: m1,m3,m5
         integer                           :: j1,j2,j3,i,ii,iii

          if(allocated(this%Cmat)) deallocate(this%Cmat)
          allocate(this%Cmat(0:this%max_order,0:this%max_order,0:this%max_order,&
                        -this%max_order:this%max_order,&
                        -this%max_order:this%max_order,&
                        -this%max_order:this%max_order))

          this%Cmat=0.0d0

          do i=0,this%max_order
          j1=this%SH(i+1)%order
           do ii=0,i
           j2=this%SH(ii+1)%order
            do iii=(i-ii),min(this%max_order,i+ii),2
            j3=this%SH(iii+1)%order

             do m1=-this%SH(i+1)%order,this%SH(i+1)%order,2
             do m3=max(-this%SH(ii+1)%order,-this%SH(iii+1)%order-m1),&
                   min(this%SH(ii+1)%order,this%SH(iii+1)%order-m1),2

              m5=m1+m3

              this%Cmat(j2,j1,j3,m3,m1,m5)=cleb(j2,m3,j1,m1,j3,m5)

             enddo
             enddo
            enddo
           enddo
          enddo

         return
         end subroutine get_CG

         subroutine get_bi(this)
         use anglib
         implicit none
         class(bispectrum_components)      :: this
         integer                           :: m1,m2,m3,m4,m5,m6
         integer                           :: j1,j2,j3,indx,i,ii,iii
         integer                           :: indx1,indx2,indx3

          indx=1

          do i=0,this%max_order
          j1=this%SH(i+1)%order
           do ii=0,i
           j2=this%SH(ii+1)%order
            do iii=(i-ii),min(this%max_order,i+ii),2
            j3=this%SH(iii+1)%order

             if (iii.lt.i) cycle             
      
             do m1=-this%SH(i+1)%order,this%SH(i+1)%order,2
             do m3=max(-this%SH(ii+1)%order,-this%SH(iii+1)%order-m1),&
                   min(this%SH(ii+1)%order,this%SH(iii+1)%order-m1),2
              m5=m1+m3

             do m2=-this%SH(i+1)%order,this%SH(i+1)%order,2
             do m4=max(-this%SH(ii+1)%order,-this%SH(iii+1)%order-m2),&
                   min(this%SH(ii+1)%order,this%SH(iii+1)%order-m2),2
              m6=m2+m4

              indx1=((this%SH(i+1)%order+m1)/2)*(this%SH(i+1)%order+1)+&
                    ((this%SH(i+1)%order+m2)/2+1)
              indx2=((this%SH(ii+1)%order+m3)/2)*(this%SH(ii+1)%order+1)+&
                    ((this%SH(ii+1)%order+m4)/2+1)
              indx3=((this%SH(iii+1)%order+m5)/2)*(this%SH(iii+1)%order+1)+&
                    ((this%SH(iii+1)%order+m6)/2+1)


               this%B(indx)=this%B(indx)+this%SH(i+1)%Y(indx1)*this%SH(ii+1)%Y(indx2)*&
                                conjg(this%SH(iii+1)%Y(indx3))*&
                                this%Cmat(j2,j1,j3,m4,m2,m6)*this%Cmat(j2,j1,j3,m3,m1,m5)

             enddo
             enddo
             enddo
             enddo
             this%B(indx)=this%B(indx)-j3-1
             indx=indx+1 
            enddo
           enddo
          enddo

         return
         end subroutine get_bi

         subroutine get_val(this,vec1,r0,weight)
         implicit none
         class(hyper_spherical_harmonic)   :: this
         double precision                  :: vec(3),vec1(3)
         double precision                  :: mm1,mm2,mm3
         double precision                  :: pi,smooth_fact,zeta,r0,weight
         double precision                  :: dist,phi,theta,theta0
         double complex                    :: jj
         double complex, allocatable       :: wig(:,:),wig1(:,:),wig2(:,:)
         integer                           :: indx,m1,m2,m3
         

          pi=acos(-1.0d0)
          jj=cmplx(0.0d0,1.0d0,8)

          vec=0.0d0
          vec(1)=vec(1)+vec1(1)**2
          vec(1)=vec(1)+vec1(2)**2
          vec(1)=vec(1)+vec1(3)**2
          vec(1)=sqrt(vec(1))

          if (vec(1).gt.r0) return

          dist=vec(1)
          vec(2)=atan2(vec1(2),vec1(1))
          vec(3)=acos(vec1(3)/vec(1))

          phi=vec(2)
          theta=vec(3)
          theta0=pi*dist/r0

          smooth_fact=0.5d0*(cos(pi*dist/r0)+1.0d0)
      
          call rot_wig(dble(this%order/2.0d0),phi,theta,-phi,wig1)
          call rot_wig(dble(this%order/2.0d0),phi,-theta,-phi,wig2)

          indx=1

          do m1=-this%order,this%order,2
           do m2=-this%order,this%order,2
            do m3=-this%order,this%order,2
             mm1=m1/2.0d0
             mm2=m2/2.0d0
             mm3=m3/2.0d0
             this%Y(indx)=this%Y(indx)+smooth_fact*&
              wig1(nint(mm1+this%order/2.0d0+1),nint(mm3+this%order/2.0d0+1))*&
              exp(-jj*mm3*theta0*2.0d0)*&  ! is the two there?
              wig2(nint(mm3+this%order/2.0d0+1),nint(mm2+this%order/2.0d0+1))
            enddo
            indx=indx+1
           enddo
          enddo
          
         return
         end subroutine get_val

         subroutine cutoff_function(dist,r0,smooth_fact)
         implicit none
         double precision       :: dist,r0,smooth_fact,pi
          pi=acos(-1.0d0)
          if (dist.le.r0)then
           smooth_fact=0.5d0*(cos(pi*dist/r0)+1.0d0)
          else
           smooth_fact=0.0d0
          endif
         return 
         end subroutine cutoff_function

         subroutine wigner(k,beta,rot)
         use anglib
         implicit none
         integer                                             :: v,t,s,mins,maxs
         double precision                                    :: beta,a,c,k,djmn
         complex(8), allocatable                             :: rot(:,:)

          if(allocated(rot)) deallocate(rot)
          allocate(rot(NINT(2*k+1),NINT(2*k+1)))
          rot=(0.0d0,0.0d0)

          do v=-NINT(2*k),NINT(2*k),2
           do t=-NINT(2*k),NINT(2*k),2

            mins=MAX( 0,NINT(t/2.-v/2.) )
            maxs=MIN( NINT(k+t/2.),NINT((k-v/2.)) )
 
            do s=mins,maxs
             a=factorial(NINT(k+t/2.0d0-s))*factorial(s)*factorial(NINT(v/2.-t/2.+s))*factorial(NINT(k-v/2.-s))
             a=((-1.0d0)**(v/2.0d0-t/2.0d0+s))/a
             c=NINT(2.0d0*k-v/2.0d0+t/2.0d0-2.0d0*s)
             a=a*((dcos(beta/2.0d0))**c)
             c=NINT(v/2.0d0-t/2.0d0+2.0d0*s)
             a=a*((dsin(beta/2.0d0))**c)
             rot(NINT(v/2.+k+1),NINT(t/2.0d0+k+1))=rot(NINT(v/2.0d0+k+1),NINT(t/2.0d0+k+1))+CMPLX(a,0.0d0,8)
            enddo

            rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*& 
            CMPLX(dsqrt(factorial(NINT(k+v/2.))*factorial(NINT(k-v/2.))*factorial(NINT(k+t/2.))*factorial(NINT(k-t/2.))),0.0d0,8)

           enddo
          enddo

         return
         end subroutine wigner

         subroutine rot_Wig(k,alpha,beta,gamma,rot)
         use anglib
         implicit none
         integer                                             :: v,t,s,mins,maxs,ialloc
         double precision                                    :: alpha,beta,gamma,a,c,k,djmn
         complex(8)                                          :: b
         complex(8), allocatable                             :: rot(:,:)


         if(allocated(rot)) deallocate(rot)
         allocate(rot(NINT(2*k+1),NINT(2*k+1)))
         rot=(0.0d0,0.0d0)

         do v=-NINT(2*k),NINT(2*k),2
          do t=-NINT(2*k),NINT(2*k),2

           mins=MAX( 0,NINT(t/2.-v/2.) )
           maxs=MIN( NINT(k+t/2.),NINT((k-v/2.)) )
 
           do s=mins,maxs
            a=factorial(NINT(k+t/2.0d0-s))*factorial(s)*factorial(NINT(v/2.-t/2.+s))*factorial(NINT(k-v/2.-s))
            a=((-1.0d0)**(v/2.0d0-t/2.0d0+s))/a
            c=NINT(2.0d0*k-v/2.0d0+t/2.0d0-2.0d0*s)
            a=a*((dcos(beta/2.0d0))**c)
            c=NINT(v/2.0d0-t/2.0d0+2.0d0*s)
            a=a*((dsin(beta/2.0d0))**c)
            rot(NINT(v/2.+k+1),NINT(t/2.0d0+k+1))=rot(NINT(v/2.0d0+k+1),NINT(t/2.0d0+k+1))+CMPLX(a,0.0d0,8)
           enddo

           rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*& 
           CMPLX(dsqrt(factorial(NINT(k+v/2.))*factorial(NINT(k-v/2.))*factorial(NINT(k+t/2.))*factorial(NINT(k-t/2.))),0.0d0,8)

          enddo
         enddo

         do v=-NINT(2*k),NINT(2*k),2
          do t=-NINT(2*k),NINT(2*k),2
           b=CMPLX(0.0d0,(-1.0d0*(t/2.0d0)*gamma),8)
           rot(NINT(v/2.+k+1),NINT(t/2.+k+1))=rot(NINT(v/2.+k+1),NINT(t/2.+k+1))*EXP(b)
           b=CMPLX(0.0d0,(-1.0d0*(v/2.0d0)*alpha),8)
           rot(NINT(v/2.0d0+k+1.0d0),NINT(t/2.0d0+k+1.0d0))=EXP(b)*rot(NINT(v/2.+k+1),NINT(t/2.+k+1))
          enddo
         enddo

         rot=transpose(rot)

        return
        end subroutine rot_Wig

        end module bispectrum_class
