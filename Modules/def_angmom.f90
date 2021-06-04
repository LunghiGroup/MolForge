!    anglib.f90: angular momentum coupling coefficients in Fortran 90
!    Copyright (C) 1998  Paul Stevenson
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
!    02110-1301,USA


module anglib
! Library of angular momentum coupling coefficient routines in fortran 90
! Paul Stevenson, p.stevenson@surrey.ac.uk

  integer, parameter :: rk = selected_real_kind(p=15)
contains

  function cleb(j1,m1,j2,m2,j,m)
    implicit none
    ! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
    ! arguments are integer and twice the true value. 

    real(rk)    :: cleb,factor,sum
    integer :: j1,m1,j2,m2,j,m,par,z,zmin,zmax

    ! some checks for validity (let's just return zero for bogus arguments)

    if (2*(j1/2)-int(2*(j1/2.0)) /= 2*(abs(m1)/2)-int(2*(abs(m1)/2.0)) .or. & 
         2*(j2/2)-int(2*(j2/2.0)) /= 2*(abs(m2)/2)-int(2*(abs(m2)/2.0)) .or. & 
         2*(j/2)-int(2*(j/2.0)) /= 2*(abs(m)/2)-int(2*(abs(m)/2.0)) .or. & 
         j1<0 .or. j2<0 .or. j<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or.& 
         abs(m)>j .or. j1+j2<j .or. abs(j1-j2)>j .or. m1+m2/=m) then
       cleb= 0.0
    else
    
       factor = 0.0
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = sqrt(factor)
       
       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum=0.0
       do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
          sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*& 
               binom((-j1+j2+j)/2,(j2+m2)/2-z)
       end do
       
       cleb = factor*sum
    end if
  end function cleb


  function sixj(a,b,c,d,e,f)
    implicit none
    integer, intent(in) :: a,b,c,d,e,f
    real(rk) :: sixj
    integer :: nlo, nhi, n
    real(rk) :: outfactors, sum, sumterm
    ! calculates a Wigner 6-j symbol. Argument a-f are integer and are
    ! twice the true value of the 6-j's arguments, in the form
    ! { a b c }
    ! { d e f }
    ! Calculated using binomial coefficients to allow for (reasonably) high
    ! arguments.

    ! First check for consistency of arguments:
    sixj=0.0
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return

    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)

    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

    sum=0.0
    do n=nlo,nhi
       sumterm = (-1)**n
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2)
       sum=sum+sumterm
    end do

    sixj = sum * outfactors

  end function sixj

  function angdelta(a,b,c)
    implicit none
    integer :: a,b,c
    real(rk)    :: angdelta, scr1
    ! calculate the function delta as defined in varshalovich et al. for
    ! use in 6-j symbol:
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta

  function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    integer :: a,b,c,d,e,f,g,h,i
    real(rk)    :: ninej, sum
    integer :: xlo, xhi
    integer :: x
    ! calculate a 9-j symbol. The arguments are given as integers twice the
    ! value of the true arguments in the form
    ! { a b c }
    ! { d e f }
    ! { g h i }

    ninej=0.0
    ! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return
    
    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)
    
    sum=0.0
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*& 
            sixj(g,h,i,x,a,d)
    end do
    ninej=sum

  end function ninej

  recursive function factorial(n) result(res)
    implicit none
    integer :: n
    real(rk) :: res

    if (n .lt. 0)then
     write(*,*) 'Factorial function crush coused by negative integer argument'
     stop
    endif

    if (n==0 .or. n==1) then
       res=1.0
    else
       res=n*factorial(n-1)
    end if
  end function factorial


  recursive function binom(n,r) result(res)
    implicit none
    integer :: n,r
    real(rk) :: res

    if(n==r .or. r==0) then
       res = 1.0
    else if (r==1) then
       res = real(n,rk)
    else
       res = real(n,rk)/real(n-r,rk)*binom(n-1,r)
    end if
  end function binom
end module anglib