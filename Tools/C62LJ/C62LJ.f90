  program main
  implicit none
  double precision      :: C6eff,C60,Reff,R0,sigma,eps,Vrel
  character(LEN=16)     :: word,label

   if (iargc().eq.0)then
    write(*,*) 'Usage: C62LJ Label C6eff'
    stop
   endif

   call getarg(1,label) 
   call getarg(2,word)
   word=trim(word)
   read(word,*) C6eff

   call GetVdWParam0(label,C60,R0)
   
   Vrel=sqrt(C6eff/C60)
   Reff=Vrel**(1.0d0/3.0d0)*R0

   sigma=0.529177*2*Reff*2**(-1.0d0/6.0d0)
   eps=627.509*C6eff/(4*(2*Reff*2**(-1.0d0/6.0d0))**6)

   write(*,*) 'eps= ',eps,' sigma= ',sigma

  return
  end


  SUBROUTINE GetVdWParam0(atom,C6,R0)
  IMPLICIT NONE
  !
  CHARACTER(LEN=3) :: atom
  double precision :: C6,alpha,R0
  !
  SELECT CASE (atom) 
  !
  CASE ('H')
  alpha=4.500000d0
  C6=6.500000d0
  R0=3.100000d0
  !
  CASE ('He')
  alpha=1.380000d0
  C6=1.460000d0
  R0=2.650000d0
  !
  CASE ('Li')
  alpha=164.200000d0
  C6=1387.000000d0
  R0=4.160000d0
  !
  CASE ('Be')
  alpha=38.000000d0
  C6=214.000000d0
  R0=4.170000d0
  !
  CASE ('B')
  alpha=21.000000d0
  C6=99.500000d0
  R0=3.890000d0
  !
  CASE ('C')
  alpha=12.000000d0
  C6=46.600000d0
  R0=3.590000d0
  !
  CASE ('N')
  alpha=7.400000d0
  C6=24.200000d0
  R0=3.340000d0
  !
  CASE ('O')
  alpha=5.400000d0
  C6=15.600000d0
  R0=3.190000d0
  !
  CASE ('F')
  alpha=3.800000d0
  C6=9.520000d0
  R0=3.040000d0
  !
  CASE ('Ne')
  alpha=2.670000d0
  C6=6.380000d0
  R0=2.910000d0
  !
  CASE ('Na')
  alpha=162.700000d0
  C6=1556.000000d0
  R0=3.730000d0
  !
  CASE ('Mg')
  alpha=71.000000d0
  C6=627.000000d0
  R0=4.270000d0
  !
  CASE ('Al')
  alpha=60.000000d0
  C6=528.000000d0
  R0=4.330000d0
  !
  CASE ('Si')
  alpha=37.000000d0
  C6=305.000000d0
  R0=4.200000d0
  !
  CASE ('P')
  alpha=25.000000d0
  C6=185.000000d0
  R0=4.010000d0
  !
  CASE ('S')
  alpha=19.600000d0
  C6=134.000000d0
  R0=3.860000d0
  !
  CASE ('Cl')
  alpha=15.000000d0
  C6=94.600000d0
  R0=3.710000d0
  !
  CASE ('Ar')
  alpha=11.100000d0
  C6=64.300000d0
  R0=3.550000d0
  !
  CASE ('K')
  alpha=292.900000d0
  C6=3897.000000d0
  R0=3.710000d0
  !
  CASE ('Ca')
  alpha=160.000000d0
  C6=2221.000000d0
  R0=4.650000d0
  !
  CASE ('Sc')
  alpha=120.000000d0
  C6=1383.000000d0
  R0=4.590000d0
  !
  CASE ('Ti')
  alpha=98.000000d0
  C6=1044.000000d0
  R0=4.510000d0
  !
  CASE ('V')
  alpha=84.000000d0
  C6=832.000000d0
  R0=4.440000d0
  !
  CASE ('Cr')
  alpha=78.000000d0
  C6=602.000000d0
  R0=3.990000d0
  !
  CASE ('Mn')
  alpha=63.000000d0
  C6=552.000000d0
  R0=3.970000d0
  !
  CASE ('Fe')
  alpha=56.000000d0
  C6=482.000000d0
  R0=4.230000d0
  !
  CASE ('Co')
  alpha=50.000000d0
  C6=408.000000d0
  R0=4.180000d0
  !
  CASE ('Ni')
  alpha=48.000000d0
  C6=373.000000d0
  R0=3.820000d0
  !
  CASE ('Cu')
  alpha=42.000000d0
  C6=253.000000d0
  R0=3.760000d0
  !
  CASE ('Zn')
  alpha=40.000000d0
  C6=284.000000d0
  R0=4.020000d0
  !
  CASE ('Ga')
  alpha=60.000000d0
  C6=498.000000d0
  R0=4.190000d0
  !
  CASE ('Ge')
  alpha=41.000000d0
  C6=354.000000d0
  R0=4.200000d0
  !
  CASE ('As')
  alpha=29.000000d0
  C6=246.000000d0
  R0=4.110000d0
  !
  CASE ('Se')
  alpha=25.000000d0
  C6=210.000000d0
  R0=4.040000d0
  !
  CASE ('Br')
  alpha=20.000000d0
  C6=162.000000d0
  R0=3.930000d0
  !
  CASE ('Kr')
  alpha=16.800000d0
  C6=129.600000d0
  R0=3.820000d0
  !
  CASE ('Rb')
  alpha=319.200000d0
  C6=4691.000000d0
  R0=3.720000d0
  !
  CASE ('Sr')
  alpha=199.000000d0
  C6=3170.000000d0
  R0=4.540000d0
  !
  CASE ('Y')
  alpha=126.7370d0
  C6=1968.580d0
  R0=4.81510d0
  !
  CASE ('Zr')
  alpha=119.97d0
  C6=1677.91d0
  R0=4.53d0
  !
  CASE ('Nb')
  alpha=101.603d0
  C6=1263.61d0
  R0=4.2365d0
  !
  CASE ('Mo')
  alpha=88.4225785d0
  C6=1028.73d0
  R0=4.099d0
  !
  CASE ('Tc')
  alpha=80.083d0
  C6=1390.87d0
  R0=4.076d0
  !
  CASE ('Ru')
  alpha=65.8950d0
  C6=609.754d0
  R0=3.99530d0
  !
  CASE ('Rh')
  alpha=56.1d0
  C6=469.0d0
  R0=3.95d0
  !
  CASE ('Pd')
  alpha=23.680000d0
  C6=157.500000d0
  R0=3.66000d0
  !
  CASE ('Ag')
  alpha=50.600000d0
  C6=339.000000d0
  R0=3.820000d0
  !
  CASE ('Cd')
  alpha=39.7d0
  C6=452.0d0
  R0=3.99d0
  !
  CASE ('In')
  alpha=70.22000d0
  C6=707.046000d0
  R0=4.23198000d0
  !
  CASE ('Sn')
  alpha=55.9500d0
  C6=587.41700d0
  R0=4.303000d0
  !
  CASE ('Sb')
  alpha=43.671970d0
  C6=459.322d0
  R0=4.2760d0
  !
  CASE ('Te')
  alpha=37.65d0
  C6=396.0d0
  R0=4.22d0
  !
  CASE ('I')
  alpha=35.000000d0
  C6=385.000000d0
  R0=4.170000d0
  !
  CASE ('Xe')
  alpha=27.300000d0
  C6=285.900000d0
  R0=4.080000d0
  !
  CASE ('Cs')
  alpha=427.12d0
  C6=6582.08d0
  R0=3.78d0
  !
  CASE ('Ba')
  alpha=275.0d0
  C6=5727.0d0
  R0=4.77d0
  !
  CASE ('Hf')
  alpha=99.52d0
  C6=1274.8d0
  R0=4.21d0
  !
  CASE ('Ta')
  alpha=82.53d0
  C6=1019.92d0
  R0=4.15d0
  !
  CASE ('W')
  alpha=71.041d0
  C6=847.93d0
  R0=4.08d0
  !
  CASE ('Re')
  alpha=63.04d0
  C6=710.2d0
  R0=4.02d0
  !
  CASE ('Os')
  alpha=55.055d0
  C6=596.67d0
  R0=3.84d0
  !
  CASE ('Ir')
  alpha=42.51d0
  C6=359.1d0
  R0=4.00d0
  !
  CASE ('Pt')
  alpha=39.68d0
  C6=347.1d0
  R0=3.92d0
  !
  CASE ('Au')
  alpha=36.5d0
  C6=298.0d0
  R0=3.86d0
  !
  CASE ('Hg')
  alpha=33.9d0
  C6=392.0d0
  R0=3.98d0
  !
  CASE ('Tl')
  alpha=69.92d0
  C6=717.44d0
  R0=3.91d0
  !
  CASE ('Pb')
  alpha=61.8d0
  C6=697.0d0
  R0=4.31d0
  !
  CASE ('Bi')
  alpha=49.02d0
  C6=571.0d0
  R0=4.32d0
  !
  CASE ('Po')
  alpha=45.013d0
  C6=530.92d0
  R0=4.097d0
  !
  CASE ('At')
  alpha=38.93d0
  C6=457.53d0
  R0=4.07d0
  !
  CASE ('Rn')
  alpha=33.54d0
  C6=390.63d0
  R0=4.23d0
  !
  CASE DEFAULT
  !
  !
  END SELECT
  !
  RETURN 
  END

