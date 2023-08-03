        module parameters_class
        use kind_class
        implicit none

        real(kind=dbl), parameter       :: PI=4.D0*DATAN(1.D0)
        real(kind=dbl), parameter       :: A_to_B= 1.8897259886d0
        real(kind=dbl), parameter       :: Har_to_Kc=627.503d0
        real(kind=dbl), parameter       :: F_conv=Har_to_Kc*A_to_B
        real(kind=dbl), parameter       :: boltz=3.11811E-06
        real(kind=dbl), parameter       :: amu_to_emass=1822.89d0

        end module parameters_class
