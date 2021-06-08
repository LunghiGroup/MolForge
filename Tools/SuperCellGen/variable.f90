MODULE variables
implicit none

character*30                                      :: input='inp',output='out',wfn_input='WFN.inp',           &
                                                     wfn_output='wfn.out',s2_output='eigen',ensemble='NVT',  &
                                                     type_dyn='CLASSICAL',kind_thermos='D_THERMO',           &
                                                     startb='MAG'

logical                                           :: &
 fulldiag=.true.,    &
 nodiag=.false.,     &
 subdiag=.false.,    &
 gsa=.false.,        &
 mix=.true.,         &
 gsa_state=.false.,  &
 SMIXING=.true. ,    &
 tennant=.true.,     &
 eso=.false.,        &
 nso=.false.,        &
 ref_D=.false.,      &
 ref_mol=.false.,     &
 ref_custom=.false., &
 b_mag=.false.,      &
 skip_diag=.false.,  &
 pGSH=.true.,        &
 pMSH=.true.,        &
 do_tao=.false.,     &
 hyst=.false.,       &
 properties=.false., &
 dynamo=.false.,     &
 start_O_rot=.false., &
 start_D_rot=.false., &
 noncolin=.false.,   &
 SPINSPIN=.false.,   &
 adapt_step=.false., &
 damp_mol=.false.  


END MODULE variables
