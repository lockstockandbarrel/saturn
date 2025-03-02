program main
use Scientific_Subroutine_Package, only : &
absnt,  acfi,    ahi,    ali,    apch,   apfs,   apll,   apmm,   arat,   array,  &
ateig,  atse,    atsg,   atsm,   auto,   avcal,  avdat,  bdtr,   besi,   besj,   &
        bound,   cadd,           ccpy,   ccut,   cdtr,   cel1,   cel2,   chisq,  &
cint,   cnp,     cnps,   convt,  corre,  cross,  cs,     csp,    csps,   csrt,   &
csum,   ctab,    ctie,   dacfi,  dahi,   dali,   dapch,  dapfs,  dapll,  darat,  &
datse,  datsg,   datsm,  dbar,   dcar,   dcel1,  dcel2,  dcla,   dcnp,   dcnps,  &
dcpy,   dcsp,    dcsps,  ddbar,  ddcar,  ddet3,  ddet5,  ddgt3,  deli1,  deli2,  &
det3,   det5,    dfmcg,  dfmfp,  dfrat,  dgelb,  dgelg,  dgels,  dgt3,   dharm,  &
dhep,   dheps,   dhpcg,  dhpcl,  discr,  djelf,  dlap,   dlaps,  dlbvp,  dlep,   &
dleps,  dlgam,   dllsq,  dmatx,  dmchb,  dmfgr,  dmfsd,  dmfss,  dminv,  dmlss,  &
dmprc,  dmtds,   dpecn,  dpecs,  dpqfb,  dprbm,  dprqd,  dqa12,  dqa16,  dqa24,  &
dqa32,  dqa4,    dqa8,   dqatr,  dqg12,  dqg16,  dqg24,  dqg32,  dqg4,   dqg8,   &
dqh16,  dqh24,   dqh32,  dqh48,  dqh64,  dqh8,   dqhe,   dqhfe,  dqhfg,  dqhg,   &
dqhse,  dqhsg,   dql12,  dql16,  dql24,  dql32,  dql4,   dql8,   dqsf,   dqtfe,  &
dqtfg,  drharm,  drkgs,  drtmi,  drtni,  drtwi,  dse13,  dse15,  dse35,  dsg13,  &
dsinv,  dtcnp,   dtcsp,  dteas,  dteul,  dthep,  dtlap,  dtlep,  eigen,  eli1,   &
eli2,   expi,    exsmo,  factr,  fmcg,   fmfp,   forif,  forit,  frat,   gauss,  &
gdata,  gelb,    gelg,   gels,   gmadd,  gmmma,  gmprd,  gmsub,  gmtra,  gtprd,  &
harm,   hep,     heps,   hpcg,   hpcl,   hsbg,   inue,   i0,     jelf,   kolm2,  &
kolmo,  krank,   lap,    laps,   lbvp,   lep,    leps,   llsq,   load,   loc,    &
madd,   mata,    mchb,   mcpy,   meanq,  mfgr,   mfsd,   mfss,   mfun,   minv,   &
misr,   mlss,    momen,  mpair,  mprc,   mprd,   mstr,   msub,   mtds,   mtra,   &
multr,  ndtr,    ndtri,  nroot,  order,  padd,   paddm,  pcla,   pcld,   pder,   &
pdiv,   pecn,    pecs,   perm,   pgcd,   phi,    pild,   pint,   pmpy,   pnorm,  &
point,  polrt,   pprcn,  pqfb,   pqsd,   prbm,   probt,  prqd,   psub,   pval,   &
pvsub,  qa10,    qa2,    qa3,    qa4,    qa5,    qa6,    qa7,    qa8,    qa9,    &
qatr,   qg10,    qg2,    qg3,    qg4,    qg5,    qg6,    qg7,    qg8,    qg9,    &
qh10,   qh2,     qh3,    qh4,    qh5,    qh6,    qh7,    qh8,    qh9,    qhfe,   &
qhfg,   qhse,    qhsg,   ql10,   ql2,    ql3,    ql4,    ql5,    ql6,    ql7,    &
ql8,    ql9,     qsf,    qtest,  qtfe,   qtfg,   radd,   randu,  rank,   rcpy,   &
rcut,   recp,    rharm,  rint,   rk1,    rk2,    rkgs,   rslmc,  rsrt,   rsum,   &
rtab,   rtie,    rtmi,   rtni,   rtwi,   sadd,   scla,   scma,   sdiv,   se13,   &
se15,   se35,    sg13,   sici,   signt,  simq,   sinv,   smirn,  smo,    smpy,   &
srank,  srate,   srma,   ssub,   stprg,  submx,  subst,  tab1,   tab2,   tally,  &
tcnp,   tcsp,    teas,   tetra,  teul,   thep,   tie,    tlap,   tlep,   tprd,   &
trace,  ttest,   twoav,  utest,  varmx,  wtest,  xcpy
implicit none
external canor, biser, besk, besy, dapmm
  print *, "hello from project Scientific-Subroutine-Package"
end program main
