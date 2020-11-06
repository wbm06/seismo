parameter(MXPOW=16,NSAC=2**MXPOW, NWIN=8192, NSTAT=4096)
parameter(NSTAT1=NSTAT*NSTAT, NSTAT2=2*NSTAT1+NSTAT)
parameter(NDAT=NSTAT1/2+NSTAT)
parameter(MXBANDS=7)     ! 6 bands plus broadband
parameter(NZ=NSTAT*MXBANDS)
