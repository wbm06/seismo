! compile: gfortran -O2 -o crosscor crosscor.f90 azidl.f /usr/local/sac/lib/sacio
! or on my Mac: g7 crosscor azidl -O2 sac bin

! in crc.inc:
! parameter(MXPOW=16,NSAC=2**MXPOW, NWIN=4096, NSTAT=1024)
! parameter(NSTAT1=NSTAT*NSTAT, NSTAT2=2*NSTAT1+NSTAT)
! parameter(NDAT=NSTAT1/2+NSTAT)
! parameter(MXBANDS=7)
! parameter(NZ=NSTAT*MXBANDS)

! Modifications January 2013:
! Allow for tbias (systematic offset in predicted arrival time)
! Read in event catalogue file name (last line in crcopt)
! Suppress ouptut of delay estimates if window < pmax
! Suppress output of delay estimates if error to large
! Fixed option kplot, it now behaves as advertized
! Fixed GMT plot error if kscale=0
! Fixed multiple writes to event catalogue file
! Added jdebug option
! Lowered spm to 0.1 for automatic determination of highest frequency

! June 2013
! added amplitudes, added arbitrary ray option 

! STILL TODO
! [3] improve noise estimates
! [4] add more complicated source option
! [5] verify degrees of freedom in computation of sigma
! [6] adapt datafile for matrix3D (now only for raydata).
! [7] add differential time measurements

! use testP.f90 to create synthetics for testing

program crosscor


! the program performs crosscorrelations among N stations
! and solves the VanDecar-Crosson equations.

! Uses GMT to produce plots during the run that allow you to refine
! the cross-correlation window (e.g. to avoid inclusion of later
! arrivals) and that can be edited for publication.
! Download GMT from http://gmt.soest.hawaii.edu/
! You need to <q>uit the GMT plot for the program to continue.

!     typical use:

!     echo 1 > incrc    [gives cluster number]
!     ls *.BHZ >> incrc	[creates input filename list]
!     vi incrc		[move low-noise signal to first place in list]
!     crosscor 

! If you intend you use clustering, program clustersac creates input
! files that can be copied to incrc (cluster001, cluster002,...)

! the various parameters can be read from screen the first time you
! run crosscor. They are saved in a file named <newcrcopt> if you do that.
! edit this file and rename it <crcopt> to be able to run crosscor
! with minimal screen input.
! If you wish to avoid all screen input set the plot option to 0, however
! be aware that visual inspection is often needed to see if the choice of
! window and tapers was a good one. So if you run on automatic, be sure
! to save the plotfiles (if you use more than one cluster the plots for
! previous clusters are overwritten if you don't!)

! where inputfile contains a list of SAC file names, one
! for each seismogram to be stacked. The SAC headers need
! to have a (theoretically determined) time pointer (A for P or S, T1 for
! PKP,... etc for other phases)
! that is used to feed the predicted travel times to crosscor.

! Make sure the first station is a low-noise station, since this is the
! one used to determine de polarity of all other seismograms.

! INPUT FILES:

! There are two input files, one that is event specific and that contains
! the SAC file names (incrc), the other that can be used for all events
! that sets the preferred parameters (crcopt)

! [1] incrc: a file with 
!     1- the cluster number (usually 1, if no clustering was applied)
!     2ff- the SAC files names of the seismograms to
!     be processed, one on each line

!     WARNING: change in input on line 8 since Oct 13, 2012, see note below
! [2] crcopt: a file with cross-correlation options with on line:
!     1- epsilon, damping parameter to center delays around IASP91
!     2- outlier level for rejection of seismograms (e.g. 10 sec)
!     3- epicentral distance range (rejects files with other distances)
!     4- source azimuth range (rejects files with other azimuth)
!     5- minimum allowable correlation coefficient (e.g. 0.8)
!     6- start and end taper (e.g. 3,3 sec)
!     7- highest period to use for filters, number of filter bands
!     8- sampling interval dtall to use for decimating SAC files to same dt

!     9-|kplot|: 0=no plot, 1=plot bb only, 2=all. If -1 or -2: do not ask 
!         for new duration after bb plot with stacked pulse
!     10- sigma0, sigmaN: lower limit for delay time uncertainty (sec), bias
!         because of small number N of stations (decreases as 1/N)
!     11- File name fr event catalogue (max 120 characters)

! An example of a crcopt file:
! .1            damping towards background model
! 5.0           outlier level (sec)
! 30.0    95.0  limits in epicentral distance
! -180  180     limits in source azimuth
! 0.85          Rmin=minimum coef crosscorrelations
! 3.00    3.00  start and end taper in s
! 32.0   6      largest passband period, number of bands
! 0.10          common sampling interval (s)
! 2             kplot (1=BB only, 0=none, 2=all)
! 0.2    5      water level for error, systematic bias decreases as 5/N
! ../Displacements/eventfile

! Before you do a first run, create an empty event file, crosscor will
! fill it in as you analyze more events

! INPUT FROM SCREEN

! 1: ident (max 16 characters)
! 2: phase name (e.g. P)
! 3: the name of the SAC variable with the time pick (e.g. A)
! 4: second phase (if differential measurement), else give None
! 4a (only if second phase was given): SAC header variable with pick
! 4b (only if second phase was given): ident for this phase
! 5: diagnostoc option (1=write file, 0=not)
! 6: tbias, tsigma, duration (see below)
! 7: phase to limit time window length (e.g. pP)
! 8 (only if kplot>0): revised window length 

! ident is used to name an output file Data.ident that can be used
! with raydata.f for subsequent finite frequency tomography
! The SAC file needs to have at least the time pick for the target phase
! (use plotpk in SAC)
! tbias is usually 0, but may be useful if you filled the SAC header with
! theoretical arrival times and you find they are all x seconds late or
! early. Sign convention: tbias is negative is the prediction is late.
! tsigma is the uncertainty in the SAC time picks. The correlation window
! is extended by tsigma to avoid that it starts too late or ends to early.
! duration is the duration of the phase in the unfiltered (broadband)
! seismogram. If you set kplot>0 in crcopt, you will be allowed to change
! this duration once you inspect a stack.


! A note on the choice of windows:
!=================================
! An earlier version used w1,w2 to define the signal window. We now define them as
! Mercerat & Nolet (GJI 2012) who use the following w1 and w2:
! pre-arrival time window w1 = sigma + taper
! post w2 = sigma + broadband duration + maximum period + taper
! where sigma is the arrival time uncertainty
! The window length for the broadband pulse is significantly reduced with respect
! to the older version. xcl() is now more correctly defined as the signal length
! before (rather than after) filtering.

include 'crc.inc'

! general purpose array
dimension d(NSAC) 

! selected seismgram windows, unfiltered and filtered
dimension s(NWIN,NSTAT),st1(NSTAT),nsn(NSTAT),stv(NSTAT)
dimension sf(5*NWIN,NSTAT),sft1(NSTAT),nsfn(NSTAT)

! data arrays
dimension ncor(NSTAT),sigma(NSTAT),polarity(NSTAT)
dimension azev(NSTAT),dist(NSTAT),stlat(NSTAT),stlon(NSTAT)
dimension stel(NSTAT)
dimension tobs(NSTAT,MXBANDS),sigt(NSTAT,MXBANDS),avcor(NSTAT,MXBANDS)
dimension aobs(NSTAT,MXBANDS),siga(NSTAT,MXBANDS)
dimension xcl(MXBANDS),rmsb(MXBANDS)

! used in inversion
dimension x(NSTAT),u(NDAT),v(NSTAT),w(NSTAT),b(NSTAT1),xa(NSTAT)
common /LSQ/ nrow,aa(NSTAT2),ia(NSTAT1),ja(NSTAT2)
dimension rhs(NSTAT1),rhsa(NSTAT1)

! event catalogue test
logical ruthere

! GMT stuff
character*18 rgb(14)

! For debugging only
common /debug/ bugfile

! Character variables
character*8 sacvar,ksacvar,fbug,bugfile
character*16 ident,ident2
character*8 phase,chan,kvar,kcmpnm,knetwk,knet(NSTAT),kcmp(NSTAT)
character*8 kstnm,code,master,wipe,kst(NSTAT),kvar2,phase2,kvar1
character*20 sacfil,gmtfil
character*30 command
character*120 eventfile,phase1file,phase2file,raydefile

data rgb/'> -W1p/120/120/120', '> -W1p/0/255/0', '> -W1p/0/0/255', &
'> -W1p/255/0/255', '> -W1p/255/0/0','> -W1p/0/0/0', &
'> -W1p/155/155/0', '> -W1p/0/255/255', '> -W1p/0/0/120',&
'> -W1p/0/120/0','> -W1p/50/120/220','> -W1p/120/0/0',&
'> -W1p/120/50/120', '> -W1p/120/120/50'/

data avcor/NZ*0./
data tobs/NZ*1.0e10/
data aobs/NZ*1.0e10/
data rmsb/MXBANDS*0./                   ! TODO (but not used elsewhere)
data pi/3.14159265/

d2r=3.14159265/180.0            ! degree to radian conversion

! the following is only for debugging
jdebug=0          ! write bug info to fort.13 if jdebug>0
bugfile='none'    ! no sac file is written by crc if bugfile='none'

print *,'Reading file <incrc> with SAC file names'

open(2,file='incrc',iostat=ios,status='old')
if(ios /= 0) then
  print *,'ERROR: Cannot open file incrc'
  print *,'This file should have all names of SAC files to be processed,'
  print *,'one on each line.'
  print *,'First line should have the cluster number (1 if all in one)'
  stop
endif

print *,'Readinf file <crcopt> with cross-correlation options'

in8=8
open(in8,file='crcopt',iostat=ios,status='old')
if(ios /= 0) then
  print *,'Cannot open file crcopt'
  print *,'Will read options from the screen instead'
  print *,'note: your options are saved in newcrcopt'
  in8=5
  close(8)
  open(8,file='newcrcopt',iostat=ios,status='old')
endif

print *,'Give run ident (max 16 characters):'
read(5,fmt='(a)') ident
open(4,file='out.crosscor.'//ident)
write(4,*) 'Run of crosscor, ident: ',ident
write(4,*)

! Open the data file later to be read by raydata
! call that copies the ray definition file to this data file.
open(11,file='Data.'//ident,iostat=ios,status='new')
if(ios /= 0) then
  print *,'A data file named ','Data.'//ident,' already exists'
  if(ident(1:4).ne.'test') then
    stop 'Remove old data file'
  else
    print *,'We assume it is a test file only and we shall overwrite it'
    open(11,file='Data.'//ident)
  endif
endif

print *,'Give phase ident (e.g. P):'
read(5,'(a)') phase
print *,'Give SAC variable  with predicted time for ',phase
print *,'(e.g. A for phase P or S):'
read(5,fmt='(a)') kvar
phase1file=trim(phase)//'Data.'//ident
print *,'If differential times, give reference phase (e.g. P), else <None>:'
read(5,'(a)') phase2
if(phase2.ne.'None') stop 'Differential times not yet implemented'
if(phase2.ne.'None') then
  print *,'Give SAC variable with predicted time for ',phase2
  read(5,fmt='(a)') kvar2
  print *,'Give ident for reference phase (e.g. xxxx for Data.xxxx):'
  read(5,'(a)') ident2
  phase2file=trim(phase2)//'Data.'//ident2
else
  phase2file='None'
endif  
print *,'Write diagnostic file (1,slow) or not (0,fast)?'
read(5,*) kdiag
if(kdiag>0) then
  open(7,file='diagnostics.crosscor.'//ident)  ! To debug data
  write(7,*) 'Run of crosscor, ident: ',ident
  write(7,*)
endif
  


! first few lines conform to data files for global tomography
! todo: adapt for matrix3D
write(11,fmt='(a)') phase1file
write(11,fmt='(a)') phase2file

! Try to find ray definition file, give preference to local definition
raydefile='absent'
! if on vif:
inquire(file='/Users/guust/data/Tables/raydefs/'//trim(phase)//'.def',exist=ruthere)
if(ruthere) raydefile='/Users/guust/data/Tables/raydefs/'//trim(phase)//'.def'
! if on my Mac:
inquire(file='/opt/linux64/kernel2.6.32/soft_sismo/GLOBALSEIS/tables/raydefs/'//trim(phase)//'.def',exist=ruthere)
if(ruthere) raydefile='/opt/linux64/kernel2.6.32/soft_sismo/GLOBALSEIS/tables/raydefs/'//trim(phase)//'.def'
! but choose local def file if available:
inquire(file=trim(phase)//'.def',exist=ruthere)
if(ruthere) raydefile=trim(phase)//'.def'
if(raydefile.eq.'absent') then
  print *,'Cannot find ray definition file: ',trim(phase)//'.def'
  stop
else
  print *,'Opening ',raydefile
  open(3,file=raydefile)
  read(3,'(a)') phase2
  read(3,'(a)') phase2
endif  

kdown=0
ktype=1
do while(kdown.ne.5)
  lastp=ktype
  read(3,*) r,kdown,ktype
  write(11,fmt='(f10.0,2i3)') r,kdown,ktype
enddo  

print *,kvar,' in SAC header has predicted time. tbias allows you to'
print *,'correct wrongly estimated time, eg because of wrong origin time.'
print *,'Sign convention: if prediction is late set time bias negative.'
print *
print *,'Give tbias, uncertainty sigma and duration of broadband pulse (in sec):'
read(5,*) tbias,tsigma,duration
if(kdiag>0) write(7,*) ,'tbias,sigma, duration=',tbias,tsigma,duration
write(4,*) ,'tbias,sigma, duration=',tbias,tsigma,duration
print *,' '

print *,'You can limit the window to exclude next phase if this is'
print *,'in the header, or type <None> to ignore:'
read(5,'(a)') kvar1
if(kdiag>0) write(7,*) 'Window limited by ',kvar1
write(4,*) 'Window limited by ',kvar1
w2base=w2       ! w2 without kvar1 limit

open(15,file='spec.'//trim(ident)//'.xy')       ! GMT file for spectra

! open(3,file='diagnostics.lsqr.'//ident)      ! test; now obsolete
open(12,file='newincrc')

print *,' '
print *,'Give epsilon (0<=eps) for damping towards IASP91 (e.g. 0.1):'
read(in8,*) epsilon
if(in8.eq.5) write(8,*) epsilon
print *,'Give maximum acceptable arrival time anomaly (e.g. 10 s):'
read(in8,*) outlier
if(in8.eq.5) write(8,*) outlier
print *,' '
print *,'Give range for epicentral distances (in degrees, eg 30,95):'
print *,' '
read(in8,*) gc1,gc2
if(in8.eq.5) write(8,*) gc1,gc2
if(kdiag>0) write(7,*) 'Epicentral distances from',gc1,' to ',gc2
write(4,*) 'Epicentral distances from',gc1,' to ',gc2
print *,'Give range for source azimuth (-180 to 180, clockwise):'
print *,' '
read(in8,*) az1,az2
if(abs(az1).gt.180.0.or.abs(az2).gt.180.0) stop 'Azimuth not in [-180,180]'
if(in8.eq.5) write(8,*) az1,az2
write(4,*) 'Source azimuth from',az1,' to ',az2
print *,'Give minimum acceptable correlation coefficient (e.g. 0.85):'
read(in8,*) cmin
if(in8.eq.5) write(8,*) cmin
if(kdiag>0) write(7,*) 'Minimum acceptable correlation coefficient:',cmin
write(4,*) 'Minimum acceptable correlation coefficient:',cmin
print *,' '
print *,'Give length of tapers iat start and end of  window (eg 3s,3s):'
read(in8,*) tap1,tap2
if(in8.eq.5) write(8,*) tap1,tap2
if(kdiag>0) write(7,*) 'tap1,tap2=',tap1,tap2,' s'
write(4,*) 'tap1,tap2=',tap1,tap2,' s'
print *,'Give largest period for filters and nr of filters (give nr=0'
print *,'to let the program decide [this is the preferred option!]:'
read(in8,*) pmax,nrbands0
if(in8.eq.5) write(8,*) pmax,nrbands0
if(nrbands0>0) then
  if(kdiag>0) write(7,*) nrbands0,' down from a period of ',pmax,' s'
else
  if(kdiag>0) write(7,*) 'Number of bands determined automatically'
  write(4,*) 'Number of bands determined automatically'
  write(6,*) 'Number of bands determined automatically'
endif  
w1=tsigma+tap1                  ! start window at Tpred-w1
w2=tsigma+tap2+duration+pmax    ! and end at Tpred+w2
print *,'Give dtall (must be 2^n of sampling period in SAC files):'
read(in8,*) dtall
if(in8.eq.5) write(8,*) dtall
if(kdiag>0) write(7,*) 'SAC files decimated to ',dtall,' s sampling interval'
write(4,*) 'SAC files decimated to ',dtall,' s sampling interval'
flowpass=0.8/dtall
print *,'Give plot option: -1=do not ask , 0=none, 1=bb only, 2=all'
read(in8,*) kplot
if(in8.eq.5) write(8,*) kplot,'          kplot'
print *,'The computed standard error can be restricted to a lower'
print *,'bound, and a regularization bias sigmaN/N an be added for'
print *,'a cluster of N stations.'
print *,'Give minimum standard error and sigmaN (eg 0.2, 5):'
read(in8,*) sigma0,sigmaN
if(in8.eq.5) write(8,*) sigma0,sigmaN,'      sigma0,sigmaN'
if(kdiag>0) write(7,*) 'sigma0,sigmaN=',sigma0,sigmaN
write(4,*) 'sigma0,sigmaN=',sigma0,sigmaN
print *,'Give event file name (max 120 char):'
read(in8,"(a)") eventfile
if(in8.eq.5) write(8,"(a)") eventfile
if(kdiag>0) write(7,*) 'eventfile=',trim(eventfile)
write(4,*) 'eventfile=',trim(eventfile)

! Read the cluster number from the seismogram list
read(2,*) kluster
write(4,*) 'Cluster=',kluster



kb=0

!========================================================
! Read the seismograms, into d, decimate to common sampling,
! window to include twice the longest period, store in s
! as broadband signals for later filtering.
! Time is redefined: time is with respect to predicted 
! arrival time of the phase of interest
!========================================================

ns=0    ! counts number of seismograms (stations)
nsnmax=0        ! largest signal length
if(kdiag>0) write(7,*)
if(kdiag>0) write(7,fmt='(2a)') 'SAC file name            beg     del    tvar',  &
  '       e      t1      t2     average'

nowrite=0
ios=0
do while(ios == 0)      ! loop over SAC files

  read(2,fmt='(a)',iostat=ios) sacfil
  if(jdebug>0) write(13,*) 'ios=',ios,', reading ',trim(sacfil)
  if(sacfil(1:4).eq.'stop') exit
  print *,'Next file ',trim(sacfil)
  if(ns > 1 .and. ios /= 0) exit          ! if last line is not 'stop'

  call rsac1(sacfil,d,nsmp,beg,del,NSAC,nerr)
  if(nerr==-803) nerr=0           ! (ignore: sac file size large)
  if(nerr /= 0) then
    print *,'SAC file read error',nerr
    if(jdebug>0) write(13,*) 'SAC file read error',nerr
    stop
  endif

  ! get SAC header variables
  kerr=0
  ! call getfhv('gcarc',gcarc,nerr)
  ! call getfhv('az',az,nerr)
  call getfhv(kvar,tvar,nerr)    
  if(nerr /= 0 .or. abs(tvar)>9999) then
    if(kdiag>0) write(7,*) trim(sacfil),': ',trim(kvar),' marker absent, dist=',&
       gcarc,' -ignored'
    cycle
  endif
  tvar=tvar+tbias                        ! correct predicted arrival time
  tvar2=99999.                           ! if no second phase limit
  if(kvar1.ne.'None') then  
    call getfhv(kvar1,tvar2,nerr)    
    if(nerr /= 0 .or. abs(tvar2)>9999) then
      tvar2=99999.
    endif
    tvar2=tvar2+tbias                    ! correct predicted arrival time
  endif  
  call getnhv('nzyear',nzyear,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getnhv('nzjday',nzjday,nerr)
  if(nerr /= 0) kerr=kerr+1
  idate=1000*mod(nzyear,100)+nzjday
  call getnhv("nzhour",nzhour,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getnhv("nzmin",nzmin,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getnhv("nzsec",nzsec,nerr)
  if(nerr /= 0) kerr=kerr+1
  iotime=10000*nzhour+100*nzmin+nzsec
  call getfhv('evla',evla,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getfhv('evlo',evlo,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getfhv('evdp',evdp,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getfhv('stlo',stlo,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getfhv('stla',stla,nerr)
  if(nerr /= 0) kerr=kerr+1
  call getfhv('stel',elev,nerr)
  if(nerr /= 0) kerr=kerr+1
  elev=0.001*elev        ! station elevation is in km
  call azidl(stla*d2r,stlo*d2r,evla*d2r,evlo*d2r,delta,dels,azis,azie)
  gcarc=delta/d2r
  az=azie/d2r           ! azimuth as seen from event (to station)
  call getkhv('kstnm',kstnm,nerr)
  call getkhv('kcmpnm',kcmpnm,nerr)
  if(nerr.ne.0) kcmpnm='UNK'
  call getkhv('knetwk',knetwk,nerr)
  if(nerr.ne.0) knetwk='UNK'
  e=beg+(nsmp-1)*del
  call clear(kstnm)
  call clear(kcmpnm)
  call clear(knetwk)
  if(kerr > 0) then
    if(kdiag>0) write(7,*) sacfil,': header incomplete, kerr=',kerr,' -ignored'
    cycle
  endif

  if(ns.eq.0) then  
    ! assign event number. For now, we assume eventcatalogue
    ! is either in /Users/guust/data/seis/ or in 
    ! /usr/local/soft_sismo/GLOBALSEIS/crceventcatalogue
    inquire(file=eventfile,exist=ruthere)
    if(ruthere) then
      open(14,file=eventfile)
    else
      print *,'Cannot find event catalogue:'
      print *,trim(eventfile)
      stop 'eventcatalogue not present'
    endif  
    ios2=0
    id=0
    it=0
    line=0
    do while(ios2.eq.0.and.nowrite.eq.0)  ! loop over events in existing catalogue
      line=line+1
      read(14,*,end=400,iostat=ios2) ievt,id,it,slat,slon,sdep
      if(line.eq.1) ios2=0   ! bricolage if first ios2 != 0 for no reason
      if(jdebug>0) write(13,*) line,ievt,id,it,' vs ',idate,iotime
      if(id.eq.idate.and.abs(it-iotime).lt.10) then
        print *,'Found an entry for this event already in the catalogue:'
        print *,'line:',line,', date: ',id,'=',idate,'i,  time:',it,'=',iotime
        nowrite=1
        if(jdebug>0) write(13,*) 'event match, nowrite=',nowrite
        goto 400
      endif  
    enddo
    400 close(14)
    if(nowrite.eq.0) ievt=line
    idate0=idate
    iotime0=iotime
  else
    if(idate.ne.idate0.or.iotime.ne.iotime0) then
      print *,'For ns=',ns,' (',trim(sacfil),') event differs from'
      print *,'date or time of first file:',idate,idate0,iotime,iotime0
      stop 'More than one event in data'
    endif  
  endif

  ! decimate all signals to dtall sampling interval

  ndecim=0
  do while(del*2**ndecim.lt.dtall-0.00001)
    ndecim=ndecim+1
  enddo
  if(abs(del*2**ndecim - dtall) > 0.0001) then
    print *,'del=',del,' but del*2**ndecim=',del*2**ndecim,' /=',dtall
    print *,'in file ',sacfil
    stop 'Fatal error because SAC file sampling rate not 2^n multiple of 10Hz'
  endif  
  !call wsac1('d1.sac',d,nsmp,beg,del,nerr)             ! debug
  call decimate(d,ndecim,nsmp,beg,del,npow,flowpass)
  !call wsac1('d2.sac',d,nsmp,beg,dtall,nerr)           ! debug

  ns=ns+1
  if(ns>NSTAT) stop 'Too many stations! Increase NSTAT in crc.inc'
  stv(ns)=tvar
  kst(ns)=trim(adjustl(kstnm))
  knet(ns)=trim(adjustl(knetwk))
  kcmp(ns)=trim(adjustl(kcmpnm))
  azev(ns)=az
  dist(ns)= gcarc
  stlat(ns)=stla
  stlon(ns)=stlo
  stel(ns)=elev

  ! all times are redefined with respect to origin time To
  ! t1 and t2 determine the window selection from this seismogram
  ! the interval (t1,t2) may still contain extraneous later arrivals 
  ! like pP that are removed by windowing more precisely later

  t1=tvar-w1
  t2=min(tvar+w2,tvar2)
  call rmean(d,nsmp,average)    ! remove mean offset of signal

  if(kdiag>0) write(7,fmt='(a20,f8.1,f8.3,4f8.1,e12.3)') sacfil,beg,del,tvar,  &
  e,t1,t2,average

  if(t1<beg) then
    if(kdiag>0) write(7,*) sacfil,': starts too late'
    if(kdiag>0) write(7,*) 'window starts at',t1,' but seismogram at',beg
    ns=ns-1
    cycle
  endif
  if(t2>beg+(nsmp-1)*del) then
    if(kdiag>0) write(7,*) sacfil,': stops too early'
    if(kdiag>0) write(7,*) 'Window ends at',t2,' but seismogram at',beg+(nsmp-1)*del
    ns=ns-1
    cycle
  endif


  zero=0.0
  call window(d,beg,del,nsmp,t1,t2,tap1,tap2)   ! window d to (t1,t2)
  xcl(1)=duration+tsigma+0.5*tap2     ! length of broadband window for t>0

  if(nsmp.gt.NWIN) then
    if(kdiag>0) write(7,*) sacfil,' has window with ',nsmp,' samples >',NWIN
    ns=ns-1
    cycle
  endif

  nsn(ns)=nsmp                  ! length of windowed seismograms
  nsnmax=max(nsmp,nsnmax)
  st1(ns)=beg-tvar              ! redefine time axis: arrival is at t=0
  do i=1,nsmp
    s(i,ns)=d(i)
  enddo  

! tests
  if(abs(del-dtall)>0.001) then
    print *,'dt=',del,' in station ',kstnm
    print *,'Expected del=',dtall
    stop 'differing sampling intervals'
  endif  
! if(nzyear /= nzyear1 .or. nzjday /= nzjday1) stop 'date differs'
  if(gcarc<gc1 .or. gcarc > gc2 ) then
    if(kdiag>0) write(7,*) sacfil,': Out of range, dist=',gcarc
    ns=ns-1
    cycle
  endif
  ! check azimuth. az,az1 and az2 are in [-180,180]
  if(az2.lt.az1) then   ! if azimuth region brackets South
    if(az<0.) az=az+360.
    if(az<az1 .or. az>az2+360. ) then
      if(kdiag>0) write(7,*) sacfil,': Out of azimuth, az=',az
      ns=ns-1
      cycle
    endif  
  else if(az<az1 .or. az > az2 ) then
    if(kdiag>0) write(7,*) sacfil,': Out of azimuth, az=',az
    ns=ns-1
    cycle
  endif

enddo   ! end of loop over SAC files

! write new line to eventcatalogue with number of stations
if(nowrite.eq.0) then
  if(jdebug>0) write(13,*) 'Adding event ',ievt,idate,iotime,nowrite
  open(14,file=eventfile,access='append')
  write(14,fmt='(i9,2i7,3f9.3)') ievt,idate,iotime,evla,evlo,evdp
  close(14)
endif

!==================================================
! STACK BROADBAND SIGNAL
!==================================================

iband=1                 ! broadband signal is first 'band'

! cross correlate broadband signal

nrow=0
kount=0
do j=1,ns
  ncor(j)=0
enddo  
if(kdiag>0) write(7,*) ns,' seismograms'
if(kdiag>0) write(7,*) 'Broadband correlations:'
if(kdiag>0) write(7,fmt='(a)') 'stat1   stat2      corcoef      tmax     delay  B    period (nominal)'

polarity(1)=1.0         ! align polarity with first station

! store window signals in sf(), do final error check
kerror=0
do i=1,ns  
  nsfn(i)=(w1+duration+tsigma+tap2)/dtall
  if(nsfn(i)<8) then
    print *,'ERROR: ',kst(i),' has only',nsfn(i),' <8 samples in window'
    kerror=kerror+1
  endif  
  sft1(i)=st1(i)
  do j=1,nsfn(i)       ! put signal in the middle
    f=1.0
    if(j*dtall>w1+duration+tsigma) f=0.5*cos(pi*(j*dtall-w1-duration-tsigma)/tap2)+0.5
    if(j*dtall>w1+duration+tsigma+tap2) f=0.
    sf(j,i)=f*s(j,i) 
  enddo  
enddo
if(kerror>0) stop 'Windowing error(s)'

do is=1,ns

  print *,'Correlating ',kst(is),' broadband, cpu=',second()

  ti=sft1(is)
  tvi=0.                ! expect arrival is at t=0

  ! add theoretical time to equations to ensure system is solvable
  ! even if clusters are disconnected, but use low weight epsilon
  nrow=nrow+1
  if(nrow > NSTAT1) stop 'Increase dimension of rhs'
  kount=kount+1
  if(kount >  NSTAT2) stop 'Increase dimenson of aa and ja'
  aa(kount)=epsilon
  rhs(nrow)=0.                
  rhsa(nrow)=0.
  ja(kount)=is
  ia(nrow)=1

  ! now correlate broadband pulse in station is with all stations js>is
  do js=is+1,ns
    tj=sft1(js)
    tvj=0.
    isgnflg=0
    taper=0.5*tap1
    call crc(sf(1,is),ti,tvi,nsfn(is),sf(1,js),tj,tvj,nsfn(js),dtall, &
      taper,outlier,tmax,ampl,rmax,isgn,isgnflg)
    if(is.eq.1) polarity(js)=isgn       ! conserve polarity of first station
    if(abs(rmax)<cmin) then
      if(kdiag>0) write(7,fmt='(2a8,f10.3,1x,a)') kst(is),kst(js),isgn*rmax,'        -'
      if(jdebug>0) write(13,*) 'js,is,rmax=',js,is,rmax,' rejected'
    else
      delay=tmax
      if(abs(delay).gt.outlier) cycle
      nrow=nrow+1
      if(nrow >  NSTAT1) stop 'Increase dimension of rhs'
      kount=kount+1
      if(kount >= NSTAT2) stop 'Increase dimenson of aa and ja'
      aa(kount)=-1.0
      ja(kount)=is
      kount=kount+1
      aa(kount)=+1.0
      ja(kount)=js
      rhs(nrow)=tmax-(tvi-tvj)        
      rhsa(nrow)=log(ampl)
      ia(nrow)=2
      ncor(is)=ncor(is)+1       ! count nr of acceptable pairs
      ncor(js)=ncor(js)+1
      avcor(is,iband)=avcor(is,iband)+rmax      ! average corr coefficient
      avcor(js,iband)=avcor(js,iband)+rmax
      if(kdiag>0) write(7,fmt='(2a8,f10.3,2f10.2,i3,f10.1)') kst(is),kst(js),isgn*rmax,tmax,delay,iband,0.
      if(jdebug>0) write(13,*) 'js,is,rmax=',js,is,rmax,' nrow=',nrow,ampl
    endif 
  enddo
enddo  

if(nrow<3) then
  if(kdiag>0) write(7,*) 'Bad data: not enough acceptable broadband correlations'
  write(4,*) 'Bad data: not enough acceptable broadband correlations'
  stop 'Bad data: not enough acceptable broadband correlations'
endif  
  
! compute average correlation coefficient in each station and find the best
ibest=0
abest=0.
do i=1,ns
  if(ncor(i)>0) avcor(i,1)=avcor(i,1)/ncor(i)
  if(avcor(i,1)>abest) then
    ibest=i
    abest=avcor(i,1)
  endif  
enddo  
write(4,fmt='(3a,f7.2)') 'Best seismogram: ',trim(kst(ibest)), &
  ', with average correlation coefficient=',abest

! solve system of equations for this broad band
print *,'Solving for broadband, nrow=',nrow
itmax=5*ns
do j=1,nrow
  u(j)=rhs(j)         ! since lsqr overwrites u, we do not use rhs itself
enddo  
call lsqr(nrow,ns,x,u,v,w,itmax,r,kdiag)      ! solves Ax=u
if(r<0.) then
  if(kdiag>0) write(7,*) 'Bad data: LSQR does not converge on the broadband correlations'
  write(4,*) 'Bad data: LSQR does not converge on the broadband correlations'
  stop 'LSQR does not converge on the broadband correlations'
endif  

! compute A*x=b to get predicted delay times for each station pair
call asol(x,b,nrow,ns)		! computes b=A*x

! estimate standard errors
kount=0
do i=1,ns
  sigma(i)=0.
enddo  
chi2=0
nchi=0
do irow=1,nrow
  if(ia(irow).eq.0) then
    cycle
  else if(ia(irow).eq.1) then
    kount=kount+1
    cycle
  else if(ia(irow).eq.2) then
    kount=kount+1
    nchi=nchi+1
    is=ja(kount)
    kount=kount+1
    js=ja(kount)
    dif2=(b(irow)-rhs(irow))**2
    sigma(is)=sigma(is)+dif2
    sigma(js)=sigma(js)+dif2
    chi2=chi2+dif2
  else
    print *,'Bug in row',irow,' ia=',ia(irow)
    stop
  endif  
enddo
do i=1,ns
  if(ncor(i).gt.3.and.polarity(i).ne.0.) then
    ! check if -2 is correct (comes from VanDecar&Crosson, 1990)
    sigma(i)=sqrt( sigma(i)/(ncor(i)-2) + (sigma0+sigmaN/ns)**2 )
  else
    sigma(i)=999.99
  endif  
  sigt(i,1)=sigma(i)
  tobs(i,1)=x(i)+stv(i)
enddo  

write(4,fmt='(//,a)') 'Broadband delays'
write(4,fmt='(a,f10.2,a,f10.2,a,f10.3,a)') 'Chi2=',chi2,', relative:', &
   chi2/nchi,' or rms=',sqrt(chi2/nchi),' sec.'
write(4,fmt='(/,2a)') 'station      T_arr    T_pred     delay     ', &
   'sigma      icor   avcor     azi   gcarc'
do j=1,ns
  if(ncor(j).gt.0) write(4,fmt='(a8,4f10.2,i10,3f8.2)') &
     kst(j),x(j)+stv(j),stv(j),x(j),sigma(j),ncor(j),avcor(j,iband), &
     azev(j),dist(j) 
enddo  

! solve same system of equations for amplitude ratios
print *,'Solving amplitudes for broadband, nrow=',nrow
itmax=5*ns
do j=1,nrow
  u(j)=rhsa(j)        ! careful - lsqr overwrites u
enddo  
! use a trick to get rid of the time regularisation
do i=1,kount
  if(abs(aa(i)-epsilon) < 0.5*epsilon) aa(i)=0.
enddo
call lsqr(nrow,ns,xa,u,v,w,itmax,r,kdiag)      ! solves Ax=u
if(r<0.) then
  if(kdiag>0) write(7,*) 'Warning: LSQR does not converge on the broadband amplitudes'
  write(4,*) 'Warning: LSQR does not converge on the broadband amplitudes'
  goto 500
endif  

! remove average, compute A*x=b to get predicted amplitudes 
! for each station pair
xav=0.
do i=1,ns
  xav=xav+xa(i)
enddo
xav=xav/ns
do i=1,ns
  xa(i)=xa(i)-xav
enddo
call asol(xa,b,nrow,ns)		! computes b=A*x

! estimate standard errors
kount=0
do i=1,ns
  sigma(i)=0.
enddo  
chi2=0
nchi=0
do irow=1,nrow
  if(ia(irow).eq.0) then
    cycle
  else if(ia(irow).eq.1) then
    kount=kount+1
    cycle
  else if(ia(irow).eq.2) then
    kount=kount+1
    nchi=nchi+1
    is=ja(kount)
    kount=kount+1
    js=ja(kount)
    dif2=(b(irow)-rhsa(irow))**2
    sigma(is)=sigma(is)+dif2
    sigma(js)=sigma(js)+dif2
    chi2=chi2+dif2
  else
    print *,'Bug in row',irow,' ia=',ia(irow)
    stop
  endif  
  if(jdebug>0) write(13,*) 'irow=',irow,', b,rhsa=',b(irow),rhsa(irow)
enddo
if(jdebug>0) write(13,*) 'chi2=',chi2
do i=1,ns
  if(ncor(i).gt.3.and.polarity(i).ne.0.) then
    ! check if -2 is correct (comes from VanDecar&Crosson, 1990)
    sigma(i)=sqrt( sigma(i)/(ncor(i)-2))
    siga(i,1)=0.1                 ! guess at error, needs improving
  else
    sigma(i)=999.99
    siga(i,1)=999.99
  endif  
  aobs(i,1)=xa(i)
enddo  

write(4,fmt='(//,a)') 'Broadband amplitudes (log)'
write(4,fmt='(a,f10.2,a,f10.3,a,f10.3,a)') 'Chi2=',chi2,', relative:', &
   chi2/nchi,' or rms=',sqrt(chi2/nchi)
write(4,fmt='(/,a)') 'station   lnA_obs     sigma      icor   avcor     azi   gcarc'
do j=1,ns
  if(ncor(j).gt.0) write(4,fmt='(a8,2f10.5,i10,3f8.2)') &
     kst(j),xa(j),sigma(j),ncor(j),avcor(j,iband),azev(j),dist(j) 
enddo  
500 continue

! align signals, stack and write GMT file
lenplot=min(99.1,max(10.1,3*pmax))        ! plot width in sec
lenhlf=lenplot/2
lenplot=2*lenhlf
tlag=-w1     ! start signal tlag s before predicted pulse
nstack=(w1+w2)/dtall+1
print *,'broadband nstack=',nstack
do i=1,nstack
  d(i)=0.
enddo  
open(9,file='gmtstack')
fourpt=4*2.54/72.               ! plot offset 4 points
lymax=max(ns*fourpt+1.0,10.01)    ! plot size (vertical)
lymax=min(28,lymax)             ! max page size 28 cm
call startgmt(9)
write(9,fmt='(a)') '# Usage: gmtstack stackx smgrsx x'
write(9,fmt='(a)') 'rm -fr $argv[1].eps '
write(9,fmt='(a,i3,"/",i2,a,i2,a,i2,a)') 'psbasemap -R',-lenhlf,  &
  int(lenhlf+duration),'/0.0/',lymax,' -JX20/',lymax,   &  
  ' -Ba5f1g5:"Delay(s)":/a5NS -P  -K > $argv[1].eps'
write(9,fmt='(a)') 'psxy -M $argv[2].xy -R -JX -W1p -K -O >> $argv[1].eps'
write(9,fmt='(a,f7.1,i3,3a)') 'echo ',-lenhlf+0.1,lymax-1,  &
  '.4 10 0 1 BL Stack band  $argv[3] | pstext -R -JX -K -O >> $argv[1].eps'
write(9,fmt='(a)') 'psxy $argv[1].xy -R -JX -W2p0/0/0 -O >> $argv[1].eps'
write(9,fmt='(a)') 'gv $argv[1].eps &'
close(9)
write(*,*) 'calling gmtstack if kplot>0'
if(kplot>0) call system('chmod a+x gmtstack')

dmax=0.
print *,'Stacking ns=',ns
do j=1,ns       ! stack stations
  if(ncor(j).le.2.or.polarity(j).eq.0.) cycle
  kdelay=(tlag+x(j)-st1(j))/dtall   ! index for start of pulse in s
  k2=min(nstack,nsn(j)-kdelay)
  k1=max(1,1-kdelay)
  do k=k1,k2
    d(k)=d(k)+polarity(j)*s(k+kdelay,j)      ! stack station j
    dmax=max(dmax,abs(d(k)))
  enddo  
  if(dmax.le.0.) then
    print *,'Bug? j,kdelay,k1,k2,polarity=',j,kdelay,k,k2,polarity(j)
  endif  
enddo  

! do a final check of polarity and correlation, but now with the
! stacked seismogram

do is=1,ns  
  if(sigt(is,1).gt.900.) cycle
  call crc(sf(1,is),ti,tvi,nsfn(is),d,tlag,0.,nstack,dtall, &
      taper,outlier,tmax,ampl,rmax,isgn,isgnflg)
  polarity(is)=isgn       
  if(isgn.eq.0.or.abs(rmax)<cmin) then
    sigt(is,1)=999.99
    siga(is,1)=999.99
    if(kdiag>0) write(7,fmt='(2a8,i3,f10.3,1x,a)') kst(is),kst(js),isgn,rmax,'   Reject'
    if(jdebug>0) write(13,*) 'js,is,rmax=',js,is,rmax,' rejected in 2nd round'
  endif  
enddo  

! call wsac1('d.sac',d,nstack,tlag,dtall,nerr)

! plot seismograms and stack, but reduce time by arrival time of phase
open(9,file='stack0.xy')
y=lymax-1.0                     ! y is in cm
call wsgmt(d,nstack,1.0,dtall,tlag,y)
close(9)
open(9,file='smgrs0.xy')
y=y-1.0
do j=1,ns
  if(ncor(j).le.2.or.polarity(j).eq.0.) cycle
  y=y-fourpt
  if(y.lt.1.0) exit             ! stop plot at bottom of the page
  kolor=mod(j,14)+1
  write(9,fmt='(a)') rgb(kolor)
  ! in reduced time first sample is at st1-delay if phase lines up
  call wsgmt(s(1,j),nsn(j),polarity(j),dtall,st1(j)-x(j),y)
enddo
close(9)
open(9,file='delay0.xy')
xmin=x(1)
xmax=xmin
do j=1,ns
  write(9,fmt='(2f8.2,f10.2)') stlon(j),stlat(j),x(j)
  xmin=min(x(j),xmin)
  xmax=max(x(j),xmax)
enddo
kscale=nint(max(abs(xmin),xmax))
kscale=max(kscale,1)
close(9)

! write GMT command file <gmtdelays>
open(9,file='gmtdelays')
open(10,file='rungmtdelays')
call startgmt(9)
write(9,fmt='(a)') '# Usage: gmtdelays delayx T Tlim (see rungmtdelays)'
write(9,fmt='(a)') 'rm -fr $argv[1].eps'
write(9,fmt='(a,f8.2)') 'set lon =',evlo
write(9,fmt='(a,f8.2)') 'set lat =',evla
write(9,fmt='(a)') 'set title = T=$argv[2]'
write(9,fmt='(a)') 'set dtlim = $argv[3]'
write(9,fmt='(a)') 'makecpt -Cpolar -D -T-${dtlim}/${dtlim}/0.2 -Z > delays.cpt'
write(9,fmt='(2a)') 'pscoast -Rg -JH${lon}/15c -B30g30:."$title": -Dc -A10000 ',  &
   '-Glightgray -Wthinnest -X2 -Y2 -P -K > $argv[1].eps'
write(9,fmt='(a)') 'echo $lon $lat | psxy -R -JH -Sa0.5 -G0 -K -O >> $argv[1].eps'   
write(9,fmt='(a)') 'psxy $argv[1].xy -M -R -JH -W1p -Cdelays.cpt -Sc0.20 -K -O >> $argv[1].eps'   
write(9,fmt='(a)') 'psscale -D0.5/3/5/0.4 -Ba1/:sec: -Cdelays.cpt -X-2 -O >> $argv[1].eps'   
write(9,fmt='(a)') 'gv $argv[1].eps &'
close(9)
call system('chmod a+x gmtdelays')

! plot broadband to screen
write(command,fmt='(a,i3)') 'gmtdelays delay0 Broadband',kscale
write(*,*) 'if kplot>0 calling ',command
if(kplot>0) call system(command)
write(10,fmt='(a)') command

if(abs(kplot) >0) then
  print *,'Plotting broadband stack'
  write(*,*) 'calling gmtstack stack0 smgrs0 0'
  call system('./gmtstack stack0 smgrs0 0')
endif  
if(kplot >0) then               ! plot GMT for broadband
  print *,'Give revised estimate for duration (now',duration,'s):'
  print *,'(0 for no change):'
  read(5,*) durnew
  if(durnew>0.) duration=durnew
  w2=tsigma+tap2+duration+pmax
  if(tpP.gt.0.) w2=min(w2,tpP-tap2)
  if(kdiag>0) write(7,*) 'Revised duration,w2=',duration,w2
endif
if(w2.le.0.) stop 'Error in time window length w2'
kplot=abs(kplot)

! compute broadband power spectrum  but do not yet write to unit 11
call powspec(d,nstack,dtall,rgb(1),pmin,peff,0,kdiag)
if(kdiag>0) write(7,*) 'pmin=',pmin,'s. Broadband spectrum peaks at ',peff,' s'
write(4,*) 'pmin=',pmin,'s. Broadband spectrum peaks at ',peff,' s'

! Determine number of frequency bands automatically so that it does not exceed pmin
! Note that nrbands0 still excludes the unfiltered delay (for which iband=1)
if(nrbands0.eq.0) nrbands0=MXBANDS-1
lowerp=0
if(jdebug>0) write(13,*) 'pmax, w2=',pmax,w2
do while(pmax>w2)
  pmax=0.5*pmax
  nrbands=nrbands-1
  if(jdebug>0) write(13,*) 'pmax, nrbands=',pmax,nrbands
  lowerp=1
enddo
if(lowerp>0) then
  if(kdiag>0) write(7,*) 'pmax lowered to ',pmax,' since w2=',w2
  write(4,*) 'pmax lowered to ',pmax,' since w2=',w2
  write(*,*) 'pmax lowered to ',pmax,' since w2=',w2
endif  
p=pmax
nrbands=0
do while(p.gt.pmin.and.nrbands.lt.nrbands0)
  p=0.5*p
  nrbands=nrbands+1
enddo
if(kdiag>0) write(7,*) 'Max period, number of frequency bands: ',pmax,nrbands
write(6,*) 'Max period, number of frequency bands: ',pmax,nrbands
write(4,*) 'Max period, number of frequency bands: ',pmax,nrbands
write(11,*) nrbands+1           ! nrbands bandpasses + broadband

! recompute broadband power spectrum to write it to unit 11
call powspec(d,nstack,dtall,rgb(1),pmin,peff,1,kdiag)

!==================================================
! DO LOOP OVER FILTERS
!==================================================

iband=2
period=pmax     ! we start with longest period

do while(iband.le.nrbands+1)              ! loop over filters

  ! filter all seismograms for band iband and store in sf

  ! WARNING: windowing has changed as of Oct 13 2012. We now
  ! shorten the original (t-w1,t+w2) interval if the period
  ! is less then w2
  
  ! set window length to pulse+/-2*period and filter
  ! unfiltered signal starts at t=-w1 and ends at t=w2
  ! filtered signals start at t=-w1-2*period, end at duration+2*period
  np=2*period/dtall
  do i=1,ns  
    do k=1,np           ! add zeros at start
      d(k)=0.
    enddo  
    do j=1,nsn(i)       ! put signal in the middle
      f=1.0
      ! taper if at end of signal
      if(j*dtall .gt. duration+period+tsigma) &
         f=0.5*cos(pi*(j*dtall-duration-period-tsigma)/tap2)+0.5
      if(j*dtall>duration+period+tsigma+tap2) f=0.
      d(np+j)=f*s(j,i)
    enddo  
    nsmp=2*np+(w1+duration+tsigma+tap2+2*period)/dtall      ! number of samples 
    do k=np+nsn(i)+1,nsmp                                  ! add zeros at end
      d(k)=0.
    enddo  
    nstack=nsmp         ! keep this for nr of samples in plot
    tlag=-2*period-w1
    sft1(i)=st1(i)-np*dtall     ! starting time
    call gaussfilter(d,nsmp,dtall,sf(1,i),nsfn(i),period)
  enddo
  xcl(iband)=duration+period+tsigma        ! true window for scatters

  ! cross correlate
  nrow=0
  kount=0
  do j=1,ns
    ncor(j)=0
  enddo  
  if(kdiag>0) write(7,fmt='(a)') 'stat1   stat2      corcoef      tmax     delay  B    period (nominal)'

  do is=1,ns

    print *,'Correlating ',kst(is),' in band ',iband,' cpu=',second()
    if(kdiag>0) write(7,*) 'Correlating ',kst(is),' in band ',iband,' cpu=',second()

    ti=sft1(is)
    tvi=0.

    ! add theoretical time to equations to ensure system is solvable
    ! even if clusters are disconnected, but use low weight epsilon
    nrow=nrow+1
    if(nrow > NSTAT1) stop 'Increase dimension of rhs'
    kount=kount+1
    if(kount >  NSTAT2) stop 'Increase dimenson of aa and ja'
    aa(kount)=epsilon   ! add damping row to the matrix (damp to IASP91)
    rhs(nrow)=0.        ! theoretical delay is time origin by definition
    rhsa(nrow)=0.
    ja(kount)=is        ! matrix column
    ia(nrow)=1          ! nr of nonzero entries in matrix row

    ! now correlate station is with all stations js>is
    do js=is+1,ns
      tj=sft1(js)
      tvj=0.
      isgnflg=0

      ! debug
      if(is.eq.0.and.iband.eq.3.and.js.eq.5) then
        bugfile=kst(js)
        print *,'Debugging ',kst(js)
      else
        bugfile='none'
      endif  

      taper=0.3*period
      call crc(sf(1,is),ti,tvi,nsfn(is),sf(1,js),tj,tvj,nsfn(js),dtall, &
        taper,outlier,tmax,ampl,rmax,isgn,isgnflg)
      if(is.eq.1) polarity(js)=isgn   
      if(abs(rmax)<cmin) then
        if(kdiag>0) write(7,fmt='(2a8,f10.3,1x,a)') kst(is),kst(js),rmax,'        -'
      else
        delay=tmax
        if(abs(delay).gt.outlier) cycle
        nrow=nrow+1
        if(nrow >  NSTAT1) stop 'Increase dimension of rhs'
        kount=kount+1
        if(kount >= NSTAT2) stop 'Increase dimenson of aa and ja'
        aa(kount)=-1.0
        ja(kount)=is
        kount=kount+1
        aa(kount)=+1.0
        ja(kount)=js
        rhs(nrow)=tmax
        rhsa(nrow)=log(ampl)
        ia(nrow)=2
        ncor(is)=ncor(is)+1       ! count nr of acceptable pairs
        ncor(js)=ncor(js)+1
        avcor(is,iband)=avcor(is,iband)+rmax    ! average corr coefficient
        avcor(js,iband)=avcor(js,iband)+rmax
        if(kdiag>0) write(7,fmt='(2a8,f10.3,2f10.2,i3,2f10.1)') kst(is),kst(js),  &
           isgn*rmax,tmax,delay,iband,period
      endif 
    enddo

  enddo  

  if(nrow<3) then
    if(kdiag>0) write(7,*) 'Not enough acceptable correlations in band',iband
    write(4,*) 'Not enough acceptable correlations in band',iband
    write(6,*) 'Not enough acceptable correlations in band',iband
    write(11,*) 0       ! empty frequency band
    iband=iband+1
    period=0.5*period
    cycle
  endif  

  do i=1,ns
    if(ncor(i)>0) avcor(i,iband)=avcor(i,iband)/ncor(i)
  enddo  

  ! solve system of equations for this frequency band
  print *,'Solving for band ',iband
  itmax=5*ns
  do j=1,nrow
    u(j)=rhs(j)         ! since lsqr overwrites u, we do not use rhs itself
  enddo  
  call lsqr(nrow,ns,x,u,v,w,itmax,r,kdiag)      ! solves Ax=u
  if(r<0.) then
    write(4,*) 'Insufficient convergence in band',iband
    write(6,*) 'Insufficient convergence in band',iband
    if(kdiag>0) write(7,*) 'Insufficient convergence in band',iband
    write(11,*) 0       ! empty frequency band
    iband=iband+1
    period=0.5*period
    cycle
  endif  

  ! compute A*x=b to get predicted delay times for each station pair
  call asol(x,b,nrow,ns)		! computes b=A*x

  ! estimate standard errors
  kount=0
  do i=1,ns
    sigma(i)=0.
  enddo  
  chi2=0
  nchi=0
  do irow=1,nrow
    if(ia(irow).eq.0) then
      cycle
    else if(ia(irow).eq.1) then
      kount=kount+1
      cycle
    else if(ia(irow).eq.2) then
      kount=kount+1
      nchi=nchi+1
      is=ja(kount)
      kount=kount+1
      js=ja(kount)
      dif2=(b(irow)-rhs(irow))**2
      sigma(is)=sigma(is)+dif2
      sigma(js)=sigma(js)+dif2
      chi2=chi2+dif2
    else
      print *,'Bug in row',irow,' ia=',ia(irow)
      stop
    endif  
  enddo
  do i=1,ns
    if(ncor(i).gt.2) then
      sigma(i)=sqrt( sigma(i)/(ncor(i)-2.0) + (sigma0+sigmaN/ns)**2 )
      siga(i,iband)=0.1
    else
      sigma(i)=999.99
      siga(i,iband)=999.99
    endif  
    sigt(i,iband)=sigma(i)
    tobs(i,iband)=x(i)+stv(i)
  enddo  

  ! compute stack
  do i=1,nstack
    d(i)=0.
  enddo  
  do j=1,ns       ! stack stations
    if(ncor(j).le.2.or.polarity(j).eq.0.) cycle
    kdelay=(tlag+x(j)-sft1(j))/dtall
    k2=min(nstack,nsfn(j)-kdelay)
    k1=max(1,1-kdelay)
    do k=k1,k2
      d(k)=d(k)+polarity(j)*sf(k+kdelay,j)      ! stack station j
    enddo  
  enddo  

  ! compute stack power spectrum and write it to unit 11
  call powspec(d,nstack,dtall,rgb(iband),pflt,peff,1,kdiag)

  write(4,fmt='(//,a,i2,a,f10.2,a,f10.2)') 'Delays in band',iband, &
     ' period=',period,' effective period=',peff
  write(4,fmt='(a,f10.2,a,f10.2,a,f10.3,a)') 'Chi2=',chi2,', relative:', &
     chi2/nchi,' or rms=',sqrt(chi2/nchi),' sec.'
  write(4,fmt='(/,2a)') 'station      T_arr    T_pred     delay     ', &
       'sigma      ncor   avcor     azi   gcarc'
  do j=1,ns
    if(ncor(j).gt.0) write(4,fmt='(a8,4f10.2,i10,3f8.2)') &
       kst(j),x(j)+stv(j),stv(j),x(j),sigma(j),ncor(j),avcor(j,iband), &
       azev(j),dist(j) 
  enddo  

  ! solve system of equations for amplitude ratios
  print *,'Solving amplitudes for band ',iband
  itmax=5*ns
  do j=1,nrow
    u(j)=rhsa(j)        ! careful - lsqr overwrites u
  enddo  
  ! use a trick to get rid of the time regularisation
  do i=1,kount
    if(abs(aa(i)-epsilon) < 0.5*epsilon) aa(i)=0.
  enddo
  call lsqr(nrow,ns,xa,u,v,w,itmax,r,kdiag)      ! solves Ax=u
  if(r<0.) then
    if(kdiag>0) write(7,*) 'Warning: LSQR does not converge on amplitudes band',iband
    write(4,*) 'Warning: LSQR does not converge on amplitudes band',iband
    do i=1,ns
      siga(i,iband)=999.99
    enddo
    goto 600
  endif  
  
  ! compute A*x=b to get predicted amplitudes for each station pair
  call asol(xa,b,nrow,ns)		! computes b=A*x
  xav=0.
  do i=1,ns
    xav=xav+xa(i)
  enddo
  xav=xav/ns
  do i=1,ns
    xa(i)=xa(i)-xav
  enddo
  
  ! estimate standard errors
  kount=0
  do i=1,ns
    sigma(i)=0.
  enddo  
  chi2=0
  nchi=0
  do irow=1,nrow
    if(ia(irow).eq.0) then
      cycle
    else if(ia(irow).eq.1) then
      kount=kount+1
      cycle
    else if(ia(irow).eq.2) then
      kount=kount+1
      nchi=nchi+1
      is=ja(kount)
      kount=kount+1
      js=ja(kount)
      dif2=(b(irow)-rhsa(irow))**2
      sigma(is)=sigma(is)+dif2
      sigma(js)=sigma(js)+dif2
      chi2=chi2+dif2
    else
      print *,'Bug in row',irow,' ia=',ia(irow)
      stop
    endif  
  enddo
  do i=1,ns
    if(ncor(i).gt.2) then
      ! check if -2 is correct (comes from VanDecar&Crosson, 1990)
      sigma(i)=sqrt( sigma(i)/(ncor(i)-2))
      siga(i,iband)=0.1           ! guess, needs improving
    else
      sigma(i)=999.99
      siga(i,iband)=999.99
    endif  
    aobs(i,iband)=xa(i)
  enddo  
  
  write(4,fmt='(//,a,i2)') 'Log Amplitudes in band ',iband
  write(4,fmt='(a,f10.2,a,f10.2,a,f10.3,a)') 'Chi2=',chi2,', relative:', &
     chi2/nchi,' or rms=',sqrt(chi2/nchi)
  write(4,fmt='(/,a)') 'station   lnA_obs     sigma      icor   avcor     azi   gcarc'
  do j=1,ns
    if(ncor(j).gt.0) write(4,fmt='(a8,2f10.5,i10,3f8.2)') &
       kst(j),xa(j),sigma(j),ncor(j),avcor(j,iband),azev(j),dist(j) 
  enddo  
  600 continue

  ! second check on sufficient correlation, now with stacked signal
  do is=1,ns  
    if(sigt(is,iband).gt.900.) cycle
    call crc(sf(1,is),ti,tvi,nsfn(is),d,tlag,0.,nstack,dtall, &
        taper,outlier,tmax,ampl,rmax,isgn,isgnflg)
    polarity(is)=isgn       
    if(isgn.eq.0.or.abs(rmax)<cmin) then
      sigt(is,iband)=999.99
      siga(is,iband)=999.99
      if(kdiag>0) write(7,fmt='(2a8,i3,f10.3,1x,a)') kst(is),kst(js),isgn,rmax,'   Reject'
      if(jdebug>0) write(13,*) 'js,is,rmax=',js,is,rmax,' rejected in 2nd round'
    endif  
  enddo  

  ! align signals and plot stack if kplot>1
  write(gmtfil,fmt='(i6)') iband
  gmtfil(1:5)='stack'
  gmtfil(7:9)='.xy'
  open(9,file=gmtfil)
  y=lymax-1.0
  call wsgmt(d,nstack,1.0,dtall,tlag,y)
  close(9)
  write(gmtfil,fmt='(i6)') iband
  gmtfil(1:5)='smgrs'
  gmtfil(7:9)='.xy'
  open(9,file=gmtfil)
  y=y-1.0
  do j=1,ns
    if(ncor(j).le.2.or.polarity(j).eq.0.) cycle
    y=y-fourpt
    if(y.lt.1.0) exit             ! stop plot at bottom of the page
    kolor=mod(j,14)+1
    write(9,fmt='(a)') rgb(kolor)
    ! plot signal with estimated time delay x(j) by shifting start time
    call wsgmt(sf(1,j),nsfn(j),polarity(j),dtall,sft1(j)-x(j),y)
  enddo
  close(9)
  write(gmtfil,fmt='(i6)') iband
  gmtfil(1:5)='delay'
  gmtfil(7:9)='.xy'
  open(9,file=gmtfil)
  do j=1,ns
    if(ncor(j)>0) write(9,fmt='(2f8.2,f10.2)') stlon(j),stlat(j),x(j)
  enddo
  close(9)
  write(command,fmt='(2a,f6.1,i3)') 'gmtdelays ',gmtfil(1:6),peff,kscale
  write(*,*) 'if kplot>1 calling ',command
  if(kplot>1) call system(command)
  write(10,fmt='(a)') command
  
  
  if(kplot >1) then               ! plot GMT for filtered signals
    print *,'calling gmtstack for band',iband,' period:',period
    if(iband.eq.1) call system('./gmtstack stack1 smgrs1 1')
    if(iband.eq.2) call system('./gmtstack stack2 smgrs2 2')
    if(iband.eq.3) call system('./gmtstack stack3 smgrs3 3')
    if(iband.eq.4) call system('./gmtstack stack4 smgrs4 4')
    if(iband.eq.5) call system('./gmtstack stack5 smgrs5 5')
    if(iband.eq.6) call system('./gmtstack stack6 smgrs6 6')
    if(iband.eq.7) call system('./gmtstack stack7 smgrs7 7')
  endif

  write(12,*) 'Band ',iband,' period=',period,' Effective period=',peff
  if(kdiag>0) write(7,fmt='(//,"  Station       N")')
  do j=1,ns
    if(kdiag>0) write(7,fmt='(a8,i8)') kst(j),ncor(j)
    if(ncor(j) > 0) write(12,fmt='(a)') kst(j)
  enddo  

  iband=iband+1
  period=0.5*period

enddo   ! end of filter loop

! write rest of data file
kunit=1                 ! assume displacement units for noise
rms0=0.                 ! TODO (but not used elsewhere)
kpole=0                 ! TODO (now: epicentral distance < 180 deg)
do k=1,ns
  nobst=0
  nobsa=0
  ! find out how many time/ampl measurements observed for this station
  do i=1,nrbands+1
    if(sigt(k,i) < 900.0) nobst=nobst+1
    if(siga(k,i) < 900.0) nobsa=nobsa+1
  enddo  
  write(11,20) idate,iotime,ievt,kluster,kst(k),knet(k),kcmp(k),evla,  &
     evlo,evdp,stlat(k),stlon(k),stel(k),nobst,nobsa,kpole
  write(11,fmt='(i2,16f10.1)') kunit,rms0,(rmsb(i),i=1,nrbands+1)
  write(11,*) nobst
  do i=1,nrbands+1
    if(sigt(k,i) < 900.0) write(11,fmt='(3f9.2,i3,f9.1)') tobs(k,i), &
             sigt(k,i),avcor(k,i),i,xcl(i)
  enddo  
  write(11,*) nobsa
  do i=1,nrbands+1
    if(siga(k,i) < 900.0) write(11,fmt='(3f9.2,i3,f9.1)') aobs(k,i), &
             siga(k,i),avcor(k,i),i,xcl(i)
  enddo  
enddo
20 format(3i8,i4,1x,a8,1x,a8,1x,a3,2f9.3,f7.1,2f9.3,f7.3,3i4)

close(11)

! Write GMT file for spectral band plotting
open(11,file='gmtspec')
call startgmt(11)
write(11,fmt='(a)') 'rm -fr spec.eps'
write(11,fmt='(5a)') 'psxy spec.'//trim(ident)//'.xy -M -H -R0.01/2.0/0/1.05 ', &
    '-JX20l/10 -Ba2f1:Hz:/a0.2:power::.power-spectra:eWnS -W1p  > spec.eps'
write(11,fmt='(a)') 'gv -resize spec.eps &'
close(11)

! plot spectra on screen
write(*,*) 'if kplot>0 calling gmtspec'
call system('chmod a+x gmtspec')
if(kplot>0) call system('./gmtspec')

end

!integer function leng(char)
!
!character char*(*)
!integer l
!
!l=len(char)
!do leng=l,1,-1
!  if(char(leng:leng) /= ' ') return
!enddo
!leng=0
!return
!end

!subroutine shift(d,n,beg,del,tshift)
!
!! shifts start of signal. adds 0's at begin or end where
!! needed to keep length nsmp
!
!include 'crc.inc'
!dimension d(NSAC)
!
!ishift=nint(tshift/del)
!if(ishift==0) return
!
!if(ishift>0) then
!  do i=1,n-ishift
!    d(i)=d(i+ishift)
!  enddo  
!  do i=n-ishift+1,n
!    d(i)=0.
!  enddo  
!else 
!  do i=n+ishift,1,-1
!    d(i-ishift)=d(i)
!  enddo
!  do i=1,-ishift-1
!    d(i)=0.
!  enddo
!endif
!
!beg=beg+ishift*del
!
!return
!end




subroutine crc(g1,b1,t1,n1,g2,b2,t2,n2,dt,taper,wt,tmax,ampl,rmax,isgn,isgnflg)

! correlates g1 with g2. Copies g1 and g2 to avoid truncating/tapering
! the originals. A maximum in correlation if searched around a delay t2-t1.

! input:
! g1, g2: times signals with start times b1, b2
! t1, t2: markers indicating start of pulses (may be rough estimates)
! dt: sampling interval
! taper: length of (cosine) taper 
! wt: half length of correlation function to be searched for delay between g1 and g2
! isgnflg: if 0, allow for sign changes to keep rmax always positive, otherwise 1.

! output:
! tmax: delay between g1 and g2
! rmax: correlation coefficient (absolute)
! isgn: -1 if sign change between g1 and g2, otherwise +1 (only if isgnflg=0)
! ampl: ratio of (integrated) amplitudes g2/g1

include 'crc.inc'
dimension g1(n1),g2(n2),cr(10*NWIN),f1(5*NWIN),f2(5*NWIN)
character*8 bugfile
common /debug/ bugfile

! debugging
if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'g1.sac',g1,n1,b1,dt,nerr)
  call wsac1(trim(bugfile)//'g2.sac',g2,n2,b2,dt,nerr)
  if(kdiag>0) write(7,*) 'crc wrote ',trim(bugfile)//'g1.sac etc'
endif  

do i=1,n1
  f1(i)=g1(i)
enddo  
do i=1,n2
  f2(i)=g2(i)
enddo  
e1=b1+(n1-1)*dt; e2=b2+(n2-1)*dt

! taper
if(n1<0) then
  print *,'Debug: ftper in crc n1,t1,t2,taper=',n1,t1,t2,taper
  print *,'b1,b2=',b1,b2
  stop
endif  
call ftaper(f1,b1,dt,n1,taper,taper)
if(n2<0) then
  print *,'Debug: ftper in crc n2,t1,t2,taper=',n2,t1,t2,taper
  print *,'b1,b2=',b1,b2
  stop
endif  
call ftaper(f2,b2,dt,n2,taper,taper)

! debugging
if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'f1.sac',f1,n1,b1,dt,nerr)
  call wsac1(trim(bugfile)//'f2.sac',f2,n2,b2,dt,nerr)
endif
tpred=t2-t1     ! predicted delay of f2 with respect to f1
call corr(f1,b1,n1,f2,b2,n2,dt,tpred,wt,tmax,ampl,rmax,isgn,isgnflg,cr,nc,tau0)

return
end

subroutine window(f,b,dt,n,t1,t2,tap1,tap2)

include 'crc.inc'
dimension f(NSAC)
character*8 bugfile
common /debug/ bugfile

! selects window between t1 and t2 from sac file, detrends this 
! input: signal f with start time b, interval dt, n samples
!        t1,t2 are the true signal start and end times, 
!        tap1 and tap2 the taper length at start and end
! output: f,b and n are modified to reflect the windowed time series

! we detrend using first/last sample in the windowed signal f to avoid 
! drastic windowing effects that might influence delay time estimates 

! debug
if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'inw.sac',f,n,b,dt,nerr)
endif  

i1=(t1-b)/dt+1
i2=(t2-b)/dt+1

! adjust t1 and t2 to actual sample times
if(i1<1) i1=1
if(i2>n) i2=n
t1=b+(i1-1)*dt
t2=b+(i2-1)*dt

if(i1>i2) stop 't outside window'
nnew=i2-i1+1
! print *,'debug window: ',t1,t2,i1,i2,nnew
if(nnew > NSAC) stop 'increase NSAC in subroutine window'
! move f(i1) to first sample etc, up to f(i2) that becomes f(nnew)
do i=1,nnew
  f(i)=f(i+i1-1)
enddo  

! remove trend
do i=1,nnew
  f(i)=f(i)-f(1)-(i-1)*(f(nnew)-f(1))/nnew
enddo  

n=nnew           ! replace n with new signal length
b=t1
if(n<0) then
  print *,'Debug: ftper in window n,t1,t2,tap1,tap2=',n,t1,t2,tap1,tap2
  print *,'i1,i2=',i1,i2
  stop
endif  
call ftaper(f,b,dt,n,tap1,tap2)         ! and taper

if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'outw.sac',f,n,b,dt,nerr)
endif  

return
end

subroutine corr(f1,b1,n1,f2,b2,n2,dt,tpr,wt,tmax,ampl,rmax,isgn,isgnflg,cr,nc,tau0)

include 'crc.inc'
character*8 bugfile
common /debug/ bugfile
dimension f1(n1),f2(n2),cr(2*NWIN)

! This routine computes the cross-correlation between signals
! f1, starting at time b1 of length n1, and f2, starting at time b2
! of length n2. Both signals have sampling interval dt secs. A
! correlation maximum is sought for a time shift in between tpr-wt
! and tpr+wt. If wt=0, tpr is ignored and the whole time window
! is searched for a max correlation.

! output: shift tmax (f2 w.r.t. f1) for which the correlation has a maximum
! absolute value rmax, isgn is the sign of the correlation at tmax.
! isgn=0 indicates an error return (tmax is then set to tpr). 
! The correlation function is in cr(1),...,cr(nc). The time shift of cr(1) is tau0.

jdebug=0

k1=1-n1
k2=n2-1
nc=k2-k1+1
tau0=b2-b1+k1*dt	! offset for cr(1)

sum1=0.; sum2=0.
!do i=1,n1
  !sum1=sum1+f1(i)**2
!enddo  
!do i=1,n2
  !sum2=sum2+f2(i)**2
!enddo  
sum1=dotpr(f1,f1,n1)
sum2=dotpr(f2,f2,n2)
snorm=sqrt(sum1*sum2)
eps=1.0e-10*max(sum1,sum2)      ! water level for amplitude ratio
ampl=sqrt(sum2/(sum1+eps))     
if(jdebug>0) write(13,*) 'n1,n2,sum1,sum2,ampl=',n1,n2,sum1,sum2,ampl

do k=k1,k2 		! loop over  offset tau
  i1=max(1,1-k)
  i2=min(n1,n2-k)
  ndot=i2-i1+1
! sum=0.
! do i=i1,i2 
!   sum=sum+f1(i)*f2(i+k)
! enddo
  sum=dotpr(f1(i1),f2(i1+k),ndot)
  sum=sum/snorm
  cr(k-k1+1)=sum
enddo

rmax=0.
imax=0
i1=1; i2=nc
if(wt>0.) then
  i1=max(1,nint((tpr-wt-tau0)/dt)+1)
  i2=min(nc,nint((tpr+wt-tau0)/dt)+2)
endif


if(wt>0.) then
  i1=max(1,nint((tpr-wt-tau0)/dt)+1)
  i2=min(nc,nint((tpr+wt-tau0)/dt)+2)
endif  
if(isgnflg.eq.1) then	! if not allowing sign change
  do i=i1,i2 
    if(cr(i)<rmax) cycle
    rmax=cr(i)
    isgn=+1
    imax=i
  enddo
endif
if(isgnflg.eq.0) then	! if allowing sign change
  do i=i1,i2 
    if(abs(cr(i))<rmax) cycle
    if(cr(i)<0.) then
      rmax=abs(cr(i))
      isgn=-1
    else 
      rmax=cr(i)
      isgn=+1
    endif
    imax=i
  enddo
endif
if(imax==i1 .or. imax==i2) then
  isgn=0
  rmax=0.
  tmax=99999.           ! error default
  ! debug: usually we get here if wt (=outlier) was underestimated, but 
  ! uncommenting the following lines allows you to inspect what happened
  ! print *,'Bug in corr? i1,i2,imax=',i1,i2,imax
  ! print *,'k1,k2,kc,tau0=',k1,k2,kc,tau0
  ! print *,'tpred,wt=',tpr,wt
  ! print *,'b1,n1=',b1,n1
  ! print *,'b2,n2=',b2,n2
  ! call wsac1('crc.sac',cr,k2-k1+1,k1*dt,dt,nerr)
  ! call wsac1('f1.sac',f1,n1,b1,dt,nerr)
  ! call wsac1('f2.sac',f2,n2,b2,dt,nerr)
  ! stop 'debug stop in corr'
  return		! return isgn=0
endif

k=imax+k1-1

! interpolate for maximum
imax=max(2,imax)
imax=min(nc-1,imax)
y0=cr(imax-1)
y1=cr(imax)
y2=cr(imax+1)

c1=0.5*(y2-y0)/dt
c2=(y2-2*y1+y0)/dt**2

tmax=tau0+(imax-1)*dt-c1/c2
rmax=y1-0.5*c1*c1/c2
rmax=abs(rmax)

return
end

real function dotpr(a,b,n)
real a(n),b(n)
dotpr=dot_product(a,b)
return
end


subroutine clear(c)
character*8 c
! removes C end-of-character symbol from station code
if(ichar(c(4:4))<32 .or. ichar(c(4:4))>126) c(4:4)=' '
if(ichar(c(5:5))<32 .or. ichar(c(5:5))>126) c(5:5)=' '
if(ichar(c(6:6))<32 .or. ichar(c(6:6))>126) c(6:6)=' '
if(ichar(c(7:7))<32 .or. ichar(c(7:7))>126) c(7:7)=' '
if(ichar(c(8:8))<32 .or. ichar(c(8:8))>126) c(8:8)=' '
return
end

subroutine rmean(d,nsmp,average)
! removes mean of signal d 
dimension d(nsmp)
sum=0.
do i=1,nsmp
  sum=sum+d(i)
enddo  
average=sum/nsmp
do i=1,nsmp
 d(i)=d(i)-average
enddo 
return
end

subroutine decimate(din,ndecim,nsmp,b,dt,npow,flowpass)

! decimate din by factor 2**ndecim using FFT, and determine npow
! while adjusting nsmp to 2**npow

! input: din,ndecim,nsmp,dt,flowpass
! output: npow, din (decimated), nsmp and dt (adjusted)

include 'crc.inc'
parameter (PI=3.141592654, TAPER=0.25)

dimension din(NSAC)
complex s(NSAC)

nsmp0=nsmp
dt0=dt			! save input


nred=2**ndecim 
dt=dt0*nred	! decimate

np=MXPOW
n2=2**np
do while(n2 >= nsmp0) 
  n2=n2/2
  np=np-1 
enddo
n2=2*n2         ! n2 is first power of 2 above nsmp0
np=np+1
if(np>NSAC) then 
  print *,'Input signal exceeds',NSAC,' samples:',nsmp0
  stop 'Fatal error, increase MXPOW'
endif
if(np < 3) then
  print *,'Decimated signal has length:',n2
  stop 'Fatal error'
endif

! fill up to power of 2
do i=nsmp0+1,n2
  din(i)=0.
enddo

if(ndecim==0) then
  npow=np
  nsmp=n2
else 
  do i=1,n2					! store in real part of s()
    s(i)=cmplx(din(i),0.)
  enddo
  call clogl(np,n2,s,+1.0,dt0)			! FFT to freq domain
  nsmp=n2/nred 					! new length
  nhalf=nsmp/2					! length of positive spectrum
  df=1.0/(nsmp0*dt0)				! frequency spacing
  i1=flowpass/df+1.
  i2=(1.0+TAPER)*flowpass/df+1.
  if(i2.le.i1) then
    print *,'np,nsmp,flowpass,df=',np,nsmp,flowpass,nsmp
    print *,'i1,i2,nhalf=',i1,i2,nhalf
    stop 'Fatal error for flowpass'
  endif
  dftaper=(i2-i1)*df
  do i=i1,i2 
    f=(i-1)*df
    x=(f-flowpass)/dftaper
    a=0.5+0.5*cos(x*PI)
    s(i)=a*s(i)
  enddo
  do i=i2+1,nhalf
    s(i)=0.
  enddo
  npow=np-ndecim					! new power of 2
  call ftinv(s,npow,nsmp,dt,din)		! FFT back to time
endif

return
end

subroutine ftinv(s,npow,nsmp,dt,r)
include 'crc.inc'
dimension r(nsmp)
complex s(nsmp)
 
! combines construction of negative spectrum, inverse fft, and
! storage in a real array. it is possible to equivalence s and
! r in the main program.
! s(1)...s(nhalf) contain positive complex spectrum of real signal
! nsmp=2*nhalf=2**npow, dt sampling interval, r destination array
nhalf=nsmp/2
call rspec(s,nhalf)
call clogl(npow,nsmp,s,-1.0,dt)
do i=1,nsmp
  r(i)=real(s(i))
enddo
return
end

subroutine rspec(s,np2)
 
! constructs negative spectrum, removes dc component of real time signal.
! input: positive spectrum in s(1)...s(np2).
! output: complete spectrum in s(1)...s(2*np2), with s(1)=s(np2+1)=0  
 
include 'crc.inc'
complex s(NSAC)
n=2*np2
n1=np2+1
s(n1)=0.
s(1)=0.
do i=1,np2
  s(np2+i)=conjg(s(np2+2-i))
enddo
return
end

      subroutine clogl(n,lx,x,zzign,dt)
! performs fft on signals with length lx=2**n and sampling interval
! of dt seconds (if in the time domain; notice that dt*df=1/2**n).
! the signal is stored in x. it may be complex.
! the spectrum is returned in x. it is almost always complex.
! a time-to-frequency transform is done with zign=+1. (conform
! the convention adopted in aki and richards - the alternative
! convention may be obtained by taking complex conjugates after
! the call to clogl).
! the normalization factor 1./twopi occurs in the frequency-to
! time transform (again aki&richards).
! normalization is such that physical dimensions are respected.
! thus, if the time signal is dimensioned in meters, the
! resulting spectral density in x is in meters/hz. for example,
! if the time signal is the unit sinc function of width dt, centered
! at t=0, the spectral density is dt for all values of the frequency.

! array locations: if x contains the spectrum, it has the spectrum
! for positive frequencies in the first 2**n/2+1 elements, such that
! x(1) is at 0 hz, x(2) at df hertz, and x(2**n/2+1) at the nyquist,
! where df=1./(2**n*dt) and the nyquist is 1./(2*dt) hz.
! the second half of x contains the spectrum for negative frequencies
! such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
! if x contains the time signal, x(1) is at time 0, x(2)
! at time dt etc.

      include 'crc.inc'
      dimension x(lx),m(25)
      complex x,wk,hold,q
      zign=zzign
      if(zign.ge.0.) then
        zign=1.
      else
        zign=-1.
      endif
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=zign*6.283185308*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q
    2 continue
      do 3 i=2,n
      ii=i
      if(k.lt.m(i)) go to 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if(k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if(k.lt.m(i)) go to 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if(zign.gt.0.) go to 9
      flx=flx*dt
      do 8 i=1,lx
    8 x(i)=x(i)/flx
      return
    9 do 10 i=1,lx
   10 x(i)=x(i)*dt
      return
      end

subroutine gaussfilter(d1,n1,dt,d2,n2,period)

! apply Gaussian filter to signal d1 of length n1
! returns signal d2 of length n2, such that n2 is at least 5*period
! period is the 'dominant' period of the bandpass filter
! dt is the sampling interval of both d1 and d2

include 'crc.inc'
character*8 fname
dimension d1(NWIN),d2(5*NWIN)
complex s(NWIN)
logical ruthere
character*8 bugfile
common /debug/ bugfile
data twopi/6.283185307/

! debugging
if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'d1.sac',d1,n1,b1,dt,nerr)
endif  

tau=period/0.627  ! See eq (42) in Hung et al., GJI 141:175, 2000
f1=tau*tau*0.5/twopi            ! tau^2/4pi
f2=f1/twopi                     ! tau^2/8pi^2

! find power of 2 > or equal to n1
np=12
nmin=max(nint(5*period/dt),n1)
do while (2**np > nmin)
  np=np-1
enddo
np=np+1
n2=2**np
nhalf=n2/2
if(n2 > 5*NWIN) then
  print *,'gaussfilter n1,period,dt=',n1,period,dt
  print *,'nmin=',nmin,'np,n2=',np,n2
  stop 'Increase dimension of d2 in gaussfilter'
endif

! store signal in complex array s
do i=1,n1
  s(i)=cmplx(d1(i),0.)
enddo
do i=n1+1,n2
  s(i)= cmplx(0.,0.)
enddo  

! step in Hz and rad/s
df=1.0/(n2*dt)
dw=twopi*df

call clogl(np,n2,s,+1.0,dt)			! FFT to freq domain

! filter
do i=1,nhalf
  w=(i-1)*dw
  w2=w*w
  s(i)=s(i)*w2*f1*exp(-w2*f2)   ! Hung 2000, eq (40) 
end do  

call ftinv(s,np,n2,dt,d2)		! FFT back to time

! debugging
if(bugfile.ne.'none') then
  call wsac1(trim(bugfile)//'d2.sac',d1,n1,b1,dt,nerr)
endif  

return
end

subroutine ftaper(f,b,dt,n,taper1,taper2)
dimension f(n)
ntaper=nint(taper1/dt)
do i=1,ntaper 
  tt=(i-1)*dt
  x=0.5-0.5*cos(3.14159265*tt/taper1)
  f(i)=x*f(i)
enddo
ntaper=nint(taper2/dt)
do i=1,ntaper 
  tt=(i-1)*dt
  x=0.5-0.5*cos(3.14159265*tt/taper2)
  f(n+1-i)=x*f(n+1-i)
enddo
return
end

subroutine lsqr(m,n,x,u,v,w,itmax,r,kdiag)
 
!   subroutine to solve the linear tomographic problem Ax=u using the 
!   lsqr algorithm.
 
!   reference: C.C.Paige and M.A.Saunders, ACM Trans.Math.Softw. 8, 43-71, 1982
!   and ACM Trans.Math.Softw. 8, 195-209, 1982. See also A.v.d.Sluis and H.v.d.
!   Vorst in: G. Nolet (ed.), Seismic Tomography, Reidel, 1987.
!
!   Input: m is the number of data (rows), n the number of unknowns (columns),
!      u contains the data (is overwritten), itmax is the number of iterations.
!   Output: x is the solution; intermediate results are written to a diskfile
!      named <tomo.int> which must not already exist before the call.
!   Scratch: arrays v(n) and w(n)
!   Subroutines: routines avpu and atupv must be supplied by the user.
!      avpu(m,n,u,v) computes u=u+A*v for given input u,v (overwrites u)
!      atupv(m,n,u,v) computes v=v+A(transpose)*u for given u,v (overwrites v)

! note to the pointers in the sparse matrix A:
! since A contains mostly zero entries, we only store the nonzero elements.
! the number of nonzero elements in row i is stored in ia(i). All these
! nonzero elements are stored one after another in a, whereby ja(j)
! identifies the column number of element a(j). kount is the total nr of
! nonzero elements.

dimension x(n),u(m),v(n),w(n)
! write(6,fmt='(" iter  cell 1",9x,"rms",9x,"rel")')
if(kdiag>0) write(7,*) 'LSQR called with n,m=',n,m
if(kdiag>0) write(7,fmt='(" iter  cell 1",9x,"rms",9x,"rel")')
do i=1,n 
  x(i)=0
  v(i)=0 
end do
! initialize
call normlz(m,u,beta)
b1=beta 
call atupv(m,n,u,v)
call normlz(n,v,alfa)
rhobar=alfa
phibar=beta; 
do i=1,n 
  w(i)=v(i) 
end do  
! write(6,fmt='(i5,3g12.4)') 0,x(1),beta,1		
if(kdiag>0) write(7,fmt='(i5,3g12.4)') 0,x(1),beta,1.0
rlast=99999.9
do iter =1,itmax 		! repeat for at most itmax iterations
   a=-alfa
   do i=1,m 
     u(i)=a*u(i) 
   end do ! bidiagonalization
   call avpu(m,n,u,v)
   call normlz(m,u,beta) 
   b=-beta
   do i=1,n 
     v(i)=b*v(i) 
   end do 
   call atupv(m,n,u,v) 
   call normlz(n,v,alfa)
   rho=sqrt(rhobar*rhobar+beta*beta)	! modified QR factorization
   c=rhobar/rho
   s=beta/rho 
   teta=s*alfa
   rhobar=-c*alfa
   phi=c*phibar
   phibar=s*phibar 
   t1=phi/rho
   t2=-teta/rho
   do i=1,n		! update solution x and storage vector w
     x(i)=t1*w(i)+x(i)
     w(i)=t2*w(i)+v(i)
   end do 
   r=phibar/b1
!  write(6,fmt='(i5,3g12.4)') iter,x(1),phibar,r
   if(kdiag>0) write(7,fmt='(i5,3g12.4)') iter,x(1),phibar,r
   rdif=abs(r-rlast)
   rlast=r
   if(rdif.lt.1.0e-6) return
end do 

! insufficient convergence
r=-1.0

return
end


subroutine normlz(n,x,s)
dimension x(n)
s=0.
do i=1,n 
  s=s+x(i)**2 
end do 
s=sqrt(s)
ss=1./s
do i=1,n 
  x(i)=x(i)*ss 
end do 
return
end

subroutine asol(sol,b,m,n)		! computes b=a*sol
dimension sol(n),b(m)
include 'crc.inc'
common /LSQ/ nrow,aa(NSTAT2),ia(NSTAT1),ja(NSTAT2)
kount=0
do i=1,m 				! work row by row
  b(i)=0.
  nc=ia(i)
  if (nc.eq.0) cycle		! nr of nonzero elements in this row
  do j=1,nc 
    kount=kount+1
    b(i)=b(i)+aa(kount)*sol(ja(kount))	! multiply nonzeroes with vector
  end do
end do
return 
end
  
subroutine avpu(m,n,u,v)		! computes u=u+a*v
dimension u(m),v(n)
include 'crc.inc'
common /LSQ/ nrow,aa(NSTAT2),ia(NSTAT1),ja(NSTAT2)
kount=0
do i=1,m 				! work row by row
  nc=ia(i)
  if (nc.eq.0) cycle
  do j=1,nc 
    kount=kount+1
    u(i)=u(i)+aa(kount)*v(ja(kount))
  end do
end do
return
end

subroutine atupv(m,n,u,v)		! computes v=v+a(transpose)*u
dimension u(m),v(n)
include 'crc.inc'
common /LSQ/ nrow,aa(NSTAT2),ia(NSTAT1),ja(NSTAT2)
kount=0
do i=1,m 				! work row by row (here too!)
  nc=ia(i)
  if (nc.eq.0) cycle
  do j=1,nc 
    kount=kount+1
    jj=ja(kount)
    v(jj)=v(jj)+aa(kount)*u(i)		! add to the right vector element
  end do
end do
return
end

subroutine wsgmt(d,n,polarity,dt,b,y)

! writes seismogram in array d to unit 9 at hight y for GMT plots
! input seismogram d of dimension n and polarity +/-1, sampling interval dt,
! start time b, plot offset y (cm).

dimension d(n)

dmax=0.
do i=1,n
  dmax=max(abs(d(i)),dmax)
enddo
if(dmax.eq.0) stop 'Error dmax=0 in wsgmt'

t=b
do i=1,n
  write(9,fmt='(f6.2,f7.2)') t,y+polarity*d(i)/dmax
  t=t+dt
enddo

return
end

subroutine powspec(d,n1,dt,rgb,pmin,peff,iwrite,kdiag)

! computes power spectrum , writes to unit 11 (data file) if iwrite=1,
! and to unit 15 (GMT plottable file)
! output: pmin, period where power is down to 0.05*maximum power
!         peff, period where power spectrum is at its maximum

! input: time series d with n1 elements, spaced dt sec

include 'crc.inc'
dimension d(NSAC)
complex s(8192)
character*18 rgb
data twopi/6.2831853/

! find power of 2 > or equal to n1
np=12
do while (2**np .ge. n1)
  np=np-1
enddo
np=np+1
nhalf=2**np      
n2=2*nhalf
if(n2 > 8192) stop 'Increase dimension of s in powspec'
np=np+1

! store signal in complex array s
do i=1,n1
  s(i)=cmplx(d(i),0.)
enddo
do i=n1+1,n2
  s(i)= cmplx(0.,0.)
enddo  

! step in Hz and rad/s
df=1.0/(n2*dt)
dw=twopi*df

call clogl(np,n2,s,+1.0,dt)			! FFT to freq domain

! find band of frequencies where amplitude is appreciable
smax=abs(s(2))
imax=2
do i=3,nhalf
  if(smax < abs(s(i))) then
    imax=i
    smax=abs(s(i))
  endif
enddo
smax2=smax*smax
peff=1.0/((imax-1)*df)

! define frequency limits (if you change, keep spm>slim)
slim=0.05*smax          ! write out until spectral amplitude below slim
spm=0.20*smax           ! sets a minimum amplitude for highest frequency band
if(kdiag>0) write(7,*) 'Powspec: smax,slim=',smax,slim,', peff=',peff

! find lowest acceptable frequency
i1=1
do while (abs(s(i1)).lt.slim.and.i1.lt.nhalf) 
  i1=i1+1
enddo
! find highest aceptable frequency, check if spectrum wide enough
i2=nhalf
do while (abs(s(i2)).lt.slim.and.i2.gt.1)
  i2=i2-1
enddo
i3=nhalf
do while (abs(s(i3)).lt.spm.and.i3.gt.1)
  i3=i3-1
enddo
if(i2-i1 < 5) stop 'Error in powspec: spectrum too narrow'
fmax=(i3-1)*df
if(kdiag>0) write(7,*) 'Powspec: fmax=',fmax,' Hz'
pmin=1.0/fmax

if(iwrite.eq.0) return

write(11,*) i2-i1+1
do i=i1,i2
  write(11,*) (i-1)*dw,abs(s(i))**2/smax2
  write(15,*) (i-1)*df,abs(s(i))**2/smax2  ! GMT file uses Hz for frequency
enddo  
write(15,fmt='(a)') rgb

return
end

subroutine startgmt(k)
! write first few lines of GMT script to unit k
write(k,fmt='(a)') '#! /bin/csh'
write(k,fmt='(a)') 'gmtset PAPER_MEDIA letter+'
write(k,fmt='(a)') 'gmtset ANOT_FONT_SIZE 10 LABEL_FONT_SIZE 10 HEADER_FONT_SIZE 14'
write(k,fmt='(a)') 'gmtset MEASURE_UNIT cm'
return
end
