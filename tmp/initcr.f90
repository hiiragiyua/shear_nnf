! ********************************************************************/
! 
!                  initializes everything 
!                    
! ********************************************************************/
subroutine initcr(iopt_ini)

  use ctes
  use statistics
  use running
  use bcs
  use temp

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer idat(50),ndat,nidat,iopt_ini
  real*8  dat(50)
  character*128 text, ctmp

  integer mybu,mzbu,i,j,k,iproc,nspsize,n1,n2
  ! --------------------------------------------------
  ! set paramter shared by 'use ctes',  see mod.f90 
  pi  = 4d0*atan(1d0)
  pi2 = 2d0*pi
  cii = dcmplx(0.d0,1.d0)
  ! set the default run options
  if ((noescru.eq.1).or.(nohist.eq.1).or.(nocf.eq.1).or.(iskip_screenout.eq.1)) then
     write(*,*) 'error, noescru, nohist, nocf, iskip_screenout should be set after calling initcr'
     stop
  end if
  ! set run output control options
  noescru=0 ! 1: for no output and skip allocation of sp, spl
  nohist =0 ! 1: skip screen output and writing cf files 
  nocf = 0 ! 1, skip writeing cf files (only for nohist = 0 )
  iskip_screenout = 0
  !
  iohre=19

  iinp=32    ! file numbers 
  ispf=35

  isn=33
  iout=31

  iocf=39
  iocfb=40
  iocfy=41   ! for local profile
  iotfb=42 ! for temperature

  filext='' ! default extension is nothing

  if (iopt_ini.eq.0) then 
     !    
     if(myid.eq.0) then
        !open(iohre,file='hre.dat',status='old')
        write(*,*)  'Note: reading hre2.dat !!!!'
        open(iohre,file='hre2.dat',status='old')

        ndat  = 0
        nidat = 0
        ! -------- mgalx, my, mgalz, ik2ki, ki2ik:  idat(1:5)
        call rtext(text,dat,idat,ndat,nidat,0,5,iohre)  
        
        ! -------- Re alp Ly/pi gam s chi CFL:  dat(1:7)
        call rtext(text,dat,idat,ndat,nidat,7,0,iohre)  
        
        ! --------  nstep, nimag, nhist, readflag, ifile:  idat(6:10)
        call rtext(text,dat,idat,ndat,nidat,0,5,iohre)  
        
        ! --------   pmesp, uprim:  dat(8:9)
        call rtext(text,dat,idat,ndat,nidat,2,0,iohre)  
        
        !
        ! --------  files      
        call rtext(filout,dat,idat,ndat,nidat,0,0,iohre)  
        call rtext(filinp,dat,idat,ndat,nidat,0,0,iohre)  
        
        close(iohre)
        
        do iproc=1,numerop-1
           
           call MPI_SEND(ndat,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
           call MPI_SEND(nidat,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
           call MPI_SEND(dat,ndat,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
           call MPI_SEND(idat,nidat,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
           
        enddo
        
     else
        
        call MPI_RECV(ndat,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
        call MPI_RECV(nidat,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
        call MPI_RECV(dat,ndat,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
        call MPI_RECV(idat,nidat,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
        
     endif

     Re    = dat(1) ! 1/nu
     alp   = dat(2)
     !Ly    = pi*dat(3)
     !bet   = dat(3)
     if (myid.eq.0) write(*,*)  ' Ly is read from hre2.dat'
     !
     Ly = dat(3)
     !bet = pi2/Ly
     if (myid.eq.0) write(*,*)  ' Ly = ', Ly
     
     ! originally dat(3) is assumed to be Ly/pi
     if ((Ly.gt.(0.636)).and.(Ly.lt.(0.637))) then 
        if(myid.eq.0) write(*,*) '(rev.583) Ly ERROR: check hre2.dat', &  
             'you need to input Ly at the column of Ly/pi', Ly
        stop
     endif
     
     ! the different Ly (say, 1.0 etc) may lead to big boundary error  
     gam   = dat(4)
     
     s     = dat(5)
     chi   = dat(6) ! rotation is not implemented
     cfl   = dat(7)
     pmesp = dat(8)
     uprim = dat(9)
     
     mgalx = idat(1)   
     my    = idat(2)   
     mgalz = idat(3)   
     blockingik2ki = idat(4) 
     blockingki2ik = idat(5)
     nstep = idat(6)
     nimag = idat(7)
     nhist = idat(8)
     readflag = idat(9)  
     ifile = idat(10)

  else if (iopt_ini.eq.1) then

     if(myid.eq.0) then
        write(*,*) 'skip reading hre2.dat'
        write(*,*) 'get params from header, set filinp'
        ctmp=''
        read(*,'(a)') ctmp
        filinp=trim(ctmp)
        ctmp=''
        write(*,*) 'set output basefile name, including output dir'
        read(*,'(a)') ctmp
        filout=trim(ctmp)
     end if 
     call  MPI_BCAST(filinp,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
     call  MPI_BCAST(filout,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 
     call get_header()
     Re=Ree
     alp=alpe
     Ly=Lye
     gam=game
     s=se
     chi=chie

     cfl=0.6 ! set as a default, later it can be overwritten in getfil(...) 
  
     pmesp=1.d0; uprim=1.d0;
     !mgalx=mgalxe
     !my=mye
     !mgalz=mgalze
     blockingik2ki = 8
     blockingki2ik = 8
     !
     nstep=2000000000
     nimag=100000
     nhist=1
     readflag=1
     ifile=0

  else 
     write(*,*) 'not implemented, initcr (iopt),iopt_ini=',iopt_ini
     stop

  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  id22=ifile

  ! --- this used to be ctes ---

  if ((mgalx.eq.0).or.(my.eq.0).or.(mgalz.eq.0)) then 
     write(*,*) 'mgalx grid error',mgalx,my,mgalz
  end if
  mx = 2*(mgalx/3) 
  mz = 2*(mgalz/3)-1
  
  mgalx1=mgalx-1
  mgalz1=mgalz-1
  mx1=mx/2-1 
  my1=my-1
  mz1=mz-1
  
  mgx=mgalx/2
  nz=(mz-1)/2
  nz1=nz
  nz2=mgalz-nz

  ! for dealiasing in y 
  mgy=my/2;
  my23 = 2*(my/3);
  ny=((my23-1)-1)/2
  ny1=ny
  ny2=my-ny

  mgalzp = mgalz/numerop+1
  myp    = my/numerop+1
  mzp    = mz/numerop+1
  mgalz1p= mgalzp-1
  nxymax = max(my*mzp,myp*mz)

  ! ---------------  compute coordinates, and modes values for FOURIER
  ! ---------------  this used to be coor
  allocate( xalp(0:mx1),xgam(0:mz1),alp2(0:mx1),gam2(0:mz1),icx(0:mz1) )
  allocate( y(0:my1))  ! where is deallocate??

  ! 
  call derivadas(my,dy) ! using 'diffy.f90' ! moved into subroutine xyz_derivadas ... rev.819 ! moved here at rev.1265
  call set_xyz_dy ! rh1, rh2, and h is set here, not in the previous derivadas, there is no problem.
  call filteras(my)
  !write(*,*) myid,'set xyz dy is done'
  !
  ! ---- set some test options ---- [moved from main.f90 (rev.561)]
  ! set end time as a big
  etime=1.d10
  dumptime=1.d10
  dtimag=etime
  timep_sh=(2.d0*pi/alp)/(abs(s)*Ly) ! set initial shearing-box period 

  ! --------------- this used to be pointers --
  allocate (jbeg(0:numerop-1),jend(0:numerop-1))
  allocate (kbeg(0:numerop-1),kend(0:numerop-1))   
  n1=my/numerop
  n2=my-numerop*n1

  jbeg(0)=0
  do i=0,n2-1
     jend(i)  = jbeg(i)+n1
     jbeg(i+1)= jend(i)+1
  enddo
  do i=n2,numerop-2
     jend(i)  = jbeg(i)+n1-1
     jbeg(i+1)= jend(i)+1
  enddo
  jend(numerop-1)=jbeg(numerop-1)+n1-1
  
  n1=mz/numerop
  n2=mz-numerop*n1

  kbeg(0)=0
  do i=0,n2-1
     kend(i)  = kbeg(i)+n1
     kbeg(i+1)= kend(i)+1
  enddo
  do i=n2,numerop-2
     kend(i)  = kbeg(i)+n1-1
     kbeg(i+1)= kend(i)+1
  enddo
  kend(numerop-1)=kbeg(numerop-1)+n1-1

  jb=jbeg(myid)
  je=jend(myid)
  kb=kbeg(myid)
  ke=kend(myid)

  if (numerop.gt.my) then
     if(myid.eq.0) write(*,*) 'Please check number of process (RANKs) < my'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     stop
  endif
  if (numerop.gt.(mz-1+1)) then
     if(myid.eq.0) write(*,*) 'Please check number of process (RANKs) < mz'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     stop
  endif
  if ((mod(numerop,5).eq.0).or.(mod(numerop,9).eq.0)) then
     if(myid.eq.0) write(*,*) 'Please check number of process (RANKs): divisor should be 2 or 3'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     stop
  endif

  buffsize = mx*max(maxval(jend-jbeg+1)*mz,maxval(kend-kbeg+1)*my)

  if (myid.eq.0) then
     write(*,*) 'pointers: myid, jb,  je,   my'
     write(*,'(x,4i10)')(i,jbeg(i),jend(i),my,i=0,numerop-1)
     write(*,*) 'pointers: myid, kb,  ke,   mz'
     write(*,'(x,4i10)')(i,kbeg(i),kend(i),mz,i=0,numerop-1)
  endif

  mmz = ke-kb+1
  mmy = je-jb+1

  !    ------------  initializes fast fourier transforms and CFDiff ----
  call cfti(mgalz)
  call rfti(mgalx) ! note: initialized again in getini_TG and so on...
                   !       clean-up is added from rev.395 
  call ffti(mgalx,mgalz) ! for threaded 2D-FFT
  
  call rfti_y(my) 
  call cfti_y(my)
  ! --------------  initialize stats -------------
  !! um vm wm up vp wp uvr uwr vwr o1m o2m o3m o1p o2p o3p eps uuv wwv vvv    
  !! 1  2  3  4  5  6  7   8   9   10  11  12  13  14  15  16  17  18  19

  !nstat = 19 
  !nstat = 20 !(rev.656)    20, Pk
  nstat = 23 !(rev.1230)    20, Pk; 21, ox*oy; 22, nuty; 23, disp_eddy
  allocate(stats(nstat,0:my1))
  allocate(stats2(nstat,0:my1))
  nacum = 0
  stats = 0.d0
  stats2 = 0.d0 ! to multiply dt
  if (itemperature.eq.1) then
     ntstat = 5 ! (rev.1708)
     allocate(tstats(ntstat,0:my1))
     tstats = 0.d0
  end if
  ! later set stime0 after reading a file...
  !--------------  allocates spectra:  DO IT !!!
  !  nspecy = 180
  !  nspecy = my/2 ! old format for channel flow 
  nspecy = my ! donot touch this

  ! --------------  initialize spectra ------------
  !! u v w uv o1 o2 o3 o1*o2
  !! 1 2 3 4  5  6  7  8
  nspec = 8
  ! t ut vt wt
  ! 1  2  3  4
  nspect= 4
  if (noescru.eq.0) then
     !allocate(sp(nspec,0:mx1,0:nz1),spl(nspec,0:mx1,0:nz1)) ! assuming symmetry
     allocate(sp(nspec,0:mx1,0:mz1),spl(nspec,0:mx1,0:mz1))
     allocate(sp2(nspec,0:mx1,0:mz1),spl2(nspec,0:mx1,0:mz1))
     nacumsp = 0
     sp = 0.0;  spl = 0.0
     sp2 = 0.0;  spl2 = 0.0
     if (itemperature.eq.1) then
        allocate(spt(nspect,0:mx1,0:mz1),splt(nspect,0:mx1,0:mz1))  
        spt = 0.d0; splt = 0.d0;
     end if

  end if
  ! ------------------ MPI Datatypes ------------------------
  call comm_setup()

!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
  allocate(shb(0:mx1),sht(0:mx1)) 
  allocate(shwkb(0:mx1),shwkt(0:mx1)) ! changed to 1d array from rev.385 

  xofft = 0.d0
  xoffb = 0.d0
  zofft = 0.d0
  zoffb = 0.d0

  xwkb = 0.d0
  xwkt = 0.d0
  zwkb = 0.d0
  zwkt = 0.d0

  shb   = dcmplx(1.d0,0.d0); sht   = dcmplx(1.d0,0.d0)
  shwkb = dcmplx(1.d0,0.d0); shwkt = dcmplx(1.d0,0.d0)
  ! --------------  write header for output -------------
  if(myid.eq.0) then

     write(*,*)
     write(*,'(a,f10.2,a8,f10.2,a8,f10.2)') '  Re(1/nu) =',Re, ' Rey =',Re*Ly*Ly*s,'Rez =',Re*Lz*Lz*s
     if (itemperature.eq.1) then
        write(*,'(a,f10.2,a8,f10.2)') '  Pr =',1.d0/Re/fkappa, ' Fr =',gbeta/abs(gbeta)/sqrt(gbeta)
     end if

     write(*,'(a7,f8.3,a8,f6.3,a10,f7.3,a9,f6.3)') 'alp =',alp,'gam =',gam,'Lx/Lz =',gam/alp,'Ly/Lz =',Ly/Lz

     write(*,'(a8,i5,a8,i5,a8,i5)') 'mgalx =',mgalx,'mgalz =',mgalz,'my =',my
     write(*,'(a8,i5,a8,i5,a8,i5)') 'mx =',mx,'mz =',mz
     write(*,*)

     write(*,'(a10,i7,a9,i6,a9,i5)') 'nstep =',nstep, & 
                                     'nimag =',nimag,'nhist =',nhist
     write(*,'(a8,f5.2)') '  CFL =',CFL
     write(*,*)

     write(*,*) ' shear-periodic b.c.'
     write(*,'(a10,f5.2,a8,f5.2,a8,i5)') '  shear = ',s,'chi = ',chi
     write(*,*)

     if (readflag.eq.0) then
        write(*,*) 'creating initial conditions, pmesp:', pmesp, &
                                             '  uprime:', uprim
     else
        write(*,'(a,a)') 'reading from:  ',trim(filinp)
     endif
     write(*,'(a,a)') '  write in :  ',trim(filout)

     write(*,'(a,G9.2)') ' HST-DNS, requests (MB/proc): ', & 
          & 8d0*(dfloat(buffsize)*dfloat(10 + 4) + dfloat(my*11) + dfloat(mx*mz*7))/1024d0/1024d0

  endif

end subroutine initcr

!   ----------  read input skipping comment 
subroutine rtext(text,dat,idat,ndat,nidat,kdat,kidat,iohre)
  implicit none
  character*128 text
  real*8  dat(*)
  integer idat(*),ndat,nidat,kdat,kidat
  integer iohre

  ! 
  do while (.true.) 
     read(iohre,'(a)') text
     if((text(1:2).eq.'cc').or.(text(1:2).eq.'CC')) then
        continue ! read next line
     else
        exit
     endif
     !if(text(1:2).ne.'CC') exit
  enddo
  
  ! 
  if (kdat.ne.0) then
     read(text,*) dat(ndat+1:ndat+kdat)
     ndat = ndat+kdat
  endif

  ! 
  if (kidat.ne.0) then
     read(text,*) idat(nidat+1:nidat+kidat)
     nidat = nidat +kidat
  endif

end subroutine rtext

subroutine set_xyz_dy

  use ctes
  use bcs
  use running
  use LES,only: iuse_LES,ifix_CsDeltag,Cles,Deltag,CsDeltag_fix
  use temp,only: fkappa

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer i,j,k,iproc

  if(myid.eq.0) write(*,*) 'update Lx,Lz,xalp,xgam,y,dx,dy,dz,cfl,update_dy'

  Lx=pi2/alp
  Lz=pi2/gam
  Lx2=pi/alp
  Lz2=pi/gam

  xalp = 0.d0; xgam = 0.d0; alp2 = 0.d0; gam2 = 0.d0
  y = 0.d0

  do k=0,nz
     xgam(k) = dcmplx(0.d0,gam*dfloat(k))
     icx(k) = k
  enddo
  
  do k=nz+1,mz1
     xgam(k) = dcmplx(0.d0 ,-gam*dfloat(mz-k))
     icx(k)  = mz-k
  enddo
  
  do i=0,mx1
     xalp(i) = cii*alp*dfloat(i)
  enddo
  alp2 = -xalp**2 ! positive
  gam2 = -xgam**2 ! positive
  
  dx  = pi2/alp/dfloat(mx)
  dz  = pi2/gam/dfloat(mz)
  dy  = Ly/dfloat(my)  ! Ly = 2.0 
  do j=0,my1
     y(j) = j*dy -0.5d0*Ly
  enddo
  !
  cflx = pi/CFL/dx         ! --- spectral
  cflz = pi/CFL/dz         ! --- spectral
  cfly = 2.5d0/CFL/dy      ! --- CFD for hepta-solver?
  !cfly = 2.d0/CFL/dy      ! --- CFD for penta-solver?
  !
  cflvv = 1.d0/re/CFL*max(pi*pi/dx/dx, &
       6.25d0/dy/dy, &
       pi*pi/dz/dz) ! --- for viscous CFL
  if (itemperature.eq.1) then
     cflkk = fkappa/CFL*max(pi*pi/dx/dx, &
          6.25d0/dy/dy, &
          pi*pi/dz/dz) ! --- for thermal diffusion CFL          
  end if

  call update_dy(dy) ! using 'diffy.f90'

  if (iuse_LES.eq.1) then
     ! this does not go into at the first call of initcr...
     ! but this is necessary for arclength updating Ly...
     ! (rev.1299)
     Deltag =(dx*dy*dz)**(1.d0/3.d0)
     if (myid.eq.0) write(*,*) ' iuse_LES, update Deltag =', Deltag

     if (ifix_CsDeltag.eq.1) then 
        Cles=CsDeltag_fix/Deltag
        if (myid.eq.0) write(*,*) ' ifix Cs*Deltag=', CsDeltag_fix, ' Cles=', Cles
     endif
  end if

end subroutine set_xyz_dy

! ====================================================================
!
!                   create initial data field 
!                                             jjs  22/8/07
! ====================================================================
subroutine getini_original(vor,phi,u00,w00,tiempo)

  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0

  ! ----------------   generate initial conditions ------------------
  !
  !         define vor <- random junk, spectrum given by sabs  (cosines)
  !                phi <- random junk, spectrum given by sabs  (sines)
  !            u00,w00 <- random junk, spectrum given by sabs  (cosines)  
  ! -----------------------------------------------------------------

  allocate ( buff1(0:2*my1+1), buff2(0:2*my1+1) ) ! for cos series ??
  buff1=0d0; buff2=0d0

  call rfti(2*my-2)
  ener = 0.d0
  bet  = pi/Ly

  iran=-123456
  
  if (myid.eq.0) then     ! the zero mode is done only once, by master
     write(*,*) ' Initial condition: original spectral profile'
     ! ------- the 00 velocities are  zero at wall to fix initial velocity jump 
     kx2 = 0.d0
     kz2 = 0.d0
     do j=1,my-2
        ky2  = (bet*j)**2
        amp  = sabs(kx2,ky2,kz2)
        ener = ener+amp**2
        buff1(2*j) = sign(amp,randu(iran)-.5d0)/ky2
        buff2(2*j) = sign(amp,randu(iran)-.5d0)/ky2
     enddo

     buff1(4:2*my-4:2) = buff1(4:2*my-4:2) - sum(buff1(0:2*my1+1:4))/size(buff1(0:2*my1+1:4))  
     buff1(2:2*my-4:2) = buff1(2:2*my-4:2) - sum(buff1(2:2*my1+1:4))/size(buff1(2:2*my1+1:4))
     buff2(4:2*my-4:2) = buff2(4:2*my-4:2) - sum(buff2(0:2*my1+1:4))/size(buff2(0:2*my1+1:4))  
     buff2(2:2*my-4:2) = buff2(2:2*my-4:2) - sum(buff2(2:2*my1+1:4))/size(buff2(2:2*my1+1:4))

     call rft(buff1,1,1,1)
     call rft(buff2,1,1,1)
     u00 = buff1(0:my1)
     w00 = buff2(0:my1)
!     write(*,*) 'my',my
 
     do iproc=1,numerop-1
        call MPI_SEND(u00,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w00,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
     enddo

  else       ! -------- this is done by all slaves

     !write(*,*) 'receive u00',myid,my
     call MPI_RECV(u00,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(w00,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     !write(*,*) 'received u00',myid

  endif

  ! -------------------  the higher harmonics -----------------------
  i1=0 ! bug fixed (rev.498)
  if (kb.eq.0) i1 = 1    ! -- skip the 00 mode
  do k= kb,ke
     kz2 = gam2(k)
     do i=i1,mx1
        kx2 = alp2(i)
        k13 = kx2+kz2
        buff1=0d0
        buff2=0d0
        do j=1,my-2                                ! phi (sines)
           ky2  = (bet*j)**2
           amp  = sabs(kx2,ky2,kz2)
           ener = ener+amp**2
           buff1(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)            
           buff2(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)
        enddo
        call rft(buff1,1,1,1)
        call rft(buff2,1,1,1)
        phi(:,i,k) = dcmplx(buff1(0:my1),buff2(0:my1))

        buff1=0d0
        buff2=0d0
        do j=0,my-2                               ! vor (cosines)
           ky2 = (bet*j)**2
           amp = sabs(kx2,ky2,kz2)
           buff1(2*j  ) = sign(amp,randu(iran)-.5d0) 
           buff2(2*j  ) = sign(amp,randu(iran)-.5d0)
        enddo
        call rft(buff1,1,1,1)
        call rft(buff2,1,1,1)
        vor(:,i,k) = dcmplx(buff1(0:my1),buff2(0:my1))! *sqrt(k13)
     enddo
     i1 = 0
  enddo

  call rfti(mgalx)
  deallocate(buff1,buff2)

  !    -------- normalise energy -----------
  amp=0.d0 ! rev.498
  call MPI_ALLREDUCE(ener,amp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  ener = uprim/sqrt(amp)

  vor  = vor*ener
  phi  = phi*ener
  u00  = u00*ener
  w00  = w00*ener

  tiempo =0d0
     
end subroutine getini_original

! ====================================================================
!
!                   create initial data field 
!                                             jjs  22/8/07
!                       modified by sekimoto
! ====================================================================
subroutine getini(vor,phi,u00,w00,tiempo)

  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0

  ! ----------------   generate initial conditions ------------------
  !
  !         define vor <- random junk, spectrum given by sabs  (cosines)
  !                phi <- random junk, spectrum given by sabs  (sines)
  !            u00,w00 <- random junk, spectrum given by sabs  (cosines)  
  !    modified for periodic in y-dir by sekimoto 2011/09/15
  ! -----------------------------------------------------------------

  allocate ( buff1(0:my1+2), buff2(0:my1+2) ) ! mod
  buff1=0d0; buff2=0d0

  ! NOTE: we need (my+2) buffer for calling rft(my,...)
  call rfti(my)
  ener = 0.d0
  !bet  = pi/Ly
  bet  = 2.0*pi/Ly

  iran=-123456
  
  if (myid.eq.0) then     ! the zero mode is done only once, by master
     write(*,*) ' Initial condition: spectral profile for periodic in y'
     ! ------- the 00 velocities are  zero at wall to fix initial velocity jump 
     kx2 = 0.d0
     kz2 = 0.d0
     do j=1,my-2
!     do j=1,my1/4
!     do j=1,3
        ky2  = (bet*j)**2
        amp  = sabs(kx2,ky2,kz2)
        ener = ener+amp**2
        !buff1(2*j) = sign(amp,randu(iran)-.5d0)/ky2
        !buff2(2*j) = sign(amp,randu(iran)-.5d0)/ky2
        buff1(j) = sign(amp,randu(iran)-.5d0)/ky2 
        buff2(j) = sign(amp,randu(iran)-.5d0)/ky2 
     enddo

!     buff1(4:2*my-4:2) = buff1(4:2*my-4:2) - sum(buff1(0:2*my1+1:4))/size(buff1(0:2*my1+1:4))  
!     buff1(2:2*my-4:2) = buff1(2:2*my-4:2) - sum(buff1(2:2*my1+1:4))/size(buff1(2:2*my1+1:4))
!     buff2(4:2*my-4:2) = buff2(4:2*my-4:2) - sum(buff2(0:2*my1+1:4))/size(buff2(0:2*my1+1:4))  
!     buff2(2:2*my-4:2) = buff2(2:2*my-4:2) - sum(buff2(2:2*my1+1:4))/size(buff2(2:2*my1+1:4))

     buff1(2:my-2) = buff1(2:my-2) - sum(buff1(0:my1+1:2))/size(buff1(0:my1+1:2))  
     buff1(1:my-2) = buff1(1:my-2) - sum(buff1(1:my1+1:2))/size(buff1(1:my1+1:2))
     buff2(2:my-2) = buff2(2:my-2) - sum(buff2(0:my1+1:2))/size(buff2(0:my1+1:2))  
     buff2(1:my-2) = buff2(1:my-2) - sum(buff2(1:my1+1:2))/size(buff2(1:my1+1:2))

     buff1(my/2:my-2)=0.d0 ! should be checked
     buff2(my/2:my-2)=0.d0

     call rft(buff1,1,1,1)
     call rft(buff2,1,1,1)
     u00 = buff1(0:my1)
     w00 = buff2(0:my1)
!     write(*,*) 'my',my
 
     do iproc=1,numerop-1
        call MPI_SEND(u00,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w00,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
     enddo

  else       ! -------- this is done by all slaves

     !write(*,*) 'receive u00',myid,my
     call MPI_RECV(u00,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(w00,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     !write(*,*) 'received u00',myid

  endif

  ! -------------------  the higher harmonics -----------------------
  i1=0
  if (kb.eq.0) i1 = 1    ! -- skip the 00 mode
  do k= kb,ke
     kz2 = gam2(k)
     do i=i1,mx1
        kx2 = alp2(i)
        k13 = kx2+kz2
        buff1=0d0
        buff2=0d0
!        do j=1,my-2                                ! phi (sines)
        do j=1,my1/2
           ky2  = (bet*j)**2
           amp  = sabs(kx2,ky2,kz2)
           ener = ener+amp**2
           !buff1(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)            
           !buff2(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)

           buff1(j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)            
           buff2(j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)
           !buff1(j+1) = amp ! *(k13+ky2)            
           !buff2(j+1) = amp ! *(k13+ky2)

        enddo
        call rft(buff1,1,1,1)
        call rft(buff2,1,1,1)
        phi(:,i,k) = dcmplx(buff1(0:my1),buff2(0:my1))

        buff1=0d0
        buff2=0d0
!        do j=0,my-2                               ! vor (cosines)
        do j=0,my1/2
           ky2 = (bet*j)**2
           amp = sabs(kx2,ky2,kz2)
           !buff1(2*j  ) = sign(amp,randu(iran)-.5d0) 
           !buff2(2*j  ) = sign(amp,randu(iran)-.5d0)

           buff1(j) = sign(amp,randu(iran)-.5d0)   
           buff2(j) = sign(amp,randu(iran)-.5d0) 
           !buff1(j) = amp   
           !buff2(j) = amp 
        enddo
        call rft(buff1,1,1,1)
        call rft(buff2,1,1,1)
        vor(:,i,k) = dcmplx(buff1(0:my1),buff2(0:my1)) !*sqrt(k13)
     enddo
     i1 = 0
  enddo

  call rfti(mgalx)
  deallocate(buff1,buff2)

  !    -------- normalise energy -----------
  amp=0.d0 ! rev.498
  call MPI_ALLREDUCE(ener,amp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  ener = uprim/sqrt(amp)

  vor  = vor*ener
  phi  = phi*ener
  u00  = u00*ener
  w00  = w00*ener

  tiempo =0d0
     
end subroutine getini

subroutine getini_tophat(vor,phi,u00,w00,chwk,tiempo,iopt)

  use ctes
  use bcs
  use LES,only :iuse_LES

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,iran1,iran2,iran3,i1,ik,kmax,irand

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 chwk(buffsize)

  complex*16 argc1, argc2, alphax, betax, ctmp1, ctmp2
  complex*16 u1, u2, u3, div
  real*8 sabs,ener,kx,ky,kz,kx2,ky2,kz2,amp,k13,bet,sabsexp, maxdiv
  real*8 randu,rand,voljac,euu,evv,eww, euu0,evv0,eww0,ephi,evor,ephi0,evor0
  real*8 sqk, krmax, kr, kd, frac, ekx, angphi,kxy,rtmp,mask
  real*8 tiempo, randa(3)

  real*8, allocatable:: efx(:)
  real*8, allocatable:: buff1(:), buff2(:)
  complex*16, allocatable:: cbuff1(:), cbuff2(:)
  complex*16, allocatable:: u00c(:), w00c(:)

  real*8 ekin
  !complex*16 v55c(0:my1)

  integer iopt
  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  !
  ! ----------------   generate initial conditions ------------------
  !
  !         define vor <- random junk, spectrum given by sabs  (cosines)
  !                phi <- random junk, spectrum given by sabs  (sines)
  !            u00,w00 <- random junk, spectrum given by sabs  (cosines)  
  !    modified for periodic in y-dir by sekimoto 2011/09/15
  !    top-hat (square-pulse) spectrum version (16<=k<=32, Roger and Moin)
  !    (1/2)*E(k) = 1 (16<=k<=32).
  !    tested by matlab/test_initial_condition.m (icase==5)
  ! BUGGY: this depends on the number of cores!! (now, rev.550)
  ! -----------------------------------------------------------------

  if (myid.eq.0) write(*,*) 'creat a field from tophat spectrum'

  ! NOTE: we need (my+2) buffer for calling rft(my,...)
  call rfti(my)
  call cfti(my)
  ! Here 'my' should be equal to magalx and mgalz...

  ener = 0.d0
  !bet  = pi/Ly
  bet  = 2.0*pi/Ly

  iran = -123456
  iran1 = -123456+myid
  iran2 = -134567+myid
  iran3 = -145678+myid
  
  if (mod(numerop,9).eq.0) then 
     if (myid.eq.0) write(*,*) 'random factor highly affects on umax for the multiple of 9'  
     !stop
  endif
  ! forced to zero to initialization 
  u00 = 0.d0
  w00 = 0.d0
  !
  allocate(efx(0:my1))
  !
  efx  = 0.d0
  !
  allocate ( buff1(0:my1+2), buff2(0:my1+2) ) ! mod
  buff1=0d0; buff2=0d0

  !write(*,*) myid, 'uprim', uprim
  ! set top hat energy spectram 
  !ekin(16:32)=1.d0*uprim;  ! Eii(k)/2 = 1.d0 (16 <= k <= 32) for rogallo
  !                         ! uprim is from hre.dat
  !                         !        uprim = (1/ss)**2 to be s=1.0
  ! now ekin is a pulse function
  !iopt=2;
  ! -------------------  the higher harmonics -----------------------
  allocate ( cbuff1(0:my1), cbuff2(0:my1) )
  allocate ( u00c(0:my1), w00c(0:my1) )
  cbuff1=0d0; cbuff2=0d0
  u00c=0d0; w00c=0d0

  !if (my.ne.mgalz) then
     ! for the Fourier transformation (cft)
     !write(*,*) 'isotropic flow can not be generated, my should be equal to mgalz'
     !stop
  !endif

  !i1=0
  !if (kb.eq.0) i1 = 1    ! -- skip the 00 mode
  
  ! for dealiasing (my=mgalx=mgalz)
  if (iuse_LES.eq.0) then
     krmax=2.d0/9.d0
     kmax=sqrt(2.)/3.d0*max(mgalx*alp,my*bet,mgalz*gam) 
  elseif (iuse_LES.eq.1) then
     krmax=1.d0/9.d0
     !krmax=sqrt(2.d0)/2.d0
     kmax=0.25*max(mgalx*alp,my*bet,mgalz*gam) 
  end if
  voljac=alp*bet*gam

  maxdiv=0.d0
  euu = 0.d0; evv = 0.d0; eww = 0.d0
  ephi = 0.d0; evor = 0.d0; 

  do k= kb,ke
     
     !write(*,*) myid, 'getting u,v,w from top hat Ek',k
  
     if (k.le.nz) then
        kz = gam*dfloat(k)
     else
        kz = -gam*dfloat(mz-k)
     end if  
     kz2 = gam2(k)

     do i=0,mx1
        !write(*,*) myid, 'getting u,v,w from top hat Ek',i,k
        kx = alp*dfloat(i)
        kx2 = alp2(i)

        k13 = kx2/alp/alp/dfloat(mgalx*mgalx) + kz2/gam/gam/dfloat(mgalz*mgalz)
        !k13 = kx2/alp/alp/dfloat(mx*mx) + kz2/gam/gam/dfloat(mz*mz)
        
        cbuff1=0d0; cbuff2=0d0
        do j=0,my1

           if (j.le.my/2+1) then
              ky = bet*dfloat(j)
           else
              ky = -bet*dfloat(my-j)
           end if
           ky2  = ky*ky
           
           kr = sqrt(kx2 + kz2 + ky2)  
           ! note: kr is kx in matlab test_initial_condition.m (icase==5)
           if (kr.gt.tiny) then
              kd = 1.d0/kr;
           else
              kd = 0.d0;
           end if
           ! mask for dealiasing
           sqk = k13 + ky2/bet/bet/(dfloat(my*my))
           if (sqk.gt.krmax) then  ! < 2/9 
              mask = 0.d0
           else
              mask = 1.d0
           end if
           !
           if (kr.gt.kmax) then
              ekx=0.d0;
           else              
              ekx=ekin(kr,uprim,iopt)*voljac
           end if
           if (ekx.lt.0.d0) then 
              write(*,*) 'error stop', i,ekx,frac,ekin(kr,uprim,iopt)
              stop
           end if
           efx(i)=mask*(1.d0/sqrt(2.d0*pi))*sqrt( ekx )*kd;

           ! random factor...
           !argc1 = randu(iran1)*pi2*cii;
           !argc2 = randu(iran2)*pi2*cii;
           !angphi= randu(iran3)*pi2;
           
           !
           !call init_random_seed(myid) ! GNU fortran
           !argc1 = rand(0)*pi2*cii;
           !argc2 = rand(0)*pi2*cii;
           !angphi= rand(0)*pi2;
           !
           call random_number(randa) ! intel fortran
           argc1 = randa(1)*pi2*cii;
           argc2 = randa(2)*pi2*cii;
           angphi= randa(3)*pi2;
           !
           kxy = sqrt(kx2+ky2)
           alphax=dcmplx(efx(i))*exp(argc1)*dcos(angphi);
           betax =dcmplx(efx(i))*exp(argc2)*dsin(angphi);
           !
           if (kxy.lt.tiny) then
              u1 = alphax/alp
              u2 = betax/bet
              u3 = 0.d0
           else
              ctmp1=alphax*kr;
              ctmp2=betax*kz;
              rtmp =kd/kxy;
              
              u1 = (ctmp1*ky+ctmp2*kx)*rtmp/alp;
              u2 = (ctmp2*ky-ctmp1*kx)*rtmp/bet;
              u3 = -betax*kxy*kd/gam; 

              ! this betax must have negative sign.
              ! (from P.52 rogallo (1981), P53 is not correct for e_3)
              ! here u1,u2,u3 are normalized by wavenumbers (rev.681) 
              euu = euu + 2.d0*u1*dconjg(u1)*alp*alp
              evv = evv + 2.d0*u2*dconjg(u2)*bet*bet
              eww = eww + 2.d0*u3*dconjg(u3)*gam*gam
             ! write(*,*) u1,u2,u3
           end if
           
           !
           ! checked with matlab 
           !if ((k.eq.5).and.(i.eq.5)) then
           !   write(33,*) dreal(u2),dimag(u2)
           !   v55c(j)= u2
           !end if
           
           ! check divergence zero
           ! note that u1,u2,u3 are normalized by wavenumbers (rev.681)
           div = cii*(alp*u1*kx + bet*u2*ky + gam*u3*kz) 
           
           maxdiv = max(abs(div),maxdiv)
           
           ! enfsym??
           ! u1, u2, u3 is normalized by wavenumvers.
           ! (rev.681)
           cbuff1(j) = -bet*u2*(kx2 + ky2 + kz2)
           cbuff2(j) = cii*(kz*u1*alp - kx*u3*gam)

           ephi = ephi + 2.d0*cbuff1(j)*dconjg(cbuff1(j))
           evor = evor + 2.d0*cbuff2(j)*dconjg(cbuff2(j))
           if ((i.eq.0).and.(k.eq.0)) then
              !if (j.le.my/2+1) then
              u00c(j)=u1; w00c(j)=u3; ! this must be symmetric-conjugate in y-dir
              !   if (j.ne.0) then
              !      u00c(my-j)=dconjg(u1); w00c(my-j)=dconjg(u3);
              !   end if
              !else
              !   u00c(j)=0.d0; w00c(j)=0.d0;                 
              !end if
           end if
        enddo

        ! OK 
        !if ((k.eq.5).and.(i.eq.5)) then
        !   !write(33,*) dreal(u2),dimag(u2)
        !   !v55c(j)=u2
        !   call cft(v55c,2,2,1,1) ! fou -> phy
        !   do j=0,my1
        !      write(34,*) dreal(v55c(j)),dimag(v55c(j))
        !   enddo
        !end if

        call cft(cbuff1,2,2,1,1) ! fou -> phy
        call cft(cbuff2,2,2,1,1) ! fou -> phy


        cbuff1(ny1+1:ny2)=0.d0 ! add more cut-off (rev.1446)
        cbuff2(ny1+1:ny2)=0.d0

        phi(:,i,k) = cbuff1(0:my1)
        vor(:,i,k) = cbuff2(0:my1)

        if ((i.eq.0).and.(k.eq.0)) then 
           !call rft(buff1,1,1,1)
           !call rft(buff2,1,1,1)
           !u00 = buff1(0:my1)
           !w00 = buff2(0:my1)
 
           phi(:,i,k) = 0.d0
           vor(:,i,k) = 0.d0           
           call cft(u00c,2,2,1,1) ! fou -> phy
           call cft(w00c,2,2,1,1) ! fou -> phy
           u00 = dreal(u00c)
           w00 = dreal(w00c)  

        end if
     enddo
  enddo

  call cfti(mgalz)
  call rfti(mgalx)

  deallocate(efx)
  deallocate(cbuff1,cbuff2)
  deallocate(u00c,w00c)
  deallocate(buff1,buff2)

  write(*,*) myid, 'div =', maxdiv
  if (maxdiv.gt.tiny) then
     write(*,*) myid, 'divergence of initial condition is big...stop'
  end if
  !    -------- normalise the energy -----------
  euu0 = 0.d0; evv0 = 0.d0; eww0 = 0.d0
  ephi0 = 0.d0; evor0 = 0.d0;

  call MPI_ALLREDUCE(euu,euu0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(evv,evv0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(eww,eww0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call MPI_ALLREDUCE(evor,evor0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(ephi,ephi0,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  !ener = uprim/sqrt(amp)
  !
  !vor  = vor*ener
  !phi  = phi*ener
  !u00  = u00*ener
  !w00  = w00*ener
  write(*,*) myid,' the initial condition: uu,vv,ww=', euu0, evv0, eww0
  write(*,*) myid, 'phi2,vor2',ephi0, evor0

  write(*,*) myid,' the initial condition is generated !', evv0/Lz/Lz, evor0/(s*s*Re*Lz*Lz)
  tiempo =0d0

  ! force to set conjugate relations for kx=0 modes
  ! since saving chikj2jik after hvhg keeps the initial perturbation forever 
  call chjik2ikj(phi,phi,chwk,chwk)  
  call chjik2ikj(vor,vor,chwk,chwk) ! ome2c
  call set_conj(vor,phi,-1)  ! rev(1158)
  call chikj2jik(phi,phi,chwk,chwk)
  call chikj2jik(vor,vor,chwk,chwk)
  write(*,*) myid,' set_conj (rev.1158) is applied in getini_tophat'     


end subroutine getini_tophat

subroutine getini_TG(vor,phi,u00,w00,tiempo)

  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  if(myid.eq.0)  write(*,*) ' Initial condition: Taylor-Green'
  ! ----------------   generate initial conditions ------------------
  !
  !         define vor <- random junk, spectrum given by sabs  (cosines)
  !                phi <- random junk, spectrum given by sabs  (sines)
  !            u00,w00 <- random junk, spectrum given by sabs  (cosines)  
  !    modified to periodic in y-dir by sekimoto 2011/09/15
  ! -----------------------------------------------------------------

  allocate ( buff1(0:my1+2), buff2(0:my1+2) ) 
  buff1=0d0; buff2=0d0

  ! NOTE: we need (my+2) buffer for calling rft(my,...)
  call rfti(my)
  ener = 0.d0
  bet  = 2*pi/Ly

  ! -------------------  the higher harmonics -----------------------
  i1 = 1  !  (rev.395) bug fixed for ke==mz1 
  !do k= kb,ke
  !do k= kb,min(ke,3)
  !do k= max(1,kb),1
     if (((kb.le.1).and.(ke.ge.1)).or.(ke.eq.mz1)) then  !  (rev.395) bug fixed for ke==mz1 
        k=1
        kz2 = gam2(k)
        !do i=i1,mx1-1
        do i=i1,1
           kx2 = alp2(i)
           k13 = kx2+kz2
           buff1=0d0
           buff2=0d0
           !        do j=1,my-2           ! phi (sines) ! cos-sin-cos
           !        do j=1,my1/2
           do j=2,2
              ky2  = (bet*j)**2
              amp  = sabs(kx2,ky2,kz2)
              ener = ener+amp**2
              !buff1(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)
              !buff2(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)
              
              buff1(j) = 0.d0 ! sin
              buff1(j+1) = uprim*0.5d0
              !buff2(j) = amp*0.5d0 ! sin
              !buff2(j) = 
              !buff1(j+1) = amp ! *(k13+ky2)            
              !buff2(j+1) = amp ! *(k13+ky2)
              
           enddo
           call rft(buff1,my+2,1,1)
           !call rft(buff2,1,1,1)
           !write(*,*) 'set Taylor-Green vortex', amp,buff1
           !phi(:,1,0) = buff1(0:my1)*dcmplx(0.5d0,0.d0)
           !phi(:,0,1) = buff1(0:my1)*dcmplx(0.5d0,0.d0)
           
           if ((kb.le.1).and.(ke.ge.1)) phi(:,1,1) = & 
                (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           if (ke.eq.mz1) phi(:,1,mz1) = &
                (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           
           buff1=0d0
           buff2=0d0
           !        do j=0,my-2                               ! vor (cosines) sin-cos-sin
           !        do j=0,my1/2                  
           do j=2,2
              ky2 = (bet*j)**2
              amp = sabs(kx2,ky2,kz2)
              !buff1(2*j  ) = sign(amp,randu(iran)-.5d0) 
              !buff2(2*j  ) = sign(amp,randu(iran)-.5d0)
              buff1(j) = uprim*0.5d0 ! real 
              buff1(j+1) = 0.d0 ! imag
              !buff1(j) = sign(amp,randu(iran)-.5d0)   
              !buff2(j) = sign(amp,randu(iran)-.5d0) 
              !buff1(j) = amp   
              !buff2(j) = amp 
           enddo
           call rft(buff1,my+2,1,1)
           !call rft(buff2,1,1,1)
           !vor(:,1,0) = -xgam(1)*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
           !vor(:,0,1) = -xgam(1)*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
           if ((kb.le.1).and.(ke.ge.1)) vor(:,1,1) = -gam/pi2*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
           if (ke.eq.mz1) vor(:,1,mz1) = -gam/pi2*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
        enddo
        i1 = 0
        !
     endif
  !enddo
  
  call rfti(mgalx)
  deallocate(buff1,buff2)

  !    -------- normalise energy -----------
  ! already normalized, ...
  !call MPI_ALLREDUCE(ener,amp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  !ener = uprim/sqrt(amp)
  !
  !vor  = vor*ener
  !phi  = phi*ener
  !u00  = u00*ener
  !w00  = w00*ener

  tiempo =0d0
     
end subroutine getini_TG

subroutine getini_seki(vor,phi,u00,w00,tiempo)

  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  if(myid.eq.0)  write(*,*) ' Initial condition: seki ini'
  ! ----------------   generate initial conditions ------------------
  !
  !         define vor <- random junk, spectrum given by sabs  (cosines)
  !                phi <- random junk, spectrum given by sabs  (sines)
  !            u00,w00 <- random junk, spectrum given by sabs  (cosines)  
  !    modified to periodic in y-dir by sekimoto 2011/09/15
  ! -----------------------------------------------------------------

  allocate ( buff1(0:my1+2), buff2(0:my1+2) ) 
  buff1=0d0; buff2=0d0

  ! NOTE: we need (my+2) buffer for calling rft(my,...)
  call rfti(my)
  ener = 0.d0
  bet  = 2*pi/Ly

  ! -------------------  the higher harmonics -----------------------
  i1 = 0 
  if (kb.eq.0) i1 = 1    ! -- skip the 00 mode
  !do k= kb,ke
  !do k= kb,min(ke,3)
  !do k= max(1,kb),1
  if ((kb.eq.0).or.(ke.eq.mz1)) then
     kz2 = gam2(k)
     !do i=i1,mx1-1
     do i=i1, 3
        kx2 = alp2(i)
        k13 = kx2+kz2
        buff1=0d0
        buff2=0d0
!        do j=1,my-2                                ! phi (sines) ! cos-sin-cos
!        do j=1,my1/2
        do j=2,2
           ky2  = (bet*j)**2
           amp  = sabs(kx2,ky2,kz2)
           ener = ener+amp**2
           !buff1(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)            
           !buff2(2*j+1) = sign(amp,randu(iran)-.5d0) ! *(k13+ky2)

           buff1(j) = 0.d0 ! sin
           buff1(j+1) = uprim*0.5d0
           !buff2(j) = amp*0.5d0 ! sin
           !buff2(j) = 
           !buff1(j+1) = amp ! *(k13+ky2)            
           !buff2(j+1) = amp ! *(k13+ky2)

        enddo
        call rft(buff1,my1+2,1,1)
        !call rft(buff2,1,1,1)
        !write(*,*) 'set Taylor-Green vortex', amp,buff1
        !phi(:,1,0) = buff1(0:my1)*dcmplx(0.5d0,0.d0)
        !phi(:,0,1) = buff1(0:my1)*dcmplx(0.5d0,0.d0)
        
        if (kb.eq.0) then 
           phi(:,i,1) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           !phi(:,i,2) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           !phi(:,i,3) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
        end if
        if (ke.eq.mz1) then 
           !phi(:,i,mz1-2) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           !phi(:,i,mz1-1) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
           phi(:,i,mz1) = (alp*alp+gam*gam+bet*bet)*buff1(0:my1)*dcmplx(0.5d0,0.d0)/(pi2*pi2)
        end if
        buff1=0d0
        buff2=0d0
!        do j=0,my-2                               ! vor (cosines) sin-cos-sin
!        do j=0,my1/2                  
        do j=2,2
           ky2 = (bet*j)**2
           amp = sabs(kx2,ky2,kz2)
           !buff1(2*j  ) = sign(amp,randu(iran)-.5d0) 
           !buff2(2*j  ) = sign(amp,randu(iran)-.5d0)
           buff1(j) = uprim*0.5d0 ! real 
           buff1(j+1) = 0.d0 ! imag
           !buff1(j) = sign(amp,randu(iran)-.5d0)   
           !buff2(j) = sign(amp,randu(iran)-.5d0) 
           !buff1(j) = amp   
           !buff2(j) = amp 
        enddo
        call rft(buff1,my1+2,1,1)
        !call rft(buff2,1,1,1)
        !vor(:,1,0) = -xgam(1)*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
        !vor(:,0,1) = -xgam(1)*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
        if (kb.eq.0) vor(:,i,1) = -gam/pi2*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
        if (ke.eq.mz1) vor(:,i,mz1) = -gam/pi2*buff1(0:my1)*dcmplx(0.d0,0.5d0) 
     enddo
     i1 = 0
  !enddo
  endif

  call rfti(mgalx)
  deallocate(buff1,buff2)

  !    -------- normalise energy -----------
  call MPI_ALLREDUCE(ener,amp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  ener = uprim/sqrt(amp)

  vor  = vor*ener
  phi  = phi*ener
  u00  = u00*ener
  w00  = w00*ener

  tiempo =0d0
     
end subroutine getini_seki


subroutine getini_streak(vor,phi,u00,w00,tiempo)

  ! assuming U = s*y + du*cos(gamma*z)  [du is constant in xyz-dir]
  !    => vor = -gamma*du*sin(gamma*z)
  !       phi = 0
  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  !real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  if(myid.eq.0)  write(*,*) ' Initial condition: streak profile: U = s*y + uprim*cos(gam*z)'
  if ((kb.le.1).and.(ke.ge.1)) then
     write(*,*) myid,'getini_streak: vor'
     vor(:,0,1) = vor(:,0,1) + dcmplx(0.d0,1.d0)*uprim*gam
  endif

  tiempo =0d0
     
end subroutine getini_streak

subroutine getini_streak_t(tb,tb00,tiempo)

  ! assuming U = s*y + du*cos(gamma*z)  [du is constant in xyz-dir]
  !    => vor = -gamma*du*sin(gamma*z)
  !       phi = 0
  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1

  complex*16 tb(0:my1,0:mx1,kb:ke)
  real*8 tb00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo
  !real*8, allocatable:: buff1(:), buff2(:)

  !      ---------------  zero everything, fitf
  !tb = 0.d0
  !tb00 = 0.d0
  if(myid.eq.0)  write(*,*) ' Initial condition: streak profile: U = s*y + uprim*cos(gam*z)'
  if ((kb.le.1).and.(ke.ge.1)) then
     write(*,*) myid,'getini_streak_t'
     tb(:,0,1) = tb(:,0,1) + dcmplx(0.d0,1.d0)*uprim*gam
     tb(:,1,1) = tb(:,1,1) + dcmplx(0.d0,1.d0)*uprim*gam
  endif
  if ((kb.le.0).and.(ke.ge.0)) then
     write(*,*) myid,'getini_streak_t'
     tb(:,1,0) = tb(:,1,0) + dcmplx(0.d0,1.d0)*uprim*gam
  endif
  if ((kb.le.2).and.(ke.ge.2)) then
     write(*,*) myid,'getini_streak_t'
     tb(:,1,2) = tb(:,1,2) + dcmplx(1.d0,1.d0)*uprim*0.01d0*gam
  endif
  if ((kb.le.3).and.(ke.ge.3)) then
     write(*,*) myid,'getini_streak_t'
     tb(:,2,3) = tb(:,2,3) + dcmplx(1.d0,1.d0)*uprim*0.001d0*gam
  endif
  tiempo =0d0
     
end subroutine getini_streak_t

subroutine getini_rdt2(vor,phi,u00,w00,tiempo)
  !
  ! assuming v = v0 (exp[i*( k0x x + k0z z + k_y y)] + exp[-i*( k0x x + k0z z + k_y y)])  
  ! (k_y0=0, please check for k_y0 =1, 2, ...)
  !  [du is constant in xyz-dir]
  !   (for the case of ky=0)
  !     => vor = 0
  !       phi = (\alpha**2 + \gamma**2) v0  
  !
  !   (for the case of ky.ne.0)
  !     => note: u0 = w0 = 0 does not satisfy the continuity
  !        set vor = 0, then
  !        -(alp2+gam2)*u0 = + kx*ky*v0 
  !        -(alp2+gam2)*w0 = + ky*kz*v0
  !        
  !     for the case of ky~=0 && kz==0, set additionally
  !         vor= - kx v0 exp[i*( k0x x + k_y y)]
  ! Note: the conjugate wave in x-dir will be also imposed, 
  !       so that the maximum velocity will be 2*uprim
  !
  ! normalize the initial condition to be sqrt(E_0) = uprim
  !
  ! not assuming conjugate mode by using cft (rev.673)
  !
  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1,k3
  integer k0(1:3) ! the index of wavenumbers
  real*8 kx,ky,kz,kzc ! the wavenumbers, k0(1:3).*[alp,bet,gam]

  !!complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo, ener0, sca0
  real*8, allocatable:: buff1(:), buff2(:), buff3(:)

  !      ---------------  zero everything ----
  vor = 0.d0; phi = 0.d0
  u00 = 0.d0; w00 = 0.d0
  if(myid.eq.0)  write(*,*) &
       ' Initial condition: one-wave for RDT check '
  if(myid.eq.0) then 
     write(*,*) 'input initial wavenumbers (index) '
     read(*,*) k0(1),k0(2),k0(3) 
     if ((k0(1).gt.mx1).or.(k0(2).gt.my/2).or.((k0(3).gt.mz1/2))) then
        write(*,*) 'wrong wave index !!!'
        stop
     end if
  endif  

  call MPI_BCAST(k0,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate ( buff1(0:2*my-1), buff2(0:2*my-1), buff3(0:2*my-1) ) 
  buff1=0d0; buff2=0d0; buff3=0d0
  ! NOTE: we need (my+2) buffer for calling rft(my,...), but not for cft

  bet  = 2*pi/Ly
  
  kx = alp*k0(1)
  ky = bet*k0(2)
  kz = gam*k0(3)
  kzc = -gam*(mz-k0(3))
  !
  ener0 = ((kx*ky)**2 + (kx*kx + kz*kz)**2 + (ky*kz)**2)/(kx*kx + kz*kz)**2
  sca0  = 1.d0/sqrt(ener0)
  !
  ! set etime to be |ky*dy|=pi  
  etime= (dfloat(my/2)*bet + ky)/kx;
  !write(*,*) myid,'k',kx,ky,kz,'xalp=',xalp(k0(1)),xgam(k0(3))
  
  if (k0(2).gt.0) then
     if (my/2.lt.k0(2)) then
        if (myid.eq.0) write(*,*) 'wrong initial wavenumber, my/2 < ', k0(2)
        stop
     end if
     !if ((k0(1).eq.0).or.(k0(3).eq.0)) then 
     if (k0(1).eq.0) then 
        if (myid.eq.0) write(*,*) 'wrong initial wavenumber'
        stop
     end if
     call cfti(my) ! set fast Fourier transformation in y-dir
     ! set fourier coefficient in y

     buff1(k0(2)*2)  = uprim    ! 
     buff2(k0(2)*2+1)= uprim*ky ! 
     buff3(k0(2)*2)  = -uprim*ky*ky 

     call cft(buff1(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! v0(y)
     call cft(buff2(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! dv0dy(y)
     call cft(buff3(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! ddv0dydy(y)
     
     call cfti(mgalz) ! reset rfti 

     ! zero-zero mode, should be zero ... for v00=0.d0

  elseif (k0(2).lt.0) then
     ! implemented for negative ky (rev.1169)
     if (my/2.lt.abs(k0(2))) then
        if (myid.eq.0) write(*,*) 'wrong initial wavenumber, my/2 < ', k0(2)
        stop
     end if
     call cfti(my) ! set fast Fourier transformation in y-dir
     ! set fourier coefficient in y

     buff1(2*my + k0(2)*2)  = uprim    ! 
     buff2(2*my + k0(2)*2 +1)= uprim*(-ky) ! 
     buff3(2*my + k0(2)*2)  = -uprim*ky*ky 

     call cft(buff1(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! v0(y)
     call cft(buff2(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! dv0dy(y)
     call cft(buff3(0),2,2*my,1,1) ! inverse FFT (dft_c2r) ! ddv0dydy(y)
     
     call cfti(mgalz) ! reset rfti 

     ! zero-zero mode, should be zero ... for v00=0.d0

  else
     buff1(0:2*my-1:2) = uprim 
     buff2 = 0.d0
     buff3 = 0.d0
  end if
     ! Note that phi and vor is (jik) (p-f-f)

  if ((kb.le.k0(3)).and.(ke.ge.k0(3))) then
     write(*,*) myid,'getini_rdt: phi',kb,ke
     ! bug fixed phi is originally complex array ... (rev.673)
     ! now imaginaly part is zero
     phi(0:2*my-1,k0(1),k0(3)) = phi(0:2*my-1,k0(1),k0(3)) + &
          sca0*buff1(0:2*my-1)*(-alp2(k0(1)) - gam2(k0(3))) &
          + sca0*buff3(0:2*my-1)
     write(*,*) myid,'set phi',k0(1),k0(3),phi(:,k0(1),k0(3))
     !
     ! for checking constant w0 
     !if (k0(3).eq.0) then
     !   vor(0:2*my-1,k0(1),k0(3)) = -kx*buff1(0:2*my-1) 
     !   ! set only imaginary part
     !   !write(*,*) myid,'set the imaginary part of vor',k0(1),k0(3),vor(:,k0(1),k0(3))
     !end if
  endif
  
  !k3 = mz-k0(3)
  !if ((kb.le.k3).and.(ke.ge.k3)) then
  !!   write(*,*) myid,'getini_rdt: phi',kb,ke
  !   phi(0:2*my-1,k0(1),mz-k0(3)) = phi(0:2*my-1,k0(1),mz-k0(3)) + &
  !        sca0*buff1(0:2*my-1)*(-alp2(k0(1)) - gam2(mz-k0(3))) &
  !        + sca0*buff3(0:2*my-1)
  !endif
  !   phi(:,k0(1),k3) = phi(:,k0(1),k3) + &
  !        buff1(0:my1)*(alp2(k0(1)) + gam2(k3)) + buff3(0:my1)
  !   vor(:,k0(1),k3) = vor(:,k0(1),k3) + &
  !        buff2(0:my1)*dcmplx(1.d0,0.d0)*kx/kzc
  !end if

  deallocate(buff1,buff2,buff3)
  tiempo =0d0
     
end subroutine getini_rdt2


subroutine getini_rdt(vor,phi,u00,w00,tiempo)
  !
  ! assuming v = v0 exp[i*( k0x x + k0z z + k_y y)] + exp[-i*( k0x x + k0z z + k_y y)]  
  ! (k_y0=0, please check for k_y0 =1, 2, ...)
  !  [du is constant in xyz-dir]
  !   (for the case of ky=0)
  !     => vor = 0
  !       phi = (\alpha**2 + \gamma**2) v0  
  !
  !   (for the case of ky.ne.0)
  !     => note: u0 = w0 = 0 does not satisfy the continuity
  !        set vor = 0, then
  !        -(alp2+gam2)*u0 = + kx*ky*v0 
  !        -(alp2+gam2)*w0 = + ky*kz*v0
  ! Note: the conjugate wave in x-dir will be also imposed, 
  !       so that the maximum velocity will be 2*uprim
  !
  ! normalize the initial condition to be sqrt(E_0) = uprim
  !
  use ctes
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,iproc,iran,i1,k3
  integer k0(1:3) ! the index of wavenumbers
  real*8 kx,ky,kz,kzc ! the wavenumbers, k0(1:3).*[alp,bet,gam]

  !!complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  real*8 phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)
  real*8 sabs,ener,kx2,ky2,kz2,amp,k13,bet,sabsexp
  real*8 randu
  real*8 tiempo, ener0, sca0
  real*8, allocatable:: buff1(:), buff2(:), buff3(:)

  !      ---------------  zero everything ----
  vor = 0.d0; phi = 0.d0
  u00 = 0.d0; w00 = 0.d0
  if(myid.eq.0)  write(*,*) &
       ' Initial condition: one-wave for RDT check '
  if(myid.eq.0) then 
     write(*,*) 'input initial wavenumbers (index) '
     read(*,*) k0(1),k0(2),k0(3) 
     if ((k0(1).gt.mx1).or.(k0(2).gt.my/2).or.((k0(3).gt.mz1/2))) then
        write(*,*) 'wrong wave index !!!'
        stop
     end if
  endif  

  call MPI_BCAST(k0,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate ( buff1(0:my1+2), buff2(0:my1+2), buff3(0:my1+2) ) 
  buff1=0d0; buff2=0d0; buff3=0d0
  ! NOTE: we need (my+2) buffer for calling rft(my,...)

  bet  = 2*pi/Ly
  
  kx = alp*k0(1)
  ky = bet*k0(2)
  kz = gam*k0(3)
  kzc = -gam*(mz-k0(3))
  !
  ener0 = ((kx*ky)**2 + (kx*kx + kz*kz)**2 + (ky*kz)**2)/(kx*kx + kz*kz)**2
  sca0  = 1.d0/sqrt(ener0)
  !
  !write(*,*) myid,'k',kx,ky,kz,'xalp=',xalp(k0(1)),xgam(k0(3))
  
  if (k0(2).ne.0) then
     if (my.le.k0(2)*2) then
        if (myid.eq.0) write(*,*) 'wrong initial wavenumber, my/2 < ', k0(2)
        stop
     end if
     !if ((k0(1).eq.0).or.(k0(3).eq.0)) then 
     if (k0(1).eq.0) then 
        if (myid.eq.0) write(*,*) 'wrong initial wavenumber'
        stop
     end if
     call rfti(my) ! set fast Fourier transformation in y-dir (my+2 ??)
     ! set fourier coefficient in y

     buff1(k0(2)*2)  = uprim    ! rft impose a conjugate mode. 
     buff2(k0(2)*2+1)= uprim*ky ! complex ! bug fixed (rev.427)
     buff3(k0(2)*2)  = -uprim*ky*ky 

     call rft(buff1(0),my+2,1,1) ! inverse FFT (dft_c2r) ! v0(y)
     call rft(buff2(0),my+2,1,1) ! inverse FFT (dft_c2r) ! dv0dy(y)
     call rft(buff3(0),my+2,1,1) ! inverse FFT (dft_c2r) ! ddv0dydy(y)
     
     !call rft(buff1,my+2,1,-1) ! forward FFT (dft_c2r) ! v0(y), OK
     !if (myid.eq.0) then
     !   write(*,*) myid,'buff1',buff1 ! rft impose a conjugate mode.      
     !   write(*,*) myid,'buff2',buff2
     !   write(*,*) myid,'buff3',buff3
     !end if
     
     call rfti(mgalx) ! reset rfti 

     ! zero-zero mode, should be zero ... for v00=0.d0

  else
     buff1 = uprim 
     buff2 = 0.d0
     buff3 = 0.d0
  end if
     ! Note that phi and vor is (jik) (p-f-f)

  if ((kb.le.k0(3)).and.(ke.ge.k0(3))) then
     write(*,*) myid,'getini_rdt: phi',kb,ke
     ! bug fixed phi is originally complex array ... (rev.673)
     ! now imaginaly part is zero
     phi(0:2*my-1:2,k0(1),k0(3)) = phi(0:2*my-1:2,k0(1),k0(3)) &
          + sca0*buff1(0:my1)*(-alp2(k0(1)) - gam2(k0(3))) &
          + sca0*buff3(0:my1)
     write(*,*) myid,'set phi',k0(1),k0(3),phi(:,k0(1),k0(3))
     !%vor(:,k0(1),k0(3)) = vor(:,k0(1),k0(3)) + &
     !     buff2(0:my1)*dcmplx(-1.d0,0.d0)*kx/kz  
     ! BUG u0=0 does not satisty the continuity 
     !vor(1:2*my-1:2,k0(1),k0(3)) = vor(1:2*my-1:2,k0(1),k0(3)) &
     !    - kx*buff1(0:my1)
  endif
  
  !k3 = mz-k0(3)
  !if ((kb.le.k3).and.(ke.ge.k3)) then
  !   write(*,*) myid,'getini_rdt: phi',kb,ke
  !   phi(:,k0(1),k3) = phi(:,k0(1),k3) + &
  !        buff1(0:my1)*(alp2(k0(1)) + gam2(k3)) + buff3(0:my1)
  !   vor(:,k0(1),k3) = vor(:,k0(1),k3) + &
  !        buff2(0:my1)*dcmplx(1.d0,0.d0)*kx/kzc
  !end if

  deallocate(buff1,buff2,buff3)
  tiempo =0d0
     
end subroutine getini_rdt

function ekin(ke,fac,iopt)
  implicit none
  real*8 ke,ekin,fac
  integer iopt
  ekin = 0.d0

  if (iopt.eq.1) then
     if ((ke.ge.16.d0).and.(ke.le.32.d0)) then 
        ekin=fac        
     end if
  elseif (iopt.eq.2) then
     if ((ke.ge.4.d0).and.(ke.le.8.d0)) then 
        ekin=fac        
     end if
  elseif (iopt.eq.3) then
     !if ((ke.ge.3.d0).and.(ke.le.4.d0)) then 
     if ((ke.ge.2.d0).and.(ke.le.3.d0)) then 
        ekin=fac        
     end if
   end if
  !write(*,*) 'ekin=',ekin, 'at ', ke, kmax

end function ekin


function sabs(kx2,ky2,kz2) 
!**********************************************************
!   isotropic spectrum with k^2 E peaking at pmesp
!**********************************************************
  use ctes
  implicit none
  real*8 kx2,ky2,kz2,kk,sabs,c0,cx,qx,q0
  parameter (q0=4., qx=5.d0/3.d0, c0=0.25*q0, cx=c0+0.25*qx) 
  !parameter (q0=4., qx=10.d0, c0=0.25*q0, cx=c0+0.25*qx) 

  !  -----  parameters for initial conditions
  ! isotropic energy spectrum E=k^q0     at k<<1
  ! isotropic energy spectrum E=k^(-qx)  at k>>1
  ! maximum of the initial E at k=pmesp -------

  kk   = (kx2+ky2+kz2)/(pmesp*1.36)**2
  sabs = kk**c0/(1+kk**cx)

end function sabs 



function sabsexp(kx2,ky2,kz2)
!**********************************************************
!   isotropic spectrum with E peaking at pmesp
!**********************************************************
  use ctes
  implicit none
  real*8 kx2,ky2,kz2,k,sabsexp

  k   = sqrt(kx2+ky2+kz2)*(1.41/pmesp)
  sabsexp = sqrt(k**4*exp(-k**2))/k   ! sqrt of energy density

end function sabsexp




! ********************************************************
!      portable random number generator
! ********************************************************
function randu(idum)
  implicit none
  integer   m,ia,ic,iff,idum,iy,j
  real*8    randu,rm

  parameter(m=714025, ia=1366, ic=150889, rm=1d0/m)
  integer   ir(97)
  data iff/0/
  save iff,iy,ir

  if (idum.lt.0. .or. iff.eq.0) then

     iff=1
     idum=mod(ic-idum,m)
     do j=1,97
        idum=mod(ia*idum+ic,m)
        ir(j) = idum
     enddo
     idum = mod(ia*idum+ic, m)
     iy = idum

  endif

  j     = 1+(97*iy)/m
  iy    = ir(j)
  randu = iy*rm
  idum  = mod(ia+idum+ic,m)
  ir(j) = idum

end function randu

SUBROUTINE init_random_seed(id)
  ! random seeding

  INTEGER :: i, n, clock, id
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  call random_seed(get=seed)   ! bug fixed (rev.1229)
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + (37+id*5) * (/ (i - 1, i = 1, n) /)

  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed

! ------------------------------------------------------------------------------- !
!  set additional options 
! DO NOT forget to BCAST
! ------------------------------------------------------------------------------- !
subroutine set_options()
  
  use ctes
  use running
  use LES
  use temp

  implicit none
  include "mpif.h"

  character*128 command
  integer nextline
  integer istat(MPI_STATUS_SIZE),ierr

  call hre_reader(ihre_error) ! read some run options from hre3.dat

  call MPI_BCAST(ihre_error,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  ! DO NOT forget to BCAST
  !
  nopt=9; nparams=13; ! here set dafault numbers of option and parameters (rev.1433)
  ! this is used throughout file dumping, readwrite and arnoldi 
  ! -----------------------------------------------------------------------------
  ! set default options 
  explicit=1; ifix_dt=0; iadd_force=0; iadd_mode=5; iadd_sym=0; iadd_damping=0;
  Deltat=0.d0; force_roll=0.d0; xforce=0.d0; zforce=0.d0; vbulk=0.d0; 
  damp_aa=0.d0; damp_up=0.d0;
  idump_mode=0; iget_cfy=0;
  iydealiasing=0;
  iuse_newton=0; 
  iread_footer=1; ! rev.1245
  inorm=0 ! the standard norm (1, for LES (mean velocity is scaled-average ..., 
          ! (2, devided by the maximum before the summation, then the best convergence for UPO2 (LB) 
          ! (3,testing DotK, double**K -order norm, this should be equivalent to inorm=0)
  ! -----------------------------------------------------------------------------
  iread_hdf5=0; iwrite_hdf5=0;
  !
  ! set default options for LES
  iuse_LES=0; idynamic=0; iadd_visc=0; Cles=0.25; ! 0.173 (Lilly) 
  ifix_CsDeltag=0; CsDeltag_fix=0.d0; ! rev.1348
  !
  itemperature=0 ! see crosst
  init_temp=0
  ! more parameters for LES dynamic
  cutoffx=0.d0;cutoffz=0.d0
  LESflt(1:3)=0.d0; ! the filter size in x,y,z
  if (ihre_error.eq.1) return
  ! ------------- go to a command to read options -------------------------------
  ! DO NOT forget to BCAST
  ! --- set time-stepping options (take care for using newton ... to set implicit)
  ! explicit ! 0, implicit; 1, explicit.
  if (myid.eq.0) then
     command='explicit'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) explicit
     end if
  end if
  call  MPI_BCAST(explicit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! iadd_mode  2 or 5 (linearly-integration)
  if (myid.eq.0) then
     command='iadd_mode'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iadd_mode
     end if
  end if
  call  MPI_BCAST(iadd_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! itemperature
  if (myid.eq.0) then
     command='itemperature'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) itemperature
        if (itemperature.eq.1) then
           read(hrefile(nextline+1),*) fkappa, bym, gbeta
           write(*,*) 'itemperature: fkappa, bym, gbeta ',fkappa, bym, gbeta
        end if
     end if
  end if
  call  MPI_BCAST(iTemperature,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(fkappa,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(bym,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(gbeta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! --- this option generate a streak-like initial temperature field
  command='init_temp'
  call switching(command,init_temp) 

  ! idump_mode
  ! 0 is original dumping mode, which uses iimag...
  ! 
  if (myid.eq.0) then
     command='idump_mode'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) idump_mode
        if (idump_mode.ge.1) then
           write(*,*) 'idump_mode=2: dumping files with a constant time interval'
           write(*,*) 'N-period, num. of intervals per 1-period'
           read(hrefile(nextline+1),*) dump2timep, dump2tint
        end if
     end if
  end if
  call  MPI_BCAST(idump_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(dump2timep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(dump2tint,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! iget_cfy
  ! =1 : save cf as a function of y... only for LES
  ! 
  if (myid.eq.0) then
     command='iget_cfy'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iget_cfy
        if (iget_cfy.eq.1) then
           write(*,*) 'iget_cfy=', iget_cfy
        end if
     end if
  end if
  call  MPI_BCAST(iget_cfy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! iadd_sym  !
  !           ! 5: shift-refrection
              ! 6: mirror-symmetry
              ! 7: rotation-refrection for varicose (UPO4) 
              ! 8: rotation-refrection for sinuous 
  if (myid.eq.0) then
     command='iadd_sym'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iadd_sym
        if ((iadd_sym.ge.5).and.(iadd_sym.le.8)) then
           read(hrefile(nextline+1),*) sym_shiftx,sym_shiftz 
           ! the fraction of Lx and Lz
           ! shifting parameter for initial condition to have the specified symmetry
           write(*,*) 'iadd_sym, initial phase shift, iadd_sym=',iadd_sym
           write(*,*) 'sym_shift(x,z)',sym_shiftx, sym_shiftz
        elseif (iadd_sym.eq.12) then
           read(hrefile(nextline+1),*) sym_shiftx,sym_shiftz 
           ! the fraction of Lx and Lz
           ! shifting parameter for initial condition to have the specified symmetry
           write(*,*) 'iadd_sym, initial phase shift, iadd_sym=',iadd_sym
           write(*,*) 'sym_shift(x,z)',sym_shiftx, sym_shiftz  
        elseif (iadd_sym.eq.-1) then
           ! for gmres_shear.f90
           if(myid.eq.0) write(*,*) 'add randum disturbance to initial vector', &
                ' to break symmetry of the initial guess' 
        elseif (iadd_sym.ne.0) then
           write(*,*) 'wrong parameter iadd_sym =', iadd_sym
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
        end if
     end if
  end if
  call  MPI_BCAST(iadd_sym,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(sym_shiftx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(sym_shiftz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! write(*,*) myid,'iadd_sym, sym_shiftx',iadd_sym, sym_shiftx, sym_shiftz 
  !
  ! iadd_force
  if (myid.eq.0) then
     command='iadd_force'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iadd_force
        if (iadd_force.ge.1) then
           read(hrefile(nextline+1),*) force_roll
        else if(iadd_force.eq.-1) then ! for shifting in x-, y-, z-dir by a force
           read(hrefile(nextline+1),*) xforce,vbulk,zforce
        end if
     end if
  end if
  call  MPI_BCAST(iadd_force,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(force_roll,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(xforce,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(vbulk,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(zforce,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !
  ! window-damping in y-dir
  if (myid.eq.0) then
     command='iadd_damping'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iadd_damping
        read(hrefile(nextline+1),*) damp_up, damp_aa
     end if
  end if
  call  MPI_BCAST(iadd_damping,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(damp_up,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(damp_aa,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if (iadd_force.ne.0) then
     if(myid.eq.0) write(*,*) ' iadd_force =', iadd_force
  end if
  if (iadd_damping.ne.0) then
     if(myid.eq.0) write(*,*) ' iadd_damping,up,aa =', iadd_damping,damp_up,damp_aa
  end if
  !
  ! --- LES options
  if (myid.eq.0) then
     command='iuse_LES'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iuse_LES
        read(hrefile(nextline+1),*) idynamic,iadd_visc,Cles
     end if
  end if
  call  MPI_BCAST(iuse_LES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(idynamic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(iadd_visc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(Cles,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if (iuse_LES.eq.1) then
     if (myid.eq.0) write(*,*) 'iuse_LES,idynamic,iadd_visc', iuse_LES,idynamic,iadd_visc
     ! grid scale Delta_g,
     !Deltag =(dx*dy*dz)**(1/3) ! BUG found by sekimoto 2014/Oct/21
     Deltag =(dx*dy*dz)**(1.d0/3.d0) ! do not forget to set again after changing dx,dy,dz...
     if (myid.eq.0) then
        if (idynamic.eq.0) then        
           write(*,*) ' ----- '
           write(*,*) ' LES Smagorinsky: [Deltag=(dx*dy*dz)**(1/3),Cles] ='
           write(*,*) ' ',Deltag,Cles
           write(*,*) ' ----- '
        elseif (idynamic.eq.1) then
           write(*,*)
           write(*,*) ' LES dynamic Smagorinsky model'
           write(*,*) ' NOT implemented ...stop' 
           stop
        end if
     endif
     !if (idynamic.eq.0) write(*,*) myid,'Cles,Deltag,idynamic=',Cles,Deltag,idynamic

  end if
 
  ! checking errors...
  if ((iuse_LES.eq.1).and.(explicit.eq.0)) then
     write(*,*) 'LES is only for explicit, stop'
     stop
  endif

  ! --- fix Cs*Deltag options
  command='ifix_CsDeltag'
  call switching(command,ifix_CsDeltag)  ! default is 0
  if (ifix_CsDeltag.eq.1) then
     if (myid.eq.0) write(*,*) 'ifix_CsDeltag '
  elseif (ifix_CsDeltag.ne.0) then
     if (myid.eq.0) write(*,*) 'ifix_CsDeltag: error ' 
     stop
  end if
  if (trim(command).eq.'ifix_CsDelta') then
     if (myid.eq.0) write(*,*) 'ifix_CsDeltag: require g in the command name'
     stop
  endif

  command='iydealiasing'
  if (myid.eq.0) then
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iydealiasing
     end if
  end if
  call  MPI_BCAST(iydealiasing,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (iydealiasing.eq.1) then
     if (myid.eq.0) write(*,*) 'iydealiasing=1 : 2/3-dealiasing only for kx=0'
  elseif (iydealiasing.eq.2) then
     if (myid.eq.0) write(*,*) 'iydealiasing=2 : remeshing-like dealiasing'
  elseif (iydealiasing.eq.3) then
     if (myid.eq.0) write(*,*) 'iydealiasing=3 : compact filter'
  elseif (iydealiasing.ne.0) then
     if (myid.eq.0) write(*,*) 'error reading option: iydealiasing' 
     stop
  end if

  ! --- reading options
  if (myid.eq.0) then
     command='iread_footer' ! set 0 for old format to read successfully, (machine-dependent)
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iread_footer
     end if
  end if
  call  MPI_BCAST(iread_footer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (iread_footer.eq.0) then
     if (myid.eq.0) write(*,*) 'skip reading footer in the field file'
  end if

  ! --- output options
  command='iskip_screenout'
  call switching(command,iskip_screenout)  ! default is 0
  if (iskip_screenout.eq.1) then
     if (myid.eq.0) write(*,*) 'skip screen output'
  elseif (iskip_screenout.ne.0) then
     if (myid.eq.0) write(*,*) 'skip screen output error: ',iskip_screenout 
     stop
  end if

  if (myid.eq.0) then
     command='inorm' ! default is 0, from rev1545
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) inorm
     end if
  end if
  call  MPI_BCAST(inorm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  command='iread_hdf5'
  call switching(command,iread_hdf5)  
  command='iwrite_hdf5'
  call switching(command,iwrite_hdf5)
  if ((iread_hdf5.le.-1).or.(iread_hdf5.ge.2)) then
     if (myid.eq.0) write(*,*) 'iread_hdf5 should be 0 or 1'
     stop
  end if
  if ((iwrite_hdf5.le.-1).or.(iwrite_hdf5.ge.2)) then
     if (myid.eq.0) write(*,*) 'iwrite_hdf5 should be 0 or 1'
     stop
  end if

end subroutine set_options

subroutine switching(command,iswitch)

  use ctes,only:myid
  use running,only:hrefile

  implicit none
  include "mpif.h"

  character*128 command
  integer istat(MPI_STATUS_SIZE),ierr,iswitch,nextline

  if (myid.eq.0) then

     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iswitch ! default is 0
     end if
  end if
  call  MPI_BCAST(iswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine switching

subroutine hre_reader(ierror)

  use ctes,only:myid
  use running,only:hrefile,hrelines
  implicit none

  ! read optional numerical parameters from hre3.dat
  ! the parameter is saved in 
  
  character*1024 file
  integer iohre,maxHreLine,il,nline,ierror

  if (myid.ne.0) return
  iohre=20
  maxHreLine=128

  file = 'hre3.dat'
  open(iohre, FILE=file, status='old', err=122)

  write(*,*) '----- reading running options from ',trim(file), ' -----'
  ierror=0
  allocate(hrefile(maxHreLine))
  !write(*,*) 'allocated'
  nline = 1
  read(iohre,100) hrefile(nline)
  !write(*,100) trim(hrefile(1))
  do while ((trim(hrefile(nline)).ne.'HREEND') & 
       & .and.(trim(hrefile(nline)).ne.'ENDHRE') & 
       & .and.(trim(hrefile(nline)).ne.'END'))
     nline = nline + 1
     read(iohre,100) hrefile(nline)
  end do

100 format(a128) ! input only, do not use for output
  hrelines = nline
  ! check
  !do il=1,nline
  !   write(*,*) trim(hrefile(il))
  !end do
  write(*,*) '----- successful read of hre3.dat -----'
  write(*,*) 
  close(iohre)

  return

122 ierror=1
  return

end subroutine hre_reader

subroutine go_command(command, nextline)
  
  use running,only:hrefile,hrelines
  implicit none

  integer i1, i2, line, nextline
  character*128 command, text
  
  i1=1
  i2=len_trim(command)

  !write(*,*) 'go command',trim(command),i1,i2

  do line = 1,hrelines
     !write(*,*) 'go command',line,'/',hrelines
     !read(hrefile(line),100) text
     text=hrefile(line)
     if (trim(text(i1:i2)).eq.trim(command)) then
        command = text
        nextline = line + 1
        return
     else if (trim(text).eq.'HREEND') then
        command = 'HREEND'
        nextline = -1
        return
     else if (trim(text).eq.'ENDHRE') then
        command = 'ENDHRE'
        nextline = -1
        return
     else if (trim(text).eq.'END') then
        command = 'END'
        nextline = -1
        return
     end if
  end do
100 format(a128)

end subroutine go_command

