
!#undef 
subroutine crosst(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy, &
     &           rhvc1,rhvc2,rhgc1,rhgc2, &
     &           tb,tb00,tbwk,ht,rhtc2,chwk)

  use ctes
  use running
  use statistics
  use bcs
  use timer
  use LES
  use temp
  use omp_lib

  implicit none
  include "mpif.h"

  !------------- Parameters for SMR R-K (A.A.Wray etc.)------------------!
  real(8)    gama(3), alpha(3), beta(3), ibeta(3), zeta(3), ca(3)
  parameter (gama =(/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /))
  parameter (alpha=(/ 29d0/96d0, -3d0/40d0,   1d0/6d0 /))
  parameter (beta =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /))
  parameter (ibeta=(/ 160d0/37d0, 24d0/5d0,   6d0     /))
  parameter (zeta =(/ -17d0/60d0, -5d0/12d0,  0d0     /))
  parameter (ca   =(/ 8d0/15d0,   2d0/3d0,    1d0 /))

  integer irun,rkstep,i,j,k,k1,ithread

  real(8)  r1,dtr1,rk,rk2,fac,fac_t
  real(8)  dtr1v,dtr1k

  complex(8) shp,shm,const
  real(8),dimension(0:my1) ::  u00,w00
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vor, phi, hv, hg
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vorwk, phiwk, dvordy

  integer iproc,istat(MPI_STATUS_SIZE),ierr
  integer icheck
  logical lex

  ! -- LES
  real*8  rhgc1(0:2*my-1,0:mx1,kb:ke),rhgc2(0:2*my-1,0:mx1,kb:ke)
  real*8  rhvc1(0:2*my-1,0:mx1,kb:ke),rhvc2(0:2*my-1,0:mx1,kb:ke)
  !real*8,dimension(:),allocatable:: txy00,tyz00

  ! -- itemperature
  real(8),dimension(0:my1) ::  tb00
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: tb, ht, tbwk

  ! -- itemperature
  real*8  rhtc2(0:2*my-1,0:mx1,kb:ke)

  real(8),dimension(buffsize) :: chwk


  real(8),dimension(:),allocatable :: buff1, buff2 ! used for addsym, too ...
  real(8),dimension(:),allocatable :: buff3 ! itemperature
  real(8),dimension(:),allocatable :: u00s, w00s ! only for addsym
  real(8),dimension(:),allocatable :: u00tmp, w00tmp ! for y-dealiasing
  real(8),dimension(:),allocatable :: u00c,w00c ! for y-dealiasing
  real(8),dimension(:),allocatable :: tb00c  ! itemperature

  real(8),dimension(:),allocatable :: rf0u, u00wk, rf0w, w00wk
  real(8),dimension(:),allocatable :: rf0tb, tb00wk ! itemperature
  real(8),dimension(:,:),allocatable :: u1r, u2r, u3r, o1r, o2r, o3r
  real(8),dimension(:,:),allocatable :: tbr  ! itemperature

  ! tmp 2d arrays
  real(8),dimension(:,:),allocatable:: tmpxzr ! for fouxz, phys2
  real(8),dimension(:),allocatable:: tmpx     ! 1d x-dir
  real(8),dimension(:),allocatable:: wk1dr,wk1dc ! for compact
  save wk1dr, wk1dc
  !$OMP THREADPRIVATE(wk1dr,wk1dc)

  character*4 ext1
  character*128 fname, fnamet
  ! ------------------------------------------


  t_derivyc = 0.d0; t_add_shear= 0.d0; t_uvw = 0.d0; t_copy=0.d0
  t_fourxz = 0.d0; t_uome=0.d0;  t_nl=0.d0; t_hvhg=0.d0;
  t_otra=0.d0; t_waiting=0.d0


  comm_time = 0d0
  commtimer = 0d0
  addsymtimer = 0d0
  transtimer= 0d0
  totaltimer= 0d0
  wtime = 0d0

  if (nohist.eq.1) then
     ihist  = 0
  else
     ihist  = 1
  end if

  icfl = 1
  iwrote = 0
  ifix_dt = 0  ! for debug

  !iadd_force = 0 ! add force (1, 1-pair-roll; 2, 2-pair-roll)
  ! do not use iadd_mode ==1 and 4 because hist is not updated for it
  iadd_mode=5 ! add_mode = 2 is original
              ! add_mode = 3 mapped by exp(ik_x s y t)
              ! add_mode = 4 mapped by exp(ik_x s y t) based on add_mode = 1
              ! [integrate the shear-term analytically]
              ! add_mode = 5 unmapped by exp(ik_x s y Delta_t) based on add_mode = 2
              ! this is equivalent with iadd_mode=3 of Siwei's code (2014/04/21)

              ! add_mode = -1 solves linearized Navier-Stokes eq.

  !iadd_damping = 0; ! window-damping...see Teramura & Toh (PRE,2014)
  ! set in main programs, which are main1.f90, newton_shear.f90, arnoldi_shear.f90, shoot_shear.f90
  !                     and some more for post-processing.

  irun  = 0  ! first time step is special in tim3rkp
  istop = 0  ! stop flag
  idump = 0

  if (itemperature.eq.0) then
     write(*,*) 'itemperature should be 1 in crosst.f90'
     stop
  endif
  !   if(myid.eq.0) write(*,*) ' A.A.Wray-type RK3'
  ! both options are checked with the RDT solution (rev.430)
  if (explicit.eq.1) then

     if((myid.eq.0).and.(nohist.eq.0)) write(*,*) ' EXPLICIT: compute viscous terms explicitly, iadd_mode =', iadd_mode
  else
     if((myid.eq.0).and.(nohist.eq.0)) write(*,*) ' IMPLICIT: compute viscous terms semi-implicitly, iadd_mode =', iadd_mode
     if ((iadd_mode.ge.3).and.(iadd_mode.le.4)) then
        if (myid.eq.0) write(*,*) 'STOP: iadd_mode=3, 4 are not implemented for implicit solver!!'
        stop
     endif

  end if
  if ((iadd_mode.eq.0).or.((iadd_mode.eq.1).or.(iadd_mode.eq.4))) then
     if (myid.eq.0) write(*,*) 'STOP: do not use iadd_mode=0,1,4, the statistics is accumulated by assuming iadd_mode=2,3,5'
     stop
  end if
  if (iadd_mode.eq.-1) then
     if (myid.eq.0) write(*,*) ' solving linearized-NS'
  end if

 if (iadd_damping.eq.1) then
     if (myid.eq.0) write(*,*) ' iadd_damping, damp_up,aa = ', iadd_damping, damp_up, damp_aa
  end if
  ! check: the statistics buffer should be initialized
  if (ihist.eq.1) then
     !if (abs(stime0 - time).gt.0.d0) then
     !   !write(*,*) 'the statistics may not be initialized..., stime0=',stime0
     !   !stime0 = time
     !   !stop
     !end if
     if ((nacum.gt.0).or.(nacumsp.gt.0)) then
        write(*,*) 'the statistics is not initialized...'
        stop
     end if
     if ((total_dt.gt.0.d0).or.(nacumsp.gt.0)) then
        write(*,*) 'the statistics is not initialized...'
        stop
     end if
     ! reset statistics
     if (nohist.eq.0) then
        nacum = 0; stats = 0.0; nacumsp = 0
        sp = 0.d0; spl = 0.d0
        stime0 = time
        stats2 = 0.d0; total_dt=0.d0; ! added (rev.852)
        sp2 = 0.d0; spl2 = 0.d0; ! added (rev.853)
     end if
  end if

  if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) ' cross: iadd_sym =',iadd_sym

  if (myid.eq.0) then
     write(ext1,'(i4.4)') id22
     fname=filout(1:index(filout,' ')-1)//'.'//ext1//'.'//'flqt'
     inquire(file=fname,exist=lex)
     if (lex) write(*,*) 'flqt file is detected, stop:',trim(fname)
  end if
  call MPI_BCAST(lex,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  if (lex) stop

  if ((myid.eq.0).and.(nohist.eq.0).and.(nocf.eq.0)) then
     write(ext1,'(i4.4)') id22
     fname=filout(1:index(filout,' ')-1)//'.'//ext1//'.cf'
     inquire(file=fname,exist=lex)
     if(lex) then
        write(*,*) 'Warning, the cf file will be overwritten !!'
        !stop
     end if
     open(iocf,file=fname,status='unknown')
     write(iocf,'(a,3(1X,F14.6))') '% (Axz,Ayz,Rez)=', gam/alp, Ly/Lz, re*Lz*Lz
     if (itemperature.eq.1) then
        write(iocf,'(a,3(1X,F14.6))') '% (Pr,Fr)=', 1.d0/re/fkappa, gbeta/abs(gbeta)/sqrt(gbeta)
     end if
     write(iocf,'(a,3(1X,I10))')   '% grids (phys)', mgalx,my,mgalz
     write(iocf,*) '% re,alp ',re,alp
     write(iocf,*) '% Ly,gam ',Ly,gam
     write(iocf,*) '% s,chi ',s,chi
     write(iocf,*) '% CFL = ',CFL
     write(iocf,'(a,3(1X,I10),2(1X,F14.6))') '% ', nstep,nimag,nhist,pmesp,uprim
     if (iadd_force.eq.1) then
        write(iocf,*) '% iadd_force=1',force_roll
     elseif (iadd_force.eq.-1) then
        write(iocf,*) '% iadd_force=-1 (x,y,z)-translations:'
        write(iocf,*) '% xf,vb,zf=',xforce,vbulk,zforce
     end if
     if (iadd_damping.ne.0) then
        write(iocf,*) '% iadd_damping',damp_up,damp_aa
     end if
     if (iuse_LES.eq.1) then
        write(iocf,*) '% iuse_LES',idynamic,iadd_visc,Cles
     end if
     if (itemperature.eq.1) then
        write(iocf,*) '% itemperature',itemperature, fkappa, bym, gbeta
     end if
     if (inonNewtonian.eq.1) then
        write(iocf,*) '% inonNewtonian',inonNewtonian,mu_0, mu_inf, lambda,nnf
     end if
     close(iocf)  ! close file

     fname=filout(1:index(filout,' ')-1)//'.'//ext1//'.cfbin'
     open(iocfb,file=fname,status='unknown',form='unformatted')

     fnamet=filout(1:index(filout,' ')-1)//'.'//ext1//'.tfbin'
     open(iotfb,file=fnamet,status='unknown',form='unformatted')

  endif

  ! --- add new buffer to save change (take care for the confusing names)---
  allocate (buff1(buffsize),buff2(buffsize))
  buff1 = 0.d0; buff2 = 0.d0;
  allocate (u00s(0:my1),w00s(0:my1))
  u00s = 0.d0; w00s = 0.d0;
  if (iydealiasing.ge.1) then
     allocate (u00tmp(0:my1+2),w00tmp(0:my1+2))
     u00tmp = 0.d0; w00tmp = 0.d0;
     allocate (u00c(0:my1*2+1),w00c(0:my1*2+1))
     u00c = 0.d0; w00c = 0.d0;

     if (iydealiasing.eq.3) then

        shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
        call filterc(u00,u00,shp,shm,1,0)
        call filterc(w00,w00,shp,shm,1,0)
        ! apply a compact filter
        do k=kb,ke
           do i=0,mx1
              shp = shwkt(i)
              shm = shwkb(i)
              call filterc(vor(0,i,k),vor(0,i,k),shp,shm,2,1)
              call filterc(phi(0,i,k),phi(0,i,k),shp,shm,2,0) ! skip _addshear
           enddo
        enddo
     end if

  end if
  !
  ! ------------------- allocate everything ------------
  allocate (u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz) )
  allocate (o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz) )
  allocate (rf0u(0:my1),u00wk(0:my1), rf0w(0:my1),w00wk(0:my1))

  ! initialization
  u1r = 0.d0; u2r = 0.d0; u3r = 0.d0;
  o1r = 0.d0; o2r = 0.d0; o3r = 0.d0;

  !-- LES
  if (iuse_LES.eq.1) then

     allocate(txy00(0:my1),tyz00(0:my1))
     txy00=0.d0; tyz00=0.d0
     allocate(nuSGS(mgalx+2,mgalz),Sij2(mgalx+2,mgalz),Tll(mgalx+2,mgalz))
     nuSGS=0.d0; Sij2=0.d0; Tll=0.d0

     allocate(ux(mgalx+2,mgalz),uy(mgalx+2,mgalz),uz(mgalx+2,mgalz))
     allocate(vx(mgalx+2,mgalz),vy(mgalx+2,mgalz),vz(mgalx+2,mgalz))
     allocate(wx(mgalx+2,mgalz),wy(mgalx+2,mgalz),wz(mgalx+2,mgalz))
     ux=0.d0;uy=0.d0;uz=0.d0;
     vx=0.d0;vy=0.d0;vz=0.d0;
     wx=0.d0;wy=0.d0;wz=0.d0;
     allocate(uxc(0:mx1,0:mz1),uyc(0:mx1,0:mz1),uzc(0:mx1,0:mz1))
     allocate(vxc(0:mx1,0:mz1),vyc(0:mx1,0:mz1),vzc(0:mx1,0:mz1))
     allocate(wxc(0:mx1,0:mz1),wyc(0:mx1,0:mz1),wzc(0:mx1,0:mz1))
     uxc=dcmplx(0.d0,0.d0);uyc=dcmplx(0.d0,0.d0);uzc=dcmplx(0.d0,0.d0)
     vxc=dcmplx(0.d0,0.d0);vyc=dcmplx(0.d0,0.d0);vzc=dcmplx(0.d0,0.d0)
     wxc=dcmplx(0.d0,0.d0);wyc=dcmplx(0.d0,0.d0);wzc=dcmplx(0.d0,0.d0)

     !if (idynamic.eq.1) then
     !   allocate(u11r(mgalx+2,mgalz),u12r(mgalx+2,mgalz),u13r(mgalx+2,mgalz))
     !   allocate(u22r(mgalx+2,mgalz),u23r(mgalx+2,mgalz),u33r(mgalx+2,mgalz))
     !   allocate(u11c(0:mx1,0:mz1),u12c(0:mx1,0:mz1),u13c(0:mx1,0:mz1))
     !   allocate(u22c(0:mx1,0:mz1),u23c(0:mx1,0:mz1),u33c(0:mx1,0:mz1))
     !   u11r=0.d0;u12r=0.d0;u13r=0.d0;u22r=0.d0;u23r=0.d0;u33r=0.d0
     !   u11c=dcmplx(0.d0,0.d0);u12c=dcmplx(0.d0,0.d0);u13c=dcmplx(0.d0,0.d0)
     !   u22c=dcmplx(0.d0,0.d0);u23c=dcmplx(0.d0,0.d0);u33c=dcmplx(0.d0,0.d0)
     !end if

  end if

    if (inonNewtonian.eq.1) then

     allocate(txy00(0:my1),tyz00(0:my1))
     txy00=0.d0; tyz00=0.d0
     allocate(mu(mgalx+2,mgalz))
     mu=0.d0;

     allocate(ux(mgalx+2,mgalz),uy(mgalx+2,mgalz),uz(mgalx+2,mgalz))
     allocate(vx(mgalx+2,mgalz),vy(mgalx+2,mgalz),vz(mgalx+2,mgalz))
     allocate(wx(mgalx+2,mgalz),wy(mgalx+2,mgalz),wz(mgalx+2,mgalz))
     ux=0.d0;uy=0.d0;uz=0.d0;
     vx=0.d0;vy=0.d0;vz=0.d0;
     wx=0.d0;wy=0.d0;wz=0.d0;
     allocate(uxc(0:mx1,0:mz1),uyc(0:mx1,0:mz1),uzc(0:mx1,0:mz1))
     allocate(vxc(0:mx1,0:mz1),vyc(0:mx1,0:mz1),vzc(0:mx1,0:mz1))
     allocate(wxc(0:mx1,0:mz1),wyc(0:mx1,0:mz1),wzc(0:mx1,0:mz1))
     uxc=dcmplx(0.d0,0.d0);uyc=dcmplx(0.d0,0.d0);uzc=dcmplx(0.d0,0.d0)
     vxc=dcmplx(0.d0,0.d0);vyc=dcmplx(0.d0,0.d0);vzc=dcmplx(0.d0,0.d0)
     wxc=dcmplx(0.d0,0.d0);wyc=dcmplx(0.d0,0.d0);wzc=dcmplx(0.d0,0.d0)
     allocate(Txx(mgalx+2,mgalz),Tyy(mgalx+2,mgalz),Tzz(mgalx+2,mgalz))
     allocate(Txy(mgalx+2,mgalz),Tzx(mgalx+2,mgalz),Tyz(mgalx+2,mgalz))
     Txx=0.d0;Tyy=0.d0;Tzz=0.d0;
     Txy=0.d0;Tyz=0.d0;Tzx=0.d0;

  end if

  if (itemperature.eq.1) then
     allocate (buff3(buffsize))
     buff3 = 0.d0;
     allocate (tbr(mgalx+2,mgalz))
     tbr=0.d0;
     allocate (rf0tb(0:my1), tb00wk(0:my1))
     rf0tb=0.d0; tb00wk=0.d0;

     allocate(uTr(mgalx+2,mgalz),vTr(mgalx+2,mgalz),wTr(mgalx+2,mgalz))
     uTr=0.d0; vTr=0.d0; wTr=0.d0;
     allocate(uTc(0:mx1,0:mz1),vTc(0:mx1,0:mz1),wTc(0:mx1,0:mz1))
     uTc=dcmplx(0.d0,0.d0);vTc=dcmplx(0.d0,0.d0);wTc=dcmplx(0.d0,0.d0)
  end if

  ! ---- allocate tmp 2D arrays ----
  allocate (tmpxzr(mgalx+2,mgalz))
  tmpxzr = 0.d0
  allocate (tmpx(mgalx+2))
  tmpx = 0.d0
  !$OMP PARALLEL
  allocate(wk1dr(my), wk1dc(2*my))
  wk1dr = 0.d0
  wk1dc = 0.d0
  !$OMP END PARALLEL
  ! force to set conjugate relations for kx=0 modes
  ! since saving chikj2jik after hvhg keeps the initial perturbation forever
  call chjik2ikj(phi,phi,chwk,chwk)
  call chjik2ikj(vor,vor,chwk,chwk) ! ome2c
  call set_conj(vor,phi,0)  ! (rev.1158)
  call chikj2jik(phi,phi,chwk,chwk)
  call chikj2jik(vor,vor,chwk,chwk)

  if (itemperature.eq.1) then
     call chjik2ikj(tb,tb,chwk,chwk)
     call set_conj_t(tb,0)
     call chikj2jik(tb,tb,chwk,chwk)
  end if
  if (iadd_sym.gt.0) then
     if ((myid.eq.0).and.(nohist.ne.1)) write(*,*) ' iadd_sym=',iadd_sym
     call add_sym(vor,phi,u00,w00,buff1,buff2,u00s,w00s,chwk,iadd_sym)
     ! 1: allocation of vors, phis
  end if

  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !

  do istep=1,nstep

     wtime=MPI_WTIME()
     if ((myid.eq.0).and.(istop.eq.0)) then
        totaltimer = totaltimer-wtime
        iter_time  = -wtime
     endif

     if (mod(istep-1,nhist) .eq. 0) then
        if (nohist.eq.0) ihist=1
     end if

     if (iydealiasing.eq.2) then
        xoffb = (xoffb/Lx - int(xoffb/Lx) )*Lx ! negative
        xofft = -xoffb ! positive
     endif
     xwkt = xofft; xwkb = xoffb
     zwkt = zofft; zwkb = zoffb

     !if ((mod(istep-1,nimag).eq.0).or.(istop.eq.1).or.(idump.eq.1)) then
     if ((istep-1).ne.0) then
     if (((mod(istep-1,nimag).eq.0).and.(istep.ne.1)).or.(istop.eq.1).or.(idump.eq.1)) then
        ! ---  write image to a file. Note: dvordy,hv are work arrays
        !
        if (noescru.eq.0) then
           if (iwrite_hdf5.eq.1) then
              if (itemperature.eq.1) then
                 call write_t_hdf5(vor,phi,tb,u00,w00,tb00,dvordy,hv,ht,chwk)
              else
                 call writehdf5(vor,phi,u00,w00,dvordy,hv,chwk)
              end if
           else
              call escru(vor,phi,u00,w00,dvordy,hv,chwk)
           end if
        end if
        ! --- do something about writing the spectra here
        if (istop .eq. 1) then
           if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) 'exit cross at time = ', time
           exit
        end if
        idump = 0
     endif
     end if
     !   ---    prepara phiwk,vorwk,u00wk,w00wk */
     !   ---    all proc. has u00 etc.?

     u00b = s*(-Ly) ! yoffb, negative
     u00t = s*Ly    ! yofft, positive
     w00b = 0.d0
     w00t = 0.d0

     if (abs(s).gt.0.001d0) then ! s shoud be 0 or 1 ?

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ithread)
        if ((istep.eq.1).and.(myid.eq.0)) write(*,*) 'I am thread', OMP_GET_THREAD_NUM(), 'of MPI &
       &                                              task', myid
        !$OMP DO SCHEDULE(STATIC)
        ! do k=0,mz1 ! each proc has full shift infomation
           do i=0,mx1
              !write(*,*) 'computing Fourier-shift parameter',xwkt,xwkb
              !shb(i,k)  = cdexp(-xalp(i)*xwkb-xgam(k)*zwkb)
              !sht(i,k)  = cdexp(-xalp(i)*xwkt-xgam(k)*zwkt)
              shb(i)  = cdexp(-xalp(i)*xwkb) ! negative shift
              sht(i)  = cdexp(-xalp(i)*xwkt) ! positive shift
           enddo
        !enddo
        !shwkt=sht; shwkb=shb; ! memcopy
        !call vcopy((mx1+1)*2,sht,shwkt)
        !call vcopy((mx1+1)*2,shb,shwkb)
        !$OMP END PARALLEL

        call dcopy((mx1+1)*2,sht,1,shwkt,1)
        call dcopy((mx1+1)*2,shb,1,shwkb,1)

     end if

     if(irun.eq.0) then   !!! ---  this is done only for the first step
        irun = 1
        !if(myid.eq.0) write(*,*)  myid,'compute v at the first step'
        ! --- calcula la v a partir de phi ---
        do k=kb,ke
           !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,rk,shp,shm) SCHEDULE(STATIC)
           do i=0,mx1
              if ((i.ne.0).or.(k.ne.0)) then
                 rk = gam2(k)+alp2(i)
                 shp = sht(i)
                 shm = shb(i)
                 call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)   ! -- v
              end if
           enddo
           !$OMP END PARALLEL DO
        enddo
        ! --- do something about writing the spectra here ---

     endif ! --- end special first step ---

     do rkstep=1,3       ! --- Runge-Kutta third order  ---
        ! --- (0) pre-computes for nonlinear terms
        !         here hg is v, hv is overwritten by dv/dy

        t1= MPI_WTIME()

        do k=kb,ke
           !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm) SCHEDULE(STATIC)
           do i=0,mx1
              shp = shwkt(i)
              shm = shwkb(i)
              ! --- d(ome_2)/ dy
              call derivyc(vor(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
              ! --- dv/dy
              call derivyc(hg(0,i,k),hv(0,i,k),shp,shm,1,2,0) ! skip _addshear
           enddo
           !$OMP END PARALLEL DO
        enddo

        t2 = MPI_WTIME()
        t_derivyc = t2-t1

        ! (1) computes nonlinear terms
        ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
        ! all arrays are  sent from y-x --> x-z
        ! phi, vor is reference variables (should be preserved)
        ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! copy phi and vor to buff1, buff2 not to destroy it. (rev.416)

        if ((iydealiasing.eq.1).or.(iydealiasing.eq.2)) then

           fac = (xwkb/Lx - int(xwkb/Lx) )*Lx/Ly
           if (fac < -Lx/Ly*0.5d0) then
              ! positively distorted
              fac=Lx/Ly+fac;
              call moving_frame(vor,phi,fac,+1)
              call dealiasing_y(vor,phi,u00,w00,u00tmp,w00tmp,u00c,w00c,fac)
              call moving_frame(vor,phi,-fac,-1)
           else
              ! negatively distorted
              call moving_frame(vor,phi,fac,-1)
              call dealiasing_y(vor,phi,u00,w00,u00tmp,w00tmp,u00c,w00c,fac)
              call moving_frame(vor,phi,-fac,1)
           end if
           !if (myid.eq.0) write(*,*)  ' apply dealiasing in y'

        elseif (iydealiasing.eq.3) then

           shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
           call filterc(u00,u00,shp,shm,1,0)
           call filterc(w00,w00,shp,shm,1,0)
           ! apply a compact filter
           do k=kb,ke
              do i=0,mx1
                 shp = shwkt(i)
                 shm = shwkb(i)
                 call filterc(vor(0,i,k),vor(0,i,k),shp,shm,2,1)
                 call filterc(phi(0,i,k),phi(0,i,k),shp,shm,2,0) ! skip _addshear
              enddo
           enddo

        end if

        buff1=0.d0; buff2=0.d0
        call dcopy(buffsize,phi,1,buff1,1)
        call dcopy(buffsize,vor,1,buff2,1)
        !
        call chjik2ikj(buff1,buff1,chwk,chwk)
        call chjik2ikj(buff2,buff2,chwk,chwk) ! ome2c
        !
        !call chjik2ikj(phi,phi,chwk,chwk)
        !call chjik2ikj(vor,vor,chwk,chwk) ! ome2c
        call set_conj(buff1,buff2,0)  ! rev(1158)
        call chjik2ikj(hv,hv,chwk,chwk)   ! dv/dy
        call chjik2ikj(hg,hg,chwk,chwk)   ! v
        call chjik2ikj(dvordy,dvordy,chwk,chwk) ! d(ome_2)/ dy

        if (itemperature.eq.1) then
           buff3=0.d0; ht=0.d0
           call dcopy(buffsize,tb,1,buff3,1)
           call chjik2ikj(buff3,buff3,chwk,chwk)
        end if
        if (iuse_LES.eq.1) then
           ! LES buffers
           rhvc1= 0.d0; rhvc2=0.d0
           rhgc1= 0.d0; rhgc2=0.d0
        end if

        if (inonNewtonian.eq.1) then
           ! LES buffers
           rhvc1= 0.d0; rhvc2=0.d0
           rhgc1= 0.d0; rhgc2=0.d0
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        t_hvhg0 = MPI_WTIME()

!        !!(rev.416) buff1:phi, buff2:vor
!        call hvhg(buff1,buff2,u00,w00,hv,hg,rf0u,rf0w,dvordy,rkstep, &
!             u1r,u2r,u3r,o1r,o2r,o3r,u1r,u2r,u3r,o1r,o2r,o3r, &
!             tmpxzr,tmpx,wk1dr)
        call hvhgt(buff1,buff2,u00,w00,hv,hg, &
             rf0u,rf0w,dvordy,rkstep, &
             u1r,u2r,u3r,o1r,o2r,o3r, &
             u1r,u2r,u3r,o1r,o2r,o3r, &
             tmpxzr,tmpx,wk1dr,rhvc1,rhvc2,rhgc1,rhgc2, &
             buff3,tb00,ht,rf0tb,tbr,tbr,rhtc2)


        t_hvhg = t_hvhg + MPI_WTIME() - t_hvhg0

        ! c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
        !  at this point: dvordy, rhv and rhg are the outs
        !  they must be trasformed from xy-xz before completion
        !  dvordy: dH1/dx + dH3/dz
        !  hg: -dH3/dx + dH1/dz
        !  hv: d^2H2/dx^2 + d^2H2/dz^2
        !  phi, vor should be recovered by shuffling, y-x <-- x-z
        ! c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c


        if (iuse_LES.eq.1) then
           call chikj2jik(rhvc1,rhvc1,chwk,chwk)
           call chikj2jik(rhvc2,rhvc2,chwk,chwk)
           call chikj2jik(rhgc1,rhgc1,chwk,chwk)
           call chikj2jik(rhgc2,rhgc2,chwk,chwk)
        end if
        if (inonNewtonian.eq.1) then
           call chikj2jik(rhvc1,rhvc1,chwk,chwk)
           call chikj2jik(rhvc2,rhvc2,chwk,chwk)
           call chikj2jik(rhgc1,rhgc1,chwk,chwk)
           call chikj2jik(rhgc2,rhgc2,chwk,chwk)
        end if
        if (itemperature.eq.1) then
           call chikj2jik(rhtc2,rhtc2,chwk,chwk)
           call chikj2jik(ht,ht,chwk,chwk)
        end if
        ! ! skip (rev.416)
        ! !call chikj2jik(phi,phi,chwk,chwk)
        ! !call chikj2jik(vor,vor,chwk,chwk)
        call chikj2jik(hv,hv,chwk,chwk)
        call chikj2jik(hg,hg,chwk,chwk)
        call chikj2jik(dvordy,dvordy,chwk,chwk)
        ! Note: now buff1(phi) and buff2(vor) is not overwritten ...
        !
        ! --- computes rhv = d^2H2/dx^2 + d^2H2/dz^2 - d(dH1/dx + dH3/dz)/dy

        t1= MPI_WTIME()


        !
        !
        !
        !

        if (iuse_LES.eq.1) then
           ! d( vor - rhvc1)/dy
           dvordy = dvordy - rhvc1
        end if
        do k=kb,ke
           !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm) SCHEDULE(STATIC)
           do i=0,mx1
              shp = shwkt(i)
              shm = shwkb(i)
              call derivyc(dvordy(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
              if (iuse_LES.eq.1) call derivyc(rhgc2(0,i,k),rhgc2(0,i,k),shp,shm,1,2,0)
              if (inonNewtonian.eq.1) call derivyc(rhgc2(0,i,k),rhgc2(0,i,k),shp,shm,1,2,0)
              if (itemperature.eq.1) then
                 ! -d(vT)/dy
                 call derivyc(rhtc2(0,i,k),rhtc2(0,i,k),shp,shm,1,2,0)
              end if
              !
              !hv(:,i,k) = hv(:,i,k) - dvordy(:,i,k)
           enddo
           !$OMP END PARALLEL DO
        enddo
        if (itemperature.eq.1) then
           shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
           call derivyc(rf0tb,rf0tb,shp,shm,1,1,0)
        end if
        if (iuse_LES.eq.1) then
           ! (d^2/dx^2+d^2/dz^2-d^2/dy^2)*(rhvc2)
           fac = -1.d0
           do k=kb,ke
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    shp = shwkt(i)
                    shm = shwkb(i)
                    rk  = -(gam2(k)+alp2(i))
                    !call derivyc(rhvc2(0,i,k) ,buff(0,i,k)  ,shp,shm,2,2,1)
                    ! BUG fixed: addshear should be called for 2nd-dev
                    !rhvc2(:,i,k)=-buff(:,i,k)-(gam2(k)+alp2(i))*rhvc2(:,i,k)
                    ! ! replaced with visc3d() with factor=-1.d0*rk
                    call visc3d(rhvc2(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',1)
                 end if
              enddo
              !$OMP END PARALLEL DO
           enddo
        end if


        t2 = MPI_WTIME()
        t_derivyc = t2-t1


        ! hv=hv+d(rhvc1)/dy+(d^2/dx^2+d^2/dz^2-d^2/dy^2)*(rhvc2)
        ! hg=hg+rhgc1+d(rhgvc2)/dy

        hv = hv - dvordy
        if (iuse_LES.eq.1) hg = hg + rhgc1  + rhgc2
        if (inonNewtonian.eq.1) hg = hg + rhgc1  + rhgc2
        if (itemperature.eq.1) ht = ht + rhtc2

        ! (1-ii) here you can add body force: ++hv, ++hg
        if (iadd_force.gt.0) then
           ! streamwise-roll force F=(Fy,Fz), div.F = 0
           !  then hv += (d^2/dz^2 + d^2/dy^2)*Fy
           !       hg += 0
           if((mod(istep-1,nhist).eq.0).and.(rkstep.eq.1)) then
              if ((ihist.eq.1).and.(myid.eq.0)) then
                 write(*,*) ' add_force',iadd_force
              end if
           end if
           call add_force(hv,iadd_force)
        end if
        if (iadd_force.eq.-1) then
           rf0u = rf0u + xforce ! inertial force relative to the moving frame
        endif
        if (iadd_damping.eq.1) then
           call add_damping(phi,vor,u00,w00,hv,hg,rf0u,rf0w,iadd_damping)
        end if
        !if (itemperature.eq.1) then
        !   ! hv += gbeta*lap.T ! add this later using visc3d for explicit
        !end if

        !---advance in time :
        r1 = ca(rkstep)*Deltat
        dtr1v = dtrv*ibeta(rkstep) ! implicit only
        dtr1k = dtrk*ibeta(rkstep) ! implicit only
        !
        if (explicit.eq.1) then
           ! --- compute  nonlinear terms and viscous terms explicitly ---
           ! (2-explicit) add  nu*lap.(phi) -> ++[hv,hg]

           t1= MPI_WTIME()

           if (iuse_LES.eq.1) then

              shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
              call derivyc(txy00,txy00,shp,shm,1,1,0)
              call derivyc(tyz00,tyz00,shp,shm,1,1,0)

              rf0u= rf0u+txy00
              rf0w= rf0w+tyz00

           endif

           if (inonNewtonian.eq.1) then
              shp = dcmplx(1.d0, 0.d0);shm = dcmplx(1.d0, 0.d0)
              call derivyc(txy00, txy00, shp, shm, 1, 1, 0)
              call derivyc(tyz00, tyz00, shp, shm, 1, 1, 0)
              rf0u = rf0u + txy00
              rf0w = rf0w + tyz00
           endif

           if (( (iadd_visc.eq.1).and.(iuse_LES.eq.1) ).or.((iuse_LES.eq.0).and.(inonNewtonian.eq.0))) then
              fac = 1.d0/Re
              do k=kb,ke          ! --- computes  lap.(vor)
                 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk) SCHEDULE(STATIC)
                 do i=0,mx1
                    if ((i.ne.0).or.(k.ne.0)) then
                       shp = shwkt(i) ! do not use the sht at t_new=t+c*dt
                       shm = shwkb(i)
                       rk = gam2(k)+alp2(i)
                       call visc3d(vor(0,i,k),hg(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',1)
                       call visc3d(phi(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',0)
                       if (itemperature.eq.1) then
                          call visc3d(tb(0,i,k),ht(0,i,k),wk1dc,shp,shm,2,rk,fkappa,'+',0)
                          ! adding buoyancy force
                          call visc3d(tb(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,gbeta,'+',0)
                       end if
                    end if
                 enddo
                 !$OMP END PARALLEL DO
              enddo
              !
              ! (*) zerozero-mode (all slaves computes 00modes (from rev.379))
              !    (2*-explicit) compute laplacian --> ++[rf0u,rf0w]
              shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
              rk=0.d0
              fac = (1.d0/Re) ! nu
              call visc3d(u00,rf0u,wk1dr,shp,shm,1,rk,fac,'+',0)
              call visc3d(w00,rf0w,wk1dr,shp,shm,1,rk,fac,'+',0)
              if (itemperature.eq.1) then
                 fac_t = fkappa
                 call visc3d(tb00,rf0tb,wk1dr,shp,shm,1,rk,fac_t,'+',0)
              end if
           end if


           t2 = MPI_WTIME()
           t_derivyc = t2-t1

           if ((iadd_mode.ge.3).and.(iadd_mode.le.4)) then
              !! multiply the right hand side by exp(i s*y(j)*alp*t)
              fac=xwkt/Ly
              call moving_frame(hg,hv,fac,1)
              if (itemperature.eq.1) then
                 write(*,*) 'use iadd_mode=5, stop'
                 stop
              end if
           endif

           ! (3-explicit) add gama(krk)*[hv,hg] --> ++[phiwk,vorwk]
           !             here hv is  R= hv + nu*lap.(vor)
           !             then, phiwk = gama_k*R_k + zeta_{k-1}*R_{k-1}
           !
           ! y*    = y(n) + gama1*h*f(y(n))
           ! y**   = y*   + gama2*h*f(y*)  + omega2*h*f(y(n))
           ! y(n+1)= y**  + gama3*h*f(y**) + omega3*h*f(y*)
           !
           if (rkstep.eq.1) then
              vorwk = gama(rkstep)*hg
              phiwk = gama(rkstep)*hv
              if (itemperature.eq.1) then
                 tbwk = gama(rkstep)*ht
              end if
           else
              vorwk = vorwk + gama(rkstep)*hg
              phiwk = phiwk + gama(rkstep)*hv
              if (itemperature.eq.1) then
                 tbwk = tbwk + gama(rkstep)*ht
              end if
           end if
           !
           !    (*) zerozero-mode
           !    (3*-explicit) add gama_k*[rf0u,rf0w] --> ++[u00wk,w00wk]
           if (rkstep.eq.1) then
              u00wk = gama(rkstep)*rf0u
              w00wk = gama(rkstep)*rf0w
              if (itemperature.eq.1) then
                 tb00wk = gama(rkstep)*rf0tb
              end if
           else
              u00wk = u00wk + gama(rkstep)*rf0u
              w00wk = w00wk + gama(rkstep)*rf0w
              if (itemperature.eq.1) then
                 tb00wk = tb00wk + gama(rkstep)*rf0tb
              end if
           end if
           !
           if ((iadd_mode.ge.3).and.(iadd_mode.le.4)) then
              fac = xwkt/Ly
              call moving_frame(vor,phi,fac,1)
           endif
           ! (4-explicit) update variables
           !
           vor = vor + Deltat*vorwk
           phi = phi + Deltat*phiwk
           if (itemperature.eq.1) then
              tb = tb + Deltat*tbwk
           end if
           !
           !    (4*-explicit) update zerozero-mode
           u00 = u00 + Deltat*u00wk
           w00 = w00 + Deltat*w00wk
           if (itemperature.eq.1) then
              tb00 = tb00 + Deltat*tb00wk
           end if
           !
           ! (5-explicit) update shear-periodic boundary factor
           if (abs(s).gt.0.01d0) then

              xwkb = xoffb + r1*u00b
              xwkt = xofft + r1*u00t
              zwkb = zoffb + r1*w00b
              zwkt = zofft + r1*w00t

              if (iydealiasing.eq.2) then
                 xwkb = (xwkb/Lx - int(xwkb/Lx) )*Lx
                 xwkt = -xwkb
              endif

              !do k=0,mz1
                 do i=0,mx1
                    !shb(i,k)  = cdexp(-xalp(i)*xwkb-xgam(k)*zwkb)
                    !sht(i,k)  = cdexp(-xalp(i)*xwkt-xgam(k)*zwkt)
                    shwkb(i)  = zexp(-xalp(i)*xwkb) ! negative shift
                    shwkt(i)  = zexp(-xalp(i)*xwkt) ! positive shift
                 enddo
              !enddo
           endif

           !
           ! (6-explicit) store [Rv,Rg]*zeta_k -> [phiwk,vorwk]
           !               here Rv and Rg are (nonlinear) + (linear)
           !         iadd_mode=5: mapped by exp( (c(i)-c(i-1))*dt*k_x*S*y ), (i=1, 2), c(0) = 0
           if (rkstep.ne.3) then
              !$OMP PARALLEL WORKSHARE
              phiwk = hv*zeta(rkstep)
              vorwk = hg*zeta(rkstep)
              !if (itemperature.eq.1) then
                 tbwk = ht*zeta(rkstep)
              !end if
              !$OMP END PARALLEL WORKSHARE
              if ((iadd_mode.eq.5).or.(iadd_mode.eq.-1)) then
                 if (rkstep.eq.1) then
                    fac = -gama(rkstep)*Deltat*s
                 elseif (rkstep.eq.2) then
                    fac = -(gama(rkstep)+zeta(rkstep-1))*Deltat*s
                 endif
                 call moving_frame(vorwk,phiwk,fac,-1)
                 if (itemperature.eq.1) then
                    call moving_frame_t(tbwk,fac,-1)
                 end if
              endif
              ! zerozero-mode
              u00wk = rf0u*zeta(rkstep)
              w00wk = rf0w*zeta(rkstep)
              if (itemperature.eq.1) then
                 tb00wk = rf0tb*zeta(rkstep)
              end if
           end if
           !
           ! un-mapping
           if ((iadd_mode.ge.3).and.(iadd_mode.le.4)) then
              fac = xwkb/Ly
              call moving_frame(vor,phi,fac,-1)
              if (itemperature.eq.1) then
                 call moving_frame_t(tb,fac,-1)
              end if
           elseif ((iadd_mode.eq.5).or.(iadd_mode.eq.-1)) then
              if (rkstep.eq.1) then
                 fac = -s*r1
              else
                 fac = -s*(ca(rkstep)-ca(rkstep-1))*Deltat
              endif
              call moving_frame(vor,phi,fac,-1)
              if (itemperature.eq.1) then
                 call moving_frame_t(tb,fac,-1)
              end if
           endif
           ! (7-explicit) solve lap.v = phi
           ! -- v
           do k=kb,ke
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,rk,shp,shm) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    rk = gam2(k)+alp2(i)
                    shp = shwkt(i)
                    shm = shwkb(i)
                    call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)
                 end if
              enddo
              !$OMP END PARALLEL DO
           enddo

        else ! semi-implicit for viscous terms
           ! --- compute nonlinear terms explicitly and viscous implicitly
           !
           !
           if (itemperature.eq.1) then
              do k=kb,ke          ! --- computes  lap.(vor),  lap.(phi)
                 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk) SCHEDULE(STATIC)
                 do i=0,mx1
                    if ((i.ne.0).or.(k.ne.0)) then
                       shp = shwkt(i) ! use the previous one before updated by (7)
                       shm = shwkb(i)
                       rk = gam2(k)+alp2(i)
                       ! adding buoyancy force
                       call visc3d(tb(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,gbeta,'+',1)
                    end if
                 end do
                 !$OMP END PARALLEL DO
              enddo
           endif
           ! (2) add gama(krk)*[hv,hg] --> ++[phiwk,vorwk]
           if (rkstep.eq.1) then
              !$OMP PARALLEL WORKSHARE
              vorwk = gama(rkstep)*hg
              phiwk = gama(rkstep)*hv
              !if (itemperature.eq.1) then
                 tbwk = gama(rkstep)*ht
              !end if
              !$OMP END PARALLEL WORKSHARE
           else
              !$OMP PARALLEL WORKSHARE
              vorwk = vorwk + gama(rkstep)*hg
              phiwk = phiwk + gama(rkstep)*hv
              !if (itemperature.eq.1) then
                 tbwk = tbwk + gama(rkstep)*ht
              !end if
              !$OMP END PARALLEL WORKSHARE
           end if
           ! (3) compute laplacian
           !          alpha(krk)*[phi,vor] --> ++[phiwk,vorwk]
           fac = alpha(rkstep)*(1.d0/Re)
           if (itemperature.eq.1) fac_t = alpha(rkstep)*fkappa
           do k=kb,ke          ! --- computes  lap.(vor),  lap.(phi)
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    shp = shwkt(i) ! use the previous one before updated by (7)
                    shm = shwkb(i)
                    rk = gam2(k)+alp2(i)
                    call visc3d(vor(0,i,k),vorwk(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',1)
                    call visc3d(phi(0,i,k),phiwk(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',0)
                    if (itemperature.eq.1) then
                       call visc3d(tb(0,i,k),tbwk(0,i,k),wk1dc,shp,shm,2,rk,fac_t,'+',0)
                    end if
                 endif
              enddo
              !$OMP END PARALLEL DO
           enddo
           ! (*) zero-zero mode
           if (myid.eq.0) then
              ! (2*) add gama(krk)*[rf0u,rf0w] --> ++[u00wk,w00wk]
              if (rkstep.eq.1) then
                 u00wk = gama(rkstep)*rf0u
                 w00wk = gama(rkstep)*rf0w
                 if (itemperature.eq.1) then
                    tb00wk = gama(rkstep)*rf0tb
                 end if
              else
                 u00wk = u00wk + gama(rkstep)*rf0u
                 w00wk = w00wk + gama(rkstep)*rf0w
                 if (itemperature.eq.1) then
                    tb00wk = tb00wk + gama(rkstep)*rf0tb
                 end if
              end if
              ! (3*) compute laplacian --> ++[u00wk,w00wk]
              shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
              rk=0.d0
              fac = alpha(rkstep)*(1.d0/Re)
              call visc3d(u00,u00wk,wk1dr,shp,shm,1,rk,fac,'+',0)
              call visc3d(w00,w00wk,wk1dr,shp,shm,1,rk,fac,'+',0)
              if (itemperature.eq.1) then
                 fac_t = alpha(rkstep)*fkappa
                 call visc3d(tb00,tb00wk,wk1dr,shp,shm,1,rk,fac_t,'+',0)
              end if
           end if

           ! (4) add (1.d0/Dt)*[phi,vor] --> ++[phiwk,vorwk]
           !$OMP PARALLEL WORKSHARE
           vorwk = vorwk + (1.d0/Deltat)*vor
           phiwk = phiwk + (1.d0/Deltat)*phi
           !if (itemperature.eq.1) then
              tbwk = tbwk + (1.d0/Deltat)*tb
           !end if
           !$OMP END PARALLEL WORKSHARE
           ! (5) driving pressure: skip
           ! [only for primitive valiables (fractional step method)]

           ! (6) multiply,  -1/(beta(krk)*nu)*[phiwk,vorwk]
           fac = -Re*ibeta(rkstep)
           if (itemperature.eq.1) fac_t = -(1.d0/fkappa)*ibeta(rkstep) !2018/11/07 fixed
           !$OMP PARALLEL WORKSHARE
           vorwk = fac*vorwk
           phiwk = fac*phiwk
           !if (itemperature.eq.1) then
              tbwk = fac_t*tbwk
           !end if
           !$OMP END PARALLEL WORKSHARE
           ! (* 4* - 6*) zero-zero mode
           if (myid.eq.0) then
              !write(*,*) 'Deltat, dtr1=', Deltat,dtr1
              u00wk = fac*(u00wk + (1.d0/Deltat)*u00)
              w00wk = fac*(w00wk + (1.d0/Deltat)*w00)
              if (itemperature.eq.1) then
                 tb00wk = fac_t*(tb00wk + (1.d0/Deltat)*tb00)
              end if
           end if

           ! (7) ------ update shear-periodic boundary -----
           if (abs(s).gt.0.5d0) then ! s shoud be 0 or 1

              xwkb = xoffb + r1*u00b
              xwkt = xofft + r1*u00t
              zwkb = zoffb + r1*w00b
              zwkt = zofft + r1*w00t

              !do k=0,mz1
                 do i=0,mx1
                    !shb(i,k)  = cdexp(-xalp(i)*xwkb-xgam(k)*zwkb)
                    !sht(i,k)  = cdexp(-xalp(i)*xwkt-xgam(k)*zwkt)
                    shwkb(i)  = zexp(-xalp(i)*xwkb) ! negative shift
                    shwkt(i)  = zexp(-xalp(i)*xwkt) ! positive shift
                 enddo
              !enddo
           endif
           !
           ! (8,9) solve helmholtz equations: phi'' - rk*phi = phiwk
           !   and store nonlinear terms: zeta(krk)*[hv, hg] -> [phiwk, vorwk]
           ! mapping r.h.s (iadd_mode=5)
           if ((iadd_mode.eq.5).or.(iadd_mode.eq.-1)) then
              if (rkstep.eq.1) then
                 fac = -s*r1
              else
                 fac = -s*(ca(rkstep)-ca(rkstep-1))*Deltat
              endif
              call moving_frame(vorwk,phiwk,fac,-1)
              if (itemperature.eq.1) then
                 call moving_frame_t(tbwk,fac,-1)
              end if
           endif

           do k=kb,ke
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk,rk2) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    rk = gam2(k)+alp2(i)
                    rk2 = dtr1v+rk
                    !
                    shp = shwkt(i)
                    shm = shwkb(i)
                    ! -- (8-i) vor viscous
                    call lapcdy(vorwk(0,i,k),rk2,vor(0,i,k),shp,shm,2)
                 endif
              enddo
              !$OMP END PARALLEL DO
           enddo
           do k=kb,ke
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk,rk2) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    rk = gam2(k)+alp2(i)
                    rk2 = dtr1v+rk
                    !
                    shp = shwkt(i)
                    shm = shwkb(i)
                    ! -- (8-ii) lap.(phi)
                    call lapcdy(phiwk(0,i,k),rk2,phi(0,i,k),shp,shm,2)
                 endif
              enddo
              !$OMP END PARALLEL DO
           enddo
           if (itemperature.eq.1) then
              do k=kb,ke
                 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk,rk2) SCHEDULE(STATIC)
                 do i=0,mx1
                    if ((i.ne.0).or.(k.ne.0)) then
                       rk = gam2(k)+alp2(i)
                       rk2 = dtr1k+rk
                       !
                       shp = shwkt(i)
                       shm = shwkb(i)
                       ! -- (8-i) temperature diffusion
                       call lapcdy(tbwk(0,i,k),rk2,tb(0,i,k),shp,shm,2)
                    endif
                 enddo
                 !$OMP END PARALLEL DO
              enddo
           end if

           ! -- (9-i)-- store zeta_k*hg --> vorwk
           ! -- (9-ii)-- store zeta(krk)*hv --> phiwk
           if (rkstep.ne.3) then
              !OMP PARALLEL WORKSHARE
              vorwk = zeta(rkstep)*hg
              phiwk = zeta(rkstep)*hv
              !if (itemperature.eq.1) then
                 tbwk = zeta(rkstep)*ht
              !end if
              !OMP END PARALLEL WORKSHARE
              if ((iadd_mode.eq.5).or.(iadd_mode.eq.-1)) then
                 if (rkstep.eq.1) then
                    fac = -gama(rkstep)*Deltat*s
                 elseif (rkstep.eq.2) then
                    fac = -(gama(rkstep)+zeta(rkstep-1))*Deltat*s
                 endif
                 call moving_frame(vorwk,phiwk,fac,-1)
                 if (itemperature.eq.1) then
                    call moving_frame_t(tbwk,fac,-1)
                 end if
              endif
           endif
           do k=kb,ke
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,shp,shm,rk,rk2) SCHEDULE(STATIC)
              do i=0,mx1
                 if ((i.ne.0).or.(k.ne.0)) then
                    rk = gam2(k)+alp2(i)
                    rk2 = dtr1v+rk
                    !
                    shp = shwkt(i)
                    shm = shwkb(i)
                    ! -- v
                    call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)
                 endif
              enddo
              !$OMP END PARALLEL DO
           enddo
           !    Now hv is free, so that we can store dv/dy to hv
           !    at the beginning of rkstep
           ! (*) zero-zero mode
           if (myid.eq.0) then
              ! (8*) solve helmholtz equations and store nonlinear terms:
              shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
              call lapcdy(u00wk,dtr1v,u00,shp,shm,1) ! real*8
              shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
              call lapcdy(w00wk,dtr1v,w00,shp,shm,1) ! real*8
              if (itemperature.eq.1) then
                 shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
                 call lapcdy(tb00wk,dtr1k,tb00,shp,shm,1) ! real*8
              end if
              !master proc. SENDING
              do iproc=1,numerop-1
                 call MPI_SEND(u00,my,MPI_REAL8,iproc,iproc, &
                      MPI_COMM_WORLD,ierr)
                 call MPI_SEND(w00,my,MPI_REAL8,iproc,iproc, &
                      MPI_COMM_WORLD,ierr)
                 if (itemperature.eq.1) then
                    call MPI_SEND(tb00,my,MPI_REAL8,iproc,iproc, &
                         MPI_COMM_WORLD,ierr)
                 end if
              enddo
              ! (9*) store nonlinear term: zeta(krk)*rf0u --> u00wk,w00wk
              if (rkstep.ne.3) then
                 u00wk = zeta(rkstep)*rf0u
                 w00wk = zeta(rkstep)*rf0w
                 if (itemperature.eq.1) then
                    tb00wk = zeta(rkstep)*rf0tb
                 end if
              end if

           else       ! -------- this is done by all slaves
              !
              !write(*,*) 'receive 00 mode',myid,my
              call MPI_RECV(u00,my,MPI_REAL8,0,MPI_ANY_TAG, &
                   MPI_COMM_WORLD,istat,ierr)
              call MPI_RECV(w00,my,MPI_REAL8,0,MPI_ANY_TAG, &
                   MPI_COMM_WORLD,istat,ierr)
              if (itemperature.eq.1) then
                 call MPI_RECV(tb00,my,MPI_REAL8,0,MPI_ANY_TAG, &
                      MPI_COMM_WORLD,istat,ierr)
              end if
           endif

        endif ! end of explicit = 0

        if (iadd_sym.gt.0) then
           !if ((myid.eq.0).and.(nohist.ne.1)) write(*,*) ' iadd_sym=',iadd_sym
           call add_sym(vor,phi,u00,w00,buff1,buff2,u00s,w00s,chwk,iadd_sym)
        end if

     enddo  !!! finalizado subpaso del RK3

     xoffb=xwkb
     xofft=xwkt
     zoffb=zwkb
     zofft=zwkt

     !sht = shwkt; shb = shwkb; ! memcopy
     !call vcopy((mx1+1)*2,shwkt,sht)
     !call vcopy((mx1+1)*2,shwkb,shb)
     call dcopy((mx1+1)*2,shwkt,1,sht,1)
     call dcopy((mx1+1)*2,shwkb,1,shb,1)
     time=time+Deltat      ! ---  update  time:

     wtime = MPI_WTIME()
     if (myid.eq.0) then

        totaltimer=totaltimer+wtime
        if ((mod(istep-1,nhist).eq.0).and.(nohist.eq.0).and.(nocf.eq.0)) then
           ! --- Time elapsed/step
           if (iskip_screenout.eq.0) then
              write(*,'(a,3(1X,E14.6))') 'Time/step:', &
                   wtime+iter_time-commtimer+comm_time,  &
                   commtimer-comm_time,wtime+iter_time
           end if
        end if
        comm_time = commtimer
     end if

  enddo   !!!!!!!  this is the istep loop

  if ((myid.eq.0).and.(nohist.eq.0)) then
     write(*,*) "Total time: ",totaltimer
     write(*,*) "Trans. time: ",transtimer
     write(*,*) "Comm. time: ",commtimer
     if(iadd_sym.gt.0) write(*,*) "Addsym. time:",addsymtimer

     write(*,*) "Trans/Total: ",transtimer/totaltimer
     write(*,*) "Comm/Total: ",commtimer/totaltimer

     write(*,*) "Aver time per step: ",totaltimer/nstep

     write(*,*) "time cost if derivyc:",t_derivyc
     write(*,*) "total of hvhg:",t_hvhg
     write(*,*) "      in hvhg (addshear, uvw, copy, fourxz):", &
          t_add_shear, t_uvw, t_copy, t_fourxz
     write(*,*) "              (uxomg, nl, otra, wait rf0uw):", &
          t_uome, t_nl, t_otra, t_waiting

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! deallocate
  deallocate(wk1dr, wk1dc)
  deallocate(tmpx)
  deallocate(tmpxzr)

  if (itemperature.eq.1) then
     deallocate(uTc,vTc,wTc)
     deallocate(uTr,vTr,wTr)
  end if

  if (iuse_LES.eq.1) then

     !if (idynamic.eq.1) then
     !
     !   deallocate(u22c,u23c,u33c)
     !   deallocate(u11c,u12c,u13c)
     !   deallocate(u22r,u23r,u33r)
     !   deallocate(u11r,u12r,u13r)
     !
     !end if

     deallocate(wxc,wyc,wzc)
     deallocate(vxc,vyc,vzc)
     deallocate(uxc,uyc,uzc)

     deallocate(wx,wy,wz)
     deallocate(vx,vy,vz)
     deallocate(ux,uy,uz)

     deallocate(nuSGS,Sij2,Tll)

     deallocate(txy00,tyz00)

  end if

  if (inonNewtonian.eq.1) then

     deallocate(txy00,tyz00)
     deallocate(wxc,wyc,wzc)
     deallocate(vxc,vyc,vzc)
     deallocate(uxc,uyc,uzc)
     deallocate(wx,wy,wz)
     deallocate(vx,vy,vz)
     deallocate(ux,uy,uz)
     deallocate(mu)
     deallocate(Txx,Tyy,Tzz)
     deallocate(Txy,Tzx,Tyz)

  end if

  deallocate(rf0u,u00wk,rf0w,w00wk)
  deallocate(o1r,o2r,o3r)
  deallocate(u1r,u2r,u3r)

  if (iydealiasing.ge.1) then
     deallocate (u00c,w00c)
     deallocate (u00tmp,w00tmp)
  endif
  deallocate(u00s,w00s)
  deallocate(buff1,buff2)
  if (itemperature.eq.1) then
     deallocate(tbr)
     deallocate(rf0tb,tb00wk)
     deallocate(buff3)
  end if

  if ((myid.eq.0).and.(nohist.eq.0)) then
     !close(iocf)  ! close file
     close(iocfb)  ! close binary file
     if (itemperature.eq.1) close(iotfb)
  end if
  ! reset ihist = 0
  ihist=0

end subroutine crosst

!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
subroutine hvhgt(phic,ome2c,u00,w00,rhvc,rhgc, &
     rf0u,rf0w,ome1c,rkstep, &
     u1r,u2r,u3r,o1r,o2r,o3r,  &
     u1c,u2c,u3c,o1c,o2c,o3c,  & ! memory is shared with u1r, etc ...
     tmpxzr,tmpx,tmpy,rhvc1,rhvc2,rhgc1,rhgc2, & ! txy00,tyz00,
     tempc,tb00,rhtc,rf0tb,tbr,tbc,rhtc2 )

  use ctes
  use statistics
  use running
  use bcs
  use timer
  use LES
  use temp
  use omp_lib

  implicit none
  include "mpif.h"

  complex(8),dimension(0:mx1,0:mz1,jb:je) :: phic, ome1c, ome2c, rhgc, rhvc
  complex(8),dimension(0:mx1,0:mz1,jb:je) :: tempc, rhtc ! itemperature
  real(8),dimension(0:my1) :: rf0u, rf0w, u00, w00
  real(8),dimension(0:my1) :: rf0tb, tb00 ! itemperature

  !------ 6 * (mgalx+2)  * mgalz buffer planes, THESE PLANES SHARE STORAGE !!
  real(8),dimension(mgalx+2,mgalz) :: u1r, u2r, u3r, o1r, o2r, o3r
  complex(8),dimension(0:mx1,0:mz1) :: u1c, u2c, u3c, o1c, o2c, o3c
  real(8),dimension(mgalx+2,mgalz) :: tbr ! itemperature
  complex(8),dimension(0:mx1,0:mz1) :: tbc ! itemperature
  !--------------------------------------------------------------------------
  real(8),dimension(mgalx+2,mgalz) :: tmpxzr
  real(8),dimension(mgalx+2) :: tmpx
  real(8),dimension(0:my1) :: tmpy

  ! --- LES buffer
  complex(8),dimension(0:mx1,0:mz1,jb:je) :: rhvc1, rhvc2, rhgc1, rhgc2
  !real(8) txy00(0:my1),tyz00(0:my1) ! moved in mod.f90
  ! --- itemperature
  complex(8),dimension(0:mx1,0:mz1,jb:je) ::  rhtc2

  integer rkstep,i,j,k,mmyr,iproc,ipo,leng,pproc
  integer istat(MPI_STATUS_SIZE),ierr
  real(8)  dt,fac
  complex(8) shp,shm, cfac
  complex(8) kkw(0:mx1)
  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
  !  at this point:
  !  rhv is dv / dy -- F-F-P
  !  rhg is v -- F-F-P
  !  phi is nabla^2(v) -- F-F-P
  !  ome1 is d(omega_2)/dy -- F-F-P
  !  ome2 is omega_2 --F-F-P
  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

  t1=MPI_WTIME()

  ! iadd_mode=1  ! using u = u' + Sy (see eq.(12) and (16) in Doc/note.pdf)
  ! iadd_mode=2  ! using u = u' + Sy, oz = oz' - S (see eq.(13) and (17))
  !
  ! iadd_mode=1 and 2 makes no difference in statistics
  ! ! both options are checked with the RDT solution (rev.430)
  ! ! I have forgotten to revert it to iadd_mode=2  (rev.440)
  ! ! take care for the statistics, o3m o3p includes mean shear
  ! ! for the case of 'iadd_mode=1'
  ! ---------------- Initialize variables outside of the y loop
  !
  !--- Note: we need special treatments for the mean flow
  shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
  call derivyc(u00,rf0u,shp,shm,1,1,0)    ! ----- ome3_0 = -du0/dy ---
  if (iadd_mode.eq.-1) then
     call derivyc(u00,tmpy,shp,shm,2,1,0)    ! ddu0/dy2 ---
  end if
  !--- Note: derivative of fluctuation of u is pure-periodic
  call derivyc(w00,rf0w,shp,shm,1,1,0)    ! ----- ome1_0 =  dw0/dy ---
  !if (itemperature.eq.1) then
  !   call derivyc(tb00,rf0tb,shp,shm,1,1,0)    ! -----  =  dtb00/dy ---
  !end if

  t2=MPI_WTIME()
  t_derivyc = t_derivyc + (t2-t1)

  if ((iadd_mode.eq.1).or.(iadd_mode.eq.4).or.(iadd_mode.eq.6)) then
     do j=0,my1
        ! u00(j) =u00(j)    ! now u00 is the mean velocity fluctuation
        rf0u(j)=rf0u(j) +s  ! --- add shear --- now rf0u is du0/dy
     end do
  !elseif ((iadd_mode.eq.2).or.(iadd_mode.eq.2)) then
  !  ! do j=0,my1
  !  !    u00(j) =u00(j)       ! now u00 is the mean velocity fluctuation
  !  !    rf0u(j)=rf0u(j)      ! now rf0u is du0/dy
  !  ! end do
  end if

  ener = 0.d0 ! initialize the history array
  ener2 = 0.d0 ! rev656, and more for LES (rev.1230)
  if ((istep.eq.1).and.(readflag.eq.0)) then
     if (abs(s).gt.0.5) then
        umaxl = 1.d0
     else
        umaxl = uprim ! for the isotropic flow simulation
     end if
  else
     umaxl = 0.d0 ! bug fixed
  endif
  if (itemperature.eq.1) then
     tbmaxl = 0.d0
     enert = 0.d0
  end if

  DO J = JB,JE       !---------------------- start operating by planes

     t1=MPI_WTIME()

     !---------------------- 00 modes, u,w,ome1,ome3
     u1c(0,0) = dcmplx(u00(j),0.d0)  ! does not include mean shear s*y
     u3c(0,0) = dcmplx(w00(j),0.d0)

     o3c(0,0) = dcmplx(-rf0u(j),0.d0)  ! include mean shear -s for the case of iadd_mode.eq.1
     o1c(0,0) = dcmplx( rf0w(j),0.d0)

     if (itemperature.eq.1) then
        call dcopy((mx1+1)*(mz1+1)*2,tempc(0,0,j),1,tbc(0,0),1)
        tbc(0,0) = tb00(j) ! does not include mean temp. d(tb99)/dy
     end if

     ! --------------------- computes non 0 modes of ome1, ome3, u1, u3
     o3c(0,1:mz1) = -ome1c(0,1:mz1,j)/xgam(1:mz1)
     o1c(0,1:mz1) = - phic(0,1:mz1,j)/xgam(1:mz1)
     u3c(0,1:mz1) = - rhvc(0,1:mz1,j)/xgam(1:mz1)
     u1c(0,1:mz1) =  ome2c(0,1:mz1,j)/xgam(1:mz1)

     kkw(0:mx1)=0.d0
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)
     !$OMP DO PRIVATE(kkw) SCHEDULE(STATIC)
     do k=0,mz1
        kkw(1:mx1)=(alp2(1:mx1)+gam2(k))
        o3c(1:mx1,k) = (ome1c(1:mx1,k,j)*xgam(k)-phic(1:mx1,k,j)*xalp(1:mx1)) &
                      /kkw(1:mx1)
        o1c(1:mx1,k) = (ome1c(1:mx1,k,j)*xalp(1:mx1)+phic(1:mx1,k,j)*xgam(k)) &
                      /kkw(1:mx1)
     enddo
     !$OMP END DO
     !$OMP DO PRIVATE(kkw) SCHEDULE(STATIC)
     do k=0,mz1
        kkw(1:mx1)=(alp2(1:mx1)+gam2(k))
        u3c(1:mx1,k) = (rhvc(1:mx1,k,j)*xgam(k)+ome2c(1:mx1,k,j)*xalp(1:mx1)) &
                      /kkw(1:mx1)
        u1c(1:mx1,k) = (rhvc(1:mx1,k,j)*xalp(1:mx1)-ome2c(1:mx1,k,j)*xgam(k)) &
                      /kkw(1:mx1)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL

     t2=MPI_WTIME()
     t_uvw = t_uvw + (t2-t1)

     ! -------------- copy v and omega_2 into their planes
     t1=MPI_WTIME()


     !$OMP PARALLEL SECTIONS
     !$OMP SECTION
     call dcopy((mx1+1)*(mz1+1)*2,rhgc(0,0,j),1,u2c(0,0),1)
     !$OMP SECTION
     call dcopy((mx1+1)*(mz1+1)*2,ome2c(0,0,j),1,o2c(0,0),1)
     !$OMP END PARALLEL SECTIONS
     u2c(0,0) =  vbulk  ! normally this should be zero...(rev934)

     t2=MPI_WTIME()
     t_copy = t_copy + (t2-t1)

     ! c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
     !  at this point:
     !  u1 is u,  u2 is v, u3 is w
     !  o1 is omega_1,  o2 is omega_2,  o3 is omega_3
     !  all variables in Fourierx -- Fourier z -- Physical y
     ! c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c

     ! ---------------- do spectra some day
     ! ---------------- only if my plane contains spectra information

     if (inonNewtonian.eq.1) then

        call dcopy((mx1+1)*(mz1+1)*2, rhvc(0,0,j),1,vyc(0,0),1)  ! dv/dy

        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              uxc(i,k) = u1c(i,k)*xalp(i)   ! du/dx
              uzc(i,k) = u1c(i,k)*cfac      ! du/dz
           end do
        end do
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              vxc(i,k) = u2c(i,k)*xalp(i)   ! dv/dx
              vzc(i,k) = u2c(i,k)*cfac      ! dv/dz
           end do
        end do
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              wxc(i,k) = u3c(i,k)*xalp(i)      ! dw/dx
              wzc(i,k) = u3c(i,k)*cfac         ! dw/dz
           end do
        end do

        !wyc = vzc + o1c     ! dw/dy=dv/dz+o1c
        !uyc = vxc - o3c     ! du/dy=dv/dx-o3c
        ! some 2d buffers are redundunt ...

        call fourxz(uxc,ux,tmpxzr,1,1)
        call fourxz(vyc,vy,tmpxzr,1,1)
        call fourxz(wzc,wz,tmpxzr,1,1)
        call fourxz(2.d0*vzc+o1c,vz,tmpxzr,1,1)
        call fourxz(wyc,wy,tmpxzr,1,1)
        call fourxz(vxc,vx,tmpxzr,1,1)
        call fourxz(2.d0*vxc-o3c,uy,tmpxzr,1,1)
        call fourxz(wxc,wx,tmpxzr,1,1)
        call fourxz(uzc+wxc,uz,tmpxzr,1,1)

        do k = 0, mz1
           do i = 0, mx1
           ! 対称テンソル D_ij の要素（中央差分/他で得た値）
              Dxx = ux(i,k)
              Dyy = vy(i,k)
              Dzz = wz(i,k)
              Dxy = 0.5d0 * (vx(i,k) + uy(i,k))
              Dzx = 0.5d0 * (wx(i,k) + uz(i,k))
              Dyz = 0.5d0 * (wy(i,k) + vz(i,k))
              ! 有効せん断速度 gamma_dot
              gamma_dot = sqrt(2.d0 * (Dxx**2 + Dyy**2 + Dzz**2 + 2.d0*(Dxy**2 + Dzx**2 + Dyz**2)))
              ! Carreau 粘度
              mu(i,k) = mu_inf + (mu_0 - mu_inf) * (1.d0 + (lambda * gamma_dot)**2)**((nnf - 1.d0) / 2.d0)
              ! せん断応力 T_ij
              Txx(i,k) = 2.d0 * mu(i,k) * Dxx
              Tyy(i,k) = 2.d0 * mu(i,k) * Dyy
              Tzz(i,k) = 2.d0 * mu(i,k) * Dzz
              Txy(i,k) = 2.d0 * mu(i,k) * Dxy
              Tzx(i,k) = 2.d0 * mu(i,k) * Dzx
              Tyz(i,k) = 2.d0 * mu(i,k) * Dyz
           end do
        end do
        ! フーリエ変換（物理空間 -> スペクトル空間）
        call fourxz(uxc, Txx, tmpxzr, -1, 1)
        call fourxz(uyc, Txy, tmpxzr, -1, 1)
        call fourxz(uzc, Tzx, tmpxzr, -1, 1)
        call fourxz(vxc, Tyy, tmpxzr, -1, 1)
        call fourxz(vyc, Tyz, tmpxzr, -1, 1)
        call fourxz(vzc, Tzz, tmpxzr, -1, 1)

        txy00(j) = Txy(0,0)
        tyz00(j) = Tyz(0,0)

        ! スペクトル空間での演算（元のコードのまま）
        do k = 0, mz1
           rhvc1(:,k,j) = alp2(:)*uxc(:,k) - 2.d0*xalp(:)*xgam(k)*uzc(:,k) + gam2(k)*vzc(:,k) - (alp2(:)+gam2(k))*vxc(:,k)
           rhvc2(:,k,j) = xalp(:)*uyc(:,k) + xgam(k)*vyc(:,k)
           rhgc1(:,k,j) = xalp(:)*xgam(k)*(uxc(:,k) - vzc(:,k)) + (-gam2(k) + alp2(:))*uzc(:,k)
           rhgc2(:,k,j) = xgam(k)*uyc(:,k) - xalp(:)*vyc(:,k)
        end do

      end if

      if (iuse_LES.eq.1) then
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              uxc(i,k) = u1c(i,k)*xalp(i)   ! du/dx
              uzc(i,k) = u1c(i,k)*cfac      ! du/dz
           end do
        end do
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              vxc(i,k) = u2c(i,k)*xalp(i)   ! dv/dx
              vzc(i,k) = u2c(i,k)*cfac      ! dv/dz
           end do
        end do
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              wxc(i,k) = u3c(i,k)*xalp(i)      ! dw/dx
              wzc(i,k) = u3c(i,k)*cfac         ! dw/dz
           end do
        end do

        !wyc = vzc + o1c     ! dw/dy=dv/dz+o1c
        !uyc = vxc - o3c     ! du/dy=dv/dx-o3c
        ! some 2d buffers are redundunt ...

        call fourxz(uxc,ux,tmpxzr,1,1)
        call fourxz(vyc,vy,tmpxzr,1,1)
        call fourxz(wzc,wz,tmpxzr,1,1)
        call fourxz(2.d0*vzc+o1c,vz,tmpxzr,1,1)
        call fourxz(wyc,wy,tmpxzr,1,1)
        call fourxz(vxc,vx,tmpxzr,1,1)
        call fourxz(2.d0*vxc-o3c,uy,tmpxzr,1,1)
        call fourxz(wxc,wx,tmpxzr,1,1)
        call fourxz(uzc+wxc,uz,tmpxzr,1,1)

        do k=1,mgalz
           do i=1,mgalx
              !SGS(i,k)=dsqrt(2.d0*ux(i,k)**2 + (uy(i,k)+1.d0+vx(i,k))**2 + &
              !               2.d0*vy(i,k)**2 + (uz(i,k)+wx(i,k))**2 +      &
              !               2.d0*wz(i,k)**2 + (vz(i,k)+wy(i,k))**2 )
              !
              Sij2(i,k) = (2.d0*ux(i,k)*ux(i,k) + (uy(i,k)+s)*(uy(i,k)+s) + &
                   &       2.d0*vy(i,k)*vy(i,k) + uz(i,k)*uz(i,k) +      &
                   &       2.d0*wz(i,k)*wz(i,k) + vz(i,k)*vz(i,k) )
           end do
        end do

        ! rev.1317 dealiasing by 2/3-rule
        wyc=0.d0
        call fourxz(wyc,Sij2,tmpxzr,-1,1)
        ! apply 2/3 rule in Fourier space
        call fourxz(wyc,Sij2,tmpxzr,1,1)

        !nu_t=(Cles*Deltag)**2|2*SS|; T=2*nu_t*S

        fac=(Cles*Cles*Deltag*Deltag)
        nuSGS=fac*dsqrt(dabs(Sij2))
        if (s.lt.0.d0) then
           nuSGS=-nuSGS
           ! rev.1426, Note that Smagorinsky model with constatn C_S is not reversible...
        end if

        ! ---------------------------------
        !  do not change this order...
        !
        wx=ux ! save ux => wx
        ux=2.d0*nuSGS*ux        ! T11
        uy=nuSGS*(uy+s)       ! T12 (uy is dv/dx-o3c)
        uz=nuSGS*uz               ! T13 (uz)
        vx=2.d0*nuSGS*vy         ! T22 (vx)
        vy=nuSGS*vz               ! T23 (vy)
        vz=2.d0*nuSGS*wz          ! T33 (vz)

        !Tll = (ux + vx + vz)/3.d0*wx
        !Tll = wx

        ! do k=1,mgalz
        !    do i=1,mgalx
        !
        !       !uy(i,k)=nuSGS(i,k)*0.5d0*(uy(i,k)+vx(i,k))   ! T12=T21
        !       !uz(i,k)=nuSGS(i,k)*0.5d0*(uz(i,k)+wx(i,k))   ! T13=T31
        !
        !       !vx(i,k)=nuSGS(i,k)*vy(i,k)                   ! T22
        !       vy(i,k)=nuSGS(i,k)*0.5d0*(vz(i,k)+wy(i,k))   ! T23=T32
        !       vz(i,k)=nuSGS(i,k)*wz(i,k)                   ! T33
        !    end do
        ! end do

        call fourxz(uxc,ux,tmpxzr,-1,1)
        call fourxz(uyc,uy,tmpxzr,-1,1)
        call fourxz(uzc,uz,tmpxzr,-1,1)
        call fourxz(vxc,vx,tmpxzr,-1,1)
        call fourxz(vyc,vy,tmpxzr,-1,1)
        call fourxz(vzc,vz,tmpxzr,-1,1)

        txy00(j)=uyc(0,0)   ! this includes mean S
        tyz00(j)=vyc(0,0)

        !
        !
        !
        !

        do k=0,mz1
           rhvc1(:,k,j) = alp2(:)*uxc(:,k)-2.d0*xalp(:)*xgam(k)*uzc(:,k)+gam2(k)*vzc(:,k)-(alp2(:)+gam2(k))*vxc(:,k)
           rhvc2(:,k,j) = xalp(:)*uyc(:,k)+xgam(k)*vyc(:,k)
           rhgc1(:,k,j) = xalp(:)*xgam(k)*(uxc(:,k)-vzc(:,k))+(-gam2(k)+alp2(:))*uzc(:,k)
           rhgc2(:,k,j) = xgam(k)*uyc(:,k)-xalp(:)*vyc(:,k)
        end do

     end if ! iuse_LES==1

     ! -- fourier statistics ! itemperature is not yet implemented...
     if (ihist.eq.1) call histt(u1r,u2r,u3r,o1r,o2r,o3r, &
                               u1c,u2c,u3c,o1c,o2c,o3c,rhvc,tbr,tbc,j,'f')

     t1=MPI_WTIME()

     ! 
     call fourxz(u1c,u1r,tmpxzr,1,1)         ! u
     call fourxz(u2c,u2r,tmpxzr,1,1)         ! v
     call fourxz(u3c,u3r,tmpxzr,1,1)         ! w
     call fourxz(o1c,o1r,tmpxzr,1,1)         !omega_1
     call fourxz(o2c,o2r,tmpxzr,1,1)         !omega_2
     call fourxz(o3c,o3r,tmpxzr,1,1)         !omega_3
     if (itemperature.eq.1) then
        call fourxz(tbc,tbr,tmpxzr,1,1)         ! t
     end if
     ! including shear (iadd_mode.eq.1, 4, 5)

     t2=MPI_WTIME()
     t_fourxz = t_fourxz + (t2-t1)

     ! -- physical statistics
     if (ihist.eq.1) call histt(u1r,u2r,u3r,o1r,o2r,o3r, &
                               u1c,u2c,u3c,o1c,o2c,o3c,rhvc,tbr,tbc,j,'p')

     if ((rkstep.eq.1).and.(icfl.eq.1)) then
        if (iuse_LES.eq.1) then
           !if (iadd_visc.eq.1) then
           max_nut  = max(max_nut, maxval(dabs(nuSGS(1:mgalx,1:mgalz))))
           !end if
        else if (inonNewtonian.eq.1) then
           max_nut = max(max_nut, maxval(dabs(mu(1:mgalx,1:mgalz))))
        else
           max_nut = 1.d0/abs(re); ! rev.1426, for negative dissipation...
        end if
        umaxl(0) = max(umaxl(0),maxval(dabs(u1r(1:mgalx,1:mgalz))))
        umaxl(1) = max(umaxl(1),maxval(dabs(u1r(1:mgalx,1:mgalz)+s*y(j))))
        umaxl(2) = max(umaxl(2),maxval(dabs(u2r(1:mgalx,1:mgalz))))
        umaxl(3) = max(umaxl(3),maxval(dabs(u3r(1:mgalx,1:mgalz))))
        if (itemperature.eq.1) then
           ! the local maximum of the intensity of temperature fulctuations
           tbmaxl = max(tbmaxl,maxval(dabs(tbr(1:mgalx,1:mgalz))))
        end if
     end if

     if (itemperature.eq.1) then
        ! computes d(uT)/dx + d(wT)/dz and vT
        ! compute uT,vT,wT
        !$OMP PARALLEL WORKSHARE
        uTr(:,:) = u1r(:,:)*tbr(:,:)
        wTr(:,:) = u3r(:,:)*tbr(:,:)
        vTr(:,:) = u2r(:,:)*tbr(:,:)
        !$OMP END PARALLEL WORKSHARE
        call fourxz(uTc,uTr,tmpxzr,-1,1)
        call fourxz(wTc,wTr,tmpxzr,-1,1)
        call fourxz(vTc,vTr,tmpxzr,-1,1)
        rhtc2(:,:,j) = -vTc  ! -vT
        rf0tb(j)     = -vTc(0,0) ! -vT_00
     end if
     ! 
     ! 
     ! 
     ! 
     ! 
     ! 

     t1=MPI_WTIME()


     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,tmpx) SCHEDULE(STATIC)
     do k=1,mgalz

        call dcopy(mgalx,u2r(1,k),1,tmpx,1)
        u2r(1:mgalx,k) = u2r(1:mgalx,k)*(o3r(1:mgalx,k)+2d0*chi)-u3r(1:mgalx,k)*o2r(1:mgalx,k)
        u3r(1:mgalx,k) = u3r(1:mgalx,k)*o1r(1:mgalx,k)-u1r(1:mgalx,k)*(o3r(1:mgalx,k)+2d0*chi)
        u1r(1:mgalx,k) = u1r(1:mgalx,k)*o2r(1:mgalx,k)-tmpx(1:mgalx)*o1r(1:mgalx,k)

     enddo
     !$OMP END PARALLEL DO

     ! --------------------- at this point
     ! --------------------- u1 : H3
     ! --------------------- u3 : H2 - 2*chi*u
     ! --------------------- u2 : H1 + 2*chi*v  ! not ! + (2*chi - ss)*v
     ! --------------------- back to F-F-T


     t2=MPI_WTIME()
     t_uome = t_uome + (t2-t1)
     !                                     ------  back to fourier
     t1=t2


     call fourxz(u1c,u1r,tmpxzr,-1,1)
     call fourxz(u2c,u2r,tmpxzr,-1,1)
     call fourxz(u3c,u3r,tmpxzr,-1,1)

     if (iadd_mode.eq.-1) then
        call fourxz(o2c,o2r,tmpxzr,-1,1)         !omega_2
     end if

     t2=MPI_WTIME()
     t_fourxz = t_fourxz + (t2-t1)

     !                                     ------  save 00 mode
     rf0u(j)=u2c(0,0)
     rf0w(j)=u1c(0,0)


     t1=MPI_WTIME()

     !
     !
      if (iadd_mode.eq.1) then
         do k=0,mz1
            o1c(:,k) =  xalp(:)*u2c(:,k) + xgam(k)*u1c(:,k)

            !u2c(:,k) = -xalp(:)*u1c(:,k) + xgam(k)*u2c(:,k)
            u2c(:,k) = -xalp(:)*(u1c(:,k) + s*y(j)*ome2c(:,k,j)) &  ! add shear
                 +xgam(k)*u2c(:,k)
            !u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k))
            u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k)) & ! note: alp2 ... is kx**2
                 -s*ome2c(:,k,j)*xgam(k)-s*y(j)*phic(:,k,j)*xalp(:) ! add shear
         enddo
         rhvc(:,:,j)=u1c  ! hv
         rhgc(:,:,j)=u2c  ! hg
         ome1c(:,:,j)=o1c ! dvordy
      else if (iadd_mode.eq.2) then

         do k=0,mz1
            !o1c(:,k) =  xalp(:)*u2c(:,k) + xgam(k)*u1c(:,k)
            ome1c(:,k,j) =  xalp(:)*u2c(:,k) + xgam(k)*u1c(:,k)

            !$$$ u2c(:,k) = -xalp(:)*u1c(:,k) + xgam(k)*u2c(:,k)
            !$$$ u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k))
            ! (rev.262)
            !u2c(:,k) = -xalp(:)*(u1c(:,k) + s*y(j)*oyc(:,k)) &  ! add shear
            !     +xgam(k)*(u2c(:,k) - s*vc(:,k))
            !u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k)) & ! note: alp2 ... is kx**2
            !     - s*y(j)*lapvc(:,k)*xalp(:) ! add shear
            ! (rev.263)
            !u2c(:,k) = -xalp(:)*(u1c(:,k) + s*y(j)*ome2c(:,k,j)) &  ! add shear
            !     +xgam(k)*(u2c(:,k) - s*rhgc(:,k,j))
            !u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k)) & ! note: alp2 ... is kx**2
            !     - s*y(j)*phic(:,k,j)*xalp(:) ! add shear
            ! (rev.264)
            rhgc(:,k,j) = -xalp(:)*(u1c(:,k) + s*y(j)*ome2c(:,k,j)) &  ! add shear
                 +xgam(k)*(u2c(:,k) - s*rhgc(:,k,j))
            rhvc(:,k,j) = -u3c(:,k)*(alp2(:)+gam2(k)) & ! note: alp2 ... is kx**2
                 - s*y(j)*phic(:,k,j)*xalp(:) ! add shear
         enddo

      elseif ((iadd_mode.eq.3).or.(iadd_mode.eq.5)) then
         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)
         if (itemperature.eq.1) then
            ! rhtc (ht) = -(d(uT)/dx + d(wT)/dz)
            !$OMP DO SCHEDULE(STATIC)
            do k=0,mz1
               rhtc(:,k,j) = -xalp(:)*uTc(:,k) - xgam(k)*wTc(:,k) - bym*rhgc(:,k,j)
            end do
            !$OMP END DO
         end if
         !$OMP DO SCHEDULE(STATIC)
         do k=0,mz1
            ome1c(:,k,j) =  xalp(:)*u2c(:,k) + xgam(k)*u1c(:,k)
            rhgc(:,k,j)  = -xalp(:)*u1c(:,k) &
                 + xgam(k)*(u2c(:,k) - s*rhgc(:,k,j))
            rhvc(:,k,j)  = -u3c(:,k)*(alp2(:)+ gam2(k))
         enddo
         !$OMP END DO
         !$OMP END PARALLEL

      elseif (iadd_mode.eq.4) then
         do k=0,mz1
            o1c(:,k) =  xalp(:)*u2c(:,k) + xgam(k)*u1c(:,k)
            u2c(:,k) = -xalp(:)*u1c(:,k) + xgam(k)*u2c(:,k)
            u1c(:,k) = -u3c(:,k)*(alp2(:)+gam2(k)) & ! note: alp2 ... is kx**2
                 -s*ome2c(:,k,j)*xgam(k) ! add shear
         enddo
         rhvc(:,:,j)=u1c  ! hv
         rhgc(:,:,j)=u2c  ! hg
         ome1c(:,:,j)=o1c ! dvordy

      elseif (iadd_mode.eq.-1) then
         ! linearized nonlinear term.
         do k=0,mz1
            ! set u1c, u2c, u3c = 0 ??
            ome1c(:,k,j) = 0.d0

            rhvc(:,k,j)  = 0.d0
            !rhvc(:,k,j) = -u00(j)*xalp(:)*rhvc(:,k,j) + tmpy(j)*xalp(:)*rhgc(:,k,j)

            rhgc(:,k,j)  = -s*xgam(k)*rhgc(:,k,j)
            !rhgc(:,k,j) = -u00(j)*xalp(:)*o2c(:,k) + xgam(k)*( (rf0u(j) - s)*rhgc(:,k,j))

         enddo

      endif

     t2=MPI_WTIME()
     t_nl = t_nl + (t2-t1)
     t1= t2

     ! --------------------- copy planes into outputs
     !rhvc(:,:,j)=u1c  ! hv
     !rhgc(:,:,j)=u2c  ! hg
     !ome1c(:,:,j)=o1c ! dvordy

     t2=MPI_WTIME()
     t_copy = t_copy + (t2-t1)

  ENDDO   ! -------------------------- finishes the y loop

  !deallocate(oyc,vc) ! for add shear

  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
  !  o1: dH1/dx + dH3/dz
  !  u2: -dH3/dx + dH1/dz
  !  u1: d^2H2/dx^2 + d^2H2/dz^2
  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

  ! -------------------- computes Deltat and write stats
  if ((rkstep.eq.1).and.(icfl.eq.1)) then
     ! update Deltat

     umax=0.d0 ! grobal value (umax(0) is fluctuation in x-dir)
     nutmax=0.d0
     call MPI_ALLREDUCE(umaxl,umax,4,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(max_nut,nutmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
     if (itemperature.eq.1) then
        tbmax=0.d0
        call MPI_ALLREDUCE(tbmaxl,tbmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
     end if

     if (ifix_dt.eq.1) then
        !if (istep.eq.1) then
           !dt = dx/14.d0 ! set smaller step for the first step
           dt = dx/4.d0 ! set smaller step for the first step
        !else
        !   dt = dx/8.d0 ! see Shumann et. al. (JFM )
        !endif
        if ((iadd_mode.ge.3).and.(iadd_mode.le.5)) then
           ! umax(0) is uprime
           CFL = dt*max(umax(0)*pi/dx,umax(2)*2.5d0/dy,umax(3)*pi/dz)
        else
           CFL = dt*max(umax(1)*pi/dx,umax(2)*2.5d0/dy,umax(3)*pi/dz)
        end if
        if (iuse_LES.eq.1) then
           if (iadd_visc.eq.1) then
              CFLv = max(abs(1.d0/re)*dt*max(pi*pi/dx/dx , &
                                          6.25d0/dy/dy, &
                                          pi*pi/dz/dz), &
                          2.d0*nutmax*dt*max(pi*pi/dx/dx , &
                                        6.25d0/dy/dy, &
                                        pi*pi/dz/dz))
           else if (iadd_visc.eq.0) then
              CFLv = 2.d0*nutmax*dt*max(pi*pi/dx/dx , &
                                   6.25d0/dy/dy, &
                                   pi*pi/dz/dz)
           end if
        else if (inonNewtonian.eq.1) then
           CFLv = 2.d0*nutmax*dt*max(pi*pi/dx/dx , &
                                     6.25d0/dy/dy, &
                                     pi*pi/dz/dz)
        else ! normal run without any model
           CFLv = (1.d0/re)*dt*max(pi*pi/dx/dx, &
                                   6.25d0/dy/dy, &
                                   pi*pi/dz/dz)
        end if
        if ((myid.eq.0).and.(iwrote.eq.0))  then
          write(*,*) ' fixed dt, = ',dt,' CFL = ', CFL, CFLv
          write(*,*) ' fixed dt, = ',dt,' CFLt = ', CFLt, CFLk
        end if
        iwrote=1
     else
        if ((iadd_mode.ge.3).and.(iadd_mode.le.5)) then
           dt = 1d0/max(umax(0)*cflx,umax(2)*cfly,umax(3)*cflz)
           ! use uprim_max
        else
           dt = 1d0/max(umax(1)*cflx,umax(2)*cfly,umax(3)*cflz)
        endif
        ! for explicit
        ! dt = 1d0/max(umax(1)*cflx,umax(2)*cfly,umax(3)*cflz,cflvv)
        if (iuse_LES.eq.1) then
           if (iadd_visc.eq.1) then
              cflvv = max(abs(1.d0/re)/CFL*max(pi*pi/dx/dx , &
                                          6.25d0/dy/dy, &
                                          pi*pi/dz/dz), &
                          nutmax/CFL*max(pi*pi/dx/dx , &
                                        6.25d0/dy/dy, &
                                        pi*pi/dz/dz))
           else
              cflvv = nutmax/CFL*max(pi*pi/dx/dx , &
                                   6.25d0/dy/dy, &
                                   pi*pi/dz/dz)
           end if
        else if (inonNewtonian.eq.1) then
           cflvv = nutmax/CFL*max(pi*pi/dx/dx , &
                                  6.25d0/dy/dy, &
                                  pi*pi/dz/dz)
        end if
        CFLv = CFL*cflvv*dt
        if (itemperature.eq.1) then
           CFLk = CFL*cflkk*dt
        end if
        if (explicit.eq.1) then
           if (CFLv.gt.CFL) then
              if (CFLk.gt.CFLv) then
                 dt = 1.d0/cflkk
                 if ( (mod(istep-1,nhist).eq.0).and.(nohist.eq.0)) then
                    if (iskip_screenout.eq.0) then
                       if (myid.eq.0) write(*,*) 'EXPLICIT: dt is reduced, CFLk <= CFL'
                    end if
                 end if
              else
                 dt = 1.d0/cflvv
                 if ( (mod(istep-1,nhist).eq.0).and.(nohist.eq.0)) then
                    if (iskip_screenout.eq.0) then
                       if (myid.eq.0) write(*,*) 'EXPLICIT: dt is reduced, CFLv <= CFL'
                    end if
                 end if
              end if
           end if
        elseif (explicit.eq.0) then
           if (CFLv.gt.2d0) then
              if (CFLk.gt.CFLv) then
                 dt = 1.d0/cflkk/CFL*6.0
                 if ( (mod(istep-1,nhist).eq.0).and.(nohist.eq.0)) then
                    if (iskip_screenout.eq.0) then
                       if (myid.eq.0) write(*,*) 'IMPLICIT: dt is reduced, CFLk <= CFL'
                    end if
                 end if
              else
                 dt = 1.d0/cflvv/CFL*6.0
                 if ( (mod(istep-1,nhist).eq.0).and.(nohist.eq.0)) then
                    if (iskip_screenout.eq.0) then
                       if (myid.eq.0) write(*,*) 'IMPLICIT: dt is reduced, CFLv <= CFL'
                    end if
                 end if
              end if
           end if
        end if
     end if

     ! for safety, pass dt for all proc.
     call MPI_ALLREDUCE(dt,Deltat,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)

     if (itemperature.eq.1) then
        dtrv=Re/Deltat ! bug fixed at rev.1724
        dtrk=(1.d0/fkappa)/Deltat ! bug fixed at rev.1724
     else
        dtr=Re/Deltat
     end if
     !if ((Deltat.lt.1.d-10).or.(umax(0)/Lz.gt.5.d0)) then
     if ((Deltat.lt.1.d-10)) then
        if (myid.eq.0) write(*,*) ' Diverged: istep, dt, time = ',istep, Deltat, time
        if (myid.eq.0) write(*,*) ' dt_x = ', 1d0/(umax(1)*cflx)
        if (myid.eq.0) write(*,*) ' dt_y = ', 1d0/(umax(2)*cfly)
        if (myid.eq.0) write(*,*) ' dt_z = ', 1d0/(umax(3)*cflz)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        istop=1 ! fixed at rev.1293
     endif

     ! idump_mode==1 for creating a pure periodic field and stop
     ! this can be used for arnoldi or newton algorithm
     if ((time+Deltat).gt.dumptime) then

        Deltat = abs(dumptime - time)
        if (itemperature.eq.1) then
           dtrv=Re/Deltat  ! bug fixed at rev.1724
           dtrk=(1.d0/fkappa)/Deltat ! bug fixed at rev.1724
        else
           dtr=Re/Deltat
        end if
        if (dumptime.gt.(etime-tiny)) istop = 1
        dumptime = dumptime + dtimag ! update next dumptime

        idump = 1
        if (myid.eq.0) write(*,*) '  dump file at next step, idump=', idump

     elseif ((time+Deltat).gt.etime) then

        Deltat = abs(etime - time)
        if (itemperature.eq.1) then
           dtrv=Re/Deltat ! bug fixed at rev.1724
           dtrk=(1.d0/fkappa)/Deltat ! bug fixed at rev.1724
        else
           dtr=Re/Deltat
        end if
        istop = 1 ! file dump + stop

     end if
  endif

  ! write screen output
  if (ihist.eq.1) call histt(u1r,u2r,u3r,o1r,o2r,o3r, &
                            u1c,u2c,u3c,o1c,o2c,o3c,rhvc,tbr,tbc,0,'w')

  t1 = MPI_WTIME()
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  t2 = MPI_WTIME()
  t_waiting=t_waiting + t2-t1
  t1 = t2

  !  ------  collect the 00 modes everywhere (not slow??)
  !  ------  [rev.406]
  !do pproc=0,numerop-1
  !   iproc=ieor(myid,pproc)
  !   if (iproc.ne.myid) then
  !      mmyr=jend(iproc)-jbeg(iproc)+1
  !      call MPI_SENDRECV(rf0u(jb),mmy,MPI_REAL8,iproc,0,rf0u(jbeg(iproc)), &
  !           mmyr,MPI_REAL8,iproc,0,MPI_COMM_WORLD,istat,ierr)
  !
  !      call MPI_SENDRECV(rf0w(jb),mmy,MPI_REAL8,iproc,0,rf0w(jbeg(iproc)), &
  !           mmyr,MPI_REAL8,iproc,0,MPI_COMM_WORLD,istat,ierr)
  !   endif
  !
  !enddo
  ! --- (rev.373) replaced by ALLREDUCE (faster than above??)
  tmpy=0.d0 ! do not forget this ! memcopy
  !tmpy(jb:je)=rf0u(jb:je)
  !call vcopy(je-jb+1,rf0u(jb),tmpy(jb))
  call dcopy(je-jb+1,rf0u(jb),1,tmpy(jb),1)
  call MPI_ALLREDUCE(tmpy,rf0u,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  tmpy=0.d0 ! do not forget this ! memcopy
  !tmpy(jb:je)=rf0w(jb:je)
  call dcopy(je-jb+1,rf0w(jb),1,tmpy(jb),1)
  call MPI_ALLREDUCE(tmpy,rf0w,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  if (iuse_LES.eq.1) then
     tmpy=0.d0
     call dcopy(je-jb+1,txy00(jb),1,tmpy(jb),1)
     call MPI_ALLREDUCE(tmpy,txy00,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     tmpy=0.d0
     call dcopy(je-jb+1,tyz00(jb),1,tmpy(jb),1)
     call MPI_ALLREDUCE(tmpy,tyz00,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end if ! iuse_LES==1
  if (inonNewtonian.eq.1) then
     tmpy=0.d0
     call dcopy(je-jb+1,txy00(jb),1,tmpy(jb),1)
     call MPI_ALLREDUCE(tmpy,txy00,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     tmpy=0.d0
     call dcopy(je-jb+1,tyz00(jb),1,tmpy(jb),1)
     call MPI_ALLREDUCE(tmpy,tyz00,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end if ! inonNewtonian==1

  if (itemperature.eq.1) then
     tmpy=0.d0
     call dcopy(je-jb+1,rf0tb(jb),1,tmpy(jb),1)
     call MPI_ALLREDUCE(tmpy,rf0tb,my,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end if


  t2 = MPI_WTIME()
  t_otra = t_otra + (t2-t1)


end subroutine hvhgt

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      hist, the diagnostics
!c      flag = 'f',  things that can be done in fourier, dv is for eps
!c      flag = 'p',                          in physical
!c      flag = 'w',  write screen output (energy etc...)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine histt(u1r,u2r,u3r,o1r,o2r,o3r,u1c,u2c,u3c,o1c,o2c,o3c,dv,&
     &           tbr,tbc,j,flag)

  use ctes
  use statistics
  use running
  use bcs
  use LES!,only:nutmax,Cles,nuSGS,Sij2,Tll,txy00, iuse_LES, iadd_visc

  implicit none
  include "mpif.h"

  integer j,i,ierr,k,kk
  real(8) up,vp,wp,o1p,o2p,o3p,o12p,uvr,uwr,vwr,eps,Pk
  real(8) tbp,utbp,vtbp,wtbp
  real(8),dimension(nner) :: uner
  real(8),dimension(nner2) :: uner2
  real(8),dimension(ntner) :: tner
  complex(8),dimension(0:mx1,0:mz1,jb:je) :: dv
  complex(8) cc
  character*1 flag
  real*8      nufac, nuty, trTij, disp_eddy

!---------- 6 * (mgalx+2) * mgalz buffer planes, THESE PLANES SHARE STORAGE !!
  real(8),dimension(mgalx+2,mgalz) :: u1r, u2r, u3r, o1r, o2r, o3r, tbr
  complex(8),dimension(0:mx1,0:mz1) :: u1c,u2c, u3c, o1c, o2c, o3c, tbc
!----------------------------------------------------------------------------
  !if ((myid.eq.0).and.(j.eq.0)) write(*,*) istep,j, 'in hist()', flag
  if (flag.eq.'f') then        !  fourier things

     !u1c(0,0) = u1c(0,0) + dcmplx(u00cut,0d0)
     !
     !if ((iadd_mode.eq.2).or.(iadd_mode.eq.3)) then
     !   o3c(0,0) = o3c(0,0) - s
     !end if
     !u1c(0,0) = u1c(0,0) + s*y(j)

     ! --- intensities ---

     up = sum(u1c(0,:)*dconjg(u1c(0,:))) &
          + 2d0*sum(u1c(1:mx1,:)*dconjg(u1c(1:mx1,:)))
     vp = sum(u2c(0,:)*dconjg(u2c(0,:))) &
          + 2d0*sum(u2c(1:mx1,:)*dconjg(u2c(1:mx1,:)))
     wp = sum(u3c(0,:)*dconjg(u3c(0,:))) &
          + 2d0*sum(u3c(1:mx1,:)*dconjg(u3c(1:mx1,:)))

     uvr = sum(u1c(0,:)*dconjg(u2c(0,:))) &
          + 2d0*sum(u1c(1:mx1,:)*dconjg(u2c(1:mx1,:)))
     uwr = sum(u1c(0,:)*dconjg(u3c(0,:))) &
          + 2d0*sum(u1c(1:mx1,:)*dconjg(u3c(1:mx1,:)))
     vwr = sum(u2c(0,:)*dconjg(u3c(0,:))) &
          + 2d0*sum(u2c(1:mx1,:)*dconjg(u3c(1:mx1,:)))

     o1p = sum(o1c(0,:)*dconjg(o1c(0,:))) &
          + 2d0*sum(o1c(1:mx1,:)*dconjg(o1c(1:mx1,:)))
     o2p = sum(o2c(0,:)*dconjg(o2c(0,:))) &
          + 2d0*sum(o2c(1:mx1,:)*dconjg(o2c(1:mx1,:)))
     o3p = sum(o3c(0,:)*dconjg(o3c(0,:))) &
          + 2d0*sum(o3c(1:mx1,:)*dconjg(o3c(1:mx1,:)))

     o12p= sum(o1c(0,:)*dconjg(o2c(0,:))) &
          + 2d0*sum(o1c(1:mx1,:)*dconjg(o2c(1:mx1,:)))

     if (itemperature.eq.1) then
        tbp = sum(tbc(0,:)*dconjg(tbc(0,:))) &
             + 2d0*sum(tbc(1:mx1,:)*dconjg(tbc(1:mx1,:)))
        utbp = sum(u1c(0,:)*dconjg(tbc(0,:))) &
             + 2d0*sum(u1c(1:mx1,:)*dconjg(tbc(1:mx1,:)))
        vtbp = sum(u2c(0,:)*dconjg(tbc(0,:))) &
             + 2d0*sum(u2c(1:mx1,:)*dconjg(tbc(1:mx1,:)))
        wtbp = sum(u3c(0,:)*dconjg(tbc(0,:))) &
             + 2d0*sum(u3c(1:mx1,:)*dconjg(tbc(1:mx1,:)))
     end if

     ! including shear for the case (iadd_mode.eq.1)

     stats(4,j) = stats(4,j) + up
     stats(5,j) = stats(5,j) + vp
     stats(6,j) = stats(6,j) + wp
     stats(7,j) = stats(7,j) + uvr
     stats(8,j) = stats(8,j) + uwr
     stats(9,j) = stats(9,j) + vwr

     stats(13,j) = stats(13,j) + o1p
     stats(14,j) = stats(14,j) + o2p
     stats(15,j) = stats(15,j) + o3p

     ! for sta4
     stats2(4,j) = stats2(4,j) + up*Deltat
     stats2(5,j) = stats2(5,j) + vp*Deltat
     stats2(6,j) = stats2(6,j) + wp*Deltat
     stats2(7,j) = stats2(7,j) + uvr*Deltat
     stats2(8,j) = stats2(8,j) + uwr*Deltat
     stats2(9,j) = stats2(9,j) + vwr*Deltat

     stats2(13,j) = stats2(13,j) + o1p*Deltat
     stats2(14,j) = stats2(14,j) + o2p*Deltat
     stats2(15,j) = stats2(15,j) + o3p*Deltat

     ener(1)=ener(1) + up
     ener(2)=ener(2) + vp
     ener(3)=ener(3) + wp

     ener(4)=ener(4) + uvr
     ener(5)=ener(5) + uwr
     ener(6)=ener(6) + vwr

     ener(7)=ener(7) + o1p
     ener(8)=ener(8) + o2p
     ener(9)=ener(9) + o3p

     if (itemperature.eq.1) then

        tstats(1,j) = tstats(1,j) + dreal(tbc(0,0))*Deltat
        tstats(2,j) = tstats(2,j) + tbp*Deltat
        tstats(3,j) = tstats(3,j) + utbp*Deltat
        tstats(4,j) = tstats(4,j) + vtbp*Deltat
        tstats(5,j) = tstats(5,j) + wtbp*Deltat

        enert(1) = enert(1) + tbp
        enert(2) = enert(2) + utbp
        enert(3) = enert(3) + vtbp
        enert(4) = enert(4) + wtbp

     end if

     ! ------------ means
     stats(1,j) = stats(1,j)+u1c(0,0)     ! -- um
     stats(2,j) = stats(2,j)+u2c(0,0)     ! -- vm
     stats(3,j) = stats(3,j)+u3c(0,0)     ! -- wm
     stats(10,j) = stats(10,j)+o1c(0,0)   ! -- o1m
     stats(11,j) = stats(11,j)+o2c(0,0)   ! -- o2m
     stats(12,j) = stats(12,j)+o3c(0,0)   ! -- o3m
     ! for sta4
     stats2(1,j) = stats2(1,j)+u1c(0,0)*Deltat      ! -- um
     stats2(2,j) = stats2(2,j)+u2c(0,0)*Deltat      ! -- vm
     stats2(3,j) = stats2(3,j)+u3c(0,0)*Deltat      ! -- wm
     stats2(10,j) = stats2(10,j)+o1c(0,0)*Deltat    ! -- o1m
     stats2(11,j) = stats2(11,j)+o2c(0,0)*Deltat    ! -- o2m
     stats2(12,j) = stats2(12,j)+o3c(0,0)*Deltat    ! -- o3m

     ! add shear effect (only for iadd_mode = 2,3,5)
     o3c(0,0) = o3c(0,0) -s
     if (iuse_LES.eq.1) then
        if (iadd_visc.eq.1) then
           nufac=(1.d0/Re)
        else
           nufac=0.d0
        end if
        nuty = sum(nuSGS(1:mgalx,1:mgalz))/dfloat(mgalx*mgalz)

        ! Total input using eddy viscosity??
        ! Pk = - (o3c(0,0))*(-uvr) + nuty*o3c(0,0)*o3c(0,0) &
        !      + (nufac)*(o3c(0,0)*o3c(0,0)) ! this is not perfect
        !
        Pk = - o3c(0,0)*(-uvr + txy00(j) )  &
             + (nufac)*(o3c(0,0)*o3c(0,0)) ! fixed at rev.1628.
     else if (inonNewtonian.eq.1) then
        Pk = - o3c(0,0)*(-uvr + txy00(j))
     else
        Pk = - (o3c(0,0))*(-uvr) + (o3c(0,0))*(o3c(0,0))*(1.d0/Re)
     end if
     ! total energy input, - <dudy>(y) <uv>(y)

     stats(20,j) =stats(20,j) + Pk
     stats(21,j) =stats(21,j) + o12p
     if (iuse_LES.eq.1) then
        stats(22,j) =stats(22,j) + nuty
     end if
     ! for sta4
     stats2(20,j)=stats2(20,j) + Pk*Deltat
     stats2(21,j)=stats2(21,j) + o12p*Deltat
     if (iuse_LES.eq.1) then
        stats2(22,j) =stats2(22,j) + nuty*Deltat
     end if

     ener2(1) = ener2(1) + Pk
     ener2(2) = ener2(2) + o12p
     if (iuse_LES.eq.1) then
        ener2(3) = ener2(3) + nuty
     end if
     ! ------------ dissipation  --- 2*<Sij*Sij>_{00}
     ! eps == sum(Sij2(1:mgalx,1:mgalz))/dfloat(mgalx*mgalz) ! checked..
     ! --------------------
     eps = sum(gam2*( u1c(0,:)*dconjg(u1c(0,:)) &
                     +u2c(0,:)*dconjg(u2c(0,:)) &
                     +u3c(0,:)*dconjg(u3c(0,:)) )) &
          +sum( dv(0,:,j)*dconjg(dv(0,:,j)) + o3c(0,:)*dconjg(o3c(0,:)) )
     do k=0,mz1
        eps = eps + 2.d0*sum( (alp2(1:mx1)+gam2(k)) &
                             *( u1c(1:mx1,k)*dconjg(u1c(1:mx1,k)) &
                               +u2c(1:mx1,k)*dconjg(u2c(1:mx1,k)) &
                               +u3c(1:mx1,k)*dconjg(u3c(1:mx1,k)) ) &
                             +dv(1:mx1,k,j)*dconjg(dv(1:mx1,k,j)))
        cc = o1c(0,k) + xgam(k)*u2c(0,k)
        eps = eps + cc*dconjg(cc)
        do i=1,mx1
           cc = o1c(i,k) + xgam(k)*u2c(i,k) ! -dw/dy
           eps = eps + 2d0*cc*dconjg(cc)
           cc = o3c(i,k) - xalp(i)*u2c(i,k) ! -du/dz
           eps = eps + 2d0*cc*dconjg(cc)
        enddo
     enddo
     stats(16,j) = stats(16,j) + eps
     ! for sta4
     stats2(16,j) = stats2(16,j) + eps*Deltat

     ener(10) = ener(10) + eps
     if (iuse_LES.eq.1) then
        ! the energy transfer rate to subgrid scale (dissipation using eddy viscosity)
        disp_eddy = sum(nuSGS(1:mgalx,1:mgalz)*Sij2(1:mgalx,1:mgalz))/dfloat(mgalx*mgalz)
        ener2(4) = ener2(4) + disp_eddy
        stats(23,j) = stats(23,j) + disp_eddy
        stats2(23,j) = stats2(23,j) + disp_eddy*Deltat
        if (nstat.lt.23) then
           write(*,*) 'check nstat, stop'
           stop
        end if
     end if

     if (nner2.lt.4) then
        write(*,*) 'set correct nner2, stop'
        stop
     end if

     o3c(0,0) = o3c(0,0) + s ! remove shear effect (only for iadd_mode = 2,3,5)
     !
     if ((noescru.eq.0).or.(nohist.eq.0)) then
        ! ------------ spectra for HSF (nspec=8)
        do kk=0,mz1
           !k=icx(kk)
           k=kk
           spl(1,0,k) = spl(1,0,k) + u1c(0,kk)*dconjg(u1c(0,kk))
           spl(2,0,k) = spl(2,0,k) + u2c(0,kk)*dconjg(u2c(0,kk))
           spl(3,0,k) = spl(3,0,k) + u3c(0,kk)*dconjg(u3c(0,kk))
           spl(4,0,k) = spl(4,0,k) + real(u1c(0,kk)*dconjg(u2c(0,kk)))
           spl(5,0,k) = spl(5,0,k) + o1c(0,kk)*dconjg(o1c(0,kk))
           spl(6,0,k) = spl(6,0,k) + o2c(0,kk)*dconjg(o2c(0,kk))
           spl(7,0,k) = spl(7,0,k) + o3c(0,kk)*dconjg(o3c(0,kk))
           spl(8,0,k) = spl(8,0,k) + real(o1c(0,kk)*dconjg(o2c(0,kk)))
           do i=1,mx1
              spl(1,i,k) = spl(1,i,k) + 2.0*u1c(i,kk)*dconjg(u1c(i,kk))
              spl(2,i,k) = spl(2,i,k) + 2.0*u2c(i,kk)*dconjg(u2c(i,kk))
              spl(3,i,k) = spl(3,i,k) + 2.0*u3c(i,kk)*dconjg(u3c(i,kk))
              spl(4,i,k) = spl(4,i,k) + 2.0*dreal(u1c(i,kk)*dconjg(u2c(i,kk))) ! energy only
              spl(5,i,k) = spl(5,i,k) + 2.0*o1c(i,kk)*dconjg(o1c(i,kk))
              spl(6,i,k) = spl(6,i,k) + 2.0*o2c(i,kk)*dconjg(o2c(i,kk))
              spl(7,i,k) = spl(7,i,k) + 2.0*o3c(i,kk)*dconjg(o3c(i,kk))
              spl(8,i,k) = spl(8,i,k) + 2.0*dreal(o1c(i,kk)*dconjg(o2c(i,kk))) ! energy only
           enddo
        enddo

        ! ------------ spectra for HSF (nspec=8) for sta4
        do kk=0,mz1
           !k=icx(kk)
           k=kk
           spl2(1,0,k) = spl2(1,0,k) + u1c(0,kk)*dconjg(u1c(0,kk))*Deltat
           spl2(2,0,k) = spl2(2,0,k) + u2c(0,kk)*dconjg(u2c(0,kk))*Deltat
           spl2(3,0,k) = spl2(3,0,k) + u3c(0,kk)*dconjg(u3c(0,kk))*Deltat
           spl2(4,0,k) = spl2(4,0,k) + real(u1c(0,kk)*dconjg(u2c(0,kk)))*Deltat
           spl2(5,0,k) = spl2(5,0,k) + o1c(0,kk)*dconjg(o1c(0,kk))*Deltat
           spl2(6,0,k) = spl2(6,0,k) + o2c(0,kk)*dconjg(o2c(0,kk))*Deltat
           spl2(7,0,k) = spl2(7,0,k) + o3c(0,kk)*dconjg(o3c(0,kk))*Deltat
           spl2(8,0,k) = spl2(8,0,k) + real(o1c(0,kk)*dconjg(o2c(0,kk)))*Deltat
           do i=1,mx1
              spl2(1,i,k) = spl2(1,i,k) + 2.0*u1c(i,kk)*dconjg(u1c(i,kk))*Deltat
              spl2(2,i,k) = spl2(2,i,k) + 2.0*u2c(i,kk)*dconjg(u2c(i,kk))*Deltat
              spl2(3,i,k) = spl2(3,i,k) + 2.0*u3c(i,kk)*dconjg(u3c(i,kk))*Deltat
              spl2(4,i,k) = spl2(4,i,k) + 2.0*dreal(u1c(i,kk)*dconjg(u2c(i,kk)))*Deltat ! energy only
              spl2(5,i,k) = spl2(5,i,k) + 2.0*o1c(i,kk)*dconjg(o1c(i,kk))*Deltat
              spl2(6,i,k) = spl2(6,i,k) + 2.0*o2c(i,kk)*dconjg(o2c(i,kk))*Deltat
              spl2(7,i,k) = spl2(7,i,k) + 2.0*o3c(i,kk)*dconjg(o3c(i,kk))*Deltat
              spl2(8,i,k) = spl2(8,i,k) + 2.0*dreal(o1c(i,kk)*dconjg(o2c(i,kk)))*Deltat ! energy only
           enddo
        enddo
        if (itemperature.eq.1) then
           do kk=0,mz1
              k=kk
              splt(1,0,k) = splt(1,0,k) + tbc(0,kk)*dconjg(tbc(0,kk))*Deltat
              splt(2,0,k) = splt(2,0,k) + real(u1c(0,kk)*dconjg(tbc(0,kk)))*Deltat
              splt(3,0,k) = splt(3,0,k) + real(u2c(0,kk)*dconjg(tbc(0,kk)))*Deltat
              splt(4,0,k) = splt(4,0,k) + real(u3c(0,kk)*dconjg(tbc(0,kk)))*Deltat
              do i=1,mx1
                 splt(1,i,k) = splt(1,i,k) + 2.0*tbc(i,kk)*dconjg(tbc(i,kk))*Deltat
                 splt(2,i,k) = splt(2,i,k) + 2.0*dreal(u1c(i,kk)*dconjg(tbc(i,kk)))*Deltat
                 splt(3,i,k) = splt(3,i,k) + 2.0*dreal(u2c(i,kk)*dconjg(tbc(i,kk)))*Deltat
                 splt(4,i,k) = splt(4,i,k) + 2.0*dreal(u3c(i,kk)*dconjg(tbc(i,kk)))*Deltat
              enddo
           enddo
        end if
        if (nspec.ne.8) then
           write(*,*) 'nspec should be set in init.f90',nspec
        end if
        if (nspect.ne.4) then
           write(*,*) 'nspect should be set in init.f90',nspect
        end if

     end if ! 

     !u1c(0,0) = u1c(0,0) - dcmplx(u00cut,0d0)
     !if ((iadd_mode.eq.2).or.(iadd_mode.eq.3)) then
     !   o3c(0,0) = o3c(0,0) + s
     !   ! remove mean shear again
     !end if
     !u1c(0,0) = u1c(0,0) - s*y(j)

  elseif (flag.eq.'p') then         ! --- things to be computed in physical

     ! u1r = u1r + u00cut
     ! Note: u00cut used be introduced for reducing the absolute value
     !       of mean flow, which contributes to CFL condition
     !
     !u1r = u1r  + s*y(j) ! rev.625 do not add for statistics

     stats(17,j) = stats(17,j) &
          + sum(u2r(1:mgalx,1:mgalz)*u1r(1:mgalx,1:mgalz)**2)   ! -- vuu
     stats(18,j) = stats(18,j) &
          + sum(u2r(1:mgalx,1:mgalz)*u3r(1:mgalx,1:mgalz)**2)   ! -- vww
     stats(19,j) = stats(19,j) + sum(u2r(1:mgalx,1:mgalz)**3)   ! -- vvv

     stats2(17,j) = stats2(17,j) &
          + sum(u2r(1:mgalx,1:mgalz)*u1r(1:mgalx,1:mgalz)**2)*Deltat   ! -- vuu
     stats2(18,j) = stats2(18,j) &
          + sum(u2r(1:mgalx,1:mgalz)*u3r(1:mgalx,1:mgalz)**2)*Deltat   ! -- vww
     stats2(19,j) = stats2(19,j) + sum(u2r(1:mgalx,1:mgalz)**3)*Deltat   ! -- vvv

     !if (iuse_LES.eq.1) then
     !   !! The trace of Tij may need to be added...but they are small
     !
     !   trTij = sum(Tll(1:mgalx,1:mgalz)*(u1r(1:mgalx,1:mgalz) + u2r(1:mgalx,1:mgalz) + u3r(1:mgalx,1:mgalz))/3.d0)/dfloat(mgalx*mgalz)
     !   Pk = Pk - trTij
     !   stats(20,j) =stats(20,j) - trTij
     !   ener2(1) = ener2(1) - trTij
     !end if

     !u1r = u1r - u00cut
     !u1r = u1r - s*y(j)

     ! moved to out of hist(...) to use icfl = 1
     !umax(1) = max(umax(1),maxval(dabs(u1r(1:mgalx,1:mgalz)+s*y(j))))
     !umax(2) = max(umax(2),maxval(dabs(u2r(1:mgalx,1:mgalz))))
     !umax(3) = max(umax(3),maxval(dabs(u3r(1:mgalx,1:mgalz))))

  else                             ! --- collect, compute dt, and write

     ! count the number of accumulated fields for the statistics
     nacumsp = nacumsp +1
     nacum   = nacum+1
     total_dt = total_dt + Deltat ! from sta4

     call MPI_ALLREDUCE(ener,uner,nner,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
     uner(1:3)=dsqrt(dabs(uner(1:3)/dfloat(my))) ! characteristic velocities...
     uner(4:6)=uner(4:6)/dfloat(my) ! bug fixed 2012/Aug/10 (rev.510)
     uner(7:9)=dsqrt(dabs(uner(7:9)/dfloat(my))) ! characteristic vorticities...
     uner(10)=uner(10)/dfloat(my)

     call MPI_ALLREDUCE(ener2,uner2,nner2,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
     uner2(1:nner2)=uner2(1:nner2)/dfloat(my) ! Pk, o1*o2, nuty, disp_eddy

     if (itemperature.eq.1) then
        call MPI_ALLREDUCE(enert,tner,ntner,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
        tner(1)=dsqrt(dabs(tner(1)/dfloat(my)))
        tner(2:4)=tner(2:4)/dfloat(my)
     end if

     if (myid.eq.0) then
        if (nocf.eq.0) then
           !write(iocf ,'(1X,F14.6,16(1X,E14.6),2(1X,F14.6),5(1X,E14.6))') &
           !     time,Deltat,uner,CFL,CFLv,umax(1:3),xoffb,xofft,umax(0),uner2
           write(iocfb) time,Deltat,uner,CFL,CFLv,umax(1:3),xoffb,xofft,umax(0),uner2, &
                &  tbmax,tner
           write(iotfb) time,Deltat,tbmax,tner
           if (iskip_screenout.eq.0) then
              write(* ,'(a,1X,I6,F14.6,16(1X,E14.6),2(1X,F14.6),5(1X,E14.6))') &
                   'UNER:',istep,time,Deltat, &
                   uner,CFL,CFLv,umax(1:3),xoffb,xofft,umax(0),uner2
              if (itemperature.eq.1) then
                 write(*,'(a,1X,I6,F14.6, 5(1X,E14.6))') &
                      'TNER:',istep,time,tbmax,tner
              end if
           end if
           if (iget_cfy.eq.1) then
              ! now only at the center plane
              !write(iocfy) time,Deltat,unerc,ucmax(1:3),xoffb,xofft,ucmax(0),uner2c
           end if
           ! (rev.547) now umax is grobal and add uprim_max
        end if
     endif

     !if ((uner(2).lt.1.d-13).and.(uner(3).lt.1.d-13)) then
     if (ISNAN(uner(1))) then
        if (myid.eq.0) write(*,*) ' NaN: crashed: '
        istop=1 !(rev.1295)
     end if
     if (uner(2)/Lz.lt.1.d-5) then ! rev.1459, laminarization criterion is vrms/(s*Lz)<1.d-5
     !if (((Deltat.gt.0.2d0).or.(uner(2).lt.1.d-8)).and.(istep.gt.10)) then
        if (myid.eq.0) write(*,*) ' Laminarized '
        !if (myid.eq.0) write(39,*) ' Laminarized '
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        istop=1 !(rev.867)
        !stop
        !call MPI_FINALIZE(ierr)
        !if (time.gt.etime) stop ! for arnoldi or newton algorithm
     endif

     ihist  = 0

  endif

end subroutine histt

