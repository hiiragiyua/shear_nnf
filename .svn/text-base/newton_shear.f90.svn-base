program newton

  use ctes
  use running
  use bcs
  use gmres
  use eig
  use precond, only:nband
  use LES

  implicit none
  include "mpif.h"

  integer i,j,k,kk,idum,info,ilin,inewton
  real(8) dp,fac,mtime,eps_j,err,res,res0,delta_backup,sxfac,szfac
  real(8),dimension(:),allocatable :: vor, phi, chwk
  real(8),dimension(:),allocatable :: u00, w00

  character*80 fname
  character*4 ext1

  integer iosol,itr,kth,knum,iloop
  logical lex
  integer iloop_stop,iconverged_newton,ihave_small_res

  integer istat(MPI_STATUS_SIZE),ierr

  ! /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)
  !
  ! version (revision number of Subversion)
  if(myid.eq.0) write(*,*) ' %--- Newton-Krylov-hookstep method for HSF by A.S., Version: rev.',irev
  !if(myid.eq.0) write(*,*) '(compute dxdt for all isol_modes, dump stat in make_mat), arclength, dre1.e-6, hre2.dat, relative in xz (isol4), isol5_6, chwk ini, ihave_small_res (fixed), close(39), reset delta, sta4, Pk, iadd_sym=5,6v (shifts_ini*[Lx Lz]),pythag,DGESVD (lwk fixed, initialization), hs_cut,cpara,set_xyz_dy, addsym fixed, sta4_fixed, addsym chwk, simply updated hookstep,ishifty, fixed dxdt, vbulk, xforce, non-zero bulk velocity,right preconditionig, prederive2 M=3, dfp=> dfp - dxp for precondition, delta is for P*dxp,ffty for u00,w00, bugfixed update_premat, fixed rfty, factorization (nband=1),itrmax=999, mpc_vorphi (vor = i*phi), dmu=0 delta=0.9*qnorm, update_hookstep is reverted(rev961), prederiv2(M=3),nband=2, iuse_hej,iuse_v,fixed xwkt, hscut in imuloop, (i4.4), single SVD, ifilter, x0mode-conj, bug fixed single solv_mu, extshiftx with vbulk (is commented), ishifty with isol3, xforce=-s*vbulk, non-integer timep, ishifty_rev, bug-fixed set_conj, timep_sh => Ly (modified), icTp sca, rescale for different gam, buff1,2 are reverted to imporse set_conj(), rev961 is reverted and merged, sca_vor sca_phi (v) is introduced for a normalization by SLz and S*sqrt(Rez),arc with iadd_damping, iuse_LES (hre3.dat), cpara128, QRorth, iuse_us, arc for Cles, meme_request_sum, arclength for alp testing (updating b.c. correctly)'
  ! (should correspond to svn revision) 
  ! --- initializes commons and things
  call initcr(0)
  !
  call set_options()
  ! explicit = 0 ! 0, implicit; 1, explicit. ==> set in hre3.dat
  !
  call set_newton_options()
  !
  iosol=98
  !  
  if(myid.eq.0) then
     write(*,*) 'set newton parameter in hre3.dat'
  end if
  !call  MPI_BCAST(para_newton,20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  !
  delta_backup=delta_ini
  !
  iget_hv=0 ! run arnoldi.f90 to get accurate eigenvalues...
  !
  if (isol_mode.eq.1) then
     if (myid.eq.0) write(*,*) & 
          &'FXP mode: Tp is now ignored'
  elseif (isol_mode.ge.3) then
     if (abs(s).lt.tiny) then
        write(*,*) 'error for timep_sh'
        stop
     endif
     !timep_sh = 2.d0*pi/alp/Ly*int(timep) ! fixed (rev.600, assuming s>0)
     timep_sh = 2.d0*pi/alp/Ly*timep ! rev.1123
     if (myid.eq.0) write(*,*) & 
          &'UPO mode: set the integration time as a multiple of the shear-period: ', & 
          &'Tp=', timep_sh
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then 
        if (myid.eq.0) write(*,*) & 
             & '  relative UPO: shiftx0, shiftz0 = ', shiftx0,shiftz0
     end if
  end if

  if (ishifty.eq.1) then     
     shifty = shifty0;
     if ( (isol_mode.eq.3).and.(iset_bulk0.eq.1) ) then
        if (myid.eq.0) write(*,*) 'option error, iset_bluk0 should be off for shifty and isol_mode=3'
        stop
     endif
     ! shifting in y-dir
     !if (isol_mode.eq.3) then
     !   !vbulk0 = -shifty0/timep_sh
     !   !vbulk = vbulk0
     !endif
     !if (myid.eq.0) write(*,*) 'ishifty=1: set initial vbulk=-sy/Tp:',vbulk
     if (myid.eq.0) write(*,*) 'ishifty=1: set initial shifty0:',shifty
  end if

  if(myid.eq.0) then 
     write(*,*) ' linear solver parameters: k_gmres', k_gmres
     write(*,'(a,3(1X,E14.6))') '  eps_gmres, eps_jacobi, eps_newton', eps_gmresk, eps_jacobi, eps_newton
  end if
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;

  !if(myid.eq.0) write(*,*) 'go to read file ...'
  if (readflag.eq.0) then
     ! generate i.c.
     !
     !-- original wall version (slowest to be turbulence)
     ! call getini_original(vor,phi,u00,w00,time)
     !
     !-- periodic in y (not checked, but faster to be turbulence than original)
     ! call getini(vor,phi,u00,w00,time)
     !
     !-- Taylor-Green (not checked)     
     ! call getini_TG(vor,phi,u00,w00,time)       
     !-- sekimoto (TG + some modes, i=2,3)
     ! call getini_seki(vor,phi,u00,w00,time)
     !-- streak profile: U=s*y+du*cos(gam*z)
     call getini_streak(vor,phi,u00,w00,time)
  else
     call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.

     !call window_y(vor,phi,u00,w00) 
     ! ==> now buggy since xofft is reset 
     ! ==> bug fixed ! set b.c. for getv (rev1395)
     !     ref) windowing like Gibson Brandt (2014,JFM)
     ! rev.1188
  endif

  
  if(myid.eq.0) write(*,*) '... file read'

  if (iadd_sym.ge.5) then

     !call  MPI_BCAST(sym_shiftx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     !call  MPI_BCAST(sym_shiftz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     write(*,'(i5,a,i5,2(1X,F14.6))') myid,': debug, iadd_sym, phase shift',iadd_sym, sym_shiftx*Lx, sym_shiftz*Lz 
     call phase_shift_fou(vor,phi,sym_shiftx*Lx,sym_shiftz*Lz)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end if

  if (iremove_singularity.eq.1) then
     ! reset boundary condition for more accurate boundary condition
     !xoffb = (xoffb/Lx - int(xoffb/Lx) )*Lx ! negative
     !xofft = -xoffb ! positive
     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1)
     fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
     call phase_shifty(vor,phi,u00,w00,(atan(phasey)/(2*pi))*Ly,fac)
     ! rev.1564, do not apply filtering ...
     !call apply_filtering(vor,phi,u00,w00,0) ! now iopt = 0 is a redundunt option.
     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1)     
     sxfac=atan(phasex)/(2d0*pi); szfac=atan(phasez)/(2d0*pi)     
     call phase_shift_fou(vor,phi,sxfac*Lx,szfac*Lz )

     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1) ! for checking zero phase 
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
     ! The extra shift in x is required due to the mean velocity U=Sy.
     ! sy= np.arctan(phasey)/(2*np.pi)*Ly
     ! sx= -sy*(Lx/Ly)*mTp
  endif
   
  ! --- setting arrays and vectors for dynamical system ---
  ! --- and allocated gmres vectors, now only by master
  ! --- allocate vv hh etc... ! 
  call gmres_ini(vor,phi,u00,w00) 
  !if(myid.eq.0) write(*,*) '... gmres allocated '
  call newton_ini(vor,phi,u00,w00,chwk)
  !
  !write(*,*) myid, 'throwing ...'
  xofft0=xofft
  xoffb0=xoffb; 
  call throwing(vor,phi,u00,w00,chwk,x0,1)

  !write(*,*) myid, 'throwing done'
  !
  if (abs(iarclength).eq.1) then
        
     ! read xm1; vor,phi,u00,w00 are used
     ! re Ly alp, etc are overwritten from the header (or footer) 
     call arclength_ini(vor,phi,u00,w00,chwk)

  end if
  ! set unknown parameters to x0(icre), etc. Note calling after arclength setting.  
  call throwing_parameter(x0,0) ! set extended unknowns for x0(icsx,icre,...)

  iloop_stop=0 ! stop option
  !
  do iloop=1,nloopmax
     ! reinitialization...
     vor=0.d0; phi=0.d0; u00=0.d0; w00=0.d0;  ! reverted rev1148
     chwk=0.d0;
     !
     vv=0.d0; hh=0.d0; gh=0.d0; rr=0.d0;
     !
     bb=0.d0; dxp=0.d0; Ax=0.d0; 
     cs=0.d0; ss=0.d0; ee=0.d0;
     yg=0.d0; 
     if (myid.eq.0) write(*,*)  ' ==== NSTEP mode, iloop =',iloop,'/',nloopmax, ' ==== '
     
     if (iarclength.eq.-1) then 
        if (myid.eq.0) then
           write(*,*) ' (testing) iarclength=-1: simple auto-modification for cpara=',trim(cpara)

           !if (trim(cpara).eq.'Rez') write(*,*) ' increase Re (use xstan for initial guess, but fixed Rez)'
           !if (trim(cpara).eq.'Ayz') write(*,*) ' increase Ly (use xstan for initial guess, but fixed Ayz)'
           !if (trim(cpara).eq.'Axz') write(*,*) ' increase Lx (use xstan for initial guess, but fixed Axz)'

        end if
     end if

     ! /*  set xf and reset vor and phi  */
     !if (isol_mode.ge.3) then ! finding orbits
     !   ! this is also used to accumulate statistiscs by nohist=0, before calling make_mat
     !   nohist=0 ! screen output and dump *.cf, accumulate stat, but not dump stat. file
     !   call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     !   nohist=1 ! reset the history output option
     !   if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
     !      call phase_shift_fou(vor,phi,shiftx,shiftz)
     !      !if (ishifty.eq.1) then
     !      !   extshiftx=s*0.5d0*timep_sh*shifty
     !      !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
     !      !endif
     !   end if
     !   call throwing(vor,phi,u00,w00,chwk,xf ,1) ! --> xf
        
     !end if
     !
     if (abs(iarclength).eq.1) then

        if (myid.eq.0) write(*,*) ' update arclength'
        if (iloop.eq.1) then
           arcl=arcl_ini
           ! rev.1467
           if (myid.eq.0) write(*,*) ' arclength parameter is adjusted only at the first trial:',arcl
        end if
        call set_arclength(vor,phi,u00,w00,chwk) ! set x0 vector 
        if (iloop.eq.1) then
           arcl=1.d0
           if (myid.eq.0) write(*,*) ' arclength parameter is recovered to be 1'
        end if

     endif !/* arclength */

     ! /* re set x0 --> vor, phi */
     time = time_ini
     !xoffb = xoffb_ini
     !xofft = xofft_ini
     xoffb=xoffb0 ! rev.1390
     xofft=xofft0
     call backthrowing(vor,phi,u00,w00,chwk,x0,1) 

     if (myid.eq.0) then 
        !write(*,*)  ' NSTEP-mode  re =', re, ' Ly=',Ly, ' alp=', alp
        if(iarclength.eq.1) write(*,*)  'NSTEP-mode, cpara = ',trim(cpara),', x0(icre) = ', x0(icre)
        ! temporarily icre significates the column of arclength parameter...
        write(*,'(a,3(1X,F14.6))')  ' NSTEP (Axz,Ayz,Rez) =', gam/alp, Ly/Lz, re*Lz*Lz
     end if
        !
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call make_mat(vor,phi,u00,w00,chwk) ! write *.cf (nohist=0)
     !     
     ! --- the newton process ---
     res=bnorm/x0norm; res0 = res;
     if(myid.eq.0) write(*,*) '==========================================='
     if(myid.eq.0) write(*,'(a,i5,4(E18.10))') ' NEWTON residual:',0,bnorm,res,x0norm,xfnorm
     if(myid.eq.0) write(*,*) '==========================================='

     !itrmax=999 ! set in hre3.dat
     iconverged_newton=0
     ihave_small_res = 0

     do inewton = 1, itrmax

        if (isnan(res)) then
           write(*,*) 'NaN is detected, stop'
           stop
        endif

        eps_j=eps_jacobi*x0norm

        if ((kth.le.1).and.(delta.lt.delta_ini*0.01)) then
           !
           ! if the previous GMRES step was only 1 and delta is shrinked too much, then
           ! the backuped field is dumped...  restart from it using small delta 
           if (myid.eq.0) write(*,*) 'recover previous delta'
           delta=delta_backup
        end if
        ! reset hook step flags
        ihookstep_equals_newtonstep=0 ! reset flag
        ihavebackup=0 
        x0_backup = x0 ! backup the present field before starting newton
        xf_backup = xf ! backup the present field before starting newton
        ! save the boundary condition (rev.1390)
        xofft0_backup = xofft0;  xoffb0_backup = xoffb0
        xofftf_backup = xofftf;  xoffbf_backup = xoffbf
        !
        delta_backup = delta
        ! solving: A*dxp = bb
        call gmres_k_shear(vor,phi,u00,w00,chwk,itr,err,kth,eps_j)
        !/* save the temporal GMRes solution (dxp, newton step) */
        !call save_vec(vor,phi,u00,w00,chwk,dxp) ! for debugging
        ! call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
        ! vor,phi is now throw-backed from x0
        !
        ! /* update the vector x0(u,v,w) = x0 + dxp                 */
        ! /* or update delta and try hookstep                       */
        !
        if (iuse_hookstep.eq.1) then
           if(myid.eq.0) write(*,*) 'updating newton by hookstep' 
           call update_hookstep(vor,phi,u00,w00,chwk,kth,eps_j)
        else
           if(myid.eq.0) write(*,*) 'updating newton by damped newton'
           call update_newton(vor,phi,u00,w00,chwk)
        end if
        !
        ! /* update xf, dxdt,b and norms */
        !    compute statistics and energy plots 
        !              ... then it will be dumped by save_vec() later
        call make_mat(vor,phi,u00,w00,chwk)

        ! /* before interupting the newton-loop iteration, save the present file */
        if (myid.eq.0) then
           ! stop option usage: touch ./pool/upo.stop  etc...  
           fname=filout(1:index(filout,' ')-1)//'.stop'
           inquire(file=fname,exist=lex)
        end if
        call MPI_BCAST(lex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(lex) then
           if (iskip_save.eq.1) then
              if (myid.eq.0) then
                 write(*,*) ' ******************************************************'
                 write(*,*) ' Stop-option is on ==> start dumping files'
                 write(*,*) ' ******************************************************'
              end if
              iskip_save=0
           end if
        end if
        !
        if (iskip_save.eq.0) then
           ! /* SKIP save the updated field  */           
           ! reset boundary condition
           time = time_ini
           !xoffb = xoffb_ini 
           !xofft = xofft_ini
           xoffb=xoffb0 ! rev.1390
           xofft=xofft0
           call save_vec(vor,phi,u00,w00,chwk,x0) ! vor,phi is now throw-backed from x0
        else
           time = time_ini           
           !xoffb = xoffb_ini 
           !xofft = xofft_ini
           xoffb=xoffb0 ! rev.1390
           xofft=xofft0
           call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
        end if
        !
        ! this update_premat should be done after make_mat, 
        ! otherwise it dumps statistics in premat..
        if (iuse_mpc.eq.1) call update_premat(dxp,dfp,chwk,1,1) 
        ! dfp is r.h.s of Newton-eq
        ! /* convergence check of newton_gmres */
        res=bnorm/x0norm 
        if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
           if (myid.eq.0) write(*,*) 'NEWTON shiftx,z: ',inewton,shiftx,shiftz 
        end if
        if (ishifty.eq.1) then
           if (myid.eq.0) write(*,*) 'NEWTON shifty: ',inewton,shifty
        end if
        if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
           if (myid.eq.0) write(*,*) 'NEWTON Ly, Tp: ',inewton, Ly, timep_sh
        end if
        if(iremove_singularity.eq.1) then
           call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1) 
           ! for checking if phases are fixed in x-, y-, z-direction 
        endif
        if(myid.eq.0) write(*,*) '==========================================='
        if(myid.eq.0) write(*,'(a,i5,4e18.10)') ' NEWTON residual:',inewton,bnorm,res,x0norm,xfnorm
        if(myid.eq.0) write(*,*) '==========================================='

        if(lex) then
           if (iskip_save.eq.0) then
              if (myid.eq.0) then
                 write(*,*) ' ******************************************************'
                 write(*,*) ' Stop-option is on ==> force to stop newton iterations '
                 write(*,*) ' Delete the file to restart: ',trim(fname)
                 write(*,*) ' ******************************************************'
              end if
              stop
           end if
        end if

        !if (iuse_mpc.eq.1) then
           ! ---- checking eigenvalues ----
           !if (myid.eq.0) call get_eigen(kth,1,'n')
        !end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !if ((inewton.gt.1).and.(res.lt.eps_newton)) then
        if (res.lt.eps_newton) then
           if (myid.eq.0) write(*,*) 'NEWTON method has converged: ',res
           iconverged_newton=1
           exit
        elseif (inewton.eq.itrmax) then
           if(myid.eq.0) write(*,*) 'NEWTON newton step reached itrmax=',itrmax
           iloop_stop = 1
        endif
        
        !if (res.lt.eps_newton*5.d0) then
        if (res.lt.eps_newton) then
           ihave_small_res=1
        end if

        if (abs(iarclength).eq.1) then
           if ((ihave_small_res.eq.1).and.(inewton.gt.10)) then
              if (myid.eq.0) then
                 write(*,*) 'NSTEP mode, more or less converged, but reduce arcl', arcl
              end if
              iconverged_newton=1
              arcl=arcl*0.8
              exit
           end if
        end if
        if (iarclength.eq.1) then
           if (iuse_LES.eq.0) then
              if ((res.gt.res0).and.(inewton.gt.1)) then
                 if (myid.eq.0) then
                    write(*,*) 'NSTEP mode, arclength assuming monotonically-decrease of residual, STOP newton'
                    write(*,*) 'NSTEP mode, the previous fields is restored, then restart by using reduced arl:',arcl
                 end if
                 
                 if (res0.lt.1.d-4) then
                    if (myid.eq.0) write(*,*) 'NSTEP mode, the error is small ... continue using x0_backup ...'
                    !BUG fixed (rev.650), re of backup was not stored...in rev.649
                    ihave_small_res = 1
                    !iloop_stop = 0
                    iloop_stop = 1
                 else
                    if (myid.eq.0) write(*,*) 'NSTEP mode, the error is big STOP ...' 
                    iloop_stop = 1 ! or continue with small arcl                
                 end if
                 ! restore previous fields ==> x0, re, Ly to save
                 x0 = x0_backup
                 xf = xf_backup
                 ! restore boundary condition (rev.1390)
                 xofft0 = xofft0_backup;  xoffb0 = xoffb0_backup
                 xofftf = xofftf_backup;  xoffbf = xoffbf_backup
                 call backthrowing_parameter(x0,0)
                 exit ! exit newton loop
              end if
           end if
        else 
           if ((res.gt.res0).and.(inewton.gt.1)) then
              if (myid.eq.0) write(*,*) ' normal mode: residual bumped-up: ', & 
                   'restore previous fields and save temporarily, and continue with smaller delta'
              ! restore previous fields ==> x0, re to save
              x0 = x0_backup 
              xf = xf_backup
              ! restore boundary condition (rev.1390)
              xofft0 = xofft0_backup;  xoffb0 = xoffb0_backup
              xofftf = xofftf_backup;  xoffbf = xoffbf_backup
              call backthrowing_parameter(x0,0) 
              if (iskip_save.eq.1) then
                 if (myid.eq.0) write(*,'(a,3(1X,F14.6))')  & 
                      ' save_vec, (Axz,Ayz,Rez) =', gam/alp, Ly/Lz, re*Lz*Lz 
                 ! reset boundary condition
                 time = time_ini
                 !xoffb = xoffb_ini 
                 !xofft = xofft_ini
                 xoffb = xoffb0 ! (rev.1390)
                 xofft = xofft0
                 call save_vec(vor,phi,u00,w00,chwk,x0)
                 if (myid.eq.0) then
                    write(ext1,'(i4.4)') id22-1 
                    fname=filout(1:index(filout,' ')-1)//'.'//ext1//'_is_tmp'
                    open(iosol,file=fname,status='unknown')
                    write(iosol,'(a,3(1X,F18.10))') '% ', gam/alp, Ly/Lz, re*Lz*Lz
                    write(iosol,'(a,i5,4(1X,F18.10))') ' NEWTON residual:',inewton,bnorm,res,x0norm,xfnorm
                    close(iosol)
                 endif
              else
                 time = time_ini
                 !xoffb = xoffb_ini 
                 !xofft = xofft_ini
                 xoffb = xoffb0 ! (rev.1390)
                 xofft = xofft0
                 call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
              end if
              ! bug fixed ... forgotten to update make_mat (2014/04/08, rev.904)
              call make_mat(vor,phi,u00,w00,chwk)

              if (iuse_mpc.eq.1) call update_premat(dxp,dfp-dxp,chwk,1,1) ! dfp is r.h.s of Newton-eq
              ! set a large nhook... (rev.1113)
              !delta = delta_backup*0.1
              !deltaMax = delta*10.0
              !if (myid.eq.0) write(*,*) '    shrink delta =',delta
           end if
        end if

        res0 = res 
        ! store present residual        

     end do
     

     ! save the updated fields at the end ... 
     if (ihave_small_res.eq.1) then

        if (myid.eq.0) write(*,*) ' ihave_small_res', res0
        if (myid.eq.0) write(*,*) '    compute .cf and sta3, then continue it'
        time = time_ini
        xoffb = xoffb0 ! (rev.1395)
        xofft = xofft0
        call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
        nohist=0 ! screen output and dump *.cf and stat. file
        !nohist=1
        call calc_dns_disk3d(vor,phi,u00,w00,chwk)
        nohist=1 ! reset the history output option
        iloop_stop = 0 ! continue

     end if

     if (abs(iarclength).eq.1) then 
        if ((iconverged_newton.eq.1).or.(ihave_small_res.eq.1)) then
           if (myid.eq.0) then 
              write(*,'(a,3(1X,F14.6))')  & 
               &  ' NSTEP save_vec, (Axz,Ayz,Rez) =', gam/alp, Ly/Lz, re*Lz*Lz 
              if(iarclength.eq.1) write(*,*)  'NSTEP save cpara = ',trim(cpara), & 
               &  ', x0(icre) = ', x0(icre)
           end if
           ! reset boundary condition
           time = time_ini
           !xoffb = xoffb_ini
           !xofft = xofft_ini
           xoffb = xoffb0 ! (rev.1390)
           xofft = xofft0
           if (iskip_save.eq.1) then
              call save_vec(vor,phi,u00,w00,chwk,x0)
              ! do not save the last one if iskip_save==0, since it saves at each
              ! vor ... are now the same with x0, in save_vec
           else
              call backthrowing(vor,phi,u00,w00,chwk,x0,1) ! fixed at rev.1611
           end if

           if (inewton.lt.5) then
              if ((iuse_LES.eq.1).and.(arcl.gt.5.d0)) then
                 if (inewton.eq.1) then
                    arcl = arcl * 1.2
                 endif
              else
                 ! increase arcl, because it is the good convergence behaviour
                 arcl = arcl * 1.2
                 !if (arcl.gt.5.0) arcl = 5.0 ! suggest max(arcl)=5.0
                 if (myid.eq.0) write(*,*) 'good convergence, increase arcl =>', arcl
              endif
           elseif (inewton.gt.10) then
              if (iuse_LES.eq.1) then
                 !rev.1348
              else
                 arcl = arcl * 0.8
              end if
              if (myid.eq.0) write(*,*) 'slow convergence, decrease arcl =>', arcl
           elseif (inewton.gt.20) then
              if (iuse_LES.eq.1) then
                 arcl = arcl * 0.8 ! rev.1348
              else
                 arcl = arcl * 0.5
              end if
              if (myid.eq.0) write(*,*) 'bad convergence, decrease arcl =>', arcl
           endif
        end if
     else
        ! reset boundary condition
        time = time_ini
        !xoffb = xoffb_ini 
        !xofft = xofft_ini
        xoffb = xoffb0 ! rev.1390)
        xofft = xofft0
        if (iskip_save.eq.1) then 
           call save_vec(vor,phi,u00,w00,chwk,x0) 
           ! do not save the last one if iskip_save==0, since it saves at each
        else
           call backthrowing(vor,phi,u00,w00,chwk,x0,1)
        end if
        ! newton iteration
     endif
       
     if (iloop_stop.ge.1) then
        if (myid.eq.0) write(*,*)  'NSTEP mode STOP!', iloop
        exit 
     else ! continue loop
        ! reset delta
        delta = delta_ini
        if(myid.eq.0) write(*,*) myid,'hookstep, re-set initial delta:', delta
     end if

     ! ---- checking eigenvalues of the pure Jacobian of Newton equation ----
     !      Note: this is different with the linear stability as eig(dudt) and 
     !            here includes eigenvalues for arclength parameters.
     write(*,*) myid,'get_eigen' 
     call get_eigen(kth,999,'v')
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  enddo ! end of NSTEP loop 

  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program newton
