program disk3d
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
  use ctes
  use running
  use bcs
  use LES
  use statistics,only: stime0
  use temp
  use hdf5
  use omp_lib

  implicit none
   include "mpif.h"

  integer i, j, idum, iniflow
  real(8),dimension(:),allocatable :: vor, phi
  real(8),dimension(:),allocatable :: u00, w00
  real(8),dimension(:),allocatable :: hv, hg, vorwk, phiwk, dvordy, chwk

  ! -- set buffers for LES --- (txy00, tzy00 are in mod.f90)
  real(8),dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2

  ! -- set buffers for itemperature --
  real(8),dimension(:),allocatable :: tb
  real(8),dimension(:),allocatable :: tb00
  real(8),dimension(:),allocatable:: ht, rhtc2, tbwk

  real(8) randu,phasex,phasez,phasey,fac
  integer istat(MPI_STATUS_SIZE),ierr,iremove_singularity,provided
  ! --- hdf5 ---
  integer h5err

  ! 
  !call MPI_INIT(ierr)

  call h5open_f(h5err)
  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,provided,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)

  ! version (revision number of Subversion)
  if(myid.eq.0) then
     write(*,*) 'HSF by A.S., Version: rev.',irev
     !write(*,*) ' optimizing hepta_solver, fftw (rfti fixed), new_change_ieor or alltoallv by siwei, two-buffer, safe pack0, new screen-output, checking with rdt, checking with isotropic flow, fixed ener, ihalf read, rogallo_ini (rand fix), idump (upo mode), umax_grobal, iadd_mode3 (moving-frame fixed), Ly = dat(3) [use hre2.dat],ihalfy (fixed),idump2 mod, ener_o3p reverted, sta4, Pk fixed, getini_rdt scaled, getini_rdt2(cft),fixe_tophat, iadd_mode=5(imp),iadd_sym=5,6(allgatherv),set_xyz_dy,iadd_sym=7(UPO4),fixed addsym,addsym8, istop=1 (laminarized), addsym chwk, vbulk, diffy 2nd(M=3), (i4.4) escru, dimension style, buff1,2 are reverted, imposed set_conj(), iadd_mode=-1 testing, igetfil_convert, idump2new, fixed interp_y, LES, hre3.dat, itemperature'
  endif
  ! (should be corresponding to svn revision)
  !
  call set_options()
  !--------------- initializes commons and things
  !call initcr(1) ! only for restart
  call initcr(0)
  !
  ! call set_options_post_initcr() ! may be need in future...see channel flow code
  !
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;
  if (itemperature.eq.1) then
     allocate(tb(buffsize), tb00(my) )
     tb=0.d0; tb00=0.d0;
  end if
  !
  if (readflag.eq.0) then
     ! generate i.c.
     if(myid.eq.0) then
        write(*,*) 'Please choose the type of i.c. to generate ...'
        write(*,*) '   [1, spectral (original wall version); 2, spectral (periodic in y); '
        write(*,*) '    3, Taylor-Green;   4, TG + some modes, i=2,3;  '
        write(*,*) '    5, streak profile; 6, RDT (input a wavemode) '
        write(*,*) '    7, tophat (16<k<32); 8, Kida-Tanaka (not yet) '
        read(*,*) iniflow
        write(*,*) ' iniflow =', iniflow
     end if
     call MPI_BCAST(iniflow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     !
     if (iniflow.eq.1) then
        !-- original wall version (slowest to become turbulence)
        call getini_original(vor,phi,u00,w00,time)
     elseif (iniflow.eq.2) then
        !-- periodic in y (not checked,
        ! but faster to be turbulence than original
        !               => it was a bug. (fixed rev.498)
        call getini(vor,phi,u00,w00,time)
     elseif (iniflow.eq.3) then
        !-- Taylor-Green (not checked, slightly shifted in z-dir??)
        call getini_TG(vor,phi,u00,w00,time)
     elseif (iniflow.eq.4) then
        !-- sekimoto (TG + some modes, i=2,3)
        call getini_seki(vor,phi,u00,w00,time)
     elseif (iniflow.eq.5) then
        !-- streak profile: U=s*y+du*cos(gam*z)
        call getini_streak(vor,phi,u00,w00,time)
     elseif (iniflow.eq.6) then
        !-- check with rdt
        !call getini_rdt(vor,phi,u00,w00,time)
        ! using rft, buggy for (k0y.ne.0)
        call getini_rdt2(vor,phi,u00,w00,time) ! (rev.673)
        !etime=Lx/Ly*10.d0 ! rev.1164
     elseif (iniflow.eq.7) then
        !-- check with rdt
        ! iopt = 1 ekin=[16:32]
        call getini_tophat(vor,phi,u00,w00,chwk,time,1)
     end if
     ! add ramdum disturbance
     !idum=-123456
     !do i=1,buffsize
     !   vor(i)=vor(i)+randu(idum)*(1.d-6)
     !end do
     !do i=1,buffsize
     !   phi(i)=phi(i)+randu(idum)*(1.d-6)
     !end do
  else

     !call getfil_convert(vor,phi,u00,w00,chwk,time)  ! read i.c in real*4
     !call getfil_yshift(vor,phi,u00,w00,chwk,time)  ! test half-Ly shift in y
     if (iread_hdf5.eq.1) then
        call getfil_t_hdf5(vor,phi,tb,u00,w00,tb00,chwk,time)  ! read i.c.
     else
        call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
     end if
        !
     if ((itemperature.eq.1).and.(init_temp.eq.1)) then
        ! add only isotropic noise
        if (myid.eq.0) write(*,*)  '     add isotropic temperature noise of uprim', uprim
        if (uprim.gt.0.5d0) then
           if (myid.eq.0) write(*,*)  '  set smaller uprim ~ 0.1, stop:'
           stop
        end if
        call getini_streak_t(tb,tb00,time) ! for pure-periodic b.c.
        !call get_kinetic_energy(vor2,phi2,u002,w002,chwk,ke_iso)
        !if (myid.eq.0) write(*,*) ' ke_iso ', ke_iso
        ! tophat of R&M does not work as a randum noise.
        !call getini_TG(vor2,phi2,u002,w002,chwk,time2) ! for pure-periodic b.c.
        fac = (xoffb_ini/Lx - int(xoffb_ini/Lx) )*Lx/Ly
        if (myid.eq.0) write(*,*)  '      the noise is distorted to satisfy B.C.', xoffb_ini,fac
        call moving_frame_t(tb,fac,-1)

        fac = uprim
        if (myid.eq.0) write(*,*)  ' ==== add iso-disturbance on temperature ==== '
        tb = fac*tb;
        tb00 =  fac*tb00;

     end if
     !call window_y(vor,phi,u00,w00)
     ! windowing like Gibson Brandt (2014,JFM)
     ! rev.1190
  endif

  time_ini=time ! (the time should be fixed to be consistent with xofft/(s*Ly)
                ! or recover by using xofft, see the following idump_mode == 2 )
  stime0=time_ini ! fixed at rev.1232

  ! set the end time as one shear-period
  !etime=time_ini + 1.d0*(2.d0*pi/alp)/(s*Ly)
  if (etime.lt.1.d9) then
     if(myid.eq.0) write(*,*) '  ------ etime =>  ', etime
  end if

  if ((iuse_LES.eq.1).and.(explicit.eq.0)) then
     write(*,*) 'LES is only for explicit, stop'
     stop
  end if
  !
  if (iadd_sym.ge.5) then
     if(myid.eq.0) then
        write(*,*) 'iadd_sym, set symmetry w.r.t, shiftx, shiftz',iadd_sym, sym_shiftx*Lx, sym_shiftz*Lz
        !read(*,*) shiftx, shiftz
     end if
     call phase_shift_fou(vor,phi,sym_shiftx*Lx,sym_shiftz*Lz)
  end if

  ! phase shift y test...
  iremove_singularity=0;
  if (iremove_singularity.eq.1) then
     ! reset boundary condition for more accurate boundary condition
     !xoffb = (xoffb/Lx - int(xoffb/Lx) )*Lx ! negative
     !xofft = -xoffb ! positive
     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1)
     fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
     call phase_shifty(vor,phi,u00,w00,atan(phasey)/(2*pi)*Ly,fac)
     ! dealiasing is required for shifting in y with an artibrary shift-periodic boundary condition...
     call apply_filtering(vor,phi,u00,w00,0) ! now iopt = 0 is a redundunt option.
     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1)
     call phase_shift_fou(vor,phi,atan(phasex)/(2*pi)*Lx,atan(phasez)/(2*pi)*Lz )

     call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,1) ! for checking zero phase
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  endif


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !
  if (idump_mode.eq.1) then
     ! for initialting isotropic flow, dumping the fields with pure-periodicity
     if(myid.eq.0) write(*,*) '  ------ dump mode for initiating isotropic --- '
     dtimag = timep_sh/dump2tint ! updated at rev.1711
     dumptime = time_ini + dtimag ! the first dump time
     !dumptime = xofft/(abs(s)*Ly) + dtimag ! fixed (rev.545)
     if(myid.eq.0) write(*,*) '  ------    dump interval, dtimag =>  ', dtimag
     nimag = nstep

  elseif (idump_mode.eq.2) then
     ! for making UPO movies, set time_interval to dump a file and etime
     !
     !
     if(myid.eq.0) write(*,*) '  ------ UPO dump mode  ----- '
     if(myid.eq.0) write(*,*) '  ------    time_ini,  =>  ', time_ini
     dtimag = timep_sh/dump2tint  ! check the time intervals and
                             ! dt should be considered for statistics...
     ! the first dump time is the next pure-periodic
     !dumptime = int(time/timep_sh)*timep_sh + timep_sh ! for UPO
     !time = (xofft/abs(s)/Ly) ! reset time stamp !
     if(myid.eq.0) write(*,*) '  ------ the time stamp is adjusted with xofft ==>  ', time
     time_ini = time
     dumptime = int(time/timep_sh)*timep_sh + timep_sh ! fixed (rev.545)
     if (abs(int(time/timep_sh)*timep_sh - time_ini).lt.1.e-9) dumptime=time+dtimag

     if(myid.eq.0) write(*,*) '  ------ the first dump time (pure-periodic),  =>  ',dumptime

     if(myid.eq.0) write(*,*) '  ------ dump interval, dtimag =>  ',dtimag
     etime = dumptime + dump2timep*timep_sh
     if(myid.eq.0) write(*,*) '  ------ etime =>', etime

     nimag = nstep

  end if

  if (iadd_damping.ne.0) then
     if(myid.eq.0) then
        write(*,*) ' iadd_damping mode, set up and aa',damp_up,damp_aa
     end if
  end if
  !if(myid.eq.0) write(*,*) '... file read'

  !  -------  allocate rest of buffers
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0;
  allocate( phiwk(buffsize),vorwk(buffsize),dvordy(buffsize) )
  phiwk = 0.d0; vorwk = 0.d0; dvordy = 0.d0;
  if (itemperature.eq.1) then
     allocate(ht(buffsize), tbwk(buffsize))
     ht=0.d0; tbwk=0.d0
  end if
  ! --- the time loop
  !if(myid.eq.0) write(*,*) myid,'just before cross ...'
  !call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,chwk)
  !
  ! allocate redundant buffers for DNS/LES
  allocate(rhvc1(buffsize),rhvc2(buffsize))
  allocate(rhgc1(buffsize),rhgc2(buffsize))
  rhvc1= 0.d0; rhvc2= 0.d0
  rhgc1= 0.d0; rhgc2= 0.d0
  if (itemperature.eq.1) then
     allocate(rhtc2(buffsize))
     rhtc2= 0.d0
  end if
  if (abs(se+s).lt.tiny) then
     if(myid.eq.0) write(*,*) 'reverse-mode, check s and se:',s,se
     vor=-vor; phi=-phi; u00=-u00; w00=-w00;
  end if
  if (itemperature.eq.1) then
     if(myid.eq.0) write(*,*) ' crosst '
     call crosst(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,&
          rhvc1,rhvc2,rhgc1,rhgc2, &
          tb,tb00,ht,tbwk,rhtc2,chwk)
  else
     call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,&
          rhvc1,rhvc2,rhgc1,rhgc2,chwk)
  end if
  if (itemperature.eq.1) then
     deallocate(ht,rhtc2)
  end if
  deallocate(hv,hg)
  deallocate(phiwk,vorwk,dvordy)
  deallocate(rhvc1,rhvc2,rhgc1,rhgc2)

  ! ------------------------------------------------------------------
  ! 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! 
  call h5close_f(h5err)
  call MPI_FINALIZE(ierr)

end program disk3d

