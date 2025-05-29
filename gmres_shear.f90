subroutine set_newton_options()
  
  ! see set_options() [initcr.f90]

  use ctes
  use running
  use gmres
  use precond, only:nband
  use LES

  implicit none
  include "mpif.h"

  character*128 command
  integer nextline
  integer istat(MPI_STATUS_SIZE),ierr

  ! DO NOT forget to BCAST
  iuse_newton=1; ! for scaled reading off...see readwrite.f90
  ! --------------------------------------------------------------------------
  ! set default values for newton and arclength method...
  iskip_save=0;
  isol_mode=3; timep=1.d0; shiftx0=0.d0; shifty0=0.d0; shiftz0=0.d0
  k_gmres=50; 
  eps_gmresk=1.d-3; eps_jacobi=1.d-6; eps_newton=1.d-6; 
  delta_ini=1.d0; hs_cut=1.d-26;
  nloopmax=1 ! the maximum trials of newton search using arclength ...
  itrmax=999 ! the maximum newton loop
  iarclength=0; arcl=1.d0; arcl_ini=1.d0; 
  ! ----- ! the default testing options !! ----- !
  ishifty=0; ! testing
  iset_bulk0=0; ! set this option on for isol_mode=4,6 (relative periodic orbit)
  iuse_mpc=0; nband=0; ! precondition, only nband=2 produce non-conj for kx=0.
  !
  iremove_singularity=0
  iuse_scaling=0  ! this mode does not work well, (rev1137)
  ifilter=0; ! filtering (zero-precondition)
  !
  iuse_compensatory_norm=0; comp_norm=2.5d0;
  iuse_QR_orth=0;
  !
  !
  iuse_ffty=0 ! check more throwing and backthrowing ... 
  !             ... now it is done for u00,w00
  iuse_us=0;

  ! ----- ! stable options ! 
  iopt_matvec=1; 
  ! eps_gmres may be changed to be 1.e-4 for central difference (iopt_matvec==2)
  iuse_hej=1; ! see Dennis & Schnabel, Sec 5.4
  iuse_v=1; ! throw v instead of lap.v
  iuse_hookstep=1
  iskip_screenout = 1
  !  
  iuse_fou_mode=2 !  
  ! this should be 2, (1 does not work, 3 is for fixing x0modes)
  ! ------------- go to a command to read options ----------------------------
  if (ihre_error.ne.0) return
  if (myid.eq.0) then
     command='isol_mode'
     !  1, a fixed point (SW); 
     !  2, Travelling wave (TW)
     !  3, periodic-orbit with Tp fixed (UPO); 
     !  4, relative periodic orbit in x-, z-dir with Tp fixed (RPO)
     !  5, periodic-orbit without fixing Tp (Ly,s,alp can be an unknown) (UPO); 
     !  6, relative periodic-orbit without fixing Tp  (RPO); 
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) isol_mode
        nextline=nextline+1
        read(hrefile(nextline),*) timep
        nextline=nextline+1
        !if (isol_mode.ge.3) then
           ! phase shift parameter (input the total shift, atan(phase)/(2*pi)*Lx )
           read(hrefile(nextline),*) shiftx0, shifty0, shiftz0
           nextline=nextline+1
        !endif
        ! dim. of Krylov subspaces
        read(hrefile(nextline),*) k_gmres
        nextline=nextline+1
        read(hrefile(nextline),*) eps_gmresk,eps_jacobi,eps_newton 
        nextline=nextline+1
        read(hrefile(nextline),*) delta_ini, hs_cut
     end if
  end if
  call  MPI_BCAST(isol_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(timep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(shiftx0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(shifty0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(shiftz0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  shiftx = shiftx0; 
  shifty = shifty0;
  shiftz = shiftz0; ! relative shift 
  call  MPI_BCAST(k_gmres,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (k_gmres.eq.0) then
     write(*,*) myid,'input format error, k_gmres',k_gmres
     stop
  end if
  !
  call  MPI_BCAST(eps_gmresk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(eps_jacobi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(eps_newton,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(delta_ini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(hs_cut,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !
  command='iskip_save'
  call switching(command,iskip_save)
  if ((iskip_save.eq.0).or.(iskip_save.eq.1)) then 
  else
     write(*,*) trim(command),': error',iskip_save
     stop
  end if
  !
  if (myid.eq.0) then
     command='iarclength'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iarclength
        nextline=nextline+1
        if (iarclength.ne.0) then
           ! for previous solution to get tangent
           ! the footer format is required for Cles and damping_up,aa continuation 
           read(hrefile(nextline),'(a)') filinm1 ! only for master
           nextline=nextline+1
           ! arclength parameter and cpara, 'Rez' 'Ayz'
           ! now testing for Ayz Axz, upr (uprim) (rev1200)
           ! iarclength=-1, => for simple modification of Rez, and so on. 
           read(hrefile(nextline),*) arcl_ini, cpara
           nextline=nextline+1
           ! num. time of trials of arclength 
           read(hrefile(nextline),*) nloopmax
        end if
     end if
  end if
  call  MPI_BCAST(iarclength,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(arcl_ini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(cpara,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(nloopmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( (nloopmax.eq.1).and.(iarclength.ne.0) ) then 
     write(*,*) 'set number of trials of arclength continuations..., stop'
     stop
  endif
  if ( iarclength.ne.0) then
     !write(*,*) myid,'iarclength =',iarclength, '; cpara=',trim(cpara)
  end if
  !
  command='itrmax'
  call switching(command,itrmax) ! default is 999
  if ((itrmax.le.0).or.(itrmax.ge.10000)) then
     write(*,*) 'itrmax error: out of 1<itrmax<10000 (4-digit extension)'
     stop
  end if
  !
  command='iopt_matvec'
  call switching(command,iopt_matvec)  ! default is 1
  if ((iopt_matvec.ne.1).and.((iopt_matvec.ne.2).and.(iopt_matvec.ne.4))) then
     write(*,*) 'iopt_matvec should be 1st- or 2nd-order or 4th-order '
     stop
  end if
  if ((iopt_matvec.eq.4).and.(iarclength.ne.0)) then
     write(*,*) 'iopt_matvec 4th-order is not implemented for dxda...stop'
     stop
  end if

  !iuse_hej=1
  command='iuse_hej'
  call switching(command,iuse_hej)  ! default is 1
  if ((iuse_hej.eq.0).or.(iuse_hej.eq.1)) then 
  else
     write(*,*) trim(command),': error',iuse_hej
     stop
  end if
  !iuse_v=1
  command='iuse_v'
  call switching(command,iuse_v) ! default is 1
  if ((iuse_v.eq.0).or.(iuse_v.eq.1)) then 
    ! write(*,*) 'iuse_v = ',iuse_v
  else
     write(*,*) trim(command),': error',iuse_v
     stop
  end if

  ! ----------  for more testing options -------------

  !iuse_us=1: throw us instead of vor, so that the unit would be the velocity
  ! vor => a part of u by converting i*kz/(kx*kx+kz*kz)*vor
  !                  note: the imaginary unit i is ignored from the scaling
  command='iuse_us'
  call switching(command,iuse_us) ! default is 0
  if ((iuse_us.eq.0).or.(iuse_us.eq.1)) then 
    ! write(*,*) 'iuse_us = ',iuse_us
  else
     write(*,*) trim(command),': error',iuse_us
     stop
  end if  
  !iuse_ffty=1
  command='iuse_ffty'
  call switching(command,iuse_ffty)
  if ((iuse_ffty.eq.0).or.(iuse_ffty.eq.1)) then     
  else
     write(*,*)  trim(command), ': error',iuse_ffty 
     stop
  end if

  !ishifty=1
  command='ishifty'
  call switching(command,ishifty)! default is 0
  if ((ishifty.eq.0).or.(ishifty.eq.1)) then     
  else
     write(*,*)  trim(command),': error',ishifty 
     stop
  end if
  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
     if(myid.eq.0) write(*,*)  'relative periodic mode, set ishifty (rev1562)= ',ishifty
     ishifty=1
  end if

  command='iset_bulk0'
  call switching(command,iset_bulk0) ! default is 0
  if (iset_bulk0.eq.1) then 
     write(*,*) 'testing use for relative periodic orbit '
  else if (iset_bulk0.ne.0) then
     write(*,*)  trim(command),': error',iset_bulk0
     stop
  end if

  if (myid.eq.0) then
     command='iuse_mpc'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iuse_mpc
        nextline=nextline+1
        if (iuse_mpc.ne.0) then
           read(hrefile(nextline),*) nband
           nextline=nextline+1
        end if
     end if
  end if
  call  MPI_BCAST(iuse_mpc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(nband,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (iuse_mpc.eq.1) then 
     if (myid.eq.0) write(*,*) ' preconditioning for newton equation,', &
          & ' nband=',nband
     if (nband.eq.1) then
        write(*,*) ' buggy mode, nband=1',nband
        stop
     elseif (nband.gt.3) then
        write(*,*) ' not implemented; nband=',nband
        stop
     end if
     
  elseif (iuse_mpc.ne.0) then
     write(*,*)  trim(command),': error',iuse_mpc,nband
     stop
  end if  


  command='iremove_singularity'
  call switching(command,iremove_singularity) ! default is 0
  if (iremove_singularity.eq.1) then 
     iset_bulk0=1;
     iuse_ffty=1
     if (myid.eq.0) write(*,*) ' for searching antisymmetric solution,', & 
          & ' phases are fixed in x,y,z'
     if (myid.eq.0) write(*,*) ' iset_bulk0 = 1 to remove singularity (phase fixed)' 
     if (myid.eq.0) write(*,*) ' iuse_ffty=1  to remove singularity (phase fixed in y)' 
  else if (iremove_singularity.ne.0) then
     write(*,*)  trim(command),': error',iremove_singularity
     stop
  end if

  command='iuse_scaling'
  call switching(command,iuse_scaling) ! default is 0
  if (iuse_scaling.eq.1) then 
     if (myid.eq.0) write(*,*) 'this mode modifies the jacobian of newton equation', & 
          & ' using diagonal-scaling'
  else if (iuse_scaling.ne.0) then
     if (myid.eq.0) write(*,*)  trim(command),': error',iuse_scaling
     stop
  end if

  command='ifilter'
  call switching(command,ifilter) ! default is 0
  if (ifilter.eq.1) then 
     if (myid.eq.0) write(*,*) 'this mode uses diagonal-scaling with zero-factor (ifilter)'
  elseif (ifilter.ne.0) then
     if (myid.eq.0) write(*,*)  trim(command),': error',ifilter
     stop
  end if

  comp_norm=1.d0
  if (myid.eq.0) then
     command='iuse_compensatory_norm'
     call go_command(command,nextline) ! DO NOT forget to BCAST
     if (nextline.ge.1) then
        read(hrefile(nextline),*) iuse_compensatory_norm
        nextline=nextline+1
        if (iuse_compensatory_norm.eq.1) then
           read(hrefile(nextline),*) comp_norm
           nextline=nextline+1
        end if
     end if
  end if
  call  MPI_BCAST(iuse_compensatory_norm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(comp_norm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  command='iuse_QR_orth'
  call switching(command,iuse_QR_orth) ! default is 0
  if (iuse_QR_orth.eq.1) then 
     if (myid.eq.0) write(*,*) 'this mode uses QR-factorization with Householder tranformations'
  elseif (iuse_QR_orth.ne.0) then
     if (myid.eq.0) write(*,*)  trim(command),': error',iuse_QR_orth
     stop
  end if

  ! ----------------------------------!
  if (iremove_singularity.eq.1) then
     ifilter=1; ! filtering (zero-precondition)
  endif
  if (iuse_scaling.eq.1) then
     ifilter=1; ! use filtering function for scaling (diagonal-scaling) !
                ! assuming iuse_v=1;
  endif

  if(myid.eq.0) then
     write(*,*) 'Optional switches: ishifty,iuse_mpc,nband,', & 
          & 'iremove_singularity,iuse_scaling,ifilter,iset_bulk0,iuse_ffty,', & 
          & 'iuse_hej,iuse_v,iuse_compensatory_norm,iuse_QR_orth,iuse_us,', & 
          & 'iydealiasing,inorm'
     write(*,'(a,15I3)')  'Optional switches: ',ishifty, & 
          & iuse_mpc,nband,iremove_singularity,iuse_scaling,ifilter, & 
          & iset_bulk0,iuse_ffty,iuse_hej,iuse_v,iuse_compensatory_norm, &
          & iuse_QR_orth,iuse_us,iydealiasing,inorm
     write(*,*) ' '
     write(*,*) ' === end of set_newton_options ==='
     write(*,*) ' '
  endif

end subroutine set_newton_options

!-------------------------------------------------------------------
!  tools for arnoldi, newton_gmres, arclength 
!--------------------------------------------------------------------
subroutine gmres_ini(vor,phi,u00,w00)

  use ctes
  use running
  use bcs
  use LES, only: iuse_LES, Cles, Clese
  use gmres
  use precond

  implicit none
  include "mpif.h"
  real(8),dimension(buffsize) :: vor,phi
  real(8),dimension(my) :: u00,w00
  real(8) sca_l(4), sca(4), mpc_ini,xmode,zmode
  real(8) mem_gmres, mem_precond, mem_eig
  real(8) mem_gmres_a, mem_precond_a, mem_eig_a
  integer ierr,index,i,k,j,iproc,leng,iband
  real(8),dimension(:),allocatable :: vel

  noescru = 1 ! not to write in matvec_dns
  nohist =  1 ! check the ener history
  nhist=1 ! the cf history is dumped at each step (rev.1373)
  ! --- save initial conditions ---- 
  time_ini=time
  etime = time + timep_sh
  xoffb_ini=xoffb
  xofft_ini=xofft
  ! save the boundary condition (used only for arclength along Axz)
  xoffb0=xoffb_ini
  xofft0=xofft_ini

  !!!write(*,*) 'save ini:',xoffb_ini,xofft_ini,xoffb,xofft
  !
  if (iuse_fou_mode.eq.3) then
     ! this is for streak instability
     ! --- save mean flow which is not thrown to x0 and xf ...
     allocate(vor_ini(buffsize),phi_ini(buffsize))
     vor_ini=0.d0; phi_ini=0.d0
     vor_ini=vor
     phi_ini=phi
     allocate(u00_ini(my),w00_ini(my))
     u00_ini=u00
     w00_ini=w00
  endif

  ! /* set normalization factor for throwing and backthrowing */
  !sca_l(3)=sum(vor*vor)
  !if (iuse_v.eq.1) then 
  !   allocate (vel(buffsize))
  !   vel(:)=0.d0
  !   call getv(phi,vel,1)
  !   sca_l(4)=sum(vel*vel)
  !   deallocate (vel)
  !else
  !   sca_l(4)=sum(phi*phi)
  !end if
  !
  !call MPI_ALLREDUCE(sca_l(3:4),sca(3:4),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  !sca_u00 = dsqrt(sum(u00*u00))
  !sca_w00 = dsqrt(sum(w00*w00))  
  !if (sca_w00.lt.tiny) then
  !   sca_w00=sca_u00
  !end if
  !sca_vor = dsqrt(sca(3))
  !sca_phi = dsqrt(sca(4))
  !
  sca_u00=s*Lz; sca_w00=s*Lz
  if (iuse_v.eq.1) then
     sca_vor=s; sca_phi=s*Lz ! here iuse_v =1 is assumed
     if (iuse_compensatory_norm.eq.1) then
        if (myid.eq.0) write(*,*) 'compensatory norm for vor and v',comp_norm
        sca_vor=s*dsqrt(s*re*Lz*Lz)/comp_norm; sca_phi=s*Lz ! here iuse_v =1 is assumed
     end if
  else
     sca_vor=s; sca_phi=1.d0  ! the scaling for phi is not implemented yet.
  endif
  
  if (iuse_scaling.eq.1) then
     if (myid.eq.0) write(*,*) 'iuse_scaling=1 mode use filter(nl) for diagonal scaling'
     ! reset off scaling factors 
     sca_u00=1.d0; sca_w00=1.d0
     sca_vor=1.d0; sca_phi=1.d0
  end if

  if (myid.eq.0) then
     write(*,*) ' sca_u00,w00,vor,phi:', sca_u00,sca_w00,sca_vor,sca_phi
     write(*,*) ' sca_time,sx,sz,re:', sca_time,sca_sx,sca_sy,sca_sz,sca_Ly,sca_re 
  endif

  if (abs(iarclength).eq.1) then
     ! ignore the parameters in hre2.dat (moved here from rev.1483)
     ! use parameters from the hearder (footer) in the field file
     !  recover the b.c. rev.1389 ! this should be the same with in the filinp. ! checked !
     xofft = xofft*alp/alpe
     xoffb = -xofft;
     xofft0=xofft ! bug fixed for updating alp, rev.1390
     xoffb0=xoffb  
     xofft_ini=xofft0
     xoffb_ini=xoffb0
     
     re = Ree; Ly = Lye; alp=alpe; 
     ! from (rev.1239)
     if (iadd_force.ne.0) then
        force_roll=force_rolle;
        xforce=xforcee;
        zforce=zforcee;
     end if
     if (iuse_LES.eq.1) then
        Cles=Clese;
     end if
     if (iadd_damping.eq.1) then
        damp_up=damp_upe
        damp_aa=damp_aae
     end if
     if (ireve.ge.1200) uprim=uprime
     call set_xyz_dy
     
     timep_sh = 2.d0*pi/alp/Ly*timep ! update timep_sh after updating Ly 
     etime = time_ini + timep_sh   
     
     if (myid.eq.0) then 
        write(*,*) ' set initial re =',re, ' (hre2.dat is ignored)'
        write(*,*) ' set initial Ly =',Ly, ' (hre2.dat is ignored)' 
        write(*,*) '            : check y0',y(0)
        write(*,*) ' set initial alp =',alp, ' (hre2.dat is ignored)'
        write(*,*) ' set initial uprim =',uprim, ' (hre2.dat is ignored)'  
     end if
  end if
  !iuse_fou_mode=2, forced to set 
  ! 1: does not work probably because of redundant 00-modes
  ! x-non-zeromodes of the result of 2 should be the same with that of 3
  if (myid.eq.0) then
     if (iuse_fou_mode.eq.1) then
        nmax=mx*mz*my*2 + my*2
     elseif (iuse_fou_mode.eq.2) then
        !nmax=mx*mz*my*2 - 2*my*2 + my*2 
        ! remove 00-modes of vor and phi, they are zero-constant
        nmax=mx*mz*my*2 - 2*my*2 + my*2 - (mz1-nz)*2*my*2
        ! remove kx0 conjugate modes, and recover those modes by set_conj in backthrowing.
     elseif (iuse_fou_mode.eq.3) then
        nmax=(mx-2)*mz*my*2 ! remove x0-modes in x-dir of vor and phi
     end if
     nmax1=nmax ! for arnoldi nmax1=nmax
  else
     !nmax=0; nmax1=0; 
     nmax=1; nmax1=1; ! bug fixed rev.1395 
  endif
   
  nf=nmax1;

  if (isol_mode.eq.3) then
     if (ishifty.eq.1) then 
        nmax=nmax + 1
        icsy=nmax
     endif
  elseif (isol_mode.eq.4) then
     ! for relative periodic orbit in (x,z)
     nmax=nmax + 2
     if (ishifty.eq.1) then 
        nmax=nmax + 1
        icsx=nmax-2
        icsy=nmax-1
        icsz=nmax
     else
        icsx=nmax-1
        icsz=nmax
     end if
  elseif (isol_mode.eq.5) then
     nmax=nmax + 1
     icTp=nmax
  elseif (isol_mode.eq.6) then
     nmax=nmax + 3
     !if (ishifty.eq.1) then 
     !   write(*,*) 'isol_mode=6 does not work with ishifty=1: stop'
     !   ! testing at rev1562
     !endif
     icsx=nmax-2
     icsz=nmax-1
     icTp=nmax
  elseif (isol_mode.gt.6) then
     write(*,*) 'not implemented isol_mode=',isol_mode
  end if

  if (iarclength.eq.1) then
     nmax=nmax + 1
     icre=nmax
  end if

  para_newton(20)=nmax
  write(*,*) myid,'set number of unknowns:', nmax

  nl=nmax ! for sequential
  ! slaves only have re,shiftx,shiftz,Tp

  ! only for sequential gmres, but DNS is parallelized
  allocate(x0(nl),xf(nl),dfp(nl),xtmp(nl),hej(nl),btmp(nl))
  x0=0.d0; xf=0.d0; dfp=0.d0; xtmp=0.d0; hej=0.d0; btmp=0.d0
  ! --- allocate vv hh etc... ! 
  allocate(vv(nl,k_gmres+1),hh(k_gmres+1,k_gmres),gh(k_gmres+1,k_gmres))
  allocate(rr(nl)); 
  vv=0.d0; hh=0.d0; gh=0.d0; rr=0.d0;
  allocate(zz(nl)); ! for storing the tmporarily vector 
  zz=0.d0; 
 
  mem_gmres=dfloat( nl*(8 + (k_gmres+1) ) + (k_gmres+1)*(k_gmres+1)*2 + 3*buffsize) 
  mem_gmres=mem_gmres + 2.d0*nl + 2.d0*(buffsize+my) + nl + 4.d0*(buffsize +my) &
                 + nl*(2.d0 + 2.d0 + 1) 
  if (iuse_QR_orth.eq.1) then
     allocate(QQ(nl,k_gmres+1))
     QQ=0.d0;
     allocate(tau(k_gmres+1))
     tau=0.d0;
     !write(*,*) myid,'allocate a big matrix'
  end if 
  mem_gmres=mem_gmres + k_gmres + nl*(k_gmres+1)
  ! add_extra_memory requests: iuse_v, make_mat, 2nd-order-arcl dxdam, throwing,
  ! iuse_ffty and ccm, ccm2, hej, xtmp and filter
  call MPI_ALLREDUCE(mem_gmres,mem_gmres_a,1,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
  if(myid.eq.0) write(*,'(i4,a,G9.2)') myid, & 
       & ': gmres requests (MB):',mem_gmres_a*8.d0/1024.d0/1024.d0

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (iuse_mpc.eq.1) then
     !write(*,*) 'allocate for preconditioner'

     ! rev1151: this mode produce non-conjugate relation for kx=0 when it is used with addsym=5

     ! /* preconditionar */
     allocate(mpcx(nl),mpcx_in(nl)) ! 1d-diagonal vector
     mem_precond=2.d0*nl
     
     mpc_ini=1.d0

     mpcx(:)=mpc_ini; mpcx_in(:)=mpc_ini

     write(*,*) myid,'preconditiner band=',nband
     allocate(u00p(0:my1),w00p(0:my1),vorp(buffsize),phip(buffsize))
     u00p=0.d0; w00p=0.d0; vorp=0.d0; phip=0.d0;     
     allocate(u00q(0:my1),w00q(0:my1),vorq(buffsize),phiq(buffsize))
     u00q=0.d0; w00q=0.d0; vorq=0.d0; phiq=0.d0;
     
     mem_precond=mem_precond + (my+buffsize)*4.d0

     if (nband.eq.0) then
        allocate(mpc_u00(1,0:my1),mpc_w00(1,0:my1),mpc_vor(1,buffsize),mpc_phi(1,buffsize))
        allocate(mpcl_u00(1,0:my1),mpcl_w00(1,0:my1),mpcl_vor(1,buffsize),mpcl_phi(1,buffsize))
        allocate(mpc_vorphi(1,buffsize),mpc_phivor(1,buffsize))
        allocate(mpcl_vorphi(1,buffsize),mpcl_phivor(1,buffsize))
     else
        allocate(mpc_u00(nband,0:my1),mpc_w00(nband,0:my1),mpc_vor(nband,buffsize),mpc_phi(nband,buffsize))
        allocate(mpcl_u00(nband,0:my1),mpcl_w00(nband,0:my1),mpcl_vor(nband,buffsize),mpcl_phi(nband,buffsize))
        allocate(mpc_vorphi(nband,buffsize),mpc_phivor(nband,buffsize))
        allocate(mpcl_vorphi(nband,buffsize),mpcl_phivor(nband,buffsize))
     end if
     mpc_u00(1,:)=mpc_ini; mpc_w00(1,:)=mpc_ini;    
     mpc_vor(1,:)=mpc_ini; mpc_phi(1,:)=mpc_ini;
     mpc_vorphi(1,:)=mpc_ini*0.d0; mpc_phivor(1,:)=mpc_ini*0.d0;     
     do iband=2,nband
        mpc_u00(iband,:)=0.d0; mpc_w00(iband,:)=0.d0; 
        mpc_vor(iband,:)=0.d0; mpc_phi(iband,:)=0.d0;
        mpc_vorphi(iband,:)=0.d0; mpc_phivor(iband,:)=0.d0;
     end do
     if (iuse_ffty.eq.1) then
        ! set constant in y-dir, in fourier space         
        mpc_u00(1,1:my1)=0.d0;
        mpc_w00(1,1:my1)=0.d0;
     end if
     mpcl_u00=mpc_u00; mpcl_w00=mpc_w00
     mpcl_vor=mpc_vor; mpcl_phi=mpc_phi
     mpcl_vorphi=mpc_vorphi; mpcl_phivor=mpc_phivor
     
     mem_precond=mem_precond + (my+buffsize)*8.d0*nband
     mem_precond=mem_precond + (buffsize)*4.d0*nband
     

     !if (iremove_singularity.eq.1) then
        if(myid.eq.0) then
           index=my*2
           do iproc=0,numerop-1
              if (iproc.ne.0) then
                 leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
                 ! Note:
                 ! all the (jend(iproc)-jbeg(iproc)+1) should be the same with 
                 ! (or smaller than) je-jb+1
              endif
              do j=jbeg(iproc),jend(iproc)
                 do k=1,mz
                    do i=1,mx 
                       
                       if (((i.eq.1).and.(k.eq.1)).or.((i.eq.2).and.(k.eq.1))) then                         
                          ! skip 00-modes of vor and phi
                          !elseif (((i.eq.4)).and.((k.eq.2).or.(k.eq.mz))) then
                       elseif ( (j.eq.my/2).and.((i.eq.4).and.(k.eq.1)) ) then
                          !mpcx(index)=0.d0 ! Imaginary part of vor(1,0,y=0)
                          ! note: index=vor,  index+1=phi 
                          index=index+2    
                       elseif ( (j.eq.my/2).and.((i.eq.2).and.(k.eq.2)) ) then
                          !mpcx(index)=0.d0 ! Imaginary part of vor(0,1,y=0)
                          index=index+2    
                       elseif ( (j.eq.my/2).and.((i.eq.2).and.(k.eq.mz)) ) then
                          !mpcx(index)=0.d0 ! Conjugate of Imaginary part of vor(0,1,y=0)
                          index=index+2    
                       else
                          index=index+2           
                       endif
                       
                    enddo
                 enddo
              enddo
              !write(*,*) 'setmpcx check, index,nl:',index,nf,nl,nmax1,nmax
           enddo
        end if
     !endif

     allocate(MTdx(nl),Mdf(nl))
     MTdx=0.d0; Mdf=0.d0
     mem_precond=mem_precond + nl*2.d0 

     if (ishifty.eq.1) mpcx(icsy)=1.d0/sca_sy
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        mpcx(icsx)=1.d0/sca_sx
        mpcx(icsz)=1.d0/sca_sz
     end if
     if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
        allocate(mpctp0(nl),mpctpf(nl))
        allocate(mpctp0_in(nl),mpctpf_in(nl))
        mpcx(icTp)=1.d0  !/sca_time
        mpcx_in=mpcx

        mpctp0(1:nf)=1.d0; mpctp0(icTp)=1.d0
        mpctpf(1:nf)=1.d0; mpctpf(icTp)=1.d0
        mpctp0_in=mpctp0; 
        mpctpf_in=mpctpf; 
        
        mem_precond=mem_precond + nl*4.d0
 
     end if
  end if

  if (ifilter.eq.1) then
     allocate(filter(nl)) ! 1d-diagonal vector
     filter=1.d0
      
     mem_precond=mem_precond + nl
     
     if (iremove_singularity.eq.1) then
        if(myid.eq.0) then
        
           index=my*2+1
           if (iuse_ffty.eq.1) then
              if (iset_bulk0.eq.1) then
                 filter(1)=0.d0; filter(my+1)=0.d0;
              end if
              filter(4)=0.d0 ! imag(u00f(1)) 
           else
              write(*,*) 'iuse_ffty should be 1 to remove singularity of translation in y-dir.'
           end if
           do iproc=0,numerop-1
              !do j=0,jend(iproc)-jbeg(iproc) ! bug fixed are rev1558
              do j=jbeg(iproc),jend(iproc)
                 do k=1,mz
                    do i=1,mx 
                       if (((i.eq.1).and.(k.eq.1)).or.((i.eq.2).and.(k.eq.1))) then
                          ! skip 00-modes of vor and phi
                          !elseif (((i.eq.4)).and.((k.eq.2).or.(k.eq.mz))) then
                       elseif ( (i.le.2).and.( k.ge.nz+2 ) ) then 
                          !filter(index:index+1)=0.d0 ! 2*2*(mz-nz+1)*my 
                          ! remove conj. imaginary part of vor(0,k=nz+1:mz1,y)
                          !!!!index=index+2 ! bug fixed at rev.1557
                          !! /* see below for checking the index */
                          !!do k=0,nz
                          !!   xgam(k) = dcmplx(0.d0,gam*dfloat(k))
                          !!   icx(k) = k
                          !!enddo
                          !!
                          !!do k=nz+1,mz1
                          !!   xgam(k) = dcmplx(0.d0 ,-gam*dfloat(mz-k))
                          !!   icx(k)  = mz-k
                          !!enddo  
                       else
                          if ( (j.eq.my/2).and.((i.eq.4).and.(k.eq.1)) ) then
                             filter(index)=0.d0 ! Imaginary part of vor(1,0,y=0)
                             ! note: index=vor,  index+1=phi                              
                          elseif ( (j.eq.my/2).and.((i.eq.2).and.(k.eq.2)) ) then
                             filter(index)=0.d0 ! Imaginary part of vor(0,1,y=0)
                          elseif ( (j.eq.my/2).and.((i.eq.2).and.(k.eq.mz)) ) then
                             filter(index)=0.d0 ! conj. of Imaginary part of vor(0,1,y=0)
                          end if
                          index=index+2           
                       endif
                    enddo
                 enddo
              enddo
              !write(*,*) 'filter check, index,nl:',index,nf,nl,nmax1,nmax
           enddo

           if (index-1.ne.nf) then
              write(*,*) 'filter index error, index,nf,nl,nmax1,nmax:',index,nf,nl,nmax1,nmax
              stop
           end if
        end if
     endif

     if ((iuse_scaling.eq.1).and.(iuse_v.eq.1)) then
        write(*,*) 'iuse_scaling is not implemented, stop'
        stop
        ! here set scaling function. the 1d array 'filter' is used 
        if(myid.eq.0) then
           index=1
           do j=0,my1
              filter(index)=(s*Lz)
              index=index+1
           enddo
           do j=0,my1
              filter(index)=(s*Lz)
              index=index+1
           enddo
           index=my*2
           do iproc=0,numerop-1
              do j=jbeg(iproc),jend(iproc)
                 do k=1,mz
                    do i=1,mx 
                       if (((i.eq.1).and.(k.eq.1)).or.((i.eq.2).and.(k.eq.1))) then
                          ! skip 00-modes of vor and phi
                          !elseif (((i.eq.4)).and.((k.eq.2).or.(k.eq.mz))) then
                       else
                         
                          if ( (i-1)/2.eq.0) then
                             xmode=1.d0;
                          else
                             if (mod(i-1,2).eq.0) then
                                xmode=(i-1)/2
                             elseif (mod(i-1,2).eq.1) then
                                xmode=(i)/2
                             end if
                          endif
                          if (icx(k-1).eq.0) then
                             zmode=1.d0
                          else
                             zmode=dfloat(icx(k-1))
                          end if
                          ! scaling for vor
                          !filter(index)=(s*Lz*dsqrt(s*re))/2.5d0 
                          !filter(index)=(s*dsqrt(s*re*Lz*Lz))  
                          filter(index)=s/zmode*(xmode*xmode+zmode*zmode) 
                          ! s*dsqrt(Rez)*2.5 ! for vorticity
                          index=index+1
                          !
                          ! scaling for v
                          filter(index)=(s*Lz) !/dsqrt(xmode*zmode)  
                          ! scale v by s*Lz, assuming iuse_v=1
                          index=index+1
                       endif
                    enddo
                 enddo
              enddo
              !write(*,*) 'filter check, index,nl:',index,nf,nl,nmax1,nmax
           enddo

        end if ! master

        !if (ishifty.eq.1) filter(icsy)= 100.d0 !1.d0/sca_sy
        !if ((isol_mode.eq.5).or.(isol_mode.eq.6)) filter(icTp)= 100.d0
        if (iarclength.eq.1) filter(icre)= sca_re
     end if
  endif
  call MPI_ALLREDUCE(mem_precond,mem_precond_a,1,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
  if(myid.eq.0) write(*,'(i4,a,G9.2)') myid, & 
       &    ': precond requests (MB):',mem_precond_a*8.d0/1024.d0/1024.d0

end subroutine gmres_ini

subroutine newton_ini(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor,phi,chwk
  real(8),dimension(my) :: u00,w00
  real(8) eps_j,mem_newton,mem_newton_a,mem_hook
  integer lwk,ierr

  ! --- check unkowns ---
  allocate(dxdt0(nf))
  mem_newton=nl
  dxdt0=0.d0;
  if (isol_mode.eq.2) then
     allocate(dx0(nf),dxf(nf))
     dx0=0.d0; dxf=0.d0;
     mem_newton=mem_newton + 2*nl
  elseif (isol_mode.eq.3) then
     allocate(dxdt(nf))
     dxdt=0.d0
     mem_newton=mem_newton + nl
  elseif (isol_mode.ge.4) then
     allocate(dx0(nf),dxf(nf),dz0(nf),dzf(nf))
     allocate(dxdt(nf))
     dx0=0.d0; dxf=0.d0; dz0=0.d0; dzf=0.d0;
     dxdt=0.d0
     mem_newton=mem_newton + 5*nl
  end if
  if (ishifty.eq.1) then
     allocate(dy0(nf),dyf(nf))
     dy0=0.d0; dyf=0.d0;
     mem_newton=mem_newton + 2*nl
  endif

  allocate(bb(nl),dxp(nl),Ax(nl))
  bb=0.d0; dxp=0.d0; Ax=0.d0; 
  mem_newton=mem_newton + 3*nl

  allocate(x0_backup(nl),xf_backup(nl)) ! xf is also backuped from rev.1115
  mem_newton=mem_newton + 2*nl

  if (abs(iarclength).eq.1) then
     allocate(dxda(nf),xstan(nf))
     allocate(xm1(nl))
     mem_newton=mem_newton + 3*nl
  endif

  ! allocate for newton-gmres
  allocate(cs(k_gmres),ss(k_gmres),ee(k_gmres+1))
  cs=0.d0; ss=0.d0; ee=0.d0
  allocate(yg(k_gmres)) ! bug fixed (rev.436)
  yg=0.d0;

  mem_newton=mem_newton + 4*(k_gmres+1);
  mem_newton=mem_newton + 3*nl; ! update_hookstep

  call MPI_ALLREDUCE(mem_newton,mem_newton_a,1,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
  if (myid.eq.0) write(*,'(i4,a,G9.2)') myid, & 
       & ': newton requests (MB):',(mem_newton_a)*8.d0/1024.d0/1024.d0

  ! memory check for hookstep
  call hookstep(k_gmres,0,chwk,2)
  call hookstep(k_gmres,0,chwk,-1)

  ! memory check for eig
  call get_eigen(k_gmres,2,'v')

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! set some paramters for updating newton steps
  if (iuse_hookstep.eq.1) then
     dampfac=2.d0
     ! for hookstep
     !delta    = 0.1d0   ! trust region
     delta    = delta_ini   ! trust region
     deltaMin = 1.d-8  ! stop if radius of trust region gets this small 
     deltaMax = delta*10.d0   ! maximum radius of trust region
     deltaFuzz= 0.01d0   ! accept steps within (1+/-deltaFuzz)*delta
     if(myid.eq.0) write(*,*) myid,'hookstep, set initial delta:', delta
     nhookstep= 50   ! this reduce delta ~> 0.01*delta 
     dlambda=0.5d0
     dlambdaMin= 0.2d0  ! minimum delta shrink rate
     dlambdaMax= 1.5d0  ! maximum delta expansion rate
     if(myid.eq.0) write(*,*) myid,'hookstep, set initial delta_rate:', dlambda

     ! for updating delta (Newton radius)
     dimprovReq=1e-3
     dimprovOk=0.10
     dimprovGood=0.75
     dimprovAcc=0.10 ! ??
  else
     dampfac=0.2d0
     damprate=1.5
     ndamp=15
  endif

end subroutine newton_ini

subroutine arclength_ini(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use bcs
  use gmres
  use LES
  !use precond

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor,phi,chwk
  real(8),dimension(my) :: u00,w00

  if (abs(iarclength).eq.1) then
     if (myid.eq.0) write(*,*) 'arclength=1: cpara =', trim(cpara)
     
     ! add error check here
     if (trim(cpara).eq.'Rez') then
        if (myid.eq.0) write(*,*) 'arclength along Rez'

     elseif (trim(cpara).eq.'Ayz') then
        if (myid.eq.0) write(*,*) 'now testing arclength for Ayz'
        if ((iuse_LES.eq.1).and.(ifix_CsDeltag.eq.0)) then
           write(*,*) 'set ifix_CsDeltag=1 in hre3.dat'
           stop
        end if

     elseif (trim(cpara).eq.'Axz') then
        if (myid.eq.0) write(*,*) 'now testing arclength for Axz'
        if ((iuse_LES.eq.1).and.(ifix_CsDeltag.eq.0)) then
           write(*,*) 'set ifix_CsDeltag=1 in hre3.dat'
           stop
        end if
        !stop
     elseif (trim(cpara).eq.'damp_up') then
        if (myid.eq.0) write(*,*) 'now testing arclength for damp_up'
        !stop
     elseif (trim(cpara).eq.'damp_aa') then
        if (myid.eq.0) write(*,*) 'now testing arclength for damp_aa'
        !stop
     elseif (trim(cpara).eq.'LES_Cs') then
        if (myid.eq.0) write(*,*) 'now testing arclength for Cs of LES'
        if ((iuse_LES.eq.1).and.(ifix_CsDeltag.eq.1)) then
           write(*,*) 'set ifix_CsDeltag=0 in hre3.dat'
           stop
        end if
     else
        if (myid.eq.0) write(*,*) ' cpara error, arclength_ini',trim(cpara)
        stop
     endif
  endif
  
  ! setting alp <= alpe is moved in gmres_ini, note that confusing xofft modification 
  ! (see also in readwrite.f90 to update xofft... rev.1483)
  !--------------------------------------------------------------------------
  if (myid.eq.0) write(*,*) 'reading previous solutions to compute xstan from:', trim(filinm1)
  filinp = filinm1     
  call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c., Ree is overwritten
  ! Note: here xofft and xoffb are overwritten!
  ! modify the b.c. (rev.1392)
  if (myid.eq.0) write(*,*) 'check xofft',xofft,xofft0
  xofft = xofft0 
  xoffb = xoffb0
  
  if (iadd_sym.ge.5) then
     !call MPI_BCAST(sym_shiftx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     !call MPI_BCAST(sym_shiftz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call phase_shift_fou(vor,phi,sym_shiftx*Lx,sym_shiftz*Lz)        
  end if
  
  call throwing(vor,phi,u00,w00,chwk,xm1,1) ! read previous solution
  rem1=Ree; Lym1=Lye; alpm1=alpe;
  if (iadd_force.ne.0) then
     force_rollm1=force_rolle;
     xforcem1=xforcee;
     zforcem1=zforcee;
  end if
  if (iuse_LES.eq.1) then
     Clesm1=Clese;
  end if
  if (iadd_damping.eq.1) then
     damp_upm1=damp_upe
     damp_aam1=damp_aae
  end if
  if (ireve.ge.1200) uprimm1=uprime
  
  ! check boundary condition
  !if (trim(cpara).eq.'Axz') then
  !   write(*,*) myid,'arclength_ini (Axz), check xofft',xofft,xoffb
  !end if

end subroutine arclength_ini


subroutine set_arclength(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use bcs
  use gmres
  use LES

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor,phi,chwk
  real(8),dimension(my) :: u00,w00
  real(8) ptan,dpara

  !
  ! update xstan 
  if (myid.eq.0) write(*,*) ' update xstan, cpara= ',trim(cpara)
  xstan(1:nf) = x0(1:nf) - xm1(1:nf)
  if (trim(cpara).eq.'Rez') then 
     ptan = re - rem1
  elseif (trim(cpara).eq.'Ayz') then
     ptan = Ly - Lym1
  elseif (trim(cpara).eq.'Axz') then
     ptan = alp - alpm1
  elseif (trim(cpara).eq.'damp_up') then
     ptan = damp_up - damp_upm1
  elseif (trim(cpara).eq.'damp_aa') then
     ptan = damp_aa - damp_aam1
     !write(*,*) myid,'damp_aa,aam1',damp_aa,damp_aam1
  elseif (trim(cpara).eq.'LES_Cs') then
     ptan = Cles - Clesm1
     !write(*,*) myid,'debug, LES_Cs',Cles,Clesm1
  else 
     write(*,*) 'cpara error, xstan, stop '
     stop
  endif
  !
  if (ifilter.eq.1) then ! rev.1129
     call scaling(xstan,zz,1)
     call norm(nf,zz,xstannorm) ! rev.(633) nl -> nf
  else
     call norm(nf,xstan,xstannorm) ! rev.(633) nl -> nf
  end if
  !
  if ((xstannorm.lt.tiny).or.(abs(ptan).lt.tiny)) then
     if (myid.eq.0) write(*,*) 'error to get tangent: tan=',xstannorm,abs(ptan)
     stop
  endif
  !xstan = xstan/xstannorm
  ! do not normalize for the better initial guess of arclength. (rev1397)

  ! save the solution to restore later, xm1
  xm1 = x0         
  if (trim(cpara).eq.'Rez') then 
     rem1 = re
  elseif (trim(cpara).eq.'Ayz') then
     Lym1 = Ly
  elseif (trim(cpara).eq.'Axz') then
     alpm1 = alp
  elseif (trim(cpara).eq.'damp_up') then
     damp_upm1 = damp_up
  elseif (trim(cpara).eq.'damp_aa') then
     damp_aam1 = damp_aa
  elseif (trim(cpara).eq.'LES_Cs') then
     Clesm1 = Cles
  else
     write(*,*) myid,'cpara error xm1'
     stop
  endif
  
  if (arcl.lt.tiny) then
     if (myid.eq.0) write(*,*) 'error,  arcl =',arcl
     stop
  endif
  
  ! update the initial unknown vector and parameter.
  if (iarclength.eq.-1) then
     ! this mode simply changes the parameter with the interval of arcl_ini
     ! which is read from hre3.dat
     x0(1:nf) = x0(1:nf) + xstan(1:nf)*arcl_ini
     dpara = ptan/abs(ptan)*arcl_ini 

  elseif (iarclength.eq.1) then
     ! arcl parameter will be adjusted automatically around 1, 
     ! only at the first step arcl will be arcl_ini read from hre3.dat
     x0(1:nf) = x0(1:nf) + xstan(1:nf)*arcl
     dpara = ptan*arcl
     ! arcl must be delta ??
  end if
  
  if (trim(cpara).eq.'Rez') then
     re = re + dpara     
  elseif (trim(cpara).eq.'Ayz') then
     Ly = Ly + dpara 
     call set_xyz_dy ! initialize also compact finite difference ...
     timep_sh = 2.d0*pi/alp/Ly*timep
     etime = time_ini + timep_sh   
  elseif (trim(cpara).eq.'Axz') then
     alp = alp + dpara
     call set_xyz_dy ! initialize also compact finite difference ...
     timep_sh = 2.d0*pi/alp/Ly*timep
     etime = time_ini + timep_sh
     ! update boundary condition (rev.1389)
     xofft = xofft*(alp-dpara)/alp
     xoffb = -xofft     
     xofft0=xofft ! bug fixed for updating alp, rev.1390
     xoffb0=xoffb  
  elseif (trim(cpara).eq.'damp_up') then
     damp_up = damp_up + dpara 
  elseif (trim(cpara).eq.'damp_aa') then
     damp_aa = damp_aa + dpara 
  elseif (trim(cpara).eq.'LES_Cs') then
     Cles = Cles + dpara
     if (ifix_CsDeltag.eq.1) then
        write(*,*) 'please set ifix_CsDeltag=0'
        stop
     end if
  else
     write(*,*) myid,'cpara error, dpara'
     stop
  endif
  !
  ! throwing numerical parameters to unknown vector ! rev(1282)
  call throwing_parameter(x0,0)
  ! check boundary condition
  !if (trim(cpara).eq.'Axz') then
  !   write(*,*) myid,'set_arclength (Axz), check xofft',xofft,xoffb
  !end if
end subroutine set_arclength


subroutine make_mat(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor,phi,chwk
  real(8),dimension(my) :: u00,w00
  real(8) eps_j, res, fac
  real(8) daa, re0, Ly0, alp0, uprim0
  real(8) randu
  integer i,idum

  real(8),dimension(:), allocatable:: dvorwk,dphiwk,du00wk,dw00wk
!  real(8),dimension(:), allocatable:: dxdam

  ! set x0, xf  
  xofft=xofft0; xoffb=xoffb0; ! rev1395 ! fixed rev.1398
  call throwing(vor,phi,u00,w00,chwk,x0, 1) ! --> x0 
  ! (we do not need this ... rev.1113)
  !   but vor,phi should be x0..., not xf  ==> reverged at rev.1150

  if ( (isol_mode.ge.4).or.((ishifty.eq.1).and.(isol_mode.eq.3)) ) then
     allocate(dvorwk(buffsize),dphiwk(buffsize))
     allocate(du00wk(my),dw00wk(my))
     du00wk=0.d0; dw00wk=0.d0;     
     !
     if (isol_mode.ge.4) then
        ! derivx --> dx0 (for shifting in x) ! only for isol_mode = 2 or 4
        dvorwk=vor; dphiwk=phi; ! just copy
        call deriv_xfou(dvorwk,dphiwk,1) ! first derivative in x
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dx0, 0) ! --> dx0
        !
        ! derivz --> dz0 (for shifting in z) ! skip now
        dvorwk=vor; dphiwk=phi; ! just copy
        call deriv_zfou(dvorwk,dphiwk,1) ! first derivative in x
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dz0, 0) ! --> dz0
        !
        ! save previous shift (not needed??)
        ! scaling for shiftx,z??
        !x0(icsx) = shiftx/sca_sx; xf(icsx) = shiftx/sca_sx
        !x0(icsz) = shiftz/sca_sz; xf(icsz) = shiftz/sca_sz
        ! moved to throwing_parameter() (rev.1239)
     end if
     if (ishifty.eq.1) then
        ! derivy --> dy0 (for shifting in y) ! only for isol_mode = 2 or 4
        dvorwk=vor; dphiwk=phi; ! just copy
        du00wk=u00; dw00wk=w00;
        ! set b.c.
        !write(*,*) 'deriv_yfou: check boundary condition stop:'
        !stop
        xofft=xofft0; xoffb=xoffb0; ! rev.1395
        call deriv_yfou(dvorwk,dphiwk,du00wk,dw00wk,1) ! first derivative in x
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dy0, 0) ! --> dy0
        !x0(icsy) = shifty/sca_sy; 
        !xf(icsy) = shifty/sca_sy;
        ! moved to throwing_parameter() (rev.1239)
     endif
     if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
        ! scaling for time period?
        !x0(icTp) = timep_sh; xf(icTp) = timep_sh !/sca_time
        ! moved to throwing_parameter() (rev.1239)
     end if
  endif
  !
  ! throwing numerical parameters to unknown vector ! rev(1239)
  call throwing_parameter(x0,0)
     
  !
  !if(myid.eq.0) write(*,*) 'compute dxdt0 by r.h.s. of NS'
  xofft=xofft0; xoffb=xoffb0; ! fixed rev.1400
  call set_dxdt(vor,phi,u00,w00,chwk,dxdt0) ! and for a fixed point
  !if (ifilter.eq.1) then
  !   call scaling(dxdt0,zz,1)
  !   call norm(nf,zz,dxdt0norm)
  !else
     call norm(nf,dxdt0,dxdt0norm)
  !end if
  !
  ! for debug rhs of NS => impossible, for shear-periodic boundary 
  ! reset boundary condition
  ! time = time_ini
  ! xoffb = xoffb_ini 
  ! xofft = xofft_ini
  !  call save_vec(vor,phi,u00,w00,chwk,dxdt0)
  !  dxdt0=0.d0
  !  call set_dxdt_dns(vor,phi,u00,w00,chwk,dxdt0) ! and for a fixed point
  !
  !  call save_vec(vor,phi,u00,w00,chwk,dxdt0)
  !  stop

  ! !!moved out of here to newton outer-outer-loop. ! ==> but reverted at (rev.1148)
  if (isol_mode.ge.3) then ! finding orbits
     ! this is also used to accumulate statistiscs by nohist=0, before calling make_mat
     nohist=0 ! screen output and dump *.cf, accumulate stat, but not dump stat. file
     xofft=xofft0; xoffb=xoffb0; ! fixed rev.1400
     call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     nohist=1 ! reset the history output option
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        call phase_shift_fou(vor,phi,shiftx,shiftz)
        if (ishifty.eq.1) then
           !   extshiftx=s*0.5d0*timep_sh*shifty
           !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
           fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
           call phase_shifty(vor,phi,u00,w00,shifty,fac)

        endif
     end if
     xofftf=xofft; xoffbf=xoffb; ! rev1395
     call throwing(vor,phi,u00,w00,chwk,xf ,1) ! --> xf
     call throwing_parameter(xf,0) ! moved here (rev.1400)
  end if

  !write(*,*) myid,'set r.h.s of Ax=bb'
  if (isol_mode.eq.1) then
     bb(1:nf) = -dxdt0(1:nf) ! set initial r.h.s of Ax=b for a fixed point
  elseif (isol_mode.ge.3) then
     bb(1:nf) = x0(1:nf)-xf(1:nf) ! set initial r.h.s of Ax=b
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        ! periodic orbit with Tp being fixed 
        bb(icsx)=0.d0
        bb(icsz)=0.d0       
     elseif ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
        ! periodic orbit with Tp being unknown
        bb(icTp)=0.d0
     end if
     if (ishifty.eq.1) bb(icsy)=0.d0
  end if

  if (iarclength.eq.1) then
     bb(icre) = 0.d0 
  endif
  !if (ifilter.eq.1) then
  !   call scaling(bb,zz,1); call norm(nf,zz,bnorm)
  !   call scaling(xf,zz,1); call norm(nf,zz,xfnorm)
  !   call scaling(x0,zz,1); call norm(nf,zz,x0norm)
  !else
     call norm(nl,bb,bnorm)
     call norm(nf,xf,xfnorm)
     call norm(nf,x0,x0norm)  
  !end if
  ! /* convergence check of newton_gmres */
  res=bnorm/x0norm 
  call norm_i(nl,x0,bb)
  if (res.lt.eps_newton) then
     if (myid.eq.0) then
        write(*,*) 'make_mat: check res= ', res
        write(*,*) '--------------------------------------------------------'
        write(*,'(a,3(1X,E18.10))') ' make_mat, norms: bb,x0,xf',bnorm,x0norm,xfnorm
        write(*,*) '--------------------------------------------------------'
     endif
     ! /* deallocate before return */
     if ( ((isol_mode.eq.4).or.(isol_mode.eq.6)).or. & 
          & ((isol_mode.eq.3).and.(ishifty.eq.1)) ) then
        deallocate(dphiwk,dvorwk)
        deallocate(du00wk,dw00wk)
     end if
     dxp=bb; ! fixed at rev.1446
     return
  end if

  ! set for the other derivatives (the sensitivities for each parameter)
  if ( ((isol_mode.eq.4).or.(isol_mode.eq.6)).or. & 
       & ((isol_mode.eq.3).and.(ishifty.eq.1)) ) then 

     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        ! derivatives in x and z should be applied to phase-shifted flow
        ! derivx --> dxf (for shifting in x) 
        dvorwk=vor; dphiwk=phi; ! just copy
        du00wk=0.d0; dw00wk=0.d0; ! dummy
        call deriv_xfou(dvorwk,dphiwk,1) ! first derivative in x
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dxf,0) ! --> dxf
        
        ! derivz --> dzf (for shifting in z) 
        dvorwk=vor; dphiwk=phi; ! just copy
        du00wk=0.d0; dw00wk=0.d0; ! dummy
        call deriv_zfou(dvorwk,dphiwk,1) ! first derivative in z
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dzf,0) ! --> dzf     
     end if
     !
     if (ishifty.eq.1) then
        ! derivy --> dyf (for shifting in y) ! only for isol_mode = 2 or 4
        dvorwk=vor; dphiwk=phi; ! just copy
        du00wk=u00; dw00wk=w00;
        call deriv_yfou(dvorwk,dphiwk,du00wk,dw00wk,1) ! first derivative in y
        call throwing(dvorwk,dphiwk,du00wk,dw00wk,chwk,dyf,0) ! --> dyf
     endif
     deallocate(dphiwk,dvorwk)
     deallocate(du00wk,dw00wk) ! dummy 
  end if

  if (isol_mode.ge.3) then 
     !if(myid.eq.0) write(*,*) 'compute dxdt by r.h.s. of NS'
     xofft=xofftf; xoffb=xoffbf; ! set b.c. 
     call set_dxdt(vor,phi,u00,w00,chwk,dxdt) ! Acutually, isol_mode=3 does not need this 
     !if (ifilter.eq.1) then
     !   call scaling(dxdt,zz,1)
     !   call norm(nf,zz,dxdtnorm)
     !else
        call norm(nf,dxdt,dxdtnorm)
     !end if
  end if

  if (iarclength.eq.1) then
     ! the sensitivity for delta re
     ! moved to set_dxda() rev.1239
     xofft=xofft0; xoffb=xoffb0; ! reset b.c.
     call set_dxda(vor,phi,u00,w00,chwk) ! vor phi are reset to be x0 inside...
  endif

  ! now vor and phi is xf...

  !set initial vector of the next newton step
  dxp=bb;
  if (iadd_sym.eq.-1) then
     if(myid.eq.0) write(*,*) 'add randum disturbance to initial vector', &
          ' to break symmetry of the initial guess' 
     ! usually this shows bad convergence property
     idum=-123456
     do i=1,nl 
        dxp(i)=dxp(i) + eps_jacobi*randu(idum); 
        ! add randum perturabation (rev.587 ==> rev.1158 )
     end do
  end if
  !
  ! --- compute A*x0 ---
  eps_j=eps_jacobi*x0norm
  ! note that vor phi, ... are buffers for calc_DNS
  call matvec_dns_disk3d(vor,phi,u00,w00,chwk,dxp,Ax, &
       eps_j,rnorm)
  !the residual of newton equation
  rr=bb-Ax
  !if (ifilter.eq.1) then
  !   call scaling(rr,zz,1)
  !   call norm(nf,zz,r0norm)    
  !elseif
     call norm(nl,rr,r0norm)
  !end if
  if(myid.eq.0) then
     write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
     write(*,'(a,3(1X,E18.10))') ' make_mat, norms: bb,xf,x0',bnorm,xfnorm,x0norm
     write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  endif
  if (bnorm.lt.tiny) then
      write(*,*) myid,'make_mat: bnorm=',bnorm
     stop
  end if

end subroutine make_mat

subroutine scaling(xx,yy,iopt)

  !use ctes
  !use running
  use gmres, only:nl,nf,isol_mode,icTp,icsx,icsy,icsz,ishifty
  use precond, only:filter

  implicit none

  real(8),dimension(nl) :: xx,yy
  integer iopt

  ! /* scaling by yy=M*dxp */
  if (iopt.eq.-1) then
     yy(1:nl)=xx(1:nl)*filter(1:nl)
  elseif (iopt.eq.1) then
     yy(1:nl)=xx(1:nl)/filter(1:nl)
  endif 
 
  if ((ishifty.eq.1).and.(isol_mode.eq.3)) then

  endif

  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then

  end if

  if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then

  end if

end subroutine scaling

subroutine set_premat_inner

  !use ctes
  !use running
  use gmres, only:isol_mode,icTp,icsx,icsy,icsz,ishifty
  use precond
  implicit none

  mpcx_in=mpcx
  mpcl_vor=mpc_vor
  mpcl_phi=mpc_phi
  mpcl_u00=mpc_u00
  mpcl_w00=mpc_w00
  mpcl_vorphi=mpc_vorphi
  mpcl_phivor=mpc_phivor

  if (ishifty.eq.1) then
     mpcsy0_in = mpcsy0    
     mpcsyf_in = mpcsyf
  endif
  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
     mpcsx0_in = mpcsx0    
     mpcsxf_in = mpcsxf
     mpcsz0_in = mpcsz0    
     mpcszf_in = mpcszf
  end if
  if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then 
     mpctp0_in = mpctp0    
     mpctpf_in = mpctpf
  end if

end subroutine set_premat_inner

subroutine prematvec(xx,yy,chwk,ct,iopt)

  use ctes
  use running
  use gmres, only:nl,nf,isol_mode,icTp,icsx,icsy,icsz,ishifty
  use precond

  implicit none
  include "mpif.h"

  ! /* iopt=1, for outer newton loop */ 
  ! /* iopt=0, for inner GMRes loop  */

  integer iopt
  character*1 ct  ! 't', multiply the transposed M
  real(8),dimension(nl) :: xx,yy
  real(8),dimension(buffsize) :: chwk
  real(8) fac,dp

  ! /* precondition by yy=M*xx */
  !write(*,*) myid,'prematvec'

  ! reset tmp phisycal field arrays
  vorp=0.d0; phip=0.d0
  u00p=0.d0; w00p=0.d0
  call backthrowing(vorp,phip,u00p,w00p,chwk,xx,-1)
  vorq=0.d0; phiq=0.d0
  u00q=0.d0; w00q=0.d0

  if (iopt.eq.1) then ! outer 
     !yy(1:nl)=xx(1:nf)*mpcx(1:nf)
    
     call factorization(vorp,phip,u00p,w00p,vorq,phiq,u00q,w00q, &
          mpc_vor,mpc_phi,mpc_u00,mpc_w00,mpc_vorphi,mpc_phivor,nband,ct)

  elseif (iopt.eq.0) then
     !yy(1:nl)=xx(1:nf)*mpcx_in(1:nf)

     call factorization(vorp,phip,u00p,w00p,vorq,phiq,u00q,w00q, &
          mpcl_vor,mpcl_phi,mpcl_u00,mpcl_w00,mpcl_vorphi,mpcl_phivor,nband,ct)

  end if
  
  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then

  end if

  if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
     !fac=xx(nl)
     !yy(icTp)=fac
     !if (ct.eq.'n') then
     !   yy(1:nf)=yy(1:nf) + fac*mpctpf(1:nf)
     !   call dotp(nf,mpctp0,xx,dp)
     !   yy(icTp)=dp
     !elseif (ct.eq.'t') then
     !   yy(1:nf)=yy(1:nf) + fac*mpctp0(1:nf)
     !   call dotp(nf,mpctpf,xx,dp)
     !   yy(icTp)=dp
     !endif
  end if

  call throwing(vorq,phiq,u00q,w00q,chwk,yy,-1)

end subroutine prematvec


subroutine factorization(vor,phi,u00,w00,vortmp,phitmp,u00tmp,w00tmp, &
          pc_vor,pc_phi,pc_u00,pc_w00,pc_vorphi,pc_phivor,nfac,ct)

  use ctes
  use running

  integer nfac,j
  character*1 ct  ! 't', multiply the transposed M

  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vor, phi
  real(8),dimension(0:my1) :: u00,w00
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vortmp, phitmp
  real(8),dimension(0:my1) :: u00tmp,w00tmp

  real(8),dimension(nfac,0:2*my-1,0:mx1,kb:ke) :: pc_vor, pc_phi
  real(8),dimension(nfac,0:my1) :: pc_u00, pc_w00

  real(8),dimension(nfac,0:2*my-1,0:mx1,kb:ke) :: pc_vorphi, pc_phivor

  if (nfac.eq.0) then
     ! only diagonal factorization
     if (ct.eq.'n') then
        do j=0,2*my-1,2
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:)
        enddo
        do j=1,2*my-1,2 
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:)
        enddo
     elseif (ct.eq.'t') then
        do j=0,2*my-1,2
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:)
        enddo
        do j=1,2*my-1,2 
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:)
        enddo
     endif
     u00tmp = pc_u00(1,:)*u00
     w00tmp = pc_w00(1,:)*w00

  elseif (nfac.eq.1) then

     if (ct.eq.'n') then
        do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) + pc_vorphi(1,j,:,:)*phi(j+1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) + pc_phivor(1,j,:,:)*vor(j+1,:,:)
        enddo
        do j=1,2*my-1,2 ! Im(oy) ~ Re(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) + pc_vorphi(1,j,:,:)*phi(j-1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) + pc_phivor(1,j,:,:)*vor(j-1,:,:)
        enddo
     elseif (ct.eq.'t') then
        do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) + pc_phivor(1,j,:,:)*phi(j+1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) + pc_vorphi(1,j,:,:)*vor(j+1,:,:)
        enddo
        do j=1,2*my-1,2 ! Im(oy) ~ Re(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) + pc_phivor(1,j,:,:)*phi(j-1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) + pc_vorphi(1,j,:,:)*vor(j-1,:,:)
        enddo
     endif
     u00tmp = pc_u00(1,:)*u00
     w00tmp = pc_w00(1,:)*w00

  elseif (nfac.eq.2) then

     if (ct.eq.'n') then
        do j=0,2*my-1,2 ! Re(oy) ~ Re(oy) + Im(oy) + Im(v) 
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) & 
                &        + pc_vor(2,j,:,:)*vor(j+1,:,:) + pc_vorphi(1,j,:,:)*phi(j+1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) & 
                &        + pc_phi(2,j,:,:)*phi(j+1,:,:) + pc_phivor(1,j,:,:)*vor(j+1,:,:)
        enddo
        do j=1,2*my-1,2 ! Im(oy) ~ Re(oy) + Im(oy) + Re(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) & 
                &        + pc_vor(2,j,:,:)*vor(j-1,:,:) + pc_vorphi(1,j,:,:)*phi(j-1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) & 
                &        + pc_phi(2,j,:,:)*phi(j-1,:,:) + pc_phivor(1,j,:,:)*vor(j-1,:,:)
        enddo
     elseif (ct.eq.'t') then
        do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) & 
                &        + pc_vor(2,j,:,:)*vor(j+1,:,:) + pc_phivor(1,j,:,:)*phi(j+1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) & 
                &        + pc_phi(2,j,:,:)*phi(j+1,:,:) + pc_vorphi(1,j,:,:)*vor(j+1,:,:)
        enddo
        do j=1,2*my-1,2 ! Im(oy) ~ Re(v)
           vortmp(j,:,:) = pc_vor(1,j,:,:)*vor(j,:,:) & 
                &        + pc_vor(2,j,:,:)*vor(j-1,:,:) + pc_phivor(1,j,:,:)*phi(j-1,:,:)
           phitmp(j,:,:) = pc_phi(1,j,:,:)*phi(j,:,:) & 
                &        + pc_phi(2,j,:,:)*phi(j-1,:,:) + pc_vorphi(1,j,:,:)*vor(j-1,:,:)
        enddo
     endif
     u00tmp = pc_u00(1,:)*u00
     w00tmp = pc_w00(1,:)*w00

  elseif (nfac.eq.3) then ! nband
     ! now ignoring the edges j=0,1 and j=2*my-2, 2*my-1
     if (ct.eq.'n') then
        do j=1,2*my-2
           vortmp(j,:,:)= pc_vor(1,j,:,:)*vor(j,:,:) + & 
                &      pc_vor(2,j,:,:)*vor(j+1,:,:) + pc_vor(3,j,:,:)*vor(j-1,:,:) + &
                &      pc_vorphi(1,j,:,:)*phi(j,:,:) + &
                &      pc_vorphi(2,j,:,:)*phi(j+1,:,:) + pc_vorphi(3,j,:,:)*phi(j-1,:,:)
        end do
        do j=1,2*my-2
           phitmp(j,:,:)= pc_phi(1,j,:,:)*phi(j,:,:) + & 
                &      pc_phi(2,j,:,:)*phi(j+1,:,:) + pc_phi(3,j,:,:)*phi(j-1,:,:) + & 
                &      pc_phivor(1,j,:,:)*vor(j,:,:) + &
                &      pc_phivor(2,j,:,:)*vor(j+1,:,:) + pc_phivor(3,j,:,:)*vor(j-1,:,:)  
        end do
        do j=1,my1-1
           u00tmp(j)= pc_u00(1,j)*u00(j) + & 
                &  pc_u00(2,j)*u00(j+1) + pc_u00(3,j)*u00(j-1) 
           w00tmp(j)= pc_w00(1,j)*w00(j) + & 
                &  pc_w00(2,j)*w00(j+1) + pc_w00(3,j)*w00(j-1)  
        end do

     elseif (ct.eq.'t') then
        do j=1,2*my-2
           vortmp(j,:,:)= pc_vor(1,j,:,:)*vor(j,:,:) + & 
                &      pc_vor(3,j,:,:)*vor(j+1,:,:) + pc_vor(2,j,:,:)*vor(j-1,:,:) + &
                &      pc_phivor(1,j,:,:)*phi(j,:,:) + &
                &      pc_phivor(3,j,:,:)*phi(j+1,:,:) + pc_phivor(2,j,:,:)*phi(j-1,:,:)
        end do
        do j=1,2*my-2
           phitmp(j,:,:)= pc_phi(1,j,:,:)*phi(j,:,:) + & 
                &      pc_phi(3,j,:,:)*phi(j+1,:,:) + pc_phi(2,j,:,:)*phi(j-1,:,:) + &
                &      pc_vorphi(1,j,:,:)*vor(j,:,:) + &
                &      pc_vorphi(3,j,:,:)*vor(j+1,:,:) + pc_vorphi(2,j,:,:)*vor(j-1,:,:)
        end do
        do j=1,my1-1
           u00tmp(j)= pc_u00(1,j)*u00(j) + & 
                &  pc_u00(3,j)*u00(j+1) + pc_u00(2,j)*u00(j-1) 
           w00tmp(j)= pc_w00(1,j)*w00(j) + & 
                &  pc_w00(3,j)*w00(j+1) + pc_w00(2,j)*w00(j-1)  
        end do
     end if
  endif

end subroutine factorization

subroutine update_mpc(vorp,phip,u00p,w00p, & 
     &                vorq,phiq,u00q,w00q, &
     &                pc_vor,pc_phi,pc_u00,pc_w00, & 
     &                pc_vorphi,pc_phivor, & 
     &                sca,nfac)

  !Mdf(vorp,...)*MTdx(vorq,...)
  ! this should be consistent with factorization...
  
  use ctes
  use running

  integer nfac,j
  real(8) sca
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vorp, phip
  real(8),dimension(0:my1) :: u00p, w00p
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vorq, phiq
  real(8),dimension(0:my1) :: u00q, w00q

  real(8),dimension(nfac,0:2*my-1,0:mx1,kb:ke) :: pc_vor, pc_phi
  real(8),dimension(nfac,0:my1) :: pc_u00, pc_w00
  
  real(8),dimension(nfac,0:2*my-1,0:mx1,kb:ke) :: pc_vorphi, pc_phivor

  ! diagonal part is done previously in mpcx

  do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
     pc_vorphi(1,j,:,:) = pc_vorphi(1,j,:,:) + vorp(j,:,:)*phiq(j+1,:,:)/sca
     pc_phivor(1,j,:,:) = pc_phivor(1,j,:,:) + phip(j,:,:)*vorq(j+1,:,:)/sca
  enddo
  do j=1,2*my-1,2 ! Im(oy) ~ Re(v)
     pc_vorphi(1,j,:,:) = pc_vorphi(1,j,:,:) + vorp(j,:,:)*phiq(j-1,:,:)/sca
     pc_phivor(1,j,:,:) = pc_phivor(1,j,:,:) + phip(j,:,:)*vorq(j-1,:,:)/sca
  enddo

  if (nfac.eq.2) then 
     do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
        pc_vor(2,j,:,:) = pc_vor(2,j,:,:) + vorp(j,:,:)*vorq(j+1,:,:)/sca
        pc_phi(2,j,:,:) = pc_phi(2,j,:,:) + phip(j,:,:)*phiq(j+1,:,:)/sca
     enddo
     do j=0,2*my-1,2 ! Re(oy) ~ Im(v)
        pc_vor(2,j,:,:) = pc_vor(2,j,:,:) + vorp(j,:,:)*vorq(j-1,:,:)/sca
        pc_phi(2,j,:,:) = pc_phi(2,j,:,:) + phip(j,:,:)*phiq(j-1,:,:)/sca
     enddo
  elseif (nfac.eq.3) then ! nband

     do j=1,2*my-2
        pc_vor(2,j,:,:)= pc_vor(2,j,:,:) + vorp(j,:,:)*vorq(j+1,:,:)/sca        
        pc_vor(3,j,:,:)= pc_vor(3,j,:,:) + vorp(j,:,:)*vorq(j-1,:,:)/sca        
     end do

     do j=1,2*my-2
        pc_phi(2,j,:,:)= pc_phi(2,j,:,:) + phip(j,:,:)*phiq(j+1,:,:)/sca        
        pc_phi(3,j,:,:)= pc_phi(3,j,:,:) + phip(j,:,:)*phiq(j-1,:,:)/sca    
     end do
     do j=1,2*my-2
        pc_vorphi(1,j,:,:)= pc_vorphi(1,j,:,:) + vorp(j,:,:)*phiq(j,:,:)/sca        
        pc_vorphi(2,j,:,:)= pc_vorphi(2,j,:,:) + vorp(j,:,:)*phiq(j+1,:,:)/sca        
        pc_vorphi(3,j,:,:)= pc_vorphi(3,j,:,:) + vorp(j,:,:)*phiq(j-1,:,:)/sca        
     end do
     do j=1,2*my-2
        pc_phivor(1,j,:,:)= pc_phivor(1,j,:,:) + phip(j,:,:)*vorq(j,:,:)/sca        
        pc_phivor(2,j,:,:)= pc_phivor(2,j,:,:) + phip(j,:,:)*vorq(j+1,:,:)/sca        
        pc_phivor(3,j,:,:)= pc_phivor(3,j,:,:) + phip(j,:,:)*vorq(j-1,:,:)/sca    
     end do
     do j=1,my1-1
        pc_u00(2,j)= pc_u00(2,j) + u00p(j)*u00q(j+1)/sca 
        pc_u00(3,j)= pc_u00(3,j) + u00p(j)*u00q(j-1)/sca 

        pc_w00(2,j)= pc_w00(2,j) + w00p(j)*w00q(j+1)/sca 
        pc_w00(3,j)= pc_w00(3,j) + w00p(j)*w00q(j-1)/sca 
     end do
  endif

end subroutine update_mpc

subroutine update_premat(dxk,dfk,chwk,iopt,isave)

  use ctes
  use running
  use precond
  use gmres, only: para_newton,nl,nf,icTp,icsx,icsy,icsz,dxdt0norm,isol_mode,ishifty,iskip_save

  implicit none
  include "mpif.h"

  integer iopt,isave
  real(8) sca,sca2
  real(8),dimension(buffsize) :: chwk
  real(8),dimension(nl) :: dfk, dxk
  character*4 ext1
  character*128 fname


  ! update preconditionar
  if (iopt.eq.1) then ! 'outer' 
     
     if(myid.eq.0) write(*,*) 'update outer preconditioner'
     mpcx = mpcx_in
     
     mpc_vor=mpcl_vor
     mpc_phi=mpcl_phi
     mpc_u00=mpcl_u00
     mpc_w00=mpcl_w00
     mpc_vorphi=mpcl_vorphi
     mpc_phivor=mpcl_phivor

     if (ishifty.eq.1) then
        mpcsy0 = mpcsy0_in
        mpcsyf = mpcsyf_in
     endif
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then 
        mpcsx0 = mpcsx0_in
        mpcsxf = mpcsxf_in
        mpcsz0 = mpcsz0_in
        mpcszf = mpcszf_in        
     endif
     if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then 
        mpctp0 = mpctp0_in
        mpctpf = mpctpf_in
     end if

  end if
  !write(*,*) myid,'update innter premat'
  ! update inner preconditioner
  call prematvec(dxk,MTdx,chwk,'t',iopt) !  zz=M'*dxk
  !if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
  !   dfk(icTp) =  dxdt0norm*dxk(icTp)
  !endif
  call dotp(nf,MTdx,dfk-dxk,sca)    ! (M'*dxk)'*dfk
  if (sca.gt.big) then 
     if (myid.eq.0) write(*,*) 'update preconditioner, sca is too big:',sca
     stop
  end if
  call prematvec(dfk-dxk,Mdf,chwk,'n',iopt) ! zz=M*dfk
  
  Mdf=dxk-Mdf;
     
  if (iopt.eq.0) then  ! 'inner'
     
     mpcx_in(1:nf)=mpcx_in(1:nf) + Mdf(1:nf)*MTdx(1:nf)/sca;

     !if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
     !   ! update by dfk(icTp)=|df/dt|*dTp, dxk=dTp
     !   mpctp0_in(1:nf)=mpctp0_in(1:nf) + Mdf(1:nf)*MTdx(icTp)/sca;
     !   mpctpf_in(1:nf)=mpctpf_in(1:nf) + Mdf(icTp)*MTdx(1:nf)/sca;
     !
     !endif
     ! reset tmp phisycal field arrays
     vorp=0.d0; phip=0.d0
     u00p=0.d0; w00p=0.d0
     ! update diagonal line
     call backthrowing(vorp,phip,u00p,w00p,chwk,mpcx_in,-1)
     mpcl_vor(1,:)=vorp
     mpcl_phi(1,:)=phip
     mpcl_u00(1,:)=u00p
     mpcl_w00(1,:)=w00p

     ! reset tmp phisycal field arrays
     vorp=0.d0; phip=0.d0
     u00p=0.d0; w00p=0.d0
     call backthrowing(vorp,phip,u00p,w00p,chwk,Mdf,-1) 
     
     ! reset tmp phisycal field arrays
     vorq=0.d0; phiq=0.d0
     u00q=0.d0; w00q=0.d0      
     call backthrowing(vorq,phiq,u00q,w00q,chwk,MTdx,-1)
     
     call update_mpc(vorp,phip,u00p,w00p,vorq,phiq,u00q,w00q, &
          mpcl_vor,mpcl_phi,mpcl_u00,mpcl_w00, & 
          mpcl_vorphi,mpcl_phivor, & 
          sca,nband)


  elseif (iopt.eq.1) then ! 'outer' 

     if (myid.eq.0) write(*,*) 'update preconditioner, sca:',sca
     mpcx(1:nf)=mpcx(1:nf) + Mdf(1:nf)*MTdx(1:nf)/sca;

     !if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
     !   ! update by dfk(icTp)=|df/dt|*dTp, dxk=dTp
     !   mpctp0(1:nf)=mpctp0(1:nf) + Mdf(1:nf)*MTdx(icTp)/sca;
     !   mpctpf(1:nf)=mpctpf(1:nf) + Mdf(icTp)*MTdx(1:nf)/sca;
     !
     !endif

     ! update diagonal line
     vorp=0.d0; phip=0.d0
     u00p=0.d0; w00p=0.d0
     call backthrowing(vorp,phip,u00p,w00p,chwk,mpcx,-1)
     mpc_vor(1,:)=vorp
     mpc_phi(1,:)=phip
     mpc_u00(1,:)=u00p
     mpc_w00(1,:)=w00p

     vorp=0.d0; phip=0.d0
     u00p=0.d0; w00p=0.d0
     call backthrowing(vorp,phip,u00p,w00p,chwk,Mdf,-1)       
     vorq=0.d0; phiq=0.d0
     u00q=0.d0; w00q=0.d0
     call backthrowing(vorq,phiq,u00q,w00q,chwk,MTdx,-1)
     
     call update_mpc(vorp,phip,u00p,w00p,vorq,phiq,u00q,w00q, &
          mpc_vor,mpc_phi,mpc_u00,mpc_w00, & 
          mpc_vorphi,mpc_phivor, & 
          sca,nband)

     if ((isave.eq.1).and.(iskip_save.eq.0)) then
        vorp=0.d0; phip=0.d0 ! tmp arrays
        u00p=0.d0; w00p=0.d0

        filext='.pre'; id22=id22-1;
        call save_vec(vorp,phip,u00p,w00p,chwk,mpcx) ! count-up id22 in escru()
        ! reset filext
        filext=''; 
        !write(*,*) 'updated preconditioner min, max(mpc)=', amin(mpc(:)),amax(mpc(:))
 
     endif
  
  endif

end subroutine update_premat 


subroutine gmres_k_shear(vor,phi,u00,w00,chwk,itr,err,kth,eps_j)

  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(my) :: u00, w00
  real(8) dp,fac,eps_j,eps,tmp,csb,err,dummy
  real(8) pythag
  real(8),dimension(:),allocatable:: wkqr
  real(8),dimension(nl-nf) :: vvend

  integer i,j,k, kk,itr, itr_k, kth, lwkqr, info, iproc
  integer istat(MPI_STATUS_SIZE),ierr

  integer iosol
  character*128 fname

  iosol=98
  ! set initial inner preconditionar
  if (iuse_mpc.eq.1)  call set_premat_inner
  
  itr_k=1
  do itr = 1, itr_k 
     ! reset gmres arrays
     vv=0.d0; gh=0.d0; hh=0.d0;
     cs=0.d0; ss=0.d0; ee=0.d0;
     yg=0.d0;
     if (iuse_QR_orth.eq.1) then
        QQ=0.d0;  tau=0.d0;
     end if
     ! iteration of gmres(k)
     !write(*,*) 'gmres outer loop is set to 1', ... 
     !                         itr_k,nl, 'ARG:',itr,err,kth,eps_j
     if (iuse_mpc.eq.1) then
        ! /* preconditioning zz=M*vv */
        call prematvec(dxp,zz,chwk,'n',1)
     else
        if (ifilter.eq.1) then
           call scaling(dxp,zz,-1)
        else
           zz(:)=dxp(:)
        endif
     endif
     call matvec_dns_disk3d(vor,phi,u00,w00,chwk,zz,Ax, & 
          eps_j,rnorm)
     ! isol_mode = 0: used only for Arnoldi method 
     !                (see makefile for building arnoldi)
     !           1: equilibrium solution (steady wave)
     !           2: Travelling-wave solution 
     !           3: periodic orbit (not relative)
     !           4: relative periodic orbit in x- and z-dir
     !           5: periodic orbit (without fixing STp)
     !           6: RPO (without fixing STp)
     ! /* Ax = b */
     fac = -1.d0
     call vecadd(nl,bb,fac,Ax,rr) ! r=b-Ax
     ! do not scale this norm
     call norm(nl,rr,r0norm) ! initial residual vector
     if (r0norm.lt.tiny) then 
        write(*,*) myid,' return: r0norm is small', r0norm
     endif
     
     if (myid.eq.0) write(*,*) 'GMRES,initial residual r0norm:',r0norm
     fac=1.d0/r0norm
     call vecscal(nl,rr,fac,vv) ! vv(:,1)=r(:)/r0
     ee(1)=r0norm     

     if(myid.eq.0) then 
        write(*,*) '---------------------------------------------'
        write(*,*) '  Gmres iteration:', itr
        write(*,*) '---------------------------------------------'
     end if


     do j=1,k_gmres  ! dimension of Krylov sub-space

        if (iuse_mpc.eq.1) then
           ! /* preconditioning zz=M*vv */
           call prematvec(vv(1,j),zz,chwk,'n',1)
        else
           if (ifilter.eq.1) then
              call scaling(vv(1,j),zz,-1)
           else
              zz(:)=vv(:,j)
           endif
        endif
        ! /* an Arnoldi process */
        call matvec_dns_disk3d(vor,phi,u00,w00,chwk,zz,Ax, & 
             eps_j,rnorm)
        if (iuse_QR_orth.eq.0) then
           vv(:,j+1)=Ax(:)
           do i=1,j
              call dotp(nl,vv(1,j+1),vv(1,i),dp) !h_ij=(vv_j+1,vv_i)
              hh(i,j)=dp; gh(i,j)=dp;
              call vecadd(nl,vv(1,j+1),-dp,vv(1,i),vv(1,j+1))
           end do
           call norm(nl,vv(1,j+1),rnorm)
           hh(j+1,j)=rnorm; gh(j+1,j)=rnorm;
           fac=1.d0/rnorm
           call vecscal(nl,vv(1,j+1),fac,vv(1,j+1))

        elseif (iuse_QR_orth.eq.1) then
           !zz=0.d0;
           QQ=vv
           QQ(:,j+1)=Ax(:)
           if (myid.eq.0) then
              ! check the optimal size
              allocate(wkqr(12*(k_gmres+1)))
              call DGEQRF(nl, j+1, QQ, nl, tau, wkqr, -1, info)
              lwkqr=wkqr(1)
              if (lwkqr.gt.12*(k_gmres+1)) then
                 deallocate(wkqr)
                 write(*,*) 'DGEQRF optimal lwk=',lwkqr
                 allocate(wkqr(lwkqr))
              end if
              call DGEQRF(nl, j+1, QQ, nl, tau, wkqr, lwkqr, info)
              ! QR factorisation with Householder transformations
              zz(1:j+1) = QQ(1:j+1,j+1)
              !zz(1) = zz(1)*dsign(1d0,QQ(1,1)) ! the QR (old)
              do kk=1,j
                 zz(kk) = zz(kk)*dsign(1d0,QQ(kk,kk)) ! fixed at rev.1387
              end do
              ! decomposition can "flip" any of the previously computed basis, 
              ! (not only the first vector) 
              ! and this information is stored on the diagonal of matrix R,
              ! that is Q after the call to DGEQRF. 
              ! Solution: if Q(1,1)=R(1,1)<0,
              ! then multiply the i-th row of R (stored in variable yy) by -1
              call DORGQR(nl, j+1, j+1, QQ, nl, tau, wkqr, lwkqr, info)
              vv(:,j+1) = QQ(:,j+1)
              if ((nl-nf).gt.0) then
                 vvend=vv((nf+1):nl,j+1)
                 do iproc=1,numerop-1
                    call MPI_SEND(vvend,nl-nf,MPI_REAL8,iproc,iproc, & 
                         &        MPI_COMM_WORLD,ierr) 
                 end do
              end if
              deallocate(wkqr)
              hh(1:j+1,j)=zz(1:j+1) 
              ! store the upper Hessian matrix for eigenvalue computation. 
           else
              if ((nl-nf).gt.0) then
                 call MPI_RECV(vvend,nl-nf,MPI_REAL8,0,MPI_ANY_TAG, & 
                      &        MPI_COMM_WORLD,istat,ierr)
                 vv((nf+1):nl,j+1) = vvend
              end if
           end if

           !write(*,*) myid,'vvend',vvend,nl-nf
           !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

           call MPI_BCAST(hh(1,j),j+1,MPI_DOUBLE_PRECISION,0, & 
                &         MPI_COMM_WORLD,ierr)
           gh(j+1,j)=hh(j+1,j)
        else
           write(*,*) 'orthogonalization error in gmres_k_shear, stop'
           stop
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! /* upper Hessenberg to upper triangle matrix */
        ! /*                 by Givens rotation matrix */
        do i=1,j-1
           tmp=hh(i,j)
           gh(i,j)  =cs(i)*tmp-ss(i)*hh(i+1,j)
           gh(i+1,j)=ss(i)*tmp+cs(i)*hh(i+1,j)
        enddo
        !csb= dsqrt(gh(j,j)*gh(j,j)+hh(j+1,j)*hh(j+1,j)) 
        csb = pythag(gh(j,j),hh(j+1,j))  
        !csb = pythag(gh(j,j),gh(j+1,j)) ! the same with above, but be sure gh(j+1,j)=hh(j+1,j)
        if(csb.gt.tiny) then
           cs(j)  = gh(j,j)/csb
           ss(j)  = -hh(j+1,j)/csb
           gh(j,j)= cs(j)*gh(j,j)-ss(j)*hh(j+1,j)
           ee(j+1)= ss(j)*ee(j)
           ee(j)  = cs(j)*ee(j)
           gh(j+1,j)=0.d0
        else
           write(*,*) myid, 'csb is small or NaN:', csb
           stop
        endif
        if (isnan(csb)) then
             write(*,*) 'NaN is detected, stop'
             stop
        endif
        
        err=ee(j+1)/r0norm
        if(myid.eq.0) write(*,*) 'GMRES,itr,step, err', itr, j, err
        if(( (j.eq.k_gmres) .or. (dabs(err).lt.eps_gmresk) ).and.(j.gt.3)) then 
           ! backward sweep
           yg(j) = ee(j)/gh(j,j)           
           do i = j-1, 1, -1
              yg(i) = ee(i)
              do kk=i+1,j
                 yg(i) = yg(i) - gh(i,kk)*yg(kk)
              enddo
              yg(i) = yg(i)/gh(i,i)              
           enddo
           
           !dxp should not be set to zero (for restart GMRES)?
           call dcopy(nl,0.d0,0,dxp,1) !rev.(990)
           do i=1,j
              ! /* dxp(:) = dxp(:)+yg(i)*vv(:,i) */
              call vecadd(nl,dxp,yg(i),vv(1,i),dxp)
           enddo
           ! 
           ! dxp is preconditioned solution !
           if (iuse_mpc.eq.1) then
              ! /* preconditioning, dxp=M*dxp */
              call prematvec(dxp,zz,chwk,'n',1) ! do not use inplace arguments
              dxp=zz
           elseif (ifilter.eq.1) then
              call scaling(dxp,dxp,-1)
           end if
           if(myid.eq.0) write(*,*) 'now we got dxp from GMRES solution at itr=',itr           
           call norm(nl,dxp,dxpnorm)
           if (iarclength.eq.1) then
              dxp(icre)=dxp(icre)/sca_re
           endif
           if(myid.eq.0) write(*,*) 'dxpnorm = ',dxpnorm
           
           if ((dabs(err).le.eps_gmresk).or.((j.eq.k_gmres).and.(itr.eq.itr_k))) then
              ! bug fixed: for the case of j==k_gmres, hookstep was not computed...
              if(myid.eq.0) write(*,*) 'CONVERGED: GMRES,itr,j,err',itr, j, err

              if (iuse_hookstep.eq.1) then
                 ! /* check |dxp| is within the trust region 'delta' */
                 if(myid.eq.0) write(*,*) ' hookstep: ', 'delta = ',delta
                 !if(dxpnorm.lt.delta) then
                 !   if(myid.eq.0) then
                 !      write(*,*) ' now hookstep is the same with dxp'
                 !      write(*,*) ' and dxp is within trust region:'
                 !      write(*,*) ' delta is set to |dxp|'
                 !   endif
                 !   delta=dxpnorm/5.d0
                 !   ihookstep_equals_newtonstep=1
                 !else
                 if(myid.eq.0) write(*,*) ' try hookstep to get dxH'
                 call hookstep(j,j,chwk,1)
                 !endif
                 if(myid.eq.0) write(*,*) &
                      '     |dxp|= ',dxpnorm,' delta= ', delta
              end if
              !
              if (isol_mode.ge.3) then
                 ! for periodic orbit
                 call dotp(nf,dxp,dxdt0,dp)
                 if(myid.eq.0) write(*,*) ' check: <dxp,dxdt(x0)>= ',dp
                 if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
                    call dotp(nf,dxp,dx0,dp)
                    if(myid.eq.0) write(*,*) ' check: <dxp,dx0>= ',dp
                    call dotp(nf,dxp,dz0,dp)
                    if(myid.eq.0) write(*,*) ' check: <dxp,dz0>= ',dp
                 endif
                 if (ishifty.eq.1) then
                    call dotp(nf,dxp,dy0,dp)
                    if(myid.eq.0) write(*,*) ' check: <dxp,dy0>= ',dp
                 end if
              end if
              kth = j
              ! write down newton increments: dxp 
              exit
           endif
           if(myid.eq.0 )write(*,*) 'gmres_k:restart,itr,j=k_gmres', &
                itr,j,k_gmres
           kth=j
           
           exit ! go out of GMRes loop

        end if
     end do
  end do

  ! write the hessenberg matrix and so on here
  if ((myid.eq.0).and.(iget_hv.eq.1)) then
     write(*,*) 'GMRes is finished, writing hh gh to .sol'!
     iosol=56
     fname=filout(1:index(filout,' ')-1)//'.sol'
     open(iosol,file=fname,status='unknown',form='unformatted')
     write(iosol) mgalx,my,mgalz,numerop
     write(iosol) time_ini,Re,alp,Ly/pi,gam,s,chi,mx,my,mz  
     write(iosol) xoffb_ini,xofft_ini,dummy,dummy
     para_newton(3)=dfloat(kth)
     para_newton(10)=dfloat(nmax)
     write(iosol) para_newton(1:10)
     ! the others are buffer.
     write(iosol) ((hh(i,k),i=1,kth+1),k=1,kth)
     !write(iosol) ((gh(i,k),i=1,kth+1),k=1,kth)
     !write(iosol) ((vv(i,k),i=1,nl),k=1,kth+1)
     !write(iosol) (bb(i),i=1,nl)
     close(iosol)
     !
  end if
  
  return

end subroutine gmres_k_shear

FUNCTION pythag(a,b)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: a,b
  REAL(8) :: pythag
  !Computes (a2 + b2 )**1/2 without destructive underflow or overflow.
  REAL(8) :: absa,absb
  absa=dabs(a)
  absb=dabs(b)
  if (absa > absb) then
     pythag=absa*dsqrt(1.0d0+(absb/absa)**2)
  else
     if (absb == 0.0) then
        pythag=0.d0
     else
        pythag=absb*dsqrt(1.0d0+(absa/absb)**2)
     end if
  end if
END FUNCTION pythag

subroutine throwing(vor,phi,u00,w00,chwk,xc,isca)

  ! the same style with escru
  ! assuming only master has xc

  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  integer i,j,k, index, iproc, ipo, leng, isca

  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vor, phi
  real(8),dimension(buffsize) :: chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(nl) :: xc
  real(8) ubulk,wbulk
  real(8) sca_u,sca_g,sca_p,xmode,zmode

  integer istat(MPI_STATUS_SIZE),ierr
  integer ptr

  real(8),dimension(:,:,:),allocatable :: wk1, wk2
  real(8),dimension(:),allocatable :: vel
  real(8),dimension(:),allocatable :: u00tmp, w00tmp

  ! scaling 
  sca_u=1.d0; sca_g=1.d0; sca_p=1.d0
  if (isca.eq.1) then
     sca_u = sca_u00
     sca_g = sca_vor
     sca_p = sca_phi
     if ((iset_bulk0.eq.1).and.(iuse_ffty.eq.0)) then
        ! set bulk velocity to be zero
        ubulk=sum(u00(0:my1))/dfloat(my)
        wbulk=sum(w00(0:my1))/dfloat(my)
        u00(:)=u00(:)-ubulk
        w00(:)=w00(:)-wbulk
     end if
  end if

  allocate (u00tmp(0:my1+2),w00tmp(0:my1+2)) ! for rfty
  u00tmp=0.d0; w00tmp=0.d0
  !write(*,*) myid,'in throwing', iuse_fou_mode, nl
  u00tmp(0:my1)=u00(0:my1); w00tmp(0:my1)=w00(0:my1)
  if ((iuse_ffty.eq.1).and.(isca.ne.-1)) then
     call rfty(u00tmp,my+2,1,-1) ! phy -> fou
     call rfty(w00tmp,my+2,1,-1) ! phy -> fou
     if (iset_bulk0.eq.1) then
        u00tmp(0:1)=0.d0; w00tmp(0:1)=0.d0;
     end if
  end if

  ! store zero-zero mode by master proc.
  index=0
  if ((iuse_fou_mode.eq.1).or.(iuse_fou_mode.eq.2)) then
     if(myid.eq.0) then
        index=index+1
        xc(index:my) = u00tmp(0:my1)/sca_u  
        index=index+my
        xc(index:2*my) = w00tmp(0:my1)/sca_u  
        index=index+my
        !write(*,*) '      ... u00, w00 thrown'
     else
        index=index+1;
     end if
  else if (iuse_fou_mode.eq.3) then
     ! does not include x0 modes
     index=index+1
  else
     write(*,*) 'Error throwing, iuse_fou_mode', iuse_fou_mode
     stop
  end if

  deallocate(u00tmp,w00tmp)

  allocate(wk1(mx,mz,jb:je),wk2(mx,mz,jb:je))
  wk1=0.d0; wk2=0.d0;
 
  chwk(:)=0.d0; 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !write(*,*) myid,'in throwing: change'
  call chjik2ikj(vor,wk1,chwk,chwk)       ! --- first  vor
  chwk(:)=0.d0; 
  if ((iuse_v.eq.1).and.(isca.ne.-1)) then
     allocate(vel(buffsize))
     vel=0.d0
     call getv(phi,vel,1)
     call chjik2ikj(vel,wk2,chwk,chwk)       ! --- v
     deallocate(vel)
  else
     call chjik2ikj(phi,wk2,chwk,chwk)       ! --- phi 
  endif

  !
  !------the master first writes its stuff,  
  !      then receives from everybody, and writes it
  if (myid.eq.0) then
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
           ! Note:
           ! all the (jend(iproc)-jbeg(iproc)+1) should be the same with 
           ! (or smaller than) je-jb+1
           call MPI_RECV(wk1,leng,MPI_REAL8,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk2,leng,MPI_REAL8,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        endif
        do j=0,jend(iproc)-jbeg(iproc)
           !
           !write(*,*) '      throwing to xc, j=',j
           !write(iout) (wk1(i,j),wk2(i,j), i=1,mx*mz)
           if (iuse_fou_mode.eq.1) then
              do k=1,mz
                 do i=1,mx
                    xc(index)=wk1(i,k,j)/sca_g
                    index=index+1
                    xc(index)=wk2(i,k,j)/sca_p
                    index=index+1
                 end do
              enddo
           elseif (iuse_fou_mode.eq.2) then
              do k=1,mz
                 do i=1,mx 
!              do k=1,4
!                 do i=1,4 
                    !if (((i.eq.1).or.(i.eq.2)).and.((k.eq.1).or.(k.eq.2))) then
                    if (((i.eq.1).and.(k.eq.1)).or.((i.eq.2).and.(k.eq.1))) then
                       ! skip 00-modes of vor and phi
                    elseif ( (i.le.2).and.( k.ge.nz+2 ) ) then 
                       ! remove conj. imaginary part of vor(0,k=nz+1:mz1,y)
                    else

                       if (iuse_us.eq.1) then
                          if ( (i-1)/2.eq.0) then
                             xmode=1.d0*alp
                          else
                             if (mod(i-1,2).eq.0) then
                                xmode=dfloat((i-1)/2)*alp
                             elseif (mod(i-1,2).eq.1) then
                                xmode=dfloat(i/2)*alp
                             end if
                          endif
                          if (icx(k-1).eq.0) then
                             zmode=1.d0*gam
                          else
                             zmode=dfloat(icx(k-1))*gam
                          end if
                          xc(index)=wk1(i,k,j)*zmode/(xmode*xmode+zmode*zmode)/(s*Lz)        
                          
                       else
                          xc(index)=wk1(i,k,j)/sca_g
                       end if                       
                       index=index+1
                       xc(index)=wk2(i,k,j)/sca_p
                       index=index+1
                    end if
                 end do
              enddo
           elseif (iuse_fou_mode.eq.3) then
              do k=1,mz
                 do i=3,mx ! skip x0-modes of vor and phi           
                    xc(index)=wk1(i,k,j)/sca_g
                    index=index+1
                    xc(index)=wk2(i,k,j)/sca_p
                    index=index+1
                 enddo
              end do
           else
              if(myid.eq.0) write(*,*) 'iuse_fou_mode error, stop'
              stop
           end if
        end do

        !write(*,*) 'throwing index,nf,nl,nmax1,nmax:',iproc,index,nf,nl,nmax1,nmax

     enddo
     if (index-1.ne.nf) then
        write(*,*) 'throwing error, index,nf,nl,nmax1,nmax:',index,nf,nl,nmax1,nmax
        stop
     end if
  else  ! --- everybody else sends things to master 

     call MPI_SEND(wk1,mx*mz*mmy,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk2,mx*mz*mmy,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  deallocate(wk1,wk2)
  chwk(:)=0.d0; ! cleaning the work buffer for change (rev.635)

end subroutine throwing

subroutine backthrowing(vor,phi,u00,w00,chwk,xc,isca)
  !
  ! the same style with getfil
  ! assuming only master has xc
  ! restore mean flow which is not included in xc
  ! 
  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr,isca

  integer i,j,k,index
  integer me(3),ntotr,iproc,mym,mxm,mzm,klen,kini1,kini2

  real(8),dimension(mx,mz,jb:je) :: vor, phi
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(buffsize) :: chwk
  real(8),dimension(nl) :: xc 
  real(8) ubulk,wbulk,sca_u,sca_g,sca_p

  real(8),dimension(:,:,:),allocatable :: wk2
  real(8),dimension(:),allocatable :: u00tmp, w00tmp

  ! ---- zero everything ----
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  ! initialize wk (rev.635)
  chwk(:)=0.d0

  me(1)=mx
  me(2)=my
  me(3)=mz

  ! scaling 
  sca_u=1.d0; sca_g=1.d0; sca_p=1.d0
  if (isca.eq.1) then
     sca_u = sca_u00
     sca_g = sca_vor
     sca_p = sca_phi
  end if

  index=0
  if ((iuse_fou_mode.eq.1).or.(iuse_fou_mode.eq.2)) then
     if (myid.eq.0) then     ! ----- this is done by the master
        ! backthrowing zero-zero mode
        !  write(*,*) 'reading 0 modes'
        mym = min(my,me(2))
        
        index=index+1
        u00(0:my1) = xc(index:my)*sca_u  
        index=index+my
        w00(0:my1) = xc(index:2*my)*sca_u 
        index=index+my
        do iproc=1,numerop-1
           call MPI_SEND(u00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
           call MPI_SEND(w00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        enddo
        
     else       ! -------- this is done by all slaves

        call MPI_RECV(u00  ,my,MPI_REAL8,0,MPI_ANY_TAG, &
             MPI_COMM_WORLD,istat,ierr)
        call MPI_RECV(w00  ,my,MPI_REAL8,0,MPI_ANY_TAG, &
             MPI_COMM_WORLD,istat,ierr)

     endif
  elseif (iuse_fou_mode.eq.3) then
     index=index+1
  end if

  ntotr=me(1)*me(3)*2
  allocate ( wk2(2,me(1),me(3)) )
  wk2 = 0.d0

  mxm=mx
  mzm=mz

  klen=(mzm+1)/2
  kini1=me(3) - klen + 1
  kini2=mz    - klen + 1

  if (myid.eq.0) then     ! ----- this is done by the master, first vor

     !write(*,*) 'reading other modes',ntotr
     
     do iproc=0,numerop-1
        do j=jbeg(iproc), jend(iproc)
           !read(iinp) wk2 ! read one fou-fou plane
           wk2=0.d0
           if (iuse_fou_mode.eq.1) then
              call vcopy(ntotr,xc(index),wk2(1,1,1)) !wk2 = xc(index:index+??)
              index=index+ntotr
           elseif (iuse_fou_mode.eq.2) then
              ! took off 00modes
              !call vcopy_00(ntotr-4,xc(index),ntotr,wk2(1,1,1),4)
              !index=index+ntotr-4 ! 00mode * (real,imag) * [2 fields]
              call vcopy_xz(xc(index),wk2(1,1,1),index,j)
           elseif (iuse_fou_mode.eq.3) then
              ! took off x0-modes
              do k=1,mz                 
                 call vcopy((mx-2)*2,xc(index),wk2(1,3,k))
                 index=index+(mx-2)*2
              end do
           end if

           if (iproc.eq.0) then  ! master stores its thing
              vor(1:mxm,1:klen,j)       =wk2(1,1:mxm,1:klen)
              vor(1:mxm,1+kini2:me(3),j)=wk2(1,1:mxm,1+kini1:me(3))
              phi(1:mxm,1:klen,j)       =wk2(2,1:mxm,1:klen)
              phi(1:mxm,1+kini2:me(3),j)=wk2(2,1:mxm,1+kini1:me(3))
           else      ! sends everything else to nodes
              call MPI_SEND(wk2,ntotr,MPI_REAL8,iproc,iproc, & 
                   MPI_COMM_WORLD,ierr)
           endif
        enddo
        !write(*,*) j,'backthrowing: index/nl',index,nl
     enddo

  else          ! -------- the slaves receive higher modes ----------

     do j=jb,je
        call MPI_RECV(wk2,ntotr,MPI_REAL8,0,MPI_ANY_TAG, & 
             MPI_COMM_WORLD,istat,ierr)
        vor(1:mxm,1:klen,j)    =wk2(1,1:mxm,1:klen)
        vor(1:mxm,1+kini2:mz,j)=wk2(1,1:mxm,1+kini1:me(3))      
        phi(1:mxm,1:klen,j)    =wk2(2,1:mxm,1:klen)
        phi(1:mxm,1+kini2:mz,j)=wk2(2,1:mxm,1+kini1:me(3))      
     enddo

  endif

  ! scaling factor (rev.601)
  if (iuse_us.eq.1) then
     vor = vor*sca_p
     phi = phi*sca_p
  else
     vor = vor*sca_g
     phi = phi*sca_p
  end if

  call set_conj(vor,phi,-1) ! put, -1 for no-report mode. 
  ! (Y-cut, plane) recover the conjugate modes of kx=0
  ! 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! reading vb vt is removed (rev260, 2012/Feb/02)
  deallocate (wk2)

  if ((iuse_ffty.eq.1).and.(isca.ne.-1)) then
     allocate( u00tmp(0:my1+2),w00tmp(0:my1+2) )
     u00tmp=0.d0; w00tmp=0.d0
     u00tmp(0:my1)=u00(0:my1)
     w00tmp(0:my1)=w00(0:my1)
     if (iset_bulk0.eq.1) then
        u00tmp(0:1)=0.d0; w00tmp(0:1)=0.d0;
     end if
     call rfty(u00tmp,my+2,1,1) ! phy <- fou
     call rfty(w00tmp,my+2,1,1) ! phy <- fou

     u00(0:my1)=u00tmp(0:my1)
     w00(0:my1)=w00tmp(0:my1)

     deallocate(u00tmp,w00tmp)
  end if
  !
  if (( (isca.eq.1).and.(iset_bulk0.eq.1) ).and.(iuse_ffty.eq.0) ) then
     ! set bulk velocity to be zero
     ubulk=sum(u00(0:my1))/dfloat(my)
     wbulk=sum(w00(0:my1))/dfloat(my)
     u00(:)=u00(:)-ubulk
     w00(:)=w00(:)-wbulk
  endif

  ! ----- before going any further exchange cuts in y with cuts in z!
  chwk(:)=0.d0; 

  call chikj2jik(vor,vor,chwk,chwk)
  call chikj2jik(phi,phi,chwk,chwk) 

  if ((iuse_v.eq.1).and.(isca.ne.-1)) then
     call getv(phi,phi,-1) ! inplace only for getting lap.(v)
  end if
  chwk(:)=0.d0; ! cleaning work buffer (rev.635)

end subroutine backthrowing

subroutine save_vec(vor,phi,u00,w00,chwk,xc)

  ! damp xc after backthrowing
  !
  ! the same style with getfil
  ! assuming only master has xc
  ! restore mean flow which is not included in xc
  ! 

  use ctes
  use running
  use gmres, only : nf,nl

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,index
  integer me(3),ntotr,iproc,mym,mxm,mzm,klen,kini1,kini2

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(my) :: u00, w00
  real(8),dimension(nl) :: xc

  real(8),dimension(:),allocatable :: vorwk, phiwk

  ! initialize work buffer for change..
  chwk=0.d0
  !if(myid.eq.0) write(*,*) myid,'save_vec ...'
  !/* save the temporal solution */
  !write(*,*) myid,'backthrowing'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  call backthrowing(vor,phi,u00,w00,chwk,xc,1)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !write(*,*) myid,'allocating '
  allocate(phiwk(buffsize),vorwk(buffsize))
  phiwk = 0.d0; vorwk = 0.d0;
  !write(*,*) myid,'escru'
  call escru(vor,phi,u00,w00,vorwk,phiwk,chwk)
  
  deallocate(phiwk,vorwk)

end subroutine save_vec


subroutine matvec_dns_disk3d(vor,phi,u00,w00,chwk,dd,cc,eps_j,cnorm)

  use ctes
  use bcs
  use running
  use gmres
  !
  !
  !
  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(nl) :: dd, cc ! nl is (unknowns-1)/numerop 
  real(8) eps,eps_j,fac,dnorm,cnorm,mtime,dp,da,sca
  ! for central difference and 4th-order difference ...
  real(8),dimension(:),allocatable:: ccm,cc2,ccm2
  
  !write(*,*) myid,'check arg.',eps_j,cnorm,iopt,isol_mode,nl
  call norm(nl,dd,dnorm) 

  !if (myid.eq.0) write(*,*) myid,'matvec dnorm: ',dnorm ! this should be unity
  if(dnorm.lt.tiny) then
     if (myid.eq.0) then
        write(*,*) myid,'matvec dnorm is tiny: ', dnorm
        write(*,*) myid,'matvec: eps will be infinity:return' 
     end if
     !return ! rev.1405
  endif

  eps=eps_j/dnorm
  if (iopt_matvec.eq.1) then ! use forward difference
     if (iuse_hej.eq.1) then
        ! see Dennis & Schnabel, P99
        xtmp=x0 + eps*dd;
        hej=xtmp-x0;
        cc=x0+hej;
     else
        call vecadd(nl,x0,eps,dd,cc) ! cc <-- x0+e_p*dd
     end if
     ! write(*,*) myid,'matvec: ... vecadd'
     xofft=xofft0; xoffb=xoffb0; ! rev1395 ! set b.c. for getv
     call backthrowing(vor,phi,u00,w00,chwk,cc,1)
     ! write(*,*) myid,'matvec: ... backthrowing'
     ! ---- restore initial base flow ---- BUG fixed: moved to after shuffling
     if (iuse_fou_mode.eq.3) then
        if(myid.eq.0) write(*,*) ' restore the baseflow from initial condition'     
        call vecadd(buffsize,vor,1.d0,vor_ini,vor)
        call vecadd(buffsize,phi,1.d0,phi_ini,phi)
     end if
     if ((isol_mode.eq.1).or.(isol_mode.eq.2)) then ! fixed point or TWs
        xofft=xofft0; xoffb=xoffb0; 
        call set_dxdt(vor,phi,u00,w00,chwk,cc)
        call vecadd(nf,cc,-1.d0,dxdt0,cc) ! cc <-- f(x0+e_p*d) - f(x0)
     else ! isol_mode==[0,3,4,5,6]
        call calc_dns_disk3d(vor,phi,u00,w00,chwk) ! cc <-- g(x0+e_p*dd;T)
        ! write(*,*) myid,'matvec: ... calc_dns_disk3d'
        if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
           !c      shiftp = shift+eps*d(n) ! do not use this
           call phase_shift_fou(vor,phi,shiftx,shiftz)
           if (ishifty.eq.1) then
           !   extshiftx=s*0.5d0*timep_sh*shifty
           !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
              fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
              call phase_shifty(vor,phi,u00,w00,shifty,fac)
           endif
        end if
        call throwing(vor,phi,u00,w00,chwk,cc,1)
        ! write(*,*) myid,'matvec: ... throwing'
        call vecadd(nf,cc,-1.d0,xf,cc) ! cc <-- g(x0+e_p*d;T) - g(x0;T)
        ! write(*,*) myid,'matvec: ... vecadd'
        if (iuse_mpc.eq.1) then
        !   !   ! update preconditioar           
           call update_premat(eps*dd,cc,chwk,0,0)
        end if
     end if

     fac=1.d0/eps
     call vecscal(nf,cc,fac,cc) ! cc <-- (g(x0+e_p*d;T) - g(x0;T))/|eps|
     ! write(*,*) myid,'matvec: ... vecscal',eps

  elseif ((iopt_matvec.eq.2).or.(iopt_matvec.eq.4)) then ! use high-order central-difference

     allocate(ccm(nl))
     ccm=0.d0
     if (iuse_hej.eq.1) then
        ! see Dennis & Schnabel, P99
        xtmp=x0 + eps*dd;
        hej =xtmp- x0;
        cc  =x0+hej; 
        xtmp=x0 - eps*dd;
        hej =xtmp- x0;
        ccm =x0+hej;
     else
        call vecadd(nf,x0,eps,dd,cc) ! cc <-- x0 + e_p*dd
        call vecadd(nf,x0,-eps,dd,ccm) ! ccm <-- x0 - e_p*dd
     endif
     ! write(*,*) myid,'matvec: ... vecadd'
     xofft=xofft0; xoffb=xoffb0; ! rev1395
     call backthrowing(vor,phi,u00,w00,chwk,cc,1)
     ! write(*,*) myid,'matvec: ... backthrowing'
     ! ---- restore initial base flow ---- BUG fixed: moved to after shuffling
     if (iuse_fou_mode.eq.3) then
        if(myid.eq.0) write(*,*) ' iuse_fou_mode=3 does not support high-order central-difference' 
        stop
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call calc_dns_disk3d(vor,phi,u00,w00,chwk) ! cc <-- f(x0+e_p*dd;T)
     ! write(*,*) myid,'matvec: ... calc_dns_disk3d'
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        !c      shiftp = shift+eps*d(n) ! do not use this
        call phase_shift_fou(vor,phi,shiftx,shiftz)
        if (ishifty.eq.1) then
           !   extshiftx=s*0.5d0*timep_sh*shifty
           !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
           fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
           call phase_shifty(vor,phi,u00,w00,shifty,fac)
        endif
     end if
     call throwing(vor,phi,u00,w00,chwk,cc,1)
     ! write(*,*) myid,'matvec: ... throwing'
    
     ! reset boundary condition (rev1395)
     xofft=xofft0; xoffb=xoffb0; ! rev1395
     call backthrowing(vor,phi,u00,w00,chwk,ccm,1)
     ! write(*,*) myid,'matvec: ... backthrowing'
     ! ---- restore initial base flow ---- BUG fixed: moved to after shuffling
     !if (iuse_fou_mode.eq.3) then
     !   if(myid.eq.0) write(*,*) ' restore the flow from initial condition'     
     !   call vecadd(buffsize,vor,1.d0,vor_ini,vor) !vor=vor+vor_ini; 
     !   call vecadd(buffsize,phi,1.d0,phi_ini,phi) !phi=phi+phi_ini
     !end if
     call calc_dns_disk3d(vor,phi,u00,w00,chwk) ! ccm <-- f(x0-e_p*dd;T)
     ! write(*,*) myid,'matvec: ... calc_dns_disk3d'
     if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
        !    shiftp = shift+eps*d(n) ! do not use this
        call phase_shift_fou(vor,phi,shiftx,shiftz)
        if (ishifty.eq.1) then
        !   extshiftx=s*0.5d0*timep_sh*shifty
        !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
           fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
           call phase_shifty(vor,phi,u00,w00,shifty,fac)
        endif
     end if
     call throwing(vor,phi,u00,w00,chwk,ccm,1)
     ! write(*,*) myid,'matvec: ... throwing'
     if (iopt_matvec.eq.2) then
        call vecadd(nl,cc,-1.d0,ccm,cc) ! cc <-- f(x0+e_p*d;T)-f(x0-e_p*dd;T)

        if (iuse_mpc.eq.1) then
        !   ! update preconditioar
           call update_premat(2.d0*eps*dd,cc,chwk,0,0)           
           !if (myid.eq.0) write(*,*) 'update_inner_preconditioner:',sca,mpc_in(icTp)
        end if

        ! write(*,*) myid,'matvec: ... vecadd'
        fac=0.5d0/eps
        call vecscal(nl,cc,fac,cc) ! cc <--(f(x0+e_p*d;T)-f(x0-e_p*dd;T))/|2*eps|
        ! write(*,*) myid,'matvec: ... vecscal',eps
        deallocate(ccm)
     elseif (iopt_matvec.eq.4) then

        write(*,*) 'I do not support forth-order stop'
        stop

        allocate(cc2(nl),ccm2(nl))
        cc2=0.d0; ccm2=0.d0
        call vecadd(nl,x0,2.d0*eps,dd,cc2) ! cc2 <-- x0 + 2e_p*dd
        call vecadd(nl,x0,-2.d0*eps,dd,ccm2) ! ccm2 <-- x0 - 2e_p*dd        
        ! cc2 <-- f(x0+ 2e_p*dd;T)
        xofft=xofft0; xoffb=xoffb0; ! reset boundary condition, rev1398
        call backthrowing(vor,phi,u00,w00,chwk,cc2,1)
        call calc_dns_disk3d(vor,phi,u00,w00,chwk) 
        if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
           call phase_shift_fou(vor,phi,shiftx,shiftz)
           if (ishifty.eq.1) then
              !   extshiftx=s*0.5d0*timep_sh*shifty
              !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
              fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
              call phase_shifty(vor,phi,u00,w00,shifty,fac)
           endif
        end if
        call throwing(vor,phi,u00,w00,chwk,cc2,1)
        ! ccm2 <-- f(x0 - 2e_p*dd;T)
        xofft=xofft0; xoffb=xoffb0; ! reset boundary condition, rev1398
        call backthrowing(vor,phi,u00,w00,chwk,ccm2,1)
        call calc_dns_disk3d(vor,phi,u00,w00,chwk) 
        if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
           call phase_shift_fou(vor,phi,shiftx,shiftz)
           if (ishifty.eq.1) then
              !   extshiftx=s*0.5d0*timep_sh*shifty              
              !   call phase_shift_fou(vor,phi,extshiftx,0.d0)«®
              fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
              call phase_shifty(vor,phi,u00,w00,shifty,fac)
           endif
        end if
        call throwing(vor,phi,u00,w00,chwk,ccm2,1)

        call vecadd(nl,cc,-1.d0,ccm,cc) ! cc <-- f(x0+e_p*d;T)-f(x0-e_p*dd;T)
        call vecadd(nl,cc2,-1.d0,ccm2,cc2) ! cc2 <-- f(x0+2e_p*d;T)-f(x0-2e_p*dd;T)        
        
        cc = (8.d0*cc - cc2)
        ! write(*,*) myid,'matvec: ... vecadd'
        fac=1d0/12d0/eps
        call vecscal(nl,cc,fac,cc) ! cc <--(   8*(f(x0+e_p*d;T)-f(x0-e_p*dd;T)) - (f(x0+2e_p*d;T)-f(x0-2e_p*dd;T))   )/|12*eps|
        ! write(*,*) myid,'matvec: ... vecscal',eps        
        deallocate(ccm,cc2,ccm2)
     endif
     
  else
     write(*,*) 'error in matvec, iopt_matvec=',iopt_matvec
     stop
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call norm(nl,cc,cnorm)

  ! 0: only for arnoldi_shear.f90
  if (isol_mode.gt.0) then 
     ! 1: equilibrium solution (steady wave, SW)
     ! 2: equilibrium solution (Travelling wave, TW)
     ! 3: periodic orbit (not relative, UPO)
     ! 4: relative periodic orbit (RPO)
     ! 5: periodic orbit without fixing Tp (UPO)
     ! 6: relative periodic orbit without fixing Tp (RPO)
     ! isol_mode.eq.1 => a fixed point (steady wave, SW)
     if (isol_mode.eq.2) then
        ! /* FINDTW, dT -> 0, this mode should be checked again, after the bug-fix of dxdt */
        !c! period dT = 0
        !dsx=d(n)
        !call vecsub(nmax1,cc,dsx,dxf,cc) ! note: d(n) is dshift 
        !c      call norm(nmax1,dd,dnorm) 
        !c-----     /* <dxp,dx0>=0 viswanath(2007) */------
        !call dotp(nmax1,dd,dx0,cc(n))
     end if
     if (isol_mode.ge.3) then
        ! for periodic orbit
        call vecadd(nf,cc,-1.d0,dd,cc) ! dd is already scaled, dd= Ds^-1 \hat{dxp}
        if (ishifty.eq.1) then
           dsy=dd(icsy)
           call vecadd(nf,cc,dsy,-dyf,cc)
           call dotp(nf,dd,dy0,cc(icsy))
        end if
     endif
     if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
        ! for relative shift
        dsx=dd(icsx)
        dsz=dd(icsz)
        ! 2014/May/23, note: phase_shift_fou is done: xf(x-sx) = exp(-sx*Dx)*xf
        !                    with positive sx in x-dir
        ! dxf is the first derivative of xf in x-dir
        call vecadd(nf,cc,dsx,-dxf,cc)         
        call vecadd(nf,cc,dsz,-dzf,cc) 
        !
        ! /* <dxp,dx0>=0 viswanath(2007) */
        call dotp(nf,dd,dx0,cc(icsx))
        call dotp(nf,dd,dz0,cc(icsz))
     endif
     if ((isol_mode.eq.5).or.(isol_mode.eq.6))  then
        ! time period 
        dTp=dd(icTp) ! *sca_time
        call vecadd(nf,cc,dTp,dxdt,cc)
        !
        ! /* <dxp,dxdt(x0)>=0 viswanath(2007) */
        call dotp(nf,dd,dxdt0,cc(icTp))
        ! write(*,*) 'matvec check index of unknown', n-1, n, cc(n-1), cc(n)
     end if
     
     if (iarclength.eq.1) then
        da=dd(icre)/sca_re
        call vecadd(nf,cc,da,dxda,cc)
        ! /* <dxp,dfa>=0 for arclength */
        call dotp(nf,dd,xstan,cc(icre))
     end if

     call norm(nl,cc,cnorm)

  endif
     
end subroutine matvec_dns_disk3d


subroutine calc_dns_disk3d(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use statistics
  use bcs
  use gmres, only: ishifty,shifty,iset_bulk0,isol_mode,iarclength,cpara, & 
       &           xofft0,xoffb0  !,xoffb_ini,xofft_ini
  use LES,only:iuse_LES,Cles,Deltag
  
  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) ::  vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8) end_time, y0, ubulk, wbulk, vbulk_save
  integer iproc, ierr
  real(8),dimension(:),allocatable :: hv, hg, vorwk, phiwk, dvordy
  ! --- set buffers for LES --- (txy00, tzy00 are in cross.f90)
  real(8),dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2

  real(8),dimension(0:numerop-1) :: chkre, chkalp, chkgam, chkLy, chks, chky0

  !  ---  allocate rest of buffers --- 
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0; 
  allocate( phiwk(buffsize),vorwk(buffsize),dvordy(buffsize) )
  phiwk = 0.d0; vorwk = 0.d0; dvordy = 0.d0;
  ! allocate redundant buffers for DNS/LES 
  allocate(rhvc1(buffsize),rhvc2(buffsize))
  allocate(rhgc1(buffsize),rhgc2(buffsize))
  rhvc1= 0.d0; rhvc2= 0.d0 
  rhgc1= 0.d0; rhgc2= 0.d0  

  !
  !if ((iarclength.eq.1).and.(trim(cpara).eq.'Axz')) then
  !   call check_boundary_condtion ! not impelemented..
  !else
  !end if
  ! reset boundary condition (this affects on Arclength method along alp...)
  !if (abs(xoffb0-xoffb).gt.tiny) then
  !   if (myid.eq.0) then
  !      
  !      write(*,*) 'checking     time xoff:',time,xoffb,xofft
  !      write(*,*) 'checking_0   time xoff:',time_ini,xoffb0,xofft0
  !      !write(*,*) 'checking_ini time xoff:',time_ini,xoffb_ini,xofft_ini
  !   end if
  !end if
  time = time_ini
  if ((abs(iarclength).eq.1).and.(trim(cpara).eq.'Axz')) then
     !xoffb = xoffb0 ! do not set ! rev.1400
     !xofft = xofft0     
  else
     ! just to safety xoffb_ini should be the same xoffb0 except for above if-statement.
     xoffb = xoffb_ini 
     xofft = xofft_ini
  end if


  if (abs(xoffb+xofft).gt.0.d0) then
     write(*,*) 'boundary condition error:STOP, xoffb,xofft',xoffb,xofft
     stop
  endif

  ! reset statistics
  if (nohist.eq.0) then ! store statistics 
     nacum = 0; stats = 0.0; nacumsp = 0
     sp = 0.d0; spl = 0.d0
     stime0 = time
     stats2 = 0.d0; total_dt=0.d0; ! added (rev.852)
     sp2 = 0.d0; spl2 = 0.d0; ! added (rev.853)

  end if

  y0=y(0)
  ! check all proc. have the same paramters [re,alp,gam,Ly]
  call MPI_GATHER(re,1,MPI_REAL8,chkre,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(alp,1,MPI_REAL8,chkalp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_GATHER(gam,1,MPI_REAL8,chkgam,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_GATHER(Ly,1,MPI_REAL8,chkLy,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(y0,1,MPI_REAL8,chky0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr) ! rev.819
  call MPI_GATHER(s,1,MPI_REAL8,chks,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  

  if (myid.eq.0) then
     do iproc=1,numerop-1
        if (abs(re-chkre(iproc)).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [re] for all proc. to be the same.'
           stop
        endif
        if (abs(alp-chkalp(iproc)).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [alp] for all proc. to be the same.'
           stop
        endif
        if (abs(gam-chkgam(iproc)).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [gam] for all proc. to be the same.'
           stop
        endif
        if (abs(Ly-chkLy(iproc)).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [Ly] for all proc. to be the same.'
           stop
        endif
        if (abs(Ly+chky0(iproc)*2.d0).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [y] for all proc. to be the same, or updated after Ly-modification'
           write(*,*) iproc,': Ly,y0',Ly,y0,chky0(iproc)
           stop
        endif
        if (abs(s-chks(iproc)).gt.tiny) then
           write(*,*) ' ERR dns: CHECK [s] for all proc. to be the same.'
           stop
        endif
     end do
  endif

  call check_parameter(damp_aa) 
  call check_parameter(damp_up)
  call check_parameter(force_roll)  
  call check_parameter(xforce)  
  call check_parameter(zforce)

  call check_parameter(vbulk)
  call check_parameter(Cles)

  if ( (iuse_LES.eq.1).and.(abs(Deltag-(dx*dy*dz)**(1.d0/3.d0) ).gt.tiny) ) then
     write(*,*) 'calc_dns (LES): Deltag is not updated in newton arclength method!! stop'
     stop
  endif

  ! set bulk velocity to be zero
  ubulk=sum(u00(0:my1))/dfloat(my)
  wbulk=sum(w00(0:my1))/dfloat(my)
  if (iset_bulk0.eq.1) then
     u00(:)=u00(:)-ubulk
     w00(:)=w00(:)-wbulk
     ubulk=sum(u00(0:my1))/dfloat(my)
     wbulk=sum(w00(0:my1))/dfloat(my)
  end if
  if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) 'calc_dns start: ubulk,wbulk',ubulk,wbulk
  ! 
  if (ishifty.eq.1) then

     !vbulk = -shifty/timep_sh
     !if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) 'calc_dns with: vbulk=',vbulk

     !if (isol_mode.eq.3) then
     !   iadd_force=-1
     !   xforce=  -s*vbulk
     !   if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) 'calc_dns with: xforce=',xforce
     !end if
  endif
  !
  !write(*,*) 'calc dns: time=',time, time_ini, xoffb_ini, xofft_ini
  !if (myid.eq.0) write(*,*) '--- start DNS ---'
  call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,&
       rhvc1,rhvc2,rhgc1,rhgc2,chwk) 

  if (ishifty.eq.1) then ! check mean velocity again if the system has an inertial force
     !if (iadd_force.ne.0) then
     ubulk=sum(u00(0:my1))/dfloat(my)
     wbulk=sum(w00(0:my1))/dfloat(my)
     if ((myid.eq.0).and.(nohist.eq.0)) write(*,*) 'calc_dns finished: ubulk,wbulk',ubulk,wbulk
     !end if
  endif

  !if(myid.eq.0) write(*,*) '--- calc_dns_disk3d finished ...'
  ! --- deallocate ---  
  deallocate(rhvc1,rhvc2,rhgc1,rhgc2)
  deallocate(phiwk,vorwk,dvordy)
  deallocate(hv,hg)

end subroutine calc_dns_disk3d


subroutine set_dxdt(vor,phi,u00,w00,chwk,dxdt)

  use ctes
  use running
  use gmres, only: nf
  use bcs

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(nf) ::  dxdt ! dimension of a field

  real(8),dimension(:),allocatable :: hv, hg, dvordy  
  real(8),dimension(:),allocatable :: rf0u, rf0w
  ! --- set buffers for LES --- (txy00, tzy00 are in cross.f90)
  real(8),dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2

  ! --- allocate rest of buffers --- 
  allocate( hv(buffsize),hg(buffsize),dvordy(buffsize) )
  hv = 0.d0; hg = 0.d0; dvordy = 0.d0;
  allocate (rf0u(0:my1), rf0w(0:my1)) 
  rf0u = 0.d0; rf0w = 0.d0
  ! --- allocate redundant buffers for DNS/LES 
  allocate(rhvc1(buffsize),rhvc2(buffsize))
  allocate(rhgc1(buffsize),rhgc2(buffsize))
  rhvc1= 0.d0; rhvc2= 0.d0 
  rhgc1= 0.d0; rhgc2= 0.d0  
  !
  ! Note: vor, phi, u00, w00 should not be broken
  !write(*,*) myid,'calc_NS_rhs'
  call calc_NS_rhs(vor,phi,u00,w00,hv,hg,rf0u,rf0w,dvordy,&
       rhvc1,rhvc2,rhgc1,rhgc2,chwk) 
  ! [hv,hg] includes nu*lap.[vor,phi]
  !write(*,*) myid,'throwing'
  call throwing(hg,hv,rf0u,rf0w,chwk,dxdt,0)
  ! BUG fixed!  hg <-> hv  2013/july/22 sekimoto
  !write(*,*) myid,'deallocate hv hg'
  deallocate(rhvc1,rhvc2,rhgc1,rhgc2)
  deallocate(rf0u,rf0w)
  deallocate(hv,hg,dvordy)

end subroutine set_dxdt


subroutine set_dxda(vor,phi,u00,w00,chwk)

  ! only for arclength=1

  use ctes
  use running
  use bcs
  use gmres
  use LES, only:Cles

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor,phi,chwk
  real(8),dimension(my) :: u00,w00

  real(8) daa, re0, Ly0, alp0, damp_up0, damp_aa0, Cles0
  real(8) xofft_tmp ! ! rev1390, xofft0 ==> xofft_tmp

  real(8),dimension(:), allocatable:: dxdam

  ! the sensitivity for delta re
  ! check the scalign option
  if ((abs(sca_u00-1.d0).gt.tiny).or.(abs(sca_vor-1.d0).gt.tiny).or.(abs(sca_phi-1.d0).gt.tiny)) then
     if (myid.eq.0) write(*,*) 'check if the scaling works or not for arclength:', sca_u00,sca_vor,sca_phi
     !stop
  endif
  time=time_ini
  xofft=xofft0; xoffb=xoffb0; ! rev1395
  call backthrowing(vor,phi,u00,w00,chwk,x0, 1) 
  daa = 1.d-6 ! use eps_jacobi
  if (trim(cpara).eq.'Rez') then
     re0 = re
     re = re + daa
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using re=',re
  elseif (trim(cpara).eq.'Ayz') then
     Ly0 = Ly
     Ly = Ly + daa
     call set_xyz_dy
     timep_sh = 2.d0*pi/alp/Ly/s*timep
     etime = time_ini + timep_sh     
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using Ly,timep_sh=',Ly,timep_sh
  elseif (trim(cpara).eq.'Axz') then
     alp0 = alp
     xofft_tmp = xofft ! save the b.c. [rev.(1389)] ! rev1390, xofft0 ==> xofft_tmp
     alp = alp + daa
     call set_xyz_dy
     timep_sh = 2.d0*pi/alp/Ly/s*timep
     etime = time_ini + timep_sh   
     ! update boundary condition (rev.1389)
     xofft = xofft_tmp*alp0/alp ! rev1390, xofft0 ==> xofft_tmp
     xoffb = -xofft
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using alp,timep_sh=',alp,timep_sh
  elseif (trim(cpara).eq.'damp_up') then
     damp_up0 = damp_up
     damp_up = damp_up + daa
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using damp_up=',damp_up
  elseif (trim(cpara).eq.'damp_aa') then
     damp_aa0 = damp_aa
     damp_aa = damp_aa + daa
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using damp_aa=',damp_aa     
  elseif (trim(cpara).eq.'LES_Cs') then
     daa=daa/10.0d0 ! rev1305 scaled by a rough factor 10
     Cles0 = Cles
     Cles  = Cles + daa 
     if (myid.eq.0) write(*,*) ' arclength: compute dxda using Cles=',Cles     
  else
     write(*,*) myid,'cpara error daa'
     stop
  end if

  nohist=1 ! reset the history output option
  call calc_dns_disk3d(vor,phi,u00,w00,chwk)
  call throwing(vor,phi,u00,w00,chwk,dxda,1)
  if (iopt_matvec.eq.1) then
     dxda(1:nf) = (dxda(1:nf) - xf(1:nf))/daa
  else
     ! 2nd order
     allocate(dxdam(nf))
     dxdam=0.d0;
     if (trim(cpara).eq.'Rez') then
        re = re0 - daa
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using re=',re
     elseif (trim(cpara).eq.'Ayz') then
        Ly = Ly0 - daa
        call set_xyz_dy
        timep_sh = 2.d0*pi/alp/Ly/s*timep
        etime = time_ini + timep_sh   
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using Ly,timep_sh=',Ly,timep_sh
     elseif (trim(cpara).eq.'Axz') then
        alp = alp0 - daa
        call set_xyz_dy
        timep_sh = 2.d0*pi/alp/Ly/s*timep
        etime = time_ini + timep_sh
        ! update boundary condition (rev.1389)
        xofft = xofft_tmp*alp0/alp
        xoffb = -xofft
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using alp,timep_sh=',alp,timep_sh
     elseif (trim(cpara).eq.'damp_up') then
        damp_up = damp_up0 - daa
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using damp_up=',damp_up
     elseif (trim(cpara).eq.'damp_aa') then
        damp_aa = damp_aa0 - daa
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using damp_aa=',damp_aa
     elseif (trim(cpara).eq.'LES_Cs') then
        Cles = Cles0 - daa
        if (myid.eq.0) write(*,*) ' arclength: compute dxda using Cles=',Cles
     else
        if (myid.eq.0) write(*,*) ' cpara error, 2nd order dxdam, stop'
        stop
     end if
     
     time=time_ini
     xofft=xofft0; xoffb=xoffb0; ! rev1395
     call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
     call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     call throwing(vor,phi,u00,w00,chwk,dxdam,1)         
     dxda(1:nf) = (dxda(1:nf) - dxdam(1:nf))/(2.d0*daa) 
     deallocate(dxdam)
  endif

  ! reset the numerical parameter (iarc_length==1)
  if (trim(cpara).eq.'Rez') then
     re = re0
  elseif (trim(cpara).eq.'Ayz') then
     Ly = Ly0
     call set_xyz_dy
     timep_sh = 2.d0*pi/alp/Ly/s*timep
     etime = time_ini + timep_sh   
  elseif (trim(cpara).eq.'Axz') then
     alp = alp0 
     xofft = xofft0 ! recover the b.c. rev.1400
     xoffb = -xofft0
     call set_xyz_dy
     timep_sh = 2.d0*pi/alp/Ly/s*timep
     etime = time_ini + timep_sh   
  elseif (trim(cpara).eq.'damp_up') then
     damp_up = damp_up0 
  elseif (trim(cpara).eq.'damp_aa') then
     damp_aa = damp_aa0 
  elseif (trim(cpara).eq.'LES_Cs') then
     Cles = Cles0
  else
     if (myid.eq.0) write(*,*) ' cpara error, set_dxda, reset, stop '
     stop
  end if

end subroutine set_dxda

subroutine set_dxdt_dns(vor,phi,u00,w00,chwk,dfdt)

  use ctes
  use running
  use bcs
  use gmres

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(nf) :: dfdt
  real(8) ddt

  real(8),dimension(:),allocatable:: hv, hg, vorwk, phiwk, dvordy
  ! -- set buffers for LES --- (txy00, tzy00 are in cross.f90)
  real(8),dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2

  write(*,*) 'set_dxdt_dns: this should be checked again' 
  stop

  xtmp=0.d0;
  call throwing(vor,phi,u00,w00,chwk,xtmp,0) 
  !  ---  allocate rest of buffers --- 
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0; 
  allocate( phiwk(buffsize),vorwk(buffsize),dvordy(buffsize) )
  phiwk = 0.d0; vorwk = 0.d0; dvordy = 0.d0;
  
  ! reset boundary condition
  time = time_ini
  if ((iarclength.eq.1).and.(trim(cpara).eq.'Axz')) then
     xoffb = xoffb0 
     xofft = xofft0     
  else
     ! just for safety, xoffb_ini should be the same xoffb0 except for above if-statement.
     xoffb = xoffb_ini 
     xofft = xofft_ini
  end if
  !
  ddt = 0.001 ! for time derivative
  etime = time_ini + ddt
  ! 
  !
  !write(*,*) 'calc dns: time=',time, time_ini, xoffb_ini, xofft_ini
  if (myid.eq.0) write(*,*) '--- start DNS ---'
  ! allocate redundant buffers for DNS/LES 
  allocate(rhvc1(buffsize),rhvc2(buffsize))
  allocate(rhgc1(buffsize),rhgc2(buffsize))
  rhvc1= 0.d0; rhvc2= 0.d0 
  rhgc1= 0.d0; rhgc2= 0.d0
  call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,&
       rhvc1,rhvc2,rhgc1,rhgc2,chwk) 
  if(myid.eq.0) write(*,*) '--- calc_dns_disk3d finished ...'

  ! reset 
  etime = time_ini + timep_sh

  ! --- deallocate ---  
  deallocate(rhvc1,rhvc2,rhgc1,rhgc2)
  deallocate(phiwk,vorwk,dvordy)
  deallocate(hv,hg)

  call throwing(vor,phi,u00,w00,chwk,dfdt,0) 
  dfdt = (dfdt - xtmp)/ddt

end subroutine set_dxdt_dns

subroutine calc_NS_rhs(vor,phi,u00,w00,hv,hg,rf0u,rf0w,dvordy,&
     rhvc1,rhvc2,rhgc1,rhgc2,chwk)

  ! please refer uvwp() in conv.f90
  ! rf0u,rf0w should be passed from above
  ! Note:vor, phi, u00, w00 should not be broken

  use ctes
  use running
  use bcs
  use LES

  implicit none
  include "mpif.h"

  integer i,j,k
  real(8) rk,rk2,rkstep,fac

  complex(8) shp,shm
  real(8),dimension(0:2*my-1,0:mx1,kb:ke) :: vor, phi, hv ,hg, dvordy
  real(8),dimension(0:my1) ::  u00, w00, rf0u, rf0w
  real(8),dimension(buffsize) :: chwk

  integer iproc,istat(MPI_STATUS_SIZE),ierr

  ! -- LES
  real*8  rhgc1(0:2*my-1,0:mx1,kb:ke),rhgc2(0:2*my-1,0:mx1,kb:ke) 
  real*8  rhvc1(0:2*my-1,0:mx1,kb:ke),rhvc2(0:2*my-1,0:mx1,kb:ke) 


  real(8),dimension(:),allocatable :: buff1, buff2
  real(8),dimension(:,:),allocatable :: u1r, u2r, u3r, o1r, o2r, o3r

  ! tmp 2d arrays
  real(8),dimension(:,:),allocatable :: tmpxzr ! for fouxz, phys2 
  real(8),dimension(:),allocatable :: tmpx ! for vectorization
  real(8),dimension(:),allocatable :: wk1dr, wk1dc ! for compact finite difference

  ! --- add new buffer to save change (take care for the confusing names)---
  allocate (buff1(buffsize),buff2(buffsize))
  buff1 = 0.d0; buff2 = 0.d0;

  ! ------------------- allocate everything ------------  
  allocate (u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz) )
  allocate (o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz) )
  !
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

  ! ---- allocate tmp 2D arrays ----
  allocate (tmpxzr(mgalx+2,mgalz))
  tmpxzr = 0.d0
  allocate (tmpx(mgalx+2))
  tmpx = 0.d0
  allocate( wk1dr(my), wk1dc(2*my))
  wk1dr = 0.d0; wk1dc = 0.d0;

  xwkt = xofft; xwkb = xoffb 
  zwkt = zofft; zwkb = zoffb 

  if (s.gt.0.5d0) then ! s shoud be 0 or 1
     do i=0,mx1
        shb(i)  = cdexp(-xalp(i)*xwkb) ! negative shift
        sht(i)  = cdexp(-xalp(i)*xwkt) ! positive shift
     enddo
  end if

  !   ---   calcula lap v a partir de phi */
  !write(*,*) myid,'lap v in NS_rhs'
  do k=kb,ke
     do i=0,mx1
        rk = gam2(k)+alp2(i)
        shp = sht(i)
        shm = shb(i)
        call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)   ! -- v
     enddo
  enddo

  !write(*,*) myid,'derivc'
  do k=kb,ke          ! --- computes 
     do i=0,mx1
        shp = sht(i)
        shm = shb(i)
        ! --- d(ome_2)/ dy
        call derivyc(vor(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
        ! --- dv/dy 
        call derivyc(hg(0,i,k),hv(0,i,k),shp,shm,1,2,0) ! skip add_shear
     enddo
  enddo

  ! copy phi and vor to buff1, buff2 not to destroy it. (rev.416)
  call dcopy(buffsize,phi,1,buff1,1)
  call dcopy(buffsize,vor,1,buff2,1)
  !
  !write(*,*) myid,'change'
  chwk(:)=0.d0
  call chjik2ikj(buff1,buff1,chwk,chwk)
  call chjik2ikj(buff2,buff2,chwk,chwk) ! ome2c
  call chjik2ikj(hv,hv,chwk,chwk)   ! dv/dy
  call chjik2ikj(hg,hg,chwk,chwk)   ! v
  call chjik2ikj(dvordy,dvordy,chwk,chwk) ! d(ome_2)/ dy

  ! LES buffers
  rhvc1= 0.d0; rhvc2=0.d0
  rhgc1= 0.d0; rhgc2=0.d0

  !(rev.416) buff1:phi, buff2:vor
  rkstep=0 ! dummy 
  !write(*,*) myid,'hvhg'
  call hvhg(buff1,buff2,u00,w00,hv,hg,rf0u,rf0w,dvordy,rkstep, & 
       u1r,u2r,u3r,o1r,o2r,o3r,u1r,u2r,u3r,o1r,o2r,o3r,      &
       tmpxzr,tmpx,wk1dr,rhvc1,rhvc2,rhgc1,rhgc2) !,txy00,tyz00)

  chwk(:)=0.d0

  if (iuse_LES.eq.1) then
     call chikj2jik(rhvc1,rhvc1,chwk,chwk)
     call chikj2jik(rhvc2,rhvc2,chwk,chwk)
     call chikj2jik(rhgc1,rhgc1,chwk,chwk)
     call chikj2jik(rhgc2,rhgc2,chwk,chwk)
  end if

  call chikj2jik(hv,hv,chwk,chwk)
  call chikj2jik(hg,hg,chwk,chwk)
  call chikj2jik(dvordy,dvordy,chwk,chwk)

  !/* rhvc1=-d^2(T11)/dx^2-d^2(T33)/dz^2-2*d^2(T13)/dxdz+(d^2/dx^2+d^2/dz^2)T22 */
  !/* rhvc2= d(T12)/dx+d(T23)/dz                                                */
  !/* rhgc1= -d^2(T33-T11)/dxdz - (d^2/dx^2-d^2/dz^2)*T13                       */
  !/* rhgc2= d(T12)/dz-d(T23)/dx                                                */
  if (iuse_LES.eq.1) then
     ! d( vor - rhvc1)/dy
     dvordy = dvordy - rhvc1
  end if
  !write(*,*) myid,'deriv of dvordy'
  do k=kb,ke
     do i=0,mx1
        shp = shwkt(i) ! should be the shift of new time step ??
        shm = shwkb(i)
        call derivyc(dvordy(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
        if (iuse_LES.eq.1) call derivyc(rhgc2(0,i,k) ,rhgc2(0,i,k) ,shp,shm,1,2,0)
        ! 
     enddo
  enddo

  if (iuse_LES.eq.1) then
     ! (d^2/dx^2+d^2/dz^2-d^2/dy^2)*(rhvc2)
     fac = -1.d0
     do k=kb,ke
        do i=0,mx1
           if ((i.ne.0).or.(k.ne.0)) then
              shp = shwkt(i) 
              shm = shwkb(i)
              rk  = -(gam2(k)+alp2(i))
              call visc3d(rhvc2(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',1)
           end if
        enddo
     enddo
  end if

  ! hv=hv*+d(rhvc1)/dy+(d^2/dx^2+d^2/dz^2-d^2/dy^2)*(rhvc2)
  ! hg=hg*+rhgc1+d(rhgvc2)/dy

  hv = hv - dvordy
  if (iuse_LES.eq.1) hg = hg + rhgc1  + rhgc2

  if (iuse_LES.eq.1) then
    
     shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
     call derivyc(txy00,txy00,shp,shm,1,1,0)
     call derivyc(tyz00,tyz00,shp,shm,1,1,0)
     
     rf0u= rf0u+txy00
     rf0w= rf0w+tyz00
     
  endif

  ! (1-ii) here you can add body force: ++hv, ++hg
  !iadd_force=0
  if (iadd_force.gt.0) then
     ! streamwise-roll force F=(Fy,Fz), div.F = 0
     !  then hv += (d^2/dz^2 + d^2/dy^2)*Fy
     !       hg += 0
     if(myid.eq.0) write(*,*) ' add_force',iadd_force
     call add_force(hv,iadd_force)
  elseif (iadd_force.eq.-1) then
     rf0u = rf0u + xforce
  endif

  !        if (iadd_damping.eq.1) then
  !           call add_damping(phi,vor,u00,w00,hv,hg,rf0u,rf0w,iadd_damping) 
  !        end if
  
  !write(*,*) myid,'visc3d'
  if (( (iadd_visc.eq.1).and.(iuse_LES.eq.1) ).or.(iuse_LES.eq.0)) then
     
     fac = (1.d0/Re)
     do k=kb,ke          ! --- computes  lap.(vor)
        do i=0,mx1
           if ((i.ne.0).or.(k.ne.0)) then
              shp = shwkt(i) ! do not use the sht at t_new=t+c*dt
              shm = shwkb(i)
              rk = gam2(k)+alp2(i)
              call visc3d(vor(0,i,k),hg(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',1)
              call visc3d(phi(0,i,k),hv(0,i,k),wk1dc,shp,shm,2,rk,fac,'+',0)
           end if
        enddo
     enddo
     !
     ! (*) zerozero-mode (all slaves computes 00modes (from rev.379))
     !    (2*-explicit) compute laplacian --> ++[rf0u,rf0w]
     shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
     rk=0.d0
     
     !write(*,*) myid,'visc3d for 00modes'   
     fac = (1.d0/Re)
     call visc3d(u00,rf0u,wk1dr,shp,shm,1,rk,fac,'+',0)
     call visc3d(w00,rf0w,wk1dr,shp,shm,1,rk,fac,'+',0)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
  end if
  !write(*,*) myid,'deallocate'   
  deallocate(wk1dr, wk1dc) ! for compact finite difference
  deallocate(tmpx)
  deallocate(tmpxzr)
  
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

  deallocate(u1r,u2r,u3r)
  deallocate(o1r,o2r,o3r)

  deallocate(buff1,buff2) ! forgotten, 2014/05/13

end subroutine calc_NS_rhs


! ---------------- vector calculus (for parallel?)-----------------------!
subroutine norm(nx,x,xnorm)

  use ctes, only:my,myid
  use running, only:inorm
  use eft

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer i, nx, imax, imin
  real(8),dimension(nx) :: x
  real(8) xnorm, xnorml,xmax,xmin
  real(8),allocatable :: y(:)
  xnorm=0.d0
  xnorml=0.d0

  if (inorm.eq.0) then
     do i=1,nx
        xnorml = xnorml + x(i)*x(i)
     end do
     xnorm = xnorml ! fixed at rev1546
     call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     xnorm = dsqrt(xnorm)

  elseif (inorm.eq.1) then
     ! this is for LES_EQ 
     if (myid.eq.0) then
        do i=1,2*my
           xnorml = xnorml + x(i)*x(i)
        end do
        xnorml=xnorml/dfloat(my) ! mean velocity is mean ! (rev.1433)

        do i=2*my+1,nx ! testing: remove mean velocities from norm (rev.1419)   
           xnorml = xnorml + x(i)*x(i)
        end do
     endif
     xnorm = xnorml
     call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE(xnorml,xnorm,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     xnorm = dsqrt(xnorm)

  elseif (inorm.eq.2) then
     ! /* check maximum point */
     ! very good convergence for test_case_upo2 with eps_gmresk = 1.e-3 
     !
     xmax=0.d0
     xmin=1.e16
     do i=1,nx
        if (xmax<abs(x(i))) then
           xmax=abs(x(i))
           imax=i
        elseif (xmin>abs(x(i))) then
           xmin=abs(x(i))
           imin=i
        endif 
     end do
     !if (myid.eq.0) write(*,*) 'norm: max', xmax, imax 
    
     do i=1,nx
        xnorml = xnorml + (x(i)/xmax)**2
     end do
     xnorm = dsqrt(xnorml)*xmax
     call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  elseif (inorm.eq.3) then
     ! very tight accurate norm, should be equivalent to those of inorm=0 and 2;
     ! and later the convergence rate is slower ...but converged 
     if (myid.eq.0) then
        call DotK(nx,x,x,xnorml) ! K>=3
        !write(*,*) '(DotK):',xnorml
     end if
     xnorm = xnorml
     call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     xnorm = dsqrt(xnorm)

  elseif (inorm.eq.4) then
     ! use only sumK for more accurate summation 
     if (myid.eq.0) then
        allocate(y(nx))
        do i=1,nx
           y(i) = x(i)*x(i)
        end do
        call SumK(nx,y,2,xnorml)
        !write(*,*) '(SumK):',xnorml, 'error=',y(nx)
        deallocate(y)
     end if
     xnorm = xnorml
     call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     xnorm = dsqrt(xnorm)
  endif


end subroutine norm


subroutine norm_i(nx,x,b)

  use ctes, only:my,myid

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer i, nx, imax, imin
  real(8),dimension(nx) :: x,b
  real(8) xnorm, xnorml,xmax,xmin,xi
  xnorm=0.d0
  xnorml=0.d0
  ! /* check maximum point */
  ! very good convergence for test_case_upo2 with eps_gmresk = 1.e-3 
  !
  xmax=0.d0
  xmin=1.e16
  do i=1,nx
     if (dabs(x(i))>0.d0) then
        xi=dabs(b(i))/dabs(x(i))
     else
        xi=0.d0
     end if
     if (xmax<xi) then
        xmax=xi
        imax=i
     elseif (xmin>xi) then
        xmin=xi
        imin=i
     endif
  end do
  if (myid.eq.0) write(*,*) 'norm_i: max', xmax, imax 

end subroutine norm_i


subroutine norm_seq(nx,x,xnorm)
 
  use eft
  
  implicit none

  integer i,nx
  real(8),dimension(nx) :: x
  real(8) xnorm

  xnorm=0.d0
  do i=1,nx
     xnorm = xnorm + x(i)*x(i)
  end do
  !call DotK(nx,x,x,xnorm) ! k>=3
  xnorm = dsqrt(xnorm)

end subroutine norm_seq

subroutine norm2(nx,x,xnorm)
  
  use ctes, only:my,myid
  use eft

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer i,nx
  real(8),dimension(nx) :: x
  real(8) xnorm,xnorml

  xnorm=0.d0
  xnorml=0.d0
  do i=1,nx
     xnorml = xnorml + x(i)*x(i)
  end do
  xnorm = xnorml
  call  MPI_BCAST(xnorm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !call MPI_ALLREDUCE(xnorml,xnorm,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine norm2

subroutine norm2_seq(nx,x,xnorm)

  use eft

  implicit none

  integer i,nx
  real(8),dimension(nx) :: x
  real(8) xnorm

  xnorm=0.d0
  do i=1,nx
     xnorm = xnorm + x(i)*x(i)
  end do

end subroutine norm2_seq


subroutine dotp(nx,x,y,dp)

  use ctes, only:my,myid
  use running, only:inorm
  use eft

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer i,nx
  real(8),dimension(nx) :: x, y
  real(8) dpl,dp
  real(8), allocatable :: z(:)
  
  ! for parallel, but now only master has a state field
  dp=0.d0
  dpl=0.d0

  if (inorm.eq.0) then
     do i=1,nx
        dpl = dpl + x(i)*y(i)
     end do

  elseif (inorm.eq.4) then
     ! use only sumK for more accurate summation 
     if (myid.eq.0) then
        allocate(z(nx))
        do i=1,nx
           z(i) = x(i)*y(i)
        end do
        call SumK(nx,z,3,dpl)
        ! write(*,*) '(SumK dpl):',dpl, 'error=',z(nx)        
        deallocate(z)
     end if

  else
     write(*,*) 'inorm =',inorm,' is not implemented for dotp'
     stop
  end if
  dp=dpl
  call  MPI_BCAST(dp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !call MPI_ALLREDUCE(dpl,dp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine dotp

!----|---------------------------------------------------------------
!         newton tools
!----|---------------------------------------------------------------
subroutine update_newton(vor,phi,u00,w00,chwk)
  ! not updated for arclength and relative periodic etc...
  use ctes
  use bcs
  use running
  use gmres

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  integer idamp
  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8) eps,eps_j,fac,dnorm,cnorm,mtime,dp
  real(8) bnorm0,bnormi
  real(8) re0,s0

  !      tfin0=tfin
  s0 = s
  shiftx0=shiftx
  if (ishifty.eq.1) then 
     shifty0=shifty
     !vbulk0=vbulk
  endif
  shiftz0=shiftz
  bnorm0=bnorm
  bnormi=big
  if (myid.eq.0) then 
     write(*,*)'----------------------------------------------'
     write(*,*)'update newton-step by damped method, bnorm0',bnorm0
  endif
  
  do idamp=1,ndamp
     dampfac=dampfac/damprate
     if(myid.eq.0) write(*,*) '   damping',idamp,'/',ndamp
     ! xf <- x0 ! for saving x0
     call dcopy(nl,x0,1,xf,1) 
     ! xf <- xf + dampfac*dxp 
     call daxpy(nl,dampfac,dxp,1,xf,1)
     if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
        shiftx = shiftx0 + dampfac*dxp(icsx)
        shiftz = shiftz0 + dampfac*dxp(icsz)
        !write(*,*) 'shiftx,z =',shiftx,shiftz 
     end if
     if (ishifty.eq.1) then 
        shifty = shifty0 + dampfac*dxp(icsy)
        !vbulk = vbulk0 - dampfac*dxp(icsy)/timep_sh
     endif
     if ((isol_mode.eq.5).or.(isol_mode.eq.6))  then
        !tfin = tfin0 + dampfac*dxp(icTp)
        !s = s0 - dampfac*dxp(icTp)*s0/timep_sh
     end if
     ! xf -> vor, phi, u00, w00
     call backthrowing(vor,phi,u00,w00,chwk,xf,1)
     if (isol_mode.eq.1) then
        xofft=xofft0; xoffb=xoffb0; 
        call set_dxdt(vor,phi,u00,w00,chwk,dxdt0)
        bb(:) = -dxdt0(:) 
     else
        call calc_dns_disk3d(vor,phi,u00,w00,chwk)
        if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
           call phase_shift_fou(vor,phi,shiftx,shiftz)
           if (ishifty.eq.1) then
           !   extshiftx=s*0.5d0*timep_sh*shifty
           !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
              fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
              call phase_shifty(vor,phi,u00,w00,shifty,fac)
           endif
        end if
        !write(*,*) 'Tfin is fixed to find a seady-wave or TWs only'
         
        ! vor,phi,u00,w00 -> xf
        call  throwing(vor,phi,u00,w00,chwk,xf,1) 
        ! b <- (x0+dampfac*dxp) - xf
        call vecadd(nl,x0,-1.d0,xf,bb)
        call daxpy(nl,dampfac,dxp,1,bb,1) ! since x0 is old one ???
     end if
     if (ishifty.eq.1) bb(icsy)=0.d0
     if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
        bb(icsx)=0.d0
        bb(icsz)=0.d0
     end if
     !if (ifilter.eq.1) then
     !   call scaling(bb,zz,1)
     !   call norm(nl,zz,bnorm)
     !else
        call norm(nl,bb,bnorm)
     !end if
     if(myid.eq.0) write(*,*) 'bnorm check',bnorm
     if(bnorm.lt.0.5*bnorm0) then
        ! great convergence 
        ! use this dxp
        if(myid.eq.0) write(*,*) '  good convergence: reset dampfac'
        dampfac=2.d0
        exit
     elseif (bnorm.ge.bnormi) then
        if(bnorm.ge.bnormi) then
           ! use previous dxp
           if(myid.eq.0) write(*,*) '  keep dampfac'
           dampfac=dampfac*damprate
        else
           if(myid.eq.0) write(*,*) '  keep the damping factor'
        endif
        exit
     elseif (bnorm.gt.bnorm0) then
        if(myid.eq.0) write(*,*) '  bad convergence: try damping newton'  
     endif
     bnormi=bnorm
  enddo
  !
  if(myid.eq.0) write(*,*) &
       'updated the newton correction by damped method fac =', dampfac
  ! 
  call daxpy(nl,dampfac,dxp,1,x0,1) ! note: x0 is now updated
  call backthrowing(vor,phi,u00,w00,chwk,x0,1) 
  !
  !   check bnorm again:DEBUG --> OK
  !call vecadd(nl,x0,-1.d0,xf,bb)
  !b(n)=0.d0
  !call norm(nl,bb,bnorm)
  !if(myid.eq.0) write(*,*) 'DEBUG in update newton; bnorm check',bnorm
  ! 
  if (ishifty.eq.1) then 
     shifty = shifty0 + dampfac*dxp(icsy)
     !vbulk = vbulk0 - dampfac*dxp(icsy)/timep_sh
  endif
  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
     shiftx = shiftx0 + dampfac*dxp(icsx)
     shiftz = shiftz0 + dampfac*dxp(icsz)
  end if
  if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
     !s = s0 - dampfac*dxp(icTp)*s0/timep_sh
  end if
  !
  !     /* set initial guess of GMRES loop at next newton step*/
  call vecscal(nl,dxp,dampfac,dxp) 
  
end subroutine update_newton
!
subroutine solv_mu(k,k_cut,ierr_hook)
  ! /*****************************************************************/
  ! /* solve 1d root finding problem to get mu from delta            */
  ! /* using newton iteration                                        */          
  !    -- see also channelflow.org/programs-octave/findorbit.cpp --
  ! /* See Dennis and Schnabel for this search-over-mu algorithm     */
  ! /*****************************************************************/
  use ctes, only: myid
  use gmres
  use svd
 
  integer i,k,imu, nstep_mu, ichannelflow,k_cut,ierr_hook,ihave_mu
  real(8)  qnorm,ph,phprim,di_mu,dr

  !dmu=0.0d0  !(rev.1148) commented
  dmu=0.001d0 ! (rev.598)
  ichannelflow=0
  ! /* finding mu by newton iteration */
  nstep_mu=200
  ierr_hook=0;
  !if(myid.eq.0) write(*,*) ' finding mu iteratively'
  ihave_mu=0
  do imu=1,nstep_mu
     if (ichannelflow.eq.0) then
        do i=1,k
           if ((hs(i)*hs(i).gt.hs_cut)) then 
              if (hs(i).gt.dmu) then
                 q(i)=tmpvec(i)/(hs(i) + dmu/hs(i))
              else
                 q(i)=hs(i)*tmpvec(i)/(hs(i)*hs(i)+dmu)
              end if
              if (i.gt.k_cut) q(i)=0.d0 ! throwing away some noisy modes.
           elseif ((hs_cut.lt.0.d0).and.(i.eq.k)) then 
              !write(*,*) 'hs_k ',i,hs(i), 'the last one is ignored'
              q(i)=0.d0 ! ignore tiny singular values
           else 
              !if((myid.eq.0).and.(ireport)) write(*,*) 'imu loop: set hs=0:',i
              q(i)=0.d0 ! ignore tiny singular values
           end if
        enddo
        call norm_seq(k,q,qnorm)
        ph=qnorm*qnorm-delta*delta
        phprim=0.d0
        do i=1,k
           phprim=phprim - 2.d0*(q(i)*q(i)/(hs(i)*hs(i)+dmu))
        enddo
        if (delta.gt.deltaMin)then
           dmu=dmu-(qnorm/delta)*(ph/phprim)
        else
           dmu=dmu-ph/phprim
        endif
        
     else ! ichannelflow==1
        ! --- copied from channelflow.org/programs-octave/findorbit.cpp ---
        !     See Dennis and Schnabel for this search-over-mu algorithm
        !     from channelflow.org
        do i=1,k
           q(i)=tmpvec(i)/(hs(i)+dmu);
        enddo
        call norm_seq(k,q,qnorm) !ndxHmu
        !dr = bnorm - bnorm0 !!??
        !dr = rnorm - r0norm !!?? 
        dr = delta
        ph=qnorm*qnorm - dr*dr;
        phprim=0.0;
        do i=1,k
           di_mu=hs(i)+dmu;
           phprim=phprim-2.d0*tmpvec(i)*tmpvec(i)/(di_mu*di_mu*di_mu);
        enddo
        dmu=dmu-(qnorm/dr)*(ph/phprim)
     endif

     if(dmu.lt.0.d0) then
        if(myid.eq.0) write(*,*) &
             &'hookstep: mu =',dmu
        dmu=0.0d0
        !delta=qnorm !(rev.1148) reverted to 0.9
        delta=qnorm*0.9
        !delta=delta*0.5
        if(myid.eq.0) write(*,*) &
             &'hookstep: mu is negative, set mu=0, keep this delta (=qnorm),', delta
!             &'hookstep: mu is negative, set mu=0 (normal GMRES step), and delta=0.9*qnorm :', delta
        if(myid.eq.0) write(*,*) &
             &'delta is applied to |M*dxp|'
        if(myid.eq.0) write(*,*) &
             'compute mu again with new delta'
        nstep_mu = nstep_mu + imu
        ph = phprim+0.1 ! dummy setting do go to imu loop again ! (rev.1148)
        !irecompute_hookstep=1
        exit
     endif
     !
     ! /* Found satisfactory value of mu and thus dxHmu and dxH s.t. */
     ! /* |dxH| < delta */
     if(abs(ph/phprim).lt.1.d-10) then
        ihave_mu=1
        if ((qnorm.lt. delta).or. &
             ((qnorm.lt.(1.0-deltaFuzz)*delta).and. &
             (qnorm.lt.(1.0+deltaFuzz)*delta))) then
           !if(myid.eq.0) write(*,*) 'hookstep: update dmu: ',dmu, &
           !     ' by this delta =',delta
           !if(myid.eq.0) write(*,*) myid,'  imu dmu delta_mu:',imu, dmu, ph/phprim
           !exit             ! exit finding mu
        else
           !if(myid.eq.0) write(*,*) 'hookstep: imu dmu  = ', imu, dmu
           !delta = delta*2.d0
           if (delta.gt.deltaMax) then
              delta = deltaMax
              if(myid.eq.0) write(*,*) 'hookstep: delta reached to max',delta 
           endif
           !if(myid.eq.0) write(*,*) 'hookstep: try wider delta',delta 
           !exit
        endif
     elseif (abs(ph/phprim).lt.1.d-10) then
        ! solve mu more accurately and 1.d-10 is acceptable like above...
        ! rev.1480
        exit
     elseif(imu.eq.nstep_mu)then
        if(myid.eq.0) then 
           write(*,*) &
             'Could not find solution of hookstep optimization eqn Phi(mu)==0'
           write(*,*) &
             'This should not be happen. It indicates an error in the algorithm.'
           write(*,*) 'hookstep:stop'
        end if
        ierr_hook=1
        !stop
     endif               !/* dmu < 0 or not*/
  enddo                  !/* imu loop */
  if (myid.eq.0) then
     write(*,*) 'hookstep report: lambda,delta=', dmu, delta
  end if
  !write(*,*) myid,'DEBUG, check qnorm and delta', qnorm, delta

end subroutine solv_mu
!
subroutine hookstep(k,k_cut,chwk,ini)

  !/***************************************************************/
  !/* locally constrained hookstep using SVD for Hessenberg matrix*/
  !/*   delta  : the trust region of Newton equation              */
  !/*            should be updated to satisfy |dxp| < delta       */
  !/* input:                                                      */
  !/*   k    : the k-th gmres step                                */
  !/*   chwk : work array for change                              */
  !/*   ini  : 1, the first trial of hookstep (to compute SVD)    */
  !/*          0, next trial using precomputed SVD,but new delta  */
  !/* output                                                      */
  !/*   yg : Vq in Viswanath(2009) paper                          */
  !/*                                                             */
  !/***************************************************************/
  use ctes
  use running
  use gmres
  use svd
  use eft

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr,ierr_hook
  integer lwk,info,i,ii,j,jj,k,imin,ini,iproc,k_cut
  integer imu,nstep_mu,ichannelflow

  real(8) time_hook,mem_hook, dp
  real(8),dimension(buffsize) :: chwk
  real(8), allocatable :: xy(:)


  ! lwk=5*(k+1) ! BUG!! 2013/Sep/30 ... !
  !    this bug affects the behavior of the newton convergence, 
  !    if the GMRES iteration is larger and larger
  !    now I can use eps_gmres=1.d-6
  ! lwk should be larger than 3*min(M,N) + max(max(M,N),4*min(M,N)*min(M,N)+3*min(M,N)+max(M,N))

  if (ini.ge.1) then
     
     time_hook=MPI_WTIME()

     allocate(vvb(k+1),tmpy(k))
     vvb=0.d0; tmpy=0.d0
     ! /* vvb=vv*b */
     do jj=1,(k+1)
        vvb(jj)=0.d0
        zz(:)=vv(:,jj)  
        call dotp(nl,zz,bb,dp)
        vvb(jj)=dp 
     enddo

     if (myid.eq.0) then
        ! bug in lapack?
        lwk=5*(k+1)*(k+1) + 10*(k+1)  ! this may be too big
        allocate(tmpvec(k+1),q(k))
        allocate(tmph(k+1,k))
        allocate(hu(k+1,k+1),hs(k+1),hvt(k,k))
        allocate(svdwk(lwk))
        mem_hook=dfloat( 2*(k+1)*numerop + 3*(k+1) + 3*(k+1)*(k+1) + lwk )/1024.d0/1024.d0*8.d0
        write(*,'(i4,a,G9.2)') myid,': hookstep requests (MB):', mem_hook
        ! initialization
        tmpvec=0.d0; q=0.d0;
        tmph=0.d0; hu=0.d0; hs=0.d0; hvt=0.d0;
        svdwk=0.d0;
        ! /* svd of hessenverg matrix, done by all proc.,   */
        ! /* but it is dangerous for a large kmax           */
        do j=1,k
           do i=1,k+1
              tmph(i,j)=hh(i,j)
           enddo
        enddo
        call DGESVD('a','a',k+1,k,tmph,k+1,hs,hu,k+1,hvt,k,svdwk,lwk,info)

        write(*,*) 'DGESVD:',info
        if (info.eq.0) then
           write(*,*) 'DGESVD:   successful exit.'
        elseif (info.lt.0) then
           write(*,*) 'DGESVD:   the i-th argument had an illegal value.'
        elseif (info.gt.0) then
           write(*,*) 'DGESVD:   DBDSDC did not converge, updating process failed.'
        end if

        deallocate(svdwk)

        if (ini.eq.1) then
           hs2max=0.d0; hs2min=1.d6; imin=0
           do i=1,k
              hs2=hs(i)*hs(i)
              if (hs2min.gt.hs2) then
                 hs2min=hs2; imin=i
              end if
           end do
           
           do i=1,k
              if ((hs(i)*hs(i).lt.hs_cut)) then 
                 write(*,*) 'hs_k ',i,hs(i), '<', dsqrt(hs_cut), 'ignored'
              elseif ((hs_cut.lt.0.d0).and.(i.eq.k)) then 
                 write(*,*) 'hs_k ',i,hs(i), ' the last one is ignored'
              else
                 write(*,*) 'hs_k ',i,hs(i) 
              end if
           end do
           !write(*,*) 'hs2_min ',imin,hs2min 
        endif
        !
        ! /*  bh =U'*vv*b, vvb=vv*b */ ! try SumK
        if (inorm.eq.4) then
           allocate(xy(k+1))
           xy=0d0
           do j=1,k+1              
              do i=1,k+1
                 xy(i)=hu(i,j)*vvb(i)
              enddo
              call SumK(k+1,xy,3,tmpvec(j))
              !write(*,*) '(SumK):',tmpvec(j), 'error=',xy(k+1)
           enddo 
           deallocate(xy)
        else
           do j=1,k+1
              tmpvec(j)=0.d0
              do i=1,k+1
                 tmpvec(j)=tmpvec(j)+hu(i,j)*vvb(i)
              enddo
           enddo
        end if
     endif

     time_hook=MPI_WTIME() - time_hook
     if(myid.eq.0) write(*,*) myid,':time_hook:',time_hook
  elseif(ini.eq.-1) then 
     ! deallocate after getting dxH...
     if (myid.eq.0) then
        deallocate(hu,hs,hvt)
        deallocate(tmph)
        deallocate(tmpvec,q)
     endif
     deallocate(vvb,tmpy)
     return
  endif ! only for ini.eq.1
  !
  if (ini.eq.2) return
  !
  if (myid.eq.0) then 
     call solv_mu(k,k_cut,ierr_hook) ! update q, mu, delta
     !
     ! /* y=hv*q; */  ! try SumK, include eft
       if (inorm.eq.4) then
           allocate(xy(k_cut))
           xy=0d0
           do i=1,k              
              do j=1,k_cut
                 xy(j)=hvt(j,i)*q(j)
              enddo
              call SumK(k_cut,xy,3,tmpy(i))
              !write(*,*) '(SumK):',tmpy(i), 'error=',xy(k_cut)
           enddo 
           deallocate(xy)
        else           
           do i=1,k
              tmpy(i)=0.d0
              do j=1,k_cut
                 tmpy(i)=tmpy(i)+hvt(j,i)*q(j)
                 !note the index, hvt is transposed V
              enddo
           enddo
        endif
  endif
  call MPI_BCAST(ierr_hook,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  if (ierr_hook.eq.1) then
     stop
  end if
  call MPI_BCAST(tmpy,k,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(dmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  !write(*,*) myid,'lambda,delta,k=',dmu,delta,k,' ... tmpy' !!,tmpy
  !
  ! /* dxp=vv*y */  ! try SumK 
  if (inorm.eq.4) then
     allocate(xy(k))
     xy=0d0
     do i=1,nl              
        do j=1,k
           xy(j)=vv(i,j)*tmpy(j)
        enddo
        call SumK(k,xy,3,dxp(i))
        !write(*,*) '(SumK):',dxp(i), 'error=',xy(k)
     enddo
     deallocate(xy)
  else
     do i=1,nl
        dxp(i)=0.d0 
        do j=1,k
           dxp(i)=dxp(i)+vv(i,j)*tmpy(j)
        enddo
     enddo
  end if
  if (iuse_mpc.eq.1) then
     ! dxp is preconditioned solution
     ! /* preconditioning, dxp=M*dxp */
     call prematvec(dxp,zz,chwk,'n',1)
     dxp=zz
     !write(*,*) myid,'hookstep check dxp(icTp)',dxp(icTp)
  elseif (ifilter.eq.1) then
     !call norm(nl,dxp,dxpnorm)
     !write(*,*) myid,'hookstep: before scaling, dxpnorm= ',dxpnorm
     call scaling(dxp,dxp,-1)
  end if
  call norm(nl,dxp,dxpnorm)
  if (iarclength.eq.1) then
     dxp(icre) = dxp(icre)/sca_re
  endif

  !write(*,*) myid,'hookstep check dxp(icTp)',dxp(icTp)
  !write(*,*) myid,'dxpnorm after scaling, = ',dxpnorm

end subroutine hookstep

subroutine update_hookstep(vor,phi,u00,w00,chwk,kth,eps_j)

!  /***********************************************************/
!  /*   hookstep update strategy is following:                */
!  /* compute (0),(1),(2),(3)                                 */
!  /* (0) x=x0+dxp                                            */
!  /* (1) actual residual ra(x0+dxp)                          */
!  /* (2) predicted residual from guadratic model of          */
!  /*     rp(x0+dxp)== (1/2)(G(x),G(x))                       */
!  /*                 + (G(x),DG dxp)                         */
!  /* (3) slope of r(x0)                                      */
!  /*     dr/dx == (r(x0+dxp)-r(x0))/|dxp|                    */
!  /*           == (G(x),DG dxp)/|dxp|                        */
!  /*  note: (a,b) is the usual inner product norm2 (self dotp)*/
!  /*                                                         */
!  /*  Then update delta using (0-3)                          */
!  /*  see <www.channelflow.org/programs-octave/findorbit.cpp>*/
!  /*   ==> not used, and now simply update                   */
!  /***********************************************************/
  use ctes
  use running
  use gmres
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00

  integer kth,improvement,irecompute_hookstep,ihook,iuse_quadratic_model,iuse_next
  integer ik,k_cut,icheck_k_cut
  real(8) eps_j
  real(8) bnorm0,rx,rH,rL,rQ,drdx,backup_rH
  real(8) re0,s0,Ly0,timep_sh0
  real(8) deltaMaxLocal,delta_new,backup_delta,Delta_rH,Delta_rL,Delta_rQ
  real(8) Delta_rH_req, Delta_rH_ok, Delta_rH_good 
  real(8) Delta_rQ_acc, Delta_rH_accP, Delta_rH_accM 

  real(8) dlambdaRequiredReduction

  real(8),dimension(:),allocatable:: Adxp,xH_backup

  !tfin0=tfin
  s0=s
  Ly0=Ly
  shiftx0=shiftx
  timep_sh0=timep_sh
  if (ishifty.eq.1) then 
     shifty0=shifty
  endif
  shiftz0=shiftz
  re0 = re
  !bnorm0=bnorm
  !if(myid.eq.0) write(*,*)' ---update newton by hookstep method bnorm0',bnorm0
  
  deltaMaxLocal=deltaMax
  ihavebackup=0
  !if (ifilter.eq.1) then
  !   call scaling(bb,zz,1)
  !   call dotp(nl,zz,zz,rx) ! = bnorm^2
  !else
     call dotp(nl,bb,bb,rx) ! = bnorm^2
  !end if
  rx=rx*0.5d0
  if(myid.eq.0) write(*,*) myid,' ---update newton by hookstep method:, rx0=', rx
     
  xtmp=0.d0; 
  allocate(Adxp(nl),xH_backup(nl))
  Adxp=0.d0; xH_backup=0.d0

  k_cut=kth

  icheck_k_cut=0
  if (icheck_k_cut.eq.1) then
     do ik=1,kth
        k_cut=kth-ik
        ! /* recompute hookstep using updated delta, only by the master proc. */
        call hookstep(kth,k_cut,chwk,0)

        call hookstep_trial(vor,phi,u00,w00,chwk)
        !
        !if (ifilter.eq.1) then
        !   call scaling(btmp,zz,1)
        !   call dotp(nl,zz,zz,rx) ! = bnorm^2
        !else
           call dotp(nl,btmp,btmp,rH) ! (G(x+dxp),G(x+dxp))
        !end if
        rH=0.5d0*rH
        if(myid.eq.0) write(*,*) k_cut,'actual residual rH=',rH
        
     enddo

     stop

  endif

  !nhookstep=0 ! debuging: skip adjusting delta 
  if (nhookstep.eq.0) then
     if(myid.eq.0) write(*,*) 'skip ihook loop, nhookstep=',nhookstep 
  else  
     call hookstep_trial(vor,phi,u00,w00,chwk)
  endif
  !
  ! !if(myid.eq.0) write(*,*) 'skip ihook loop, nhookstep=',nhookstep 
  ! !nhookstep=0
  !
  iuse_next=0;
  do ihook=1,nhookstep

     !! checking hookstep fields
     !  xtmp = x0 + dxp
     !  xoffb=xoffb0 ! rev.1390
     !  xofft=xofft0
     !  call save_vec(vor,phi,u00,w00,chwk,x0)       
     !stop
     if (ihook.gt.1)then
        ! /* recompute hookstep using updated delta */
        if(myid.eq.0) write(*,*) myid,'    hookstep loop, ihook =',ihook
        call hookstep(kth,k_cut,chwk,0)    
        call hookstep_trial(vor,phi,u00,w00,chwk)

     endif
     !write(*,*) myid,'    hookstep loop,  dotp ,ihook=',ihook
     !if (ifilter.eq.1) then
     !   call scaling(btmp,zz,1)
     !   call dotp(nl,zz,zz,rH) ! = bnorm^2
     !else
        call dotp(nl,btmp,btmp,rH) ! (G(x+dxp),G(x+dxp))
     !endif
     rH=0.5d0*rH
     if(myid.eq.0) write(*,*) myid,'actual residual rH=',rH

     Delta_rH=rH-rx

     if (Delta_rH.lt.0.d0) then
        if(myid.eq.0) write(*,*) ' good delta: increase delta by 10% '
        if(dmu.le.0.d0) iuse_next=1 ! normal GMRes step
        if (iuse_next.eq.1) then
           if(myid.eq.0) write(*,*) myid,'i simply use this hookstep '
           irecompute_hookstep=0
        else
           if(myid.eq.0) write(*,*) myid,'i keep this step, try hookstep again'
           irecompute_hookstep=1
           ihavebackup=1; rx=rH; !(rev.1005)
           backup_delta=delta;
           xH_backup=xtmp; ! backup temporarily computed, xf 
        endif
        delta=delta*1.1d0
        if (iuse_mpc.eq.0) irecompute_hookstep=0 ! rev961-mode 
     else
        if (ihavebackup.eq.1) then ! (rev.990)
           if(myid.eq.0) write(*,*) ' i use hookstep with the backuped delta, '
           delta=backup_delta
           call hookstep(kth,k_cut,chwk,0)
           xtmp=xH_backup; 
           irecompute_hookstep=0
        else
           if(myid.eq.0) write(*,*) ' shrink delta to be 0.8*delta, try hookstep again'
           Delta = Delta*0.8d0
           irecompute_hookstep=1
           iuse_next=1;! older rev961-mode 
        end if
     endif

     iuse_quadratic_model=0 
     ! ! scaling is not applied for dotp and norm and matvec...
     ! copied from channelflow.org ... should be checked for its usage...
     if (iuse_quadratic_model.eq.1) then
        !    /* (2) predicted residual from guadratic model of     */
        !    /*     rp(x0+dxp)== (1/2)(G(x),G(x))                  */
        !    /*                 + (G(x),DG dxp)                    */
        !    /*    [DG dxP = 'Adxp', (G(x),DG dxp)='Delta_rL']     */
        !    /*  ref.)                                             */
        !    /*    rL is predicted residual from linear model      */
        !    /*       rL == 1/2 (G(x),G(x)) + (G(x), DG dx)        */
        !    /*    rQ is predicted residual from quadratic model   */
        !    /*       rQ == 1/2 |G(x) + DG dx|^2                  */
        call matvec_dns_disk3d(vor,phi,u00,w00,chwk,dxp,Adxp, &
             eps_j,rnorm)
        !     /* linear model */
        !         call dotp(n,b,b,rx)
        !         rx=rx*0.5d0
        call dotp(nl,bb,Adxp,Delta_rL)
        rL=rx+Delta_rL
        !     /* quadratic model */
        !     /*       rQ == 1/2 |G(x) + DG dx|^2                        */
        !call vecnul(n,btmp)
        btmp=0.d0
        call vecadd(nl,bb,1.d0,Adxp,btmp)
        call dotp(nl,btmp,btmp,rQ)
        rQ=rQ*0.5d0
        Delta_rQ=rQ-rx   
        !     
        !     /* (3) slope of r(x0)                                      */
        !     /*     dr/dx == (r(x0+dxp)-r(x0))/|dxp|                    */
        !     /*           == (G(x),DG dxp)/|dxp|                        */
        drdx=rL/dxpnorm
        !     
        !     /* Try to minimize a quadratic model of the residual  */
        !     /* r(x+dx) = r(x) + r' |dx| + 1/2 r'' |dx|^2          */ 
        dlambda = -0.5d0*drdx*dxpnorm/(rH - rx - drdx*dxpnorm)
        if(myid.eq.0) then
           write(*,*) &
                'hookstep: dxpnorm delta rx rH rL rQ drdx delta_rate;', &
                dxpnorm,delta,rx,rH,rL,rQ,drdx,dlambda
        end if
        if(rL.ge.rx) then
           if(myid.eq.0) then
              write(*,*) 'error : local linear model of residual', &
                   'is increasing, indicating', &
                   'that the solution to the Newton equations', &
                   'is inaccurate... exit' 
              write(*,*) 'please try to use linealized mat-vec prod.'
           end if
           exit
        endif
     
        !  ---/* update delta */----
        if(myid.eq.0) then
           write(*,*) &
                'Analysing what to do with current', &
                'Newton/hookstep and trust region'
        end if
        !     
        !    /* Compare the actual reduction in residual                */
        !    /* to the quadratic model of residual.                     */
        !    /* How well the model matches the actual reduction         */ 
        !    /* will determine, later on, how and                       */
        !    /* when the radius of the trust region should be changed.  */
        !       dimprovReq=1e-3
        !       dimprovOk=0.10
        !       dimprovGood=0.75
        !       dimprovAcc=0.10 !
        Delta_rH_req  = dimprovReq*drdx*dxpnorm 
        ! the minimum acceptable change in residual
        
        Delta_rH_ok   = dimprovOk*Delta_rQ 
        ! acceptable, but reduce trust region for next newton step
        
        Delta_rH_good = dimprovGood*Delta_rQ 
        ! acceptable, keep same trust region in next newton step
        
        Delta_rQ_acc  = dabs(dimprovAcc*Delta_rQ) 
        ! for accurate models, increase trust region and recompute hookstep
        
        Delta_rH_accP = Delta_rQ + Delta_rQ_acc 
        ! upper bound for accurate predictions
        
        Delta_rH_accM = Delta_rQ - Delta_rQ_acc 
        ! lower bound for accurate predictions
        
        !$$$         write(*,*) 'DEBUG, Deltas',Delta_rH,Delta_rH_req,Delta_rH_ok,
        !$$$     $        Delta_rH_good,Delta_rQ_acc,Delta_rH_accP,Delta_rH_accM 
        !     /* Place improvement in contiguous spectrum:*/
        !     /* improvement: (1)Unacceptable > (2)Poor > (3)Ok > (4)Good  */
        if (Delta_rH .gt. Delta_rH_req) then
           improvement = 1 ! Unacceptable;  
           !// not even a tiny fraction of linear rediction
        elseif(Delta_rH .gt. Delta_rH_ok) then
           improvement = 2 ! Poor;
           !// worse than small fraction of quadratic prediction
        elseif((Delta_rH .gt. Delta_rH_ok) &
             .and.(Delta_rH .gt. Delta_rH_good)) then
           improvement = 3 ! Ok;
        else 
           improvement = 4 ! Good;  
           !// not much worse or better than large fraction of prediction
           if( (Delta_rH_accM .le. Delta_rH) &
                .and.(Delta_rH .le.Delta_rH_accP) ) then
              improvement = 5  ! Accurate;
              !// close to quadratic prediction
           elseif(Delta_rH .lt. Delta_rL) then
              improvement = 6  ! NegaCurve;
              !// negative curvature in r(|s|) => try bigger step
           endif! /* option of Good improvement (4) */ 
        endif !/* improvement spectrum */
        if (myid.eq.0) then
           write(*,*) myid,'improvement=',improvement
           write(*,*) myid,'Deltas: ', Delta_rH_req, '>',Delta_rH_ok, '>', &
                Delta_rH_good, &
                Delta_rH_accP,Delta_rH_accM 
        end if
        ! /* CASE 1: UNACCEPTABLE IMPROVEMENT */
        ! /* Improvement is so bad (relative to gradient at current position) */
        ! /* that we'll in all likelihood get better results                  */ 
        ! /* in a smaller trust region. */
        if(improvement.eq.1) then
           if(myid.eq.0) write(*,*) 'improvement is unacceptable (case1):'
           if((ihavebackup.eq.1).and.(rH.gt.backup_rH))then
              if(myid.eq.0) write(*,*) 'case1: use backup step'
              rH=backup_rH
              delta=backup_delta
              call hookstep(kth,kth,chwk,0)
              irecompute_hookstep=0
           else
              if(myid.eq.0) write(*,*) myid,'case1: i do not have backup step'
              ! /* Reduce trust region by minimizing local quadratic model */
              ! /*        and recompute hookstep.                          */
              deltaMaxLocal = delta
              dlambdaRequiredReduction=0.5 ! reduce at least by this rate
              call adjustLambda(dlambda, dlambdaMin, &
                   dlambdaRequiredReduction,myid)
              if(myid.eq.0) write(*,*)  &
                   'case1: reduce delta by delta_rate=',dlambda
              call adjustDelta(delta, dlambda, deltaMin, deltaMax,myid)
              if (delta.gt.dxpnorm)then
                 if(myid.eq.0) then
                    write(*,*) &
                         'case1: that delta is still bigger than the Newton step.'
                    write(*,*) &
                         'case1: Reducing delta to half the length of the Newton step.'
                 end if
                 delta = 0.5*dxpnorm;
              endif
              irecompute_hookstep=1 ! or use backuped step which is acceptable
           endif
           !
           ! /* CASE 2: EXCELLENT IMPROVEMENT AND ROOM TO GROW */
        elseif(((improvement.eq.5).or.(improvement.eq.6)).and. &
             ((ihavebackup.eq.0).or.(backup_rH.gt.rH)).and. &
             (delta.lt.deltaMax).and. &
             (ihookstep_equals_newtonstep.eq.0)) then
           delta_new=delta
           if(myid.eq.0) write(*,*) 'case2: adjustDelta to: ',delta_new 
           call adjustDelta(delta_new, dlambdaMax, deltaMin, deltaMax,myid)
           if(myid.eq.0) write(*,*) 'case2: increase delta and recompute hookstep.'
           if (delta_new .lt. deltaMaxLocal)then
              write(*,*) myid,'case2: backup current step'
              ! backup this step: 
              backup_rH=rH
              backup_delta=delta
              ihavebackup=1
              irecompute_hookstep = 1
              delta=delta_new
           else
              if(myid.eq.0)write(*,*) 'Stop adjusting trust region radius delta'
              ! and take step because the new delta reached a local limit:  
              ! new_delta >= deltaMaxLocal 
              delta = deltaMaxLocal
              if(myid.eq.0)write(*,*) 'case2: delta set to local limit', delta
              irecompute_hookstep = 0
           endif
           ! /* CASE 3: MODERATE IMPROVEMENT, NO ROOM TO GROW, OR BACKUP IS BETTER */
        else
           irecompute_hookstep = 0
           if ((ihavebackup.eq.1).and.(rH.gt.backup_rH))then
              ! if backup step is better, take backup step instead of current step
              delta=backup_delta
              rH=backup_rH
              call hookstep(kth,kth,chwk,0)
              if (myid.eq.0) then
                 write(*,*) 'case3: backuped step is better'
                 write(*,*) 'case3: use backuped step'
              end if
           else
              if(ihavebackup.eq.1) then
                 if (myid.eq.0) write(*,*) 'take backuped step'
              elseif(improvement.eq.2 ) then ! improvement == Poor
                 if (myid.eq.0) write(*,*) 'case2: take current step'
                 call adjustLambda(dlambda, dlambdaMin, 1.0d0,myid)
                 call adjustDelta(delta, dlambda, deltaMin, deltaMax,myid)
                 if (myid.eq.0) write(*,*) 'case3: set delta and lambda;', &
                      delta,dlambda
              elseif((improvement.eq.3).or. &
                   (ihookstep_equals_newtonstep.eq.0).or. &
                   (delta .gt. deltaMax)) then 
                 ! improvement == OK
                 ! Take current step and leave delta unchanged, 
                 ! for the following reasons:
                 if (myid.eq.0) write(*,*) 'case3: take current step and keep delta '
              else !improvement == Good, Accurate, or NegaCurve 
                 !and no restriction on increasing apply
                 if (myid.eq.0) then 
                    write(*,*) 'case3: take step and increase delta ' 
                    write(*,*) 'case3: if there is room.'
                 end if
                 call adjustLambda(dlambda, 1.0d0, dlambdaMax,myid)
                 call adjustDelta(delta, dlambda, deltaMin, deltaMax,myid)
              endif
           endif ! /* ihavebackup */
        endif !/* end of cases*/
     end if ! iuse_quadratic_mode ==1

     
     if(irecompute_hookstep.eq.1)then
        if(myid.eq.0) then
           write(*,*) '--- recompute hookstep, using new delta = ',delta
        end if
        continue;  ! recompute hookstep
     else
        !if(myid.eq.0) write(*,*) 'got current step:lambda,delta=',dmu,delta
        ! here delta is already adjusted for next step
        if(myid.eq.0) write(*,*) 'got current step:|dxp|=',dxpnorm
        exit ! go out of hookstep loop
     endif
  enddo
  call hookstep(kth,k_cut,chwk,-1) ! deallocate only
 
  ! now finally use this dxp as a newton correction 
  x0 = x0 + dxp
  !!! xf = xtmp  (rev.1151) xf is computed again later in matvec .
  call backthrowing(vor,phi,u00,w00,chwk,x0,1)
  call backthrowing_parameter(x0,1)
  ! update boundary condition
  xofft0=xofft; xoffb0=xoffb;  ! rev.1399

  ! here only reporting (rev.1151)
  if ((isol_mode.eq.4).or.(isol_mode.eq.6)) then
     if (myid.eq.0) then 
        write(*,*)  'updated_hookstep: shiftx,z =',shiftx,shiftz
        write(*,*)  'updated_hookstep: dxp(icsx,z) =',dxp(icsx),dxp(icsz) 
     end if
  end if
  if (ishifty.eq.1) then 
     if (myid.eq.0) write(*,*)  'updated_hookstep: dxp(icsy) =', dxp(icsy) 
  endif
  if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
     !if (myid.eq.0) write(*,*) 'updated_hookstep: shear =',s 
     !if (myid.eq.0) write(*,*) 'updated_hookstep: re,Rez =',re,re*Lz*Lz*s
     if (myid.eq.0) write(*,*) 'updated_hookstep: Ly =', Ly     
     if (myid.eq.0) write(*,*) 'updated_hookstep: Tp =', timep_sh     
  endif

  if (iarclength.eq.1) then
          
     if (myid.eq.0) then 

        write(*,*) 'updated_hookstep: cpara =',trim(cpara), ' x0(icre)=', x0(icre) 

        if (trim(cpara).eq.'Rez') write(*,*) 'updated_hookstep: re, Rez =',re, re*Lz*Lz*s 
        if (trim(cpara).eq.'Ayz') write(*,*) 'updated_hookstep: Ly, Ayz =',Ly, Ly/Lz 
        if (trim(cpara).eq.'Axz') write(*,*) 'updated_hookstep: alp, Axz =',alp, gam/alp 

     end if

  endif

  ! do not forget to syncronize re and shiftx, shiftz, tfin (rev.594)
  ! see => dotp 

  if (nhookstep.ne.0) deallocate(Adxp,xH_backup)
  
end subroutine update_hookstep

subroutine hookstep_trial(vor,phi,u00,w00,chwk)
!  /***********************************************************/
!  /*   hookstep trial : copied from update_hookstep()        */
!  /***********************************************************/
  use ctes
  use bcs
  use running
  use gmres

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr
  
  real(8),dimension(buffsize) :: vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00
  real(8) fac
  !real(8) re0,s0,Ly0,timep_sh0

  ! /* (0) x=x0+dxp */
  xtmp = x0 + dxp
  ! update boundary condition
  xofft=xofft0; xoffb=xoffb0; ! this is tmp,  rev.1398 
  call backthrowing(vor,phi,u00,w00,chwk,xtmp,1)
  call backthrowing_parameter(xtmp,0)
  !
  ! /* (1) actual residual ra(x0+dxp) [='btmp'] */
  if ((isol_mode.eq.1).or.(isol_mode.eq.2)) then ! fixed point or TWs
     ! use temporarily updated xofft and xoffb 
     call set_dxdt(vor,phi,u00,w00,chwk,xtmp)
     btmp(:) = -dxdt0(:)
  elseif (isol_mode.ge.3) then 
     !write(*,*) myid,'    hookstep loop, backthrowing,ihook=',ihook
     nohist=1 ! nohist=0 => dump *.cf, and accumulate stat
     call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     nohist=1
     if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
        call phase_shift_fou(vor,phi,shiftx,shiftz)
        if (ishifty.eq.1) then
           !   extshiftx=s*0.5d0*timep_sh*shifty
           !   call phase_shift_fou(vor,phi,extshiftx,0.d0)
           fac = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
           call phase_shifty(vor,phi,u00,w00,shifty,fac)       
        endif
     endif
     call throwing(vor,phi,u00,w00,chwk,xtmp,1) 
     dfp = xtmp - xf  ! for updating preconditioner (iuse_mpc==1)
     ! keep xf to be the old data
     btmp = x0+dxp-xtmp ! because x0 is kept to be the old one 

  end if

  !btmp for extended lines for shiftx etc... should be zero...        
  if (ishifty.eq.1) btmp(icsy) = 0.d0
  if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
     btmp(icsx) = 0.d0           
     btmp(icsz) = 0.d0
  end if
  if ((isol_mode.eq.5).or.(isol_mode.eq.6))  then
     btmp(icTp) = 0.d0
  end if
  if (iarclength.eq.1) then
     btmp(icre) = 0.d0
  endif

end subroutine hookstep_trial

subroutine throwing_parameter(xvec,ireport)

  use ctes
  use running
  use gmres, only: nl,nf,shiftx,shifty,shiftz,timep,ishifty,isol_mode,iarclength, &
       &           icsx,icsy,icsz,icTp,icre,cpara, & 
       &  sca_u00,sca_w00, sca_vor, sca_phi, sca_sx, sca_sy,sca_sz, sca_time, sca_Ly, sca_re 
  use LES, only:Cles

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr,ireport

  real(8),dimension(nl) :: xvec
 
  if ( (isol_mode.ge.4).or.((ishifty.eq.1).and.(isol_mode.eq.3)) ) then
     if (isol_mode.ge.4) then 
        xvec(icsx) = shiftx/sca_sx; 
        xvec(icsz) = shiftz/sca_sz;
     endif
     if (ishifty.eq.1) then
        xvec(icsy) = shifty/sca_sy; 
     endif
     if ((isol_mode.eq.5).or.(isol_mode.eq.6)) then
        ! scaling for time period?
        xvec(icTp) = timep_sh; !/sca_time
     end if
  endif

  if (iarclength.eq.1) then

     if (trim(cpara).eq.'Rez') then
        ! scaling for re?
        xvec(icre) = re
     elseif (trim(cpara).eq.'Ayz') then
        xvec(icre) = Ly 
     elseif (trim(cpara).eq.'Axz') then
        xvec(icre) = alp 
     elseif (trim(cpara).eq.'damp_up') then 
        xvec(icre) = damp_up
     elseif (trim(cpara).eq.'damp_aa') then 
        xvec(icre) = damp_aa
     elseif (trim(cpara).eq.'LES_Cs') then 
        xvec(icre) = Cles
     else
        write(*,*) myid,'cpara error'
     end if

  end if

end subroutine throwing_parameter

subroutine backthrowing_parameter(xvec,ireport)

  use ctes
  use bcs, only:xofft,xoffb
  use running
  use gmres, only: nl,nf,shiftx,shifty,shiftz,timep,ishifty,isol_mode,iarclength, & 
       &           icsx,icsy,icsz,icTp,icre,cpara,sca_re,xofft0,xoffb0
  use LES, only:Cles

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr,ireport

  real(8) alp0,xofft_tmp
  real(8),dimension(nl) :: xvec
  
  if (ishifty.eq.1) then 
     shifty = xvec(icsy) 
     if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'shifty =',shifty
  endif
  if ((isol_mode.eq.4).or.(isol_mode.eq.6))  then
     shiftx = xvec(icsx) 
     shiftz = xvec(icsz) 
     if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) myid,'shiftx,shiftz =',shiftx,shiftz 
  end if
  if ((isol_mode.eq.5).or.(isol_mode.eq.6))  then
     !s = s0 - (s0/timep_sh)*dxp(icTp)
     ! this means changing Rez, but fixing Rez means
     ! that nothing has changed by using dxp(icTp)

     !re = re0 - (re0/timep_sh)*dxp(icTp) ! keep Rez constant

     ! (rev.600) there are only Axz,Ayz contributing S*Tb.
     ! now changing Ly using dxp(icTp)
     !Ly = Ly0 - (Ly0/timep_sh)*dxp(icTp)*sca_time
     !Ly = Ly0 - s*Ly0/(2.d0*pi/alp)*dxp(icTp) !*sca_time
     !timep_sh = 2.d0*pi/alp/Ly/s*timep ! assuming s>0    

     !timep_sh = timep_sh0 + dxp(icTp)  
     timep_sh = xvec(icTp)
     Ly = timep*(2.d0*pi/alp)/(s*timep_sh)
     call set_xyz_dy
     !timep_sh = 2.d0*pi/alp/Ly/s*timep ! assuming s>0
     etime = time_ini + timep_sh ! bug fixed (rev.600) 
     
     if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'Ly, Ayz, Tp =',Ly, Ly/Lz, timep_sh

  end if
  if (iarclength.eq.1) then
     if(trim(cpara).eq.'Rez') then
        re = xvec(icre) !*sca_re do not scale here
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 're, Rez =',re, re*Lz*Lz*s 
     elseif(trim(cpara).eq.'Ayz') then
        Ly = xvec(icre) 
        call set_xyz_dy
        timep_sh = 2.d0*pi/alp/Ly/s*timep
        etime = time_ini + timep_sh   
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'Ly, Ayz, Tp =',Ly, Ly/Lz, timep_sh 
     elseif(trim(cpara).eq.'Axz') then
        alp0 = alp
        xofft_tmp = xofft ! rev1390, xofft0 ==> xofft_tmp
        alp = xvec(icre)
        ! update boundary condition
        xofft =  xofft_tmp*alp0/alp
        xoffb = -xofft
        !xofft0=xofft ! (rev.1390) ! this should not be modified here... (rev.1399)
        !xoffb0=xoffb
        call set_xyz_dy
        timep_sh = 2.d0*pi/alp/Ly/s*timep
        etime = time_ini + timep_sh   
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'alp, Axz, Tp =',alp, gam/alp, timep_sh 

     elseif(trim(cpara).eq.'damp_up') then
        damp_up = xvec(icre)
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'damp_up =',damp_up
     elseif(trim(cpara).eq.'damp_aa') then
        damp_aa = xvec(icre)
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'damp_aa =',damp_aa 
     elseif (trim(cpara).eq.'LES_Cs') then 
        Cles = xvec(icre)  
        if ((myid.eq.0).and.(ireport.eq.1)) write(*,*) 'Cles =',Cles 
     endif
  endif

end subroutine backthrowing_parameter



subroutine get_eigen(kth,ireport,cv)

  !/***************************************************************/
  !/* compute eigenvalues and vector of Hessenberg matrix .       */
  !/* input:                                                      */
  !/*   kth  : the k-th gmres step                                */
  !/*   ireport: 1, screen output                                 */
  !/*   ireport: 2, for memory check                              */
  !/*   ireport: 999, to dump eigen functions (testing)           */
  !/*   cv: this is passed to the option for dgeev                */
  !/*       cv='n', only eigenvalues                              */
  !/*       cv='v', compute eigen functions using                 */ 
  !/*               the right eigenvector will be computed        */
  !/* output                                                      */
  !/*                                                             */
  !/***************************************************************/
  use ctes
  use running
  use gmres
  use eig
  use precond
  use LES,only:iuse_LES,idynamic,iadd_visc,Cles

  implicit none
  include "mpif.h"

  integer kth,i,k,ireport,info,iflqt
  character*1 cv
  character*128 fname
  character*4 ext1
  real(8) dummy,mem_eig,time_eig
  real(8),dimension(:,:),allocatable:: hvtmp

  iflqt=99

  if(myid.gt.0) then 
     write(*,*) 'eig is done only by master'
     return
  endif
  time_eig=MPI_WTIME()
  ! restore Hessenberg matrix
  allocate(hvtmp(kth,kth))
  hvtmp=hh(1:kth,1:kth)
  mem_eig=kth*kth   
  !
  ! eig module
  lwkeig=4*kth ! not optimal value
  allocate(Dnr(kth),Dni(kth),wkeig(lwkeig)) ! in the module 'eig'
  Dnr=0.d0; Dni=0.d0; wkeig=0.d0
  mem_eig=mem_eig + kth*6.d0    
  
  allocate(VR(kth,kth))
  VR=0.d0
  mem_eig=mem_eig + kth*kth    
  ! check optimal work array size
  call dgeev('n',cv,kth,hvtmp,kth,Dnr,Dni,dummy,1, &
       VR,kth,wkeig,-1,info)  
  lwkeig=wkeig(1)
  write(*,'(i4,a,i8)')  myid,': eig: optimal lwk=',lwkeig
  mem_eig=mem_eig + lwkeig + nl*kth    
  
  write(*,'(i4,a,G9.2)') myid,': eig requests (MB):',(mem_eig)*8.d0/1024.d0/1024.d0
  deallocate(wkeig) 
  if (ireport.eq.2) then 
     deallocate(VR,Dnr,Dni,hvtmp)
     return
  endif
  !
  allocate(wkeig(lwkeig))
  wkeig=0.d0
  ! Solve the eigenvalue problem (need: lapack) ![Vn,Dn]=eig(hh);
  call dgeev('n',cv,kth,hvtmp,kth,Dnr,Dni,dummy,1, &
       VR,kth,wkeig,lwkeig,info)

  ! (rev.1137) commented because the code does not go to finalize ... 
  !if (cv.eq.'v') then
  !   allocate(eigf(nl,kth))
  !   eigf=0.d0
  !   do k=1,kth     
  !      do i=1,nl
  !         eigf(i,k)=sum(vv(i,:)*VR(:,k)) ! eigen function (a flow field) 
  !         ! note that vv(nmax1,kth+1)
  !      end do
  !   end do
  !endif

  if ((ireport.eq.1).and.(myid.eq.0)) then
     write(*,*) 'eig: Floquet multiplier of preconditioned Jacobian, [M-I]*P dyp'
     do k=1,kth
        write(*,*) 'eig: (',Dnr(k),Dni(k),')'
     end do
  elseif ((ireport.eq.999).and.((iuse_scaling.eq.0).and.(iuse_mpc.eq.0))) then 
     ! at the end of newton-step in iloop_newton
     do k=1,kth
        if (isol_mode.ge.3) then
           Dnr(k) = Dnr(k)+1.d0
        end if
        !if (myid.eq.0) write(*,*) 'B: Dnr Dni',Dnr(k), Dni(k)
     end do
     allocate(Dnc(kth))     
     do k=1,kth
        Dnc(k)=dcmplx(Dnr(k),Dni(k))
        if (myid.eq.0) write(*,*) 'A: sigma ',log(Dnc(k))/(etime-time_ini)
     end do
     deallocate(Dnc)
     if (myid.eq.0) then
        write(ext1,'(i4.4)') id22-1 ! fixed ...
        fname=filout(1:index(filout,' ')-1)//'.'//ext1//'.flqt'
        open(iflqt,file=fname,status='unknown')
        write(iflqt,'(a,4(1X,F14.6))') '% ', gam/alp, Ly/Lz, re*Lz*Lz, (etime-time_ini)
        write(iflqt,'(a,3(1X,I10))') '% grids (phys) ', mgalx,my,mgalz 
        write(iflqt,*) '% re,alp ',re,alp
        write(iflqt,*) '% Ly,gam ',Ly,gam
        write(iflqt,*) '% s,chi ',s,chi
        !write(iflqt,'(a,3(1X,I10),3(1X,F14.6))') '% ', nstep,nimag,nhist,pmesp,uprim,CFL
        if (iadd_force.eq.1) then
           write(iflqt,*) '% iadd_force=1',force_roll
        elseif (iadd_force.eq.-1) then
           write(iflqt,*) '% iadd_force=-1 (x,y,z)-translations:'
           write(iflqt,*) '% xf,vb,zf=',xforce,vbulk,zforce
        end if
        if (iadd_damping.ne.0) then
           write(iflqt,*) '% iadd_damping',damp_up,damp_aa 
        end if
        if (iuse_LES.eq.1) then
           write(iflqt,*) '% iuse_LES',idynamic,iadd_visc,Cles
        end if
        !
        do k=1,kth
           write(iflqt,'(2(1X,F14.6))') Dnr(k), Dni(k)
        enddo
        close(iflqt)
     endif
  endif
  !if (cv.eq.'v') deallocate(eigf)
  deallocate(wkeig,VR,hvtmp)           
  deallocate(Dnr,Dni)

  time_eig=MPI_WTIME() - time_eig
  write(*,*) myid,'time_eig',time_eig

end subroutine get_eigen

!-------------------------------------------------------------------
subroutine adjustLambda(dlam,dlamMin,dlamMax,myid)

  implicit none

  real*8 dlam,dlamMin,dlamMax
  integer myid

  if(dlam.gt.dlamMax) then
     if(myid.eq.0) write(*,*) 'dlambda is maximum'
     dlam= dlamMax
  elseif(dlam.lt.dlamMin) then
     if(myid.eq.0) write(*,*) 'dlambda is minimum'
     dlam= dlamMin
  else
     dlam=dlam
  endif
  if(myid.eq.0) write(*,*) 'delta_rate is set to: ',dlam

end subroutine adjustLambda
!--------------------------------------------------------------------
subroutine adjustDelta(dlt,dlt_rate,dltMin,dltMax,myid)

  implicit none

  real*8 dlt,dlt_rate,dltMin,dltMax
  integer myid

  dlt=dlt*dlt_rate
  if(dlt.gt.dltMax) then
     if(myid.eq.0) write(*,*) 'delta is maximum'
     dlt=dltMax
  elseif(dlt.lt.dltMin) then
     if(myid.eq.0)    write(*,*) 'delta is minimum'
     dlt=dltMin
  endif
  if(myid.eq.0) write(*,*) 'newton trust region is set to: ',dlt
  return
end subroutine adjustDelta

!-------------------------------------------------------------------!
! vector tools
!-------------------------------------------------------------------!
subroutine vecadd(nx,x,fac,y,z)

  ! z <- x+fac*y 

  implicit none
  !include "mpif.h"

  integer i,nx
  real(8),dimension(nx) :: x, y, z
  real(8) fac

  z(:) = x(:) + fac*y(:) 

end subroutine vecadd

subroutine vecscal(nx,x,fac,y)

  ! y <- fac*x

  implicit none
  !include "mpif.h"

  integer i,nx
  real(8),dimension(nx) :: x, y
  real(8) fac

  y(:) = fac*x(:) 

end subroutine vecscal

subroutine vcopy_00(nx,x,ny,y,yoffset)

  ! y <- x

  implicit none
  !include "mpif.h"

  integer i,nx,ny,yoffset
  real(8),dimension(nx) :: x, y

  if (nx.ne.(ny-yoffset)) then
     write(*,*) 'vcopy_00 error, check nx=ny-yoffset',nx,ny,yoffset
  endif
  y(yoffset+1:ny) = x(1:nx) 

end subroutine vcopy_00

subroutine vcopy_xz(x,y,index,j)

  ! y <- x
  use ctes, only : mx,my,mz,icx,nz,alp,gam
  use gmres, only : iuse_us
  implicit none
  !include "mpif.h"

  integer i,k,index,j
  real(8) xmode,zmode
  real(8),dimension(index:*) :: x
  real(8),dimension(2,mx,mz) :: y
  do k=1,mz
     do i=1,mx 
        if (((i.eq.1).and.(k.eq.1)).or.((i.eq.2).and.(k.eq.1))) then
           ! skip 00-modes of vor and phi
        elseif ( (i.le.2).and.( k.ge.nz+2 ) ) then 
           !filter(index:index+1)=0.d0 ! 2*2*(mz-nz+1)*my 
           ! remove conj. imaginary part of vor(0,k=nz+1:mz1,y)
        else

           if (iuse_us.eq.1) then
              if ( (i-1)/2.eq.0) then
                 xmode=1.d0*alp
              else
                 if (mod(i-1,2).eq.0) then
                    xmode=(i-1)/2*alp
                 elseif (mod(i-1,2).eq.1) then
                    xmode=(i)/2*alp
                 end if
              endif
              if (icx(k-1).eq.0) then
                 zmode=1.d0*gam
              else
                 zmode=dfloat(icx(k-1))*gam
              end if
              y(1,i,k)=x(index)/zmode*(xmode*xmode+zmode*zmode)        
           else
              y(1,i,k)=x(index)
           end if
           index=index+1
           y(2,i,k)=x(index)
           index=index+1
        endif
     enddo
  enddo


end subroutine vcopy_xz

subroutine vcopy(nx,x,y)

  ! y <- x

  implicit none

  integer i,nx
  real*8 x(nx),y(nx)

  y(1:nx) = x(1:nx) 

end subroutine vcopy


