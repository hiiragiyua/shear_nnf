program shooting

  use ctes
  use running
  use bcs
  ! do not use statistics
!  use gmres
!  use eig

  implicit none
  include "mpif.h"

  integer iproc,itags,newtag,imess,i,j,k,kk,idum,info
  real*8  val,dp,fac,mtime,eps_j,err,dummy
  real*8, allocatable:: vor(:),phi(:),chwk(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: vorwk(:),phiwk(:) ! for escru
  real*8, allocatable:: u00wk(:),w00wk(:)

  real*8, allocatable:: vor1(:),phi1(:)
  real*8, allocatable:: u001(:),w001(:)

  real*8, allocatable:: vor2(:),phi2(:)
  real*8, allocatable:: u002(:),w002(:)

  ! for shooting goal
  real*8, allocatable:: vorg(:),phig(:)
  real*8, allocatable:: u00g(:),w00g(:)

  ! reading parameter
  real*8 para_shoot(1:10),shtmu,shtmu_base,timep,rem1, time1, time2, timeg
  !real*8 kener(1:3),kener2(1:3),kener_fin(1:3)
  real*8 ke1, ke2, ke_iso, ke_fin, ke12r, res, res0
  character*128 fname,filinmg,filout_base
  character*3 extl
 

  integer ishoot_mode,nloop,iloop,igo

  integer istat(MPI_STATUS_SIZE),ierr


  ! /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)


  ! version (revision number of Subversion)
  if(myid.eq.0) write(*,*) 'SHOOTING method for HSF by A.S., Version: rev.',irev

!       '(beta, compute dxdt for all isol_modes, hookstep fixed, dump stat in make_mat), arclength, dre1.e-6, readable-arcl, stop criteria of newton: 1.e-6, timp fixed, hre2.dat, eigf, deltaMax=1, ihv_0, iadd_sym=-1, relative in xz (isol4), isol5_6,dmu=0.1,sca_vor, nosca, auto-arc-1, chwk ini, auto-arc1, ihave_small_res (fixed), close(39), reset delta, improving arcl, arc-1 new, ihist_reset,fixed ihave_small (arc-1), sta4, Pk, iadd5(exp,imp), shooting, iso noise(shoot2), fixed_fac (iso), fixed Ree, iadd_sym=5,6, fixed fac, error check getfil, addsym fixed,addsym8, ishoot10, addsym chwk, fixed random seed , set_conj, fixed argument for getini_tophat'
  ! (should correspond to svn revision) 
  !--------------- initializes commons and things
  call initcr(1)
  call set_options()
  !
  iskip_screenout = 1
  nloop=1 ! default 1-loop
  para_shoot=0.d0
  !  
  if(myid.eq.0) then
     write(*,*) 'set shooting parameter ', & 
          ' ... [ ]'
     read(*,'(a)') filinm1
     read(*,*) para_shoot(1) ! shooting parameter 0 < lam < 1 (for laminar-turbulent boundary) 
     read(*,*) para_shoot(2) ! time shear-period 
     para_shoot(3)=iadd_sym
     read(*,*) para_shoot(4) ! num. trial  
     read(*,*) para_shoot(5) ! shooting mode
     if (int(para_shoot(5)).eq.1) then
        write(*,*) 'set base  perturbation (shoot parameter will be (shtmu + base)'
        read(*,*) para_shoot(6) ! base perturbation
     endif
     if (int(para_shoot(5)).eq.10) then
        read(*,'(a)') filinmg
     end if
  end if

  call  MPI_BCAST(para_shoot,10,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

  shtmu = para_shoot(1) ! shooting parameter
  shtmu_base = para_shoot(6) ! shooting parameter
  !
  timep = para_shoot(2) ! time period for periodic orbit 
  !
  nloop = int(para_shoot(4))
  ishoot_mode = int(para_shoot(5))
  !
  !timep_sh = 2.d0*pi/alp/Ly*timep ! fixed (rev.600, assuming s>0)
  if (myid.eq.0) write(*,*) & 
       & 'set the integration time as a multiple of the shear-period: ', & 
       & 'Tp=', timep_sh
  !
  noescru=0 ! 1: for no output and skip allocation of sp, spl
  nohist =0
  !

  !
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;

  !write(*,*) 'allocating vor1'
  allocate(vor1(buffsize), phi1(buffsize), u001(my), w001(my) )
  vor1 = 0.d0; phi1 = 0.d0; u001 = 0.d0; w001 = 0.d0;

  if (readflag.eq.0) then
     if (myid.eq.0) write(*,*) 'error: shooting needs initial condition'
     stop
  else
     call getfil(vor1,phi1,u001,w001,chwk,time1)  ! read i.c.
     if (iadd_sym.ge.5) then
        write(*,*) 'iadd_sym, phase shift for vor1 phi1',iadd_sym, sym_shiftx*Lx, sym_shiftz*Lz 
        call phase_shift_fou(vor1,phi1,sym_shiftx*Lx,sym_shiftz*Lz)
     end if
     rem1=Ree
     call get_kinetic_energy(vor1,phi1,u001,w001,chwk,ke1)
     if (myid.eq.0) write(*,*) ' ke1 (upo)', ke1 
  endif

  ! error check to read the UPOs.
  if (abs(Lye-Ly).gt.tiny) then 
     write(*,*)' getfil: Ly is different, STOP!!',Lye
     stop
  end if
  if (abs(Ree-re).gt.tiny) then
     write(*,*)' getfil: re is different, STOP!!',ree
     stop
  end if
  if (abs(alpe-alp).gt.tiny) then
     write(*,*)' getfil: alp is different, STOP!!',alpe
     stop
  end if
  if (abs(game-gam).gt.tiny) then
     write(*,*)' getfil: gam is different, STOP!!',game
     stop
  end if

  time=time1
  time_ini=time
  etime = time + timep_sh*timep
  xoffb_ini=xoffb
  xofft_ini=xofft

  if (idump_mode.eq.2) then
     ! for making UPO movies, set time_interval to dump a file and etime 
     !
     if(myid.eq.0) write(*,*) '  ------ UPO dump mode  ----- '
     if(myid.eq.0) write(*,*) '  ------    time_ini,  =>  ', time_ini
     dtimag = timep_sh/dump2tint  ! check the time intervals and 
                             ! dt should be considered for statistics... 
     time = (xofft_ini/abs(s)/Ly)
     if(myid.eq.0) write(*,*) '  ------ the time stump is adjusted with xofft ==>  ', time
     
     time_ini = time
     !dumptime = int(time/timep_sh)*timep_sh + timep_sh ! fixed (rev.545)
     !
     nimag = nstep

  end if


  ! set run options

  if(myid.eq.0) write(*,*) '... file1 read'
  !
  if (myid.eq.0) write(*,*) ' reading a unstable direction:', trim(filinm1)
  !write(*,*) 'allocating vor2'
  allocate(vor2(buffsize), phi2(buffsize), u002(my), w002(my) )
  vor2 = 0.d0; phi2 = 0.d0; u002 = 0.d0; w002 = 0.d0;
  if ((ishoot_mode.eq.1).or.(ishoot_mode.eq.2)) then
     
     filinp = filinm1     
     call getfil(vor2,phi2,u002,w002,chwk,time2)  ! read i.c., Ree is overwritten
     if (iadd_sym.ge.5) then
        
        write(*,*) 'iadd_sym, phase shift for vor2 phi2',iadd_sym,sym_shiftx*Lx,sym_shiftz*Lz 
        call phase_shift_fou(vor2,phi2,sym_shiftx*Lx,sym_shiftz*Lz)
        
     end if

     !call energy(vor2,phi2,u002,w002,ener2)
     call get_kinetic_energy(vor2,phi2,u002,w002,chwk,ke2)
     if (myid.eq.0) write(*,*) ' ke2 (eigf)', ke2 
     
     ke12r=ke1/ke2
  end if
  ! error check to read the UPOs.
  if (abs(Lye-Ly).gt.tiny) then 
     write(*,*)' getfil2: Ly is different, STOP!!',Lye
     stop
  end if
  if (abs(Ree-re).gt.tiny) then
     write(*,*)' getfil2: re is different, STOP!!',ree
     stop
  end if
  if (abs(alpe-alp).gt.tiny) then
     write(*,*)' getfil2: alp is different, STOP!!',alpe
     stop
  end if
  if (abs(game-gam).gt.tiny) then
     write(*,*)' getfil2: gam is different, STOP!!',game
     stop
  end if


  if (abs(re-rem1).gt.tiny) then
     igo=0
     if (myid.eq.0) write(*,*) 're is different ...re(hre.dat)=',re,' <= re', rem1
     if (myid.eq.0) write(*,*) ' ... 1: OK, 0: Stop'
     if (myid.eq.0) read(*,*) igo
     call  MPI_BCAST(igo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     if (igo.eq.0) then
        stop
     end if
  end if
  if(myid.eq.0) write(*,*) '... file2 read'

  
  if (ishoot_mode.eq.10) then
     if (myid.eq.0) write(*,*) ' set the goal UPO:', trim(filinmg)
     !write(*,*) 'allocating vorg'
     allocate(vorg(buffsize), phig(buffsize), u00g(my), w00g(my) )
     vorg = 0.d0; phig = 0.d0; u00g = 0.d0; w00g = 0.d0;
     filinp = filinmg  
     call getfil(vorg,phig,u00g,w00g,chwk,timeg)  ! read i.c., Ree is overwritten     

     if (iadd_sym.ge.5) then
        write(*,*) 'iadd_sym, phase shift for vorg phig',iadd_sym, sym_shiftx*Lx, sym_shiftz*Lz 
        call phase_shift_fou(vorg,phig,sym_shiftx*Lx,sym_shiftz*Lz)
     end if
  end if
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  filout_base=trim(filout)
  !
  do iloop=1,nloop
     ! reinitialization...
     vor=0.d0; phi=0.d0; u00=0.d0; w00=0.d0; chwk=0.d0;
     id22=0; ! set the index of the first file
     !
     if (myid.eq.0) write(*,*)  ' ==== shooting, iloop =',iloop,'/',nloop, ' ==== '
     !
     write(extl,'(i3.3)') iloop 
     filout=trim(filout_base)//'_'//extl  
  
     if (idump_mode.eq.2) then
        ! for making UPO movies, set time_interval to dump a file and etime 
        !
        if(myid.eq.0) write(*,*) '  ------ UPO dump mode  ----- '
        if(myid.eq.0) write(*,*) '  ------    time_ini,  =>  ', time_ini
        ! reset time to set the first dumptime
        time = time_ini
        dumptime = int(time/timep_sh)*timep_sh + timep_sh ! fixed (rev.545)

        if (abs(int(time/timep_sh)*timep_sh - time_ini).lt.1.e-6) dumptime=time+dtimag
        
        if(myid.eq.0) write(*,*) '  ------ the first dump time (pure-periodic),  =>  ',dumptime        
        
     end if     

     !
     if (ishoot_mode.eq.1) then
        ! add fac*eigf + UPO-LB, or -UB (real part)
        fac = (shtmu_base + shtmu*dfloat(iloop-1))
        if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
        fac = fac*ke12r ! rev.818 (reverted, but the output is different)
        vor=vor1+fac*vor2; phi=phi1+fac*phi2
        u00=u001+fac*u002; w00=w001+fac*w002

     elseif (ishoot_mode.eq.2) then
        ! add randum to above, fac*eigf + UPO-LB, to get probability of the separatorix
        !fac = shtmu
        !if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu =', fac
        if (iloop.eq.1) then
           ! add the primary unstable direction
           fac=shtmu
           if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode 
           if (myid.eq.0) write(*,*)  '     mu, fac =', shtmu, fac
           fac=fac*ke12r
           vor1=vor1+fac*vor2; phi1=phi1+fac*phi2
           u001=u001+fac*u002; w001=w001+fac*w002
           
        end if
        ! add white noise with the order of 0.1*fac
        !uprim=0.01
        
        if (myid.eq.0) write(*,*)  '     add isotropic noise of uprim', uprim
        ! top hat option: 1, R&M [8 16]; 2, [4 32]
        ! reinitialize vor2 ...
        !vor2=0.d0; phi2=0.d0; u002=0.d0; w002=0.d0;
        call getini_tophat(vor2,phi2,u002,w002,chwk,time2,2) ! for pure-periodic b.c.
        ! ! Note: check this tophat is really randum for each iloop !!!
        !call energy(vor2,phi2,u002,w002,ener2)
        call get_kinetic_energy(vor2,phi2,u002,w002,chwk,ke_iso)
        if (myid.eq.0) write(*,*) ' ke_iso ', ke_iso 
        ! tophat of R&M does not work as a randum noise.
        !call getini_TG(vor2,phi2,u002,w002,chwk,time2) ! for pure-periodic b.c.
        fac = (xoffb_ini/Lx - int(xoffb_ini/Lx) )*Lx/Ly
        if (myid.eq.0) write(*,*)  '      the noise is distorted to satisfy B.C.', xoffb_ini,fac
        !if (myid.eq.0) write(*,*)  '          ener2 =', ener2
        call moving_frame(vor2,phi2,fac,-1)    

        fac=0.01
        if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
        fac=fac*ke2/ke_iso
        vor=vor1+fac*vor2; phi=phi1+fac*phi2;
        u00=u001+fac*u002; w00=w001+fac*w002; 
        
     elseif (ishoot_mode.eq.3) then
        
        if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu =', shtmu

        ! add only isotropic noise
        uprim=0.1
        
        if (myid.eq.0) write(*,*)  '     add isotropic noise of uprim', uprim
        ! top hat option: 1, R&M [8 16]; 2, [4 32]
        ! reinitialize vor2 ...
        vor2=0.d0; phi2=0.d0; u002=0.d0; w002=0.d0;
        call getini_tophat(vor2,phi2,u002,w002,chwk,time2,2) ! for pure-periodic b.c.
        !call energy(vor2,phi2,u002,w002,ener2)
        call get_kinetic_energy(vor2,phi2,u002,w002,chwk,ke_iso)
        if (myid.eq.0) write(*,*) ' ke_iso ', ke_iso 
        ! tophat of R&M does not work as a randum noise.
        !call getini_TG(vor2,phi2,u002,w002,chwk,time2) ! for pure-periodic b.c.
        fac = (xoffb_ini/Lx - int(xoffb_ini/Lx) )*Lx/Ly
        if (myid.eq.0) write(*,*)  '      the noise is distorted to satisfy B.C.', xoffb_ini,fac
        call moving_frame(vor2,phi2,fac,-1)              

        fac = shtmu
        if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
        fac = fac*ke1/ke_iso ! BUG fixe rev861
        vor = vor1 + fac*vor2; phi = phi1 + fac*phi2;
        u00 = u001 + fac*u002; w00 = w001 + fac*w002; 

     elseif (ishoot_mode.eq.4) then               
        
        if (iloop.eq.1) then
           if (myid.eq.0) write(*,*)  ' only at first, create an isotropic noise', uprim
           ! add only isotropic noise
           uprim=0.1
           ! top hat option: 1, R&M [8 16]; 2, [4 32]
           ! reinitialize vor2 ...
           vor2=0.d0; phi2=0.d0; u002=0.d0; w002=0.d0;
           call getini_tophat(vor2,phi2,u002,w002,chwk,time2,2) ! only for pure-periodic b.c.
           !call energy(vor2,phi2,u002,w002,ener2)
           call get_kinetic_energy(vor2,phi2,u002,w002,chwk,ke_iso)
           if (myid.eq.0) write(*,*) ' ke_iso ', ke_iso 
           ! tophat of R&M does not work as a randum noise.
           !call getini_TG(vor2,phi2,u002,w002,chwk,time2) ! for pure-periodic b.c.
           fac = (xoffb_ini/Lx - int(xoffb_ini/Lx) )*Lx/Ly
           if (myid.eq.0) write(*,*)  '      the noise is distorted to satisfy B.C.', xoffb_ini,fac
           call moving_frame(vor2,phi2,fac,-1)              
        end if

        fac = shtmu*float(iloop)
        if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu*float(iloop), fac
        fac = fac*ke1/ke_iso
        vor = vor1 + fac*vor2; phi = phi1 + fac*phi2;
        u00 = u001 + fac*u002; w00 = w001 + fac*w002; 
     
     elseif (ishoot_mode.eq.10) then
        ! newton method based on ishoot==1 
        if (iloop.eq.1) then
           fac = shtmu ! eps_newton
           if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
           fac = fac*ke12r 
           vor=vor1+fac*vor2; phi=phi1+fac*phi2
           u00=u001+fac*u002; w00=w001+fac*w002
        elseif (iloop.eq.2) then
           fac = shtmu + 1.d-6 ! eps_newton
           if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
           fac = fac*ke12r 
           vor=vor1+fac*vor2; phi=phi1+fac*phi2
           u00=u001+fac*u002; w00=w001+fac*w002
        elseif (iloop.ge.3) then
           fac = shtmu  - res0/((res-res0)/1.d-6) ! eps_newton
           if (myid.eq.0) write(*,*)  ' ====   mode =', ishoot_mode, 'mu, fac=', shtmu, fac
           fac = fac*ke12r 
           vor=vor1+fac*vor2; phi=phi1+fac*phi2
           u00=u001+fac*u002; w00=w001+fac*w002
        end if
     else
        write(*,*) 'not implemented ishoot_mode =', ishoot_mode
        stop
     end if
     ! --- the shooting process ---
     nohist=0 ! screen output and dump *.cf and stat. file
     !nohist=1
     call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     allocate( phiwk(buffsize),vorwk(buffsize),u00wk(my),w00wk(my) )
     phiwk = 0.d0; vorwk = 0.d0; 
     vorwk = vor-vor1; phiwk = phi-phi1;
     u00wk = u00-u001; w00wk = w00-w001;
     
     call get_kinetic_energy(vor,phi,u00,w00,chwk,ke_fin)
     !if (myid.eq.0) write(*,*)  'growth_rate ', log(sqrt(sum(ke_fin))/ke1)/timep_sh
      
     if (ishoot_mode.eq.10) then
        phiwk = 0.d0; vorwk = 0.d0; 
        vorwk = vor-vorg; phiwk = phi-phig;
        u00wk = u00-u00g; w00wk = w00-w00g;
        ! compute L_2 norm of the difference of the goal states vorg,phig
        res0=res; ! store prefious residual 
        call get_kinetic_energy(vorwk,phiwk,u00wk,w00wk,chwk,res)
        if(myid.eq.0) write(*,*) iloop,':ishoot=10, residual;', res/ke1
        ! check convergence
        if (res/ke1.lt.1.e-6) then
           if(myid.eq.0) write(*,*) iloop,':ishoot=10, converged;', res/ke1
           exit
        end if
     end if
     deallocate(vorwk,phiwk,u00wk,w00wk)

  enddo ! end of NSTEP loop 

  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)   !! used to be statement 200

end program shooting

subroutine calc_dns_disk3d(vor,phi,u00,w00,chwk)

  ! the same with in gmres_shear.f90, but without setting some gmres-newton options...
  
  use ctes
  use running
  use statistics
  use bcs
  use LES,only:iuse_LES,Cles,Deltag

  implicit none
  include "mpif.h"

  real(8),dimension(buffsize) ::  vor, phi, chwk
  real(8),dimension(0:my1) :: u00, w00

  real*8  end_time,y0
  integer iproc,ierr
  real*8, allocatable:: hv(:),hg(:),vorwk(:),phiwk(:),dvordy(:)

  ! -- set buffers for LES --- (txy00, tzy00 are in cross.f90)
  real(8),dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2
  real*8 chkre(0:numerop-1),chkalp(0:numerop-1), & 
       chkgam(0:numerop-1),chkLy(0:numerop-1),chks(0:numerop-1), &
       chky0(0:numerop-1)

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
  ! reset boundary condition
  time = time_ini
  xoffb = xoffb_ini 
  xofft = xofft_ini
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
  !
  !if ( (iuse_LES.eq.1).and.(abs(Deltag-(dx*dy*dz)**(1.d0/3.d0) ).gt.tiny) ) then
  !   write(*,*) 'calc_dns (LES): Deltag is not updated in newton arclength method!! stop'
  !   stop
  !endif
  !write(*,*) 'calc dns: time=',time, time_ini, xoffb_ini, xofft_ini
  if (myid.eq.0) write(*,*) '--- start DNS ---'
  call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy, &
       &     rhvc1,rhvc2,rhgc1,rhgc2,chwk) 
  if(myid.eq.0) write(*,*) '--- calc_dns_disk3d finished ...'

  ! --- deallocate ---  
  deallocate(rhvc1,rhvc2,rhgc1,rhgc2)
  deallocate(phiwk,vorwk,dvordy)
  deallocate(hv,hg)

end subroutine calc_dns_disk3d

subroutine get_kinetic_energy(vor,phi,u00,w00,chwk,kine)

  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"

  real*8  u00(0:my1),w00(0:my1)
  real*8  vor(buffsize),phi(buffsize),chwk(buffsize)
  real*8  end_time
  real*8 up,vp,wp,kine
  real*8 uner(nner),uner2(nner2)
  integer iproc,ierr
  real*8, allocatable:: hv(:),hg(:),dvordy(:)

  real*8 chkre(0:numerop-1),chkalp(0:numerop-1), & 
       chkgam(0:numerop-1),chkLy(0:numerop-1),chks(0:numerop-1)

  !  ---  allocate rest of buffers --- 
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0; 
  allocate( dvordy(buffsize) )
  dvordy = 0.d0;

  call calc_energy(vor,phi,u00,w00,hv,hg,dvordy,chwk)
  ! --- deallocate ---  
  deallocate(dvordy)
  deallocate(hv,hg)

  call MPI_ALLREDUCE(ener,uner,nner,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
  uner(1:3)=dsqrt(dabs(uner(1:3)/my)) ! characteristic velocities...
  uner(4:6)=uner(4:6)/my ! bug fixed 2012/Aug/10 (rev.510)
  uner(7:9)=dsqrt(dabs(uner(7:9)/my)) ! characteristic vorticities...
  uner(10)=uner(10)/my

  call MPI_ALLREDUCE(ener2,uner2,nner2,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
  uner2(1)=uner2(1)/my ! Pk

  if(myid.eq.0) then
     write(* ,'(a,14(1X,E14.6))') & 
          'CALC ENERGY:', &
          uner(1:10),uner2(1:nner2) 
  endif

  up=uner(1);
  vp=uner(2);
  wp=uner(3);

  kine = 0.5*(up + vp + wp)

end subroutine get_kinetic_energy

subroutine calc_energy(vor,phi,u00,w00,hv,hg,dvordy,chwk) 
  use ctes
  use running
  use statistics 
  use bcs
  
  implicit none
  include "mpif.h"

  integer i,j,k
  real*8  rk,rk2

  complex*16 shp,shm
  real*8  u00(0:my1),w00(0:my1)
  real*8  phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  real*8  hv(0:2*my-1,0:mx1,kb:ke) 
  real*8  hg(0:2*my-1,0:mx1,kb:ke)
  real*8  dvordy(0:2*my-1,0:mx1,kb:ke)
  real*8  chwk(buffsize)

  integer iproc,istat(MPI_STATUS_SIZE),ierr

  real*8, allocatable:: buff1(:), buff2(:)

  real*8, allocatable:: rf0u(:),rf0w(:)
  real*8, allocatable:: u1r(:,:),u2r(:,:),u3r(:,:),o1r(:,:),o2r(:,:),o3r(:,:)
  !real*8, allocatable:: divr(:,:)

  ! tmp 2d arrays
  real*8, allocatable:: tmpxzr(:,:) ! for fouxz, phys2 
  real*8, allocatable:: tmpx(:)
  real*8, allocatable:: wk1dr(:),wk1dc(:)
  character fname*128

  ! --- add new buffer to save change (take care for the confusing names)---
  allocate (buff1(buffsize),buff2(buffsize))
  buff1 = 0.d0; buff2 = 0.d0;

  ! ------------------- allocate everything ------------  
  allocate (u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz) )
  allocate (o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz) )
  !allocate (divr(mgalx+2,mgalz))
  allocate (rf0u(0:my1), rf0w(0:my1))
  u1r = 0.d0; u2r = 0.d0; u3r = 0.d0;
  o1r = 0.d0; o2r = 0.d0; o3r = 0.d0;
  ! ---- allocate tmp 2D arrays ----
  allocate (tmpxzr(mgalx+2,mgalz))
  !allocate (oyc(0:mx1,0:mz1),vc(0:mx1,0:mz1),lapvc(0:mx1,0:mz1))
  tmpxzr = 0.d0
  allocate (tmpx(mgalx+2))
  tmpx = 0.d0
  allocate( wk1dr(my), wk1dc(2*my))
  wk1dr = 0.d0; wk1dc = 0.d0;

  xwkt = xofft; xwkb = xoffb 
  zwkt = zofft; zwkb = zoffb 
  
  if (s.gt.0.5d0) then ! s shoud be 0 or 1
     
     !do k=0,mz1
     do i=0,mx1
        shb(i)  = cdexp(-xalp(i)*xwkb) ! negative shift
        sht(i)  = cdexp(-xalp(i)*xwkt) ! positive shift
     enddo
     !enddo

  end if

  !   ---   calcula la v a partir de phi */
  do k=kb,ke
     do i=0,mx1
        if ((i.ne.0).or.(k.ne.0)) then
           rk = gam2(k)+alp2(i)
           shp = sht(i)
           shm = shb(i)
           call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)   ! -- v
        end if
     enddo
  enddo
  
  u00b = s*(-Ly) ! negative
  u00t = s*Ly ! positive
  w00b = 0.d0
  w00t = 0.d0

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
  call chjik2ikj(buff1,buff1,chwk,chwk)
  call chjik2ikj(buff2,buff2,chwk,chwk) ! ome2c
  !call chjik2ikj(phi,phi,chwk,chwk)
  !call chjik2ikj(vor,vor,chwk,chwk)
  call chjik2ikj(hv,hv,chwk,chwk)
  call chjik2ikj(hg,hg,chwk,chwk)
  call chjik2ikj(dvordy,dvordy,chwk,chwk)

  ! allocate tmp 3D arrays
  !allocate(ur(mgalx+2,mgalz,jb:je),vr(mgalx+2,mgalz,jb:je),wr(mgalx+2,mgalz,jb:je))
  !allocate(oxr(mgalx+2,mgalz,jb:je),oyr(mgalx+2,mgalz,jb:je),ozr(mgalx+2,mgalz,jb:je))
  !ur = 0.d0; vr = 0.d0; wr = 0.d0
  !oxr = 0.d0; oyr = 0.d0; ozr = 0.d0
  
  !write(*,*) 'hvhg_uvwp'
  !(rev.416) buff1:phi, buff2:vor
  call hvhg_uvwp_fou(buff1,buff2,u00,w00,hv,hg,rf0u,rf0w,dvordy,   & 
       u1r,u2r,u3r,o1r,o2r,o3r,u1r,u2r,u3r,o1r,o2r,o3r, &
       tmpxzr,tmpx,wk1dr)

  deallocate(wk1dr,wk1dc)
  deallocate(tmpx)
  deallocate(tmpxzr)
  deallocate(rf0u,rf0w)
  deallocate(u1r,u2r,u3r,o1r,o2r,o3r)
  deallocate(buff1, buff2)

end subroutine calc_energy

subroutine hvhg_uvwp_fou(phic,ome2c,u00,w00,rhvc,rhgc,rf0u,rf0w,ome1c, &
     u1r,u2r,u3r,o1r,o2r,o3r,  &
     u1c,u2c,u3c,o1c,o2c,o3c,tmpxzr,tmpx,tmpy)

  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"
  
  complex*16 phic(0:mx1,0:mz1,jb:je),  & 
             ome1c(0:mx1,0:mz1,jb:je),ome2c(0:mx1,0:mz1,jb:je), &
             rhgc (0:mx1,0:mz1,jb:je),rhvc (0:mx1,0:mz1,jb:je)

  real*8 rf0u(0:my1),rf0w(0:my1),u00(0:my1),w00(0:my1),u00s(0:my1)

  !------ 6 * (mgalx+2) * mgalz buffer planes, THESE PLANES SHARE STORAGE !!
  real*8 u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz),  &
         o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz)
  complex*16 u1c(0:mx1,0:mz1),u2c(0:mx1,0:mz1),u3c(0:mx1,0:mz1),  &
             o1c(0:mx1,0:mz1),o2c(0:mx1,0:mz1),o3c(0:mx1,0:mz1)
  ! -------------------------------------------------------------------------
  real*8 divr(mgalx+2,mgalz)
  complex*16 divc(0:mx1,0:mz1)
  real*8 tmpxzr(mgalx+2,mgalz)
  real*8 tmpx(mgalx+2)
  real*8 tmpy(0:my1)

  real*8      up,vp,wp,o1p,o2p,o3p,uvr,uwr,vwr,eps,Pk
  complex*16  cc
  ! -------------------------------------------------------------------------
  integer i,j,k,iproc
  integer istat(MPI_STATUS_SIZE),ierr

  real*8 aa, v00, div_max, sytime, diff_xoff
  complex*16  shp,shm
  complex*16  sydt(0:mx1,0:my1),sydb(0:mx1,0:my1)

  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
  !  at this point:
  !  rhv is dv / dy -- F-F-P
  !  rhg is v -- F-F-P
  !  phi is nabla^2(v) -- F-F-P
  !  ome1 is d(omega_2)/dy -- F-F-P
  !  ome2 is omega_2 --F-F-P
  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

  ! ---------------- Initialize variables outside of the y loop
  !
  !--- Note: we need special treatments for the mean flow
  shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
  call derivyc(u00,rf0u,shp,shm,1,1,0)    ! ----- ome3_0 = -du0/dy ---
  !--- Note: derivative of fluctuation of u is pure-periodic
  call derivyc(w00,rf0w,shp,shm,1,1,0)    ! ----- ome1_0 =  dw0/dy ---
 
  do j=0,my1
     u00(j)=u00(j) !+s*y(j)  ! --- add shear --- now u00 is the mean shear flow
     rf0u(j)=rf0u(j) !+s  ! --- add shear --- now rf0u is du0/dy
  end do

 ! umaxl = 0d0
 ! umax = 0d0
  div_max=0.d0

  ! Take care ... check Ly ...
  sytime = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly


  DO J = JB,JE       !---------------------- start operating by planes
     !---------------------- 00 modes, u,w,ome1,ome3
     u1c(0,0) =  u00(j)
     u3c(0,0) =  w00(j)
     o3c(0,0) = -rf0u(j)
     o1c(0,0) =  rf0w(j)
     
     ! --------------------- computes non 0 modes of ome1, ome3, u1, u3 

     o3c(0,1:mz1) = -ome1c(0,1:mz1,j)/xgam(1:mz1)
     o1c(0,1:mz1) = - phic(0,1:mz1,j)/xgam(1:mz1)
     u3c(0,1:mz1) = - rhvc(0,1:mz1,j)/xgam(1:mz1)
     u1c(0,1:mz1) =  ome2c(0,1:mz1,j)/xgam(1:mz1)

     do k=0,mz1
        o3c(1:mx1,k) = (ome1c(1:mx1,k,j)*xgam(k)-phic(1:mx1,k,j)*xalp(1:mx1)) &
                       /(alp2(1:mx1)+gam2(k))
        o1c(1:mx1,k) = (ome1c(1:mx1,k,j)*xalp(1:mx1)+phic(1:mx1,k,j)*xgam(k)) &
                       /(alp2(1:mx1)+gam2(k))
        u3c(1:mx1,k) = (rhvc(1:mx1,k,j)*xgam(k)+ome2c(1:mx1,k,j)*xalp(1:mx1)) &
                       /(alp2(1:mx1)+gam2(k))
        u1c(1:mx1,k) = (rhvc(1:mx1,k,j)*xalp(1:mx1)-ome2c(1:mx1,k,j)*xgam(k)) &
                       /(alp2(1:mx1)+gam2(k))
     enddo

     ! -------------- copy v and omega_2 into their planes
     u2c = rhgc(:,:,j)
     o2c = ome2c(:,:,j)  

     ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
     !  at this point:
     !  u1 is u,  u2 is v, u3 is w
     !  o1 is omega_1,  o2 is omega_2,  o3 is omega_3
     !  all variables in Fourierx -- Fourier z -- Physical y
     ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
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

     ener(1)=ener(1) + up
     ener(2)=ener(2) + vp
     ener(3)=ener(3) + wp

     ener(4)=ener(4) + uvr
     ener(5)=ener(5) + uwr
     ener(6)=ener(6) + vwr

     ener(7)=ener(7) + o1p
     ener(8)=ener(8) + o2p
     ener(9)=ener(9) + o3p
     ! -----------------  do spectra some day
     ! ---------------- only if my plane contains spectra information

     ! add shear effect (only for iadd_mode = 2,3)
     o3c(0,0) = o3c(0,0) -s

     Pk = - (o3c(0,0))*(-uvr) + (o3c(0,0))*(o3c(0,0))*(1/Re) 

     ener2(1) = ener2(1) + Pk

     eps = sum(gam2*( u1c(0,:)*dconjg(u1c(0,:)) & 
                     +u2c(0,:)*dconjg(u2c(0,:)) & 
                     +u3c(0,:)*dconjg(u3c(0,:)) )) & 
          +sum( rhvc(0,:,j)*dconjg(rhvc(0,:,j)) + o3c(0,:)*dconjg(o3c(0,:)) )
     do k=0,mz1
        eps = eps + 2.d0*sum( (alp2(1:mx1)+gam2(k)) &
                             *( u1c(1:mx1,k)*dconjg(u1c(1:mx1,k)) & 
                               +u2c(1:mx1,k)*dconjg(u2c(1:mx1,k)) &
                               +u3c(1:mx1,k)*dconjg(u3c(1:mx1,k)) ) & 
                             +rhvc(1:mx1,k,j)*dconjg(rhvc(1:mx1,k,j)))
        cc = o1c(0,k) + xgam(k)*u2c(0,k)
        eps = eps + cc*dconjg(cc)
        do i=1,mx1
           cc = o1c(i,k) + xgam(k)*u2c(i,k)
           eps = eps + 2d0*cc*dconjg(cc)
           cc = o3c(i,k) - xalp(i)*u2c(i,k)
           eps = eps + 2d0*cc*dconjg(cc)
        enddo
     enddo
     ener(10) = ener(10) + eps
     o3c(0,0) = o3c(0,0) + s ! remove shear effect (only for iadd_mode = 2,3)
     ! /* check divergence in this plane */
     call compute_div(u1c,u3c,rhvc(0,0,j),divc) ! rhvc = dv/dy

     call fourxz(divc,divr,tmpxzr,1,1)       ! divergence
     do k=1,mgalz 
        do i=1,mgalx
           if (div_max.lt.abs(divr(i,k)))then
              div_max=abs(divr(i,k))              
           end if
        end do
     end do
     ! umax should be syncronized

  enddo
  
  !call MPI_ALLREDUCE(umaxl,umax,4,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  !if (myid.eq.0) write(*,*) ' umax check: ', umax
  if (myid.eq.0) write(*,*) 'Divergence check: ', div_max

end subroutine hvhg_uvwp_fou
