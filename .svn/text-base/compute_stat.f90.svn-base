program compute_stat

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

  ! reading parameter
  real*8 para(1:10),timep, rem1, time1, time2
  !real*8 kener(1:3),kener2(1:3),kener_fin(1:3)
  real*8 ke1, ke2, ke_iso, ke_fin, ke12r
  character*80 fname,filinm1

  integer ishoot_mode,nloop,iloop,igo,indx,lex

  integer istat(MPI_STATUS_SIZE),ierr

  ! /*   initializes everything    */
  call MPI_INIT(ierr)
  ! MPI_COMM_WORLD, MPI_COMM_SELF
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)

  ! version (revision number of Subversion)
  if(myid.eq.0) write(*,*) ' compute sta4_fixed, addsym chwk, rev.882'
  ! (should correspond to svn revision) 
  !--------------- initializes commons and things
  call initcr(0)
  explicit = 0 ! 0, implicit; 1, explicit.
  !
  idump_mode=0 ! normal dumping files mode
  nloop=1 ! default 1-loop
  para=0.d0
  !  
  if(myid.eq.0) then
     write(*,*) 'set UPO '
     read(*,*) para(1) ! time shear-period 
     read(*,*) para(2) ! iadd_sym  (iget_hv=1) rev589
  end if

  call MPI_BCAST(para,10,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  
  iadd_force=0
  if (iadd_force.ne.0) then
     if(myid.eq.0) write(*,*) ' iadd_force =', iadd_force
  end if
  !
  timep = para(1) ! time period for periodic orbit 
  !
  iadd_sym=int(para(2)) ! > 1: add symmetry
  !
  !nloop = int(para(2))
  !
  !
  noescru=1 ! 1: for no output and skip allocation of sp, spl
  nohist =0
  nocf = 1
  !
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  !vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;

  ! set run options
  !
  do while(1)
     if (myid.eq.0) then
        read(*,'(a)') fname
        !write(*,*) 'flist', fname
        indx=index(fname,'.flqt')-1
        filinp =  trim(fname(1:indx))
        write(*,*) 'recompute sta4, reading:',trim(filinp)
        inquire(file=filinp,exist=lex)
        if(lex.eq.0) then
           write(*,*) 'The upo file does not exist !!'
           stop
        end if

        !id22=int(filinp(end-2:end))
        read(filinp(indx-2:indx),*) id22
        write(*,*) ' i get, id22 = ', id22
        !filout = './pool/'//fname( (index(fname,"/",back=.true.)+1):indx) 
        !filout = './pool/'//fname( (index(fname,"/",back=.true.)+1):(index(filinp,".",back=.true.)-1)) 
        filout = fname(1:(index(filinp,".",back=.true.)-1)) 
        write(*,*) 'recompute sta4, basename to write:',trim(filout)
       
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ! reinitialization...
     vor=0.d0; phi=0.d0; u00=0.d0; w00=0.d0; chwk=0.d0;
     !
     if (myid.eq.0) write(*,*)  ' ====  iloop =',iloop, ' ==== '
     !
     call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
     time_ini=time
     xoffb_ini=xoffb
     xofft_ini=xofft
     ! reset re and Ly and the others
     re=Ree
     Ly=Lye
     alp=alpe
     gam=game
     call set_xyz_dy
     timep_sh = 2.d0*pi/alp/Ly*int(timep) ! fixed (rev.600, assuming s>0)
     if (myid.eq.0) write(*,*) & 
          &'set the integration time as a multiple of the shear-period: ', & 
          &'Tp=', timep_sh
     etime = time + timep_sh
     !
     if (myid.eq.0)  write(*,'(i3,a,3(1X,f14.6))') myid,' (Axz,Ayz,Rez) = ', Lx/Lz, Ly/Lz, re*Lz*Lz
     ! 
     write(*,*) 'nocf',nocf,nohist,noescru
     call calc_dns_disk3d(vor,phi,u00,w00,chwk)
     call write_stat
     ! 
  enddo ! end of NSTEP loop 
  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program compute_stat

subroutine calc_dns_disk3d(vor,phi,u00,w00,chwk)

  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"

  real*8  vor(0:2*my-1,0:mx1,kb:ke),phi(0:2*my-1,0:mx1,kb:ke)
  real*8  u00(0:my1),w00(0:my1)
  real*8  chwk(buffsize)
  real*8  end_time,y0
  integer iproc,ierr
  real*8, allocatable:: hv(:),hg(:),vorwk(:),phiwk(:),dvordy(:)
  real*8, allocatable:: wk1dr(:),wk1dc(:)

  real*8 chkre(0:numerop-1),chkalp(0:numerop-1), & 
       chkgam(0:numerop-1),chkLy(0:numerop-1),chks(0:numerop-1), &
       chky0(0:numerop-1)

  !  ---  allocate rest of buffers --- 
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0; 
  allocate( phiwk(buffsize),vorwk(buffsize),dvordy(buffsize) )
  phiwk = 0.d0; vorwk = 0.d0; dvordy = 0.d0;
  
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

  !
  !write(*,*) 'calc dns: time=',time, time_ini, xoffb_ini, xofft_ini
  if (myid.eq.0) write(*,*) '--- start DNS ---'
  call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,chwk) 
  if(myid.eq.0) write(*,*) '--- calc_dns_disk3d finished ...'

  ! --- deallocate ---  
  deallocate(phiwk,vorwk,dvordy)
  deallocate(hv,hg)

end subroutine calc_dns_disk3d
