module gmres
  !-------------- gmres paremeters -------------------------------------|
  integer k_gmres,nmax,nmax1,nf,nl,iuse_fou_mode,isol_mode,iget_hv,iskip_save
  ! iadd_sym is moved to running module (mod.f90)
  integer iarclength,iopt_matvec,ishifty,iuse_mpc,ifilter,iremove_singularity,& 
    &     iuse_scaling,iget_phase,iset_bulk0,iuse_ffty,iuse_hej,iuse_v, &
    &     iuse_compensatory_norm,iuse_QR_orth,iuse_us
  real(8)  arcl,arcl_ini,comp_norm
  integer nloopmax, itrmax
  integer ictp, icsx, icsy, icsz, icre

  character*128 cpara ! setting for arclength variable, Rez is implemented ...
  !                   and [Ayz, Axz, upr] are testing...from (rev.1200)

  ! --- gmres vector and arrays ---
  real(8),dimension(:),allocatable :: bb, dxp, Ax ! A is monodoromy matrix
  real(8) bnorm, dxpnorm, Axnorm
  real(8),dimension(:,:),allocatable :: vv, QQ, hh, gh
  real(8),dimension(:),allocatable :: vvb, tau
  
  real*8, allocatable:: cs(:),ss(:),ee(:),yg(:) 
  real*8, allocatable:: zz(:),mpc(:),mpc_in(:),rr(:)

  ! they are not changed while running

  real*8, allocatable:: x0(:),xf(:),dx0(:),dxf(:),dz0(:),dzf(:),dxdt0(:),dxdt(:),dfp(:)
  real*8  x0norm,xfnorm,dx0norm,dxfnorm,dxdt0norm,dxdtnorm,r0norm,rnorm
  real*8, allocatable:: dy0(:), dyf(:)
  
  ! save the boundary condition for arclength method along alp (Axz) (rev.1390)
  real*8 xofft0,xoffb0,xofftf,xoffbf
  real*8 xofft0_backup,xoffb0_backup,xofftf_backup,xoffbf_backup
  
  real*8, allocatable:: dxda(:), xstan(:), xm1(:), x0_backup(:), xf_backup(:), hej(:), xtmp(:), btmp(:)
  real*8  dxdanorm, xstannorm,  re_backup, Ly_backup
  real*8  rem1, Lym1, alpm1, uprimm1 ! testing for arclength (rev.1200)
  real*8  force_rollm1, xforcem1, zforcem1, Clesm1, damp_upm1, damp_aam1 ! testing for arclength (rev.1239)
  real*8 dsx0,dsx,dsy,dsy0,dsz0,dsz,dTp0,dTp,dss ! for relative periodic
  real*8 shiftx0,shiftx,shiftz0,shiftz,tfin0,tfin,shear,shear0 
  real*8 shifty0,shifty,shifty_ini,vbulk0,vbulk_ini ! rev935 (ishifty==1)
  real*8 extshiftx ! rev.1107 ==  shifty*(s*Tp/2)
  real*8 phasex,phasez ! constant phase in x and z (rev.979)
  real*8 phasey ! (rev.1101)

  !real*8 xoffb_ini, xofft_ini  ! move to the module 'bcs'
  ! for relative periodic
  !
  ! reading parameter
  real*8 para_arnoldi(1:20)
  real*8 para_newton(1:20)
 
  real*8 eps_gmresk,eps_sol,eps_jacobi,eps_newton,timep
  real*8 sca_u00, sca_w00, sca_vor, sca_phi, sca_sx, sca_sy,sca_sz, sca_time, sca_Ly, sca_re, sca_alp 
  !
  ! scaling factors
  !parameter(sca_vor=1.d0)
  !parameter(sca_phi=1.d0)    ! rev(951)
  parameter(sca_Ly=1.d0) 
  parameter(sca_time=1.d0)
  !parameter(sca_time=100.d0)
  parameter(sca_sx=1.d0) 
  parameter(sca_sy=1.d0) 
  parameter(sca_sz=1.d0) 
  !parameter(sca_re=100.d0) 
  parameter(sca_re=1.d0) 
  !parameter(sca_re=1.d-2) ! testing for a curved line 
  parameter(sca_alp=0.01d0) 
  ! ----- damped Newton parameters ---------------------------------|
  real*8 dampfac,damprate
  integer ndamp
  ! ----- hookstep paremeters ---------------------------------|
  integer iuse_hookstep
  integer nhookstep,ihookstep_equals_newtonstep,ihavebackup 
  real*8  delta,dmu,delta_ini,hs_cut !,dr
  real*8  deltaFuzz,deltaMin,deltaMax,dlambda,dlambdaMin,dlambdaMax
  real*8  dimprovReq,dimprovOk,dimprovGood,dimprovAcc

  ! for saving initial stats ! for iuse_fou_mode==3 (streak instability)
  real(8),dimension(:),allocatable :: vor_ini, phi_ini 
  real(8),dimension(:),allocatable :: u00_ini, w00_ini

end module gmres

module svd

  real*8 hs2max, hs2min, hs2
  
!  real*8, allocatable:: tmpvec(:),qq(:) ! QQ is used above
  real*8, allocatable:: tmpvec(:),q(:) !

  real*8, allocatable:: tmpy(:)
  real*8, allocatable:: tmph(:,:)
  real*8, allocatable:: hu(:,:), hs(:), hvt(:,:), svdwk(:)

end module svd

module precond

  real(8), allocatable:: filter(:)

  real(8),dimension(:),allocatable:: mpcx,mpctp0,mpctpf,mpcsx0,mpcsy0,mpcsz0, & 
       &                    mpcsxf,mpcsyf,mpcszf
  real(8),dimension(:),allocatable:: mpcx_in,mpctp0_in,mpctpf_in, &
       &                    mpcsx0_in,mpcsy0_in,mpcsz0_in, & 
       &                    mpcsxf_in,mpcsyf_in,mpcszf_in
  real(8), allocatable:: MTdx(:), Mdf(:)

  integer nband
  real(8), allocatable:: mpc_u00(:,:),mpc_w00(:,:),mpcl_u00(:,:),mpcl_w00(:,:)
  real(8), allocatable:: mpc_vor(:,:),mpc_phi(:,:),mpcl_vor(:,:),mpcl_phi(:,:)
  real(8), allocatable:: mpc_vorphi(:,:),mpc_phivor(:,:),mpcl_vorphi(:,:),mpcl_phivor(:,:)
  ! preconditioners  equivalent with mpcx... 
  real(8), allocatable:: vorp(:),phip(:),u00p(:),w00p(:) ! tmp arrays for a field
  real(8), allocatable:: vorq(:),phiq(:),u00q(:),w00q(:) ! tmp arrays for a field

end module precond

module eig
  !
  ! for eigen value problems
  ! using dgeev in lapack
  !
  real*8, allocatable:: Dnr(:),Dni(:)
  complex*16, allocatable:: Dnc(:)
  real*8, allocatable:: VR(:,:) ! right eigen vector
  real*8, allocatable:: wkeig(:)
  real*8, allocatable:: eigf(:,:)
  integer lwkeig

end module eig
