program arnoldi

  use ctes
  use running
  use bcs
  use gmres
  use eig
  use LES

  implicit none
  include "mpif.h"

  integer i,j,k,kk,idum,info,ilin,iopt,iproc
  real*8  val,dp,fac,mtime,eps_j,eps
  real*8, allocatable:: vor(:),phi(:),chwk(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: vorwk(:),phiwk(:) ! for escru
  real*8, allocatable:: hvtmp(:,:)

  real*8, allocatable:: buff1(:), buff2(:) ! used for addsym, too ...
  real*8, allocatable:: u00s(:),w00s(:) ! only for addsym

  real(8),dimension(:),allocatable:: wkqr

  real*8  randu,dummy
  character*80 fname
  character*4 ext1
  integer iolin, ioarn, knum,lwkqr

  integer istat(MPI_STATUS_SIZE),ierr

  ! /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)

  ! version (revision number of Subversion)
  if(myid.eq.0) write(*,*) 'Arnoldi method for HSF by A.S., Version: rev.',irev
  !if(myid.eq.0) write(*,*) 'use [hre2.dat], timep, dump_eigef, isolmode0, check convergence, 
  !error check getfil, addsym fixed, addsym8,iso_rr, addsym chwk, addsym set_conj'
  ! (should be svn revision) 
  !--------------- initializes commons and things
  call initcr(1) ! iopt=1, read the header of a field file.

  call set_options()
  !call set_arnoldi_options()

  isol_mode=0;
  iarclength=0;
  idump_mode=0 ! normal dumping files mode
  !

  ishifty=0; iset_bulk0=0; 
  iuse_mpc=0; 

  iremove_singularity=0
  iuse_scaling=0  ! this mode does not work well, (rev1137)
  ifilter=0;
  iuse_hookstep=1
  iskip_screenout = 1
  !
  ! these are set in set_options_newton() usually, but not in arnold_shear.f90, set here
  iuse_ffty=0; 
  iopt_matvec=1
  iuse_hej=1;
  !
  iuse_v=0;
  iuse_us=0;
  iuse_QR_orth=1; ! rev.1448, QR-decomposition improves the accuracy of orthogonarisation
  !
  if (iuse_LES.eq.1) then
     if(myid.eq.0) write(*,*) 'iuse_LES, iuse_v,iuse_us are on'
     iuse_v=1; iuse_us=1; 
  
  endif
  !
  ! --------------- allocates buffers

  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;
  !
  sca_u00 = 1.d0
  sca_w00 = 1.d0
  sca_vor = 1.d0
  sca_phi = 1.d0
  !  
  if(myid.eq.0) then
     write(*,*) 'set arnoldi parameter ...', & 
          '[uprim(for streak ini.), m*Ts, krylov_dim,', &
          ' eps_jacobi]'
     read(*,*) para_arnoldi(1) ! uprim is used for initial isotropic disturbance 
     read(*,*) para_arnoldi(2) ! time period n,  n*Lx/(SLy)
     read(*,*) para_arnoldi(3) ! number of Krylov subspace
     read(*,*) para_arnoldi(4) ! eps_jacobi
   !  read(*,*) para_arnoldi(5) ! iopt=1: only for streak instability
   !                            ! iopt=2: for linear stability of UPOs
  end if
  call  MPI_BCAST(para_arnoldi,20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  !if(myid.eq.0) write(*,*) 'go to read file ...'
  !iadd_force=0
  !iadd_damping=0
  if (iadd_force.ne.0) then
     if(myid.eq.0) write(*,*) ' iadd_force =', iadd_force
  end if
  if (iadd_damping.ne.0) then
     if(myid.eq.0) write(*,*) ' iadd_damping '
  end if
  uprim = para_arnoldi(1)
  timep = para_arnoldi(2)
  k_gmres = int(para_arnoldi(3))
  eps_jacobi = para_arnoldi(4)
  ! 
  !iopt=int(para_arnoldi(5))
  iopt=2

  iget_hv=0;
  if (iopt.eq.1) then
     write(*,*) 'iopt=1, isol_mode, iuse_fou_mode should be checked...they are old format, stop'
     stop
     !iuse_fou_mode=int(para_arnoldi(5)) 
     !iuse_fou_mode=3 (do not count x0modes) !
     isol_mode = 0 ! only for arnoldi
     timep_sh = timep
  elseif (iopt.eq.2) then
     ! for analysing UPO
     iuse_fou_mode=2   
     ! this should be 2, (1 does not work, 3 is for fixing x0modes)
     isol_mode = 0 ! only for arnoldi
     timep_sh = 2.d0*pi/alp/Ly*timep
  endif
  !
  if (myid.eq.0) write(*,*) & 
       'isol_mode=', isol_mode, 'iuse_fou_mode=', iuse_fou_mode
  if (myid.eq.0) write(*,*) & 
       &'UPO mode: set the integration time as a multiple of the shear-period: ', & 
       &'Tp=', timep_sh
  if (iopt.eq.1) then
     ! generate i.c.
     !-- streak profile: U=s*y+du*cos(gam*z)
     call getini_streak(vor,phi,u00,w00,time)
  else
     call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
  endif

  call arnoldi_ini()

  if(myid.eq.0) write(*,*) '... file read'

  if (iadd_sym.ge.5) then

     write(*,*) 'iadd_sym, phase shift',iadd_sym, sym_shiftx*Lx, sym_shiftz*Lz 
     call phase_shift_fou(vor,phi,sym_shiftx*Lx,sym_shiftz*Lz)
          
  end if

  ! --- setting arrays and vectors for dynamical system ---
  ! and allocated gmres vectors only by master
  ! --- allocate vv hh etc... ! 
  call gmres_ini(vor,phi,u00,w00) ! set_conj (rev.1155~)
  if(myid.eq.0) write(*,*) '... gmres allocated '
  ! set x0, xf as in make_mat
  call throwing(vor,phi,u00,w00,chwk,x0,1) ! --> x0 (base flow) ! BUG fixed (rev.1289) with option 1
  call calc_dns_disk3d(vor,phi,u00,w00,chwk) ! set_conj (rev.1155~)
  ! mtime=etime-time_ini
  ! call moving_frame(vor,phi,mtime) !???
  call throwing(vor,phi,u00,w00,chwk,xf,1) ! --> xf (time-integrated base flow)
  !
  ! set initial random vector of disturbance
  if(myid.eq.0) write(*,*) 'set random disturbance'
  idum=-123456
  do i=1,nl 
     rr(i)=randu(idum)
  end do

  ! top hat option: 1, R&M [8 16]; 2, [4 32]
  !call getini_tophat(vor,phi,u00,w00,chwk,time,2+iuse_LES) ! for pure-periodic b.c.
  fac = (xoffb_ini/Lx - int(xoffb_ini/Lx) )*Lx/Ly
  if (abs(fac).gt.0.5d0*Lx/Ly) then
     fac=(Lx/Ly+fac)
     call moving_frame(vor,phi,fac,+1)  
  else
     call moving_frame(vor,phi,fac,-1)    
  end if
  if (myid.eq.0) write(*,*)  '      the noise is distorted to satisfy B.C.', xoffb_ini,fac


  if (iadd_sym.ge.5) then

     !call backthrowing(vor,phi,u00,w00,chwk,rr,1)
     ! --- add new buffer to save change (take care for the confusing names)---
     allocate (buff1(buffsize),buff2(buffsize))
     buff1 = 0.d0; buff2 = 0.d0;
     allocate (u00s(0:my1),w00s(0:my1))
     u00s = 0.d0; w00s = 0.d0;  
     if (myid.eq.0) write(*,*) 'addsym to rr, iadd_sym=',iadd_sym
     call add_sym(vor,phi,u00,w00,buff1,buff2,u00s,w00s,chwk,iadd_sym) ! set_conj (rev.1158)
     rr=0.d0 ! reset rr to remove spurias nonsymmetric elements
  
     !call apply_filtering(vor,phi,u00,w00,3) ! apply compact 2/3-like filter 
     call throwing(vor,phi,u00,w00,chwk,rr,1)
     
     deallocate(u00s,w00s)
     deallocate(buff1,buff2)
     ! note: vor and phi are destroyed
  else
     call throwing(vor,phi,u00,w00,chwk,rr,1)
  end if
  !
  !rr=x0-xf ! for testing.
  !rr=x0
  !rr(1)=randu(idum)
  ! /* compute x0,xf norm */
  call norm(nl,xf,xfnorm)
  call norm(nl,x0,x0norm)  
  call norm(nl,rr,r0norm)  
  ! --- the arnoldi process ---
  if(myid.eq.0) write(*,*) 'Arnoldi, initial norm r0,x0,xfnorm:',r0norm,x0norm,xfnorm
  if(r0norm.lt.tiny) stop
  fac=1.d0/r0norm
  call vecscal(nl,rr,fac,vv) ! vv(:,1)=r(:)/r0

  do j=1,k_gmres  ! dimension of Krylov sub-space
     if(myid.eq.0) write(*,*) '-----------------------------------------------'
     if(myid.eq.0) write(*,*) '   Arnoldi iteration:', j
     if(myid.eq.0) write(*,*) '-----------------------------------------------'

     eps_j=eps_jacobi
     if (eps_j.gt.0.1) then
        if(myid.eq.0) then 
           write(*,*) ' Arnoldi: checking linearity'
           fname=filout(1:index(filout,' ')-1)//'.lincheck'
           open(iolin,file=fname,status='unknown')
        end if
        eps_j=1.d0
        do ilin=1,13 ! check linearity of matvec           
           call matvec_dns_disk3d(vor,phi,u00,w00,chwk,vv(1,j),vv(1,j+1), &
                eps_j,rnorm) ! iorder = 1, isol_mode = 0 (only for arnoldi)
           if(myid.eq.0) write(iolin,*) eps_j, rnorm 
           eps_j=eps_j/(10.d0)
        enddo ! check linearity
        if(myid.eq.0) then 
           close(iolin)
           write(*,*) ' Arnoldi: check linearity (eps_jacobi, norm) in', trim(fname)
        endif
        stop 
     else ! /* check-linearity */
        call matvec_dns_disk3d(vor,phi,u00,w00,chwk,vv(1,j),xtmp, & 
             eps_j,rnorm) 
     end if 

     if (iuse_QR_orth.eq.0) then
        vv(:,j+1)=xtmp(:)
        do i=1,j
           call dotp(nl,vv(1,j+1),vv(1,i),dp) !h_ij=(vv_j+1,vv_i)
           hh(i,j)=dp
           call vecadd(nl,vv(1,j+1),-dp,vv(1,i),vv(1,j+1))
        end do
        call norm(nl,vv(1,j+1),rnorm)
        hh(j+1,j)=rnorm
        fac=1.d0/rnorm
        call vecscal(nl,vv(1,j+1),fac,vv(1,j+1))

     elseif (iuse_QR_orth.eq.1) then
        
        QQ=vv
        QQ(:,j+1)=xtmp(:)
        if (myid.eq.0) then
           ! check the optimal size
           allocate(wkqr(20*(k_gmres+1)))
           call DGEQRF(nl, j+1, QQ, nl, tau, wkqr, -1, info)
           lwkqr=wkqr(1)
           if (lwkqr.gt.20*(k_gmres+1)) then
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
           deallocate(wkqr)
           hh(1:j+1,j)=zz(1:j+1) 
           ! store the upper Hessian matrix for eigenvalue computation. 
        end if        
        call MPI_BCAST(hh(1,j),j+1,MPI_DOUBLE_PRECISION,0, & 
             &         MPI_COMM_WORLD,ierr)

     else
        write(*,*) 'orthogonalization error in gmres_k_shear, stop'
        stop
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     !if ((myid.eq.0).and.(mod(j,1).eq.0)) then
     !   write(*,*) 'Arnoldi writing hh(56) for checking convergence '
     !   write(ext1,'(i4.4)') j
     !   ioarn=56
     !   fname=filout(1:index(filout,' ')-1)//'_kth_'//ext1//'.arn'
     !   open(ioarn,file=fname,status='unknown',form='unformatted')
     !   write(ioarn) mgalx,my,mgalz,numerop
     !   write(ioarn) time_ini,Re,alp,Ly/pi,gam,s,chi,mx,my,mz  
     !   write(ioarn) xoffb_ini,xofft_ini,dummy,dummy
     !   para_arnoldi(2)=timep_sh 
     !   para_newton(3)=dfloat(j) ! rev656
     !   para_arnoldi(6)=dfloat(nmax)
     !   write(ioarn) para_arnoldi(1:20)
     !   ![dU, timep_sh, k_gmres, eps_jacobi, iuse_fou_mode] and nmax. 
     !   ! the others are buffer.
     !   write(ioarn) ((hh(i,k),i=1,k_gmres+1),k=1,k_gmres)
     !   ! skip writing the orthogonal matrix
     !   !write(ioarn) ((vv(i,k),i=1,nl),k=1,k_gmres+1)
     !   close(ioarn)
     !   !
     !end if
     ! skip, since Arnoldi
     !c        /* upper Hessenberg to upper triangle matrix */
     !c        /*                 by Givens rotation matrix */
  end do


  if (myid.eq.0) then
     write(*,*) 'Arnoldi finished, writing hh(56) vv(57) for test'
     ioarn=56
     fname=filout(1:index(filout,' ')-1)//'.arn'
     open(ioarn,file=fname,status='unknown',form='unformatted')
     write(ioarn) mgalx,my,mgalz,numerop
     write(ioarn) time_ini,Re,alp,Ly/pi,gam,s,chi,mx,my,mz  
     write(ioarn) xoffb_ini,xofft_ini,dummy,dummy
     para_arnoldi(2)=timep_sh 
     para_arnoldi(3)=k_gmres 
     para_arnoldi(6)=dfloat(nmax)
     write(ioarn) para_arnoldi(1:10) ! now keep 10
     ![dU, timep_sh, k_gmres, eps_jacobi, iuse_fou_mode] and nmax. 
     ! the others are buffer.
     write(ioarn) ((hh(i,k),i=1,k_gmres+1),k=1,k_gmres)
     ! skip writing the orthogonal matrix
     !write(ioarn) ((vv(i,k),i=1,nl),k=1,k_gmres+1)

     ! ---- rev. 1433  -----
     write(ioarn) CFL, uprim, vbulk, xforce ! rev934
     write(ioarn) irev,nopt,nparams
     write(ioarn) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
          &      iuse_LES,iadd_visc,idynamic ! rev1200 
     write(ioarn) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(1:3)
     !
     close(ioarn)
     !
  end if
  !
  ! restore Hessenberg matrix
  allocate(hvtmp(k_gmres,k_gmres))
  hvtmp=hh(1:k_gmres,1:k_gmres)
  !
  ! eig module
  lwkeig=4*k_gmres ! not optimal value
  allocate(Dnr(k_gmres),Dni(k_gmres),wkeig(lwkeig))
  allocate(VR(k_gmres,k_gmres))
  allocate(eigf(nl,k_gmres))
  ! check optimal work array size
  call dgeev('n','v',k_gmres,hvtmp,k_gmres,Dnr,Dni,dummy,1, &
       VR,k_gmres,wkeig,-1,info)  
  ! write(*,*) myid,'dgeev optimal work array length', wkeig(1)
  lwkeig=wkeig(1)
  deallocate(wkeig)
  !
  allocate(wkeig(lwkeig))
  ! Solve the eigenvalue problem (need: lapack) ![Vn,Dn]=eig(hh);
  call dgeev('n','v',k_gmres,hvtmp,k_gmres,Dnr,Dni,dummy,1, &
       VR,k_gmres,wkeig,lwkeig,info)
  !*          If the j-th eigenvalue is real, then v(j) = VR(:,j),
  !*          the j-th column of VR.
  !*          If the j-th and (j+1)-st eigenvalues form a complex
  !*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
  !*          v(j+1) = VR(:,j) - i*VR(:,j+1).
  !
  ! write(*,*) myid,'dgeev info', info
  ! matrix vector product
  do k=1,k_gmres     
     do i=1,nl
        eigf(i,k)=sum(vv(i,1:k_gmres)*VR(:,k)) ! eigen function (a flow field) 
        ! note that vv(nmax1,k_gmres+1)
     end do
  end do

  do k=1,k_gmres
     if ((isol_mode.eq.3).or.(isol_mode.eq.4)) then
        Dnr(k) = Dnr(k)+1.d0
     end if
     if (myid.eq.0)  write(*,*) 'B: Dnr Dni',Dnr(k),Dni(k)
  end do
  allocate(Dnc(k_gmres))     
  do k=1,k_gmres
     Dnc(k)=dcmplx(Dnr(k),Dni(k))
     if (myid.eq.0)  write(*,*) 'A: sigma ',log(Dnc(k))/(etime-time_ini)
  end do
  deallocate(Dnc)  
  !
  ! ---- matlab memo ----
  !vvm1=vv(1:nmax1,1:kk); % Q'_n
  !
  !% check eigen values
  !
  !lamb=diag(Dn)        %  du/dt   = f_{NS} (u) = (df/du)*u =A*u 
  !lama=log(lamb)/etime %  u(t+Tp) = \int_t^{t+Tp} f_{NS} (u) dt = g(u)
  !                     %             = (dg/du)*u = B*u
  !
  !f1=vvm1*Vn(:,iv); % eigen function
  !
  ! writing eigen function corresponding to max eigen value
  allocate(vorwk(buffsize), phiwk(buffsize))  
  if(myid.eq.0) then 
     write(*,*) &
       'writing eigen function corresponding to all the eigen values'
     write(*,*) &
       '        note: filout and id22 are overwritten'
  end if
  filout = filout(1:index(filout,' ')-1)//'_eigf'
  do k=1,k_gmres
     if ((Dnr(k)**2+Dni(k)**2).gt.0.99d0) then
        ! set the number of unstable modes, including conjugate 
        knum=k
     end if
  end do

  do k=1,knum
     id22=k
     if (myid.eq.0)  write(*,*) &
          'save eigen function of B:Dnr Dni',Dnr(k),Dni(k)
     call backthrowing(vor,phi,u00,w00,chwk,eigf(1,k),1) ! bug fixed at rev1295
     call escru(vor,phi,u00,w00,vorwk,phiwk,chwk)
  end do
  deallocate(eigf)


  deallocate(vor, phi, chwk, u00, w00)
  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program arnoldi

subroutine arnoldi_ini()

  use ctes
  use running
  use bcs
  ! error check to read the UPOs.
  if (abs(Lye-Ly).gt.tiny) then 
     write(*,*) myid,' arnoldi: Ly is different, STOP!!',Lye,Ly
     stop
  end if
  if (abs(Ree-re).gt.tiny) then
     write(*,*) myid,' arnoldi: re is different, STOP!!',Ree,re
     stop
  end if
  if (abs(alpe-alp).gt.tiny) then
     write(*,*) myid,' arnoldi: alp is different, STOP!!',alpe,alp
     stop
  end if
  if (abs(game-gam).gt.tiny) then
     write(*,*) myid,' arnoldi: gam is different, STOP!!',game,gam
     stop
  end if

  ! ignores some newton options...

end subroutine arnoldi_ini
