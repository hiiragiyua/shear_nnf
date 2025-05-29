!=======================================================================
!       PAQUETE DE SUBRUTINAS PARA DIFERENCIAS FINITAS COMPACTAS
!       O.F, nov 2002
!       S.H. may 2005
!       J.J. july 2007,  uniform grid for disk
!       A.S. August 2011, uniform grid for shear-periodic (disk) 
!                         using only central part.
!=======================================================================
module yderiv
  real*8, allocatable:: pr1(:), rh1(:) ! (constant) B u' = A u
  real*8, allocatable:: pr2(:), rh2(:) ! (constant) B u''= A u
  real*8, allocatable:: lapmat(:)      !  (A+ B k2*I) u = B Rhs

  !7*my real*8  + 3*my complex*16 for Cholesky decomposition
  ! c1r => first derivatives, c2r => second derivatives, chr => helmholtz
  !real*8, allocatable:: c1r(:,:),c2r(:,:),chr(:,:) 
  !complex*16, allocatable:: c1i(:,:),c2i(:,:),chi(:,:) 

  real*8, allocatable:: a1(:),b1(:),d1(:)
  !complex*16, allocatable:: a1(:),b1(:),d1(:)
  complex*16, allocatable:: f1(:),g1(:) 
  complex*16, allocatable:: cf1(:),cg1(:)

  real*8, allocatable:: a2(:),b2(:),d2(:)
  !complex*16, allocatable:: a2(:),b2(:),d2(:)
  complex*16, allocatable:: f2(:),g2(:) 
  complex*16, allocatable:: cf2(:),cg2(:)

  real*8, allocatable:: ah(:),bh(:),ch(:),dh(:)
  !complex*16, allocatable:: ah(:),bh(:),ch(:),dh(:)
  complex*16, allocatable:: eh(:),fh(:),gh(:) 
  complex*16, allocatable:: ceh(:),cfh(:),cgh(:)

  complex*16, allocatable:: wk1(:), wk2(:), ur(:)
  complex*16, allocatable:: uc(:), duc(:), pc(:)
  real*8  h
  integer N1,N2, my

end module yderiv

!**********************************************************************
! prederiv1:
!
! Computes the coefficients for the compact finite differences
!     for the first derivatives of bandwidth 2*Nb+1
!     with shear-periodic boundary condition. 
!
! References:  /shear/trunk/Doc/maple/test.ms   
!     and      /shear/trunk/Doc/diffyc.tex   
!      by A.S.
!
! Scheme is: B f' = A f
!            sum_{j=-Nb}^Na b_ij f'_{i+j} = 
!            sum_{j=-Na}^Nb a_ij f_{i+j}
!      
! Although this is optimized for Na=4, Nb=2
!
!*********************************************************************
subroutine prederiv1
  use yderiv
  implicit none

  integer N, M

  integer j,k,kk,kp
  real*8  pi
  complex*16 zrr,zii

  ! ----------- define overall parameters ----
  N=4    ! semi-band width for r.h.s
  M=2    ! semi-band width for l.h.s
  ! number of spectral points ! not implemented.
  pi=4d0*atan(1d0)
  !write(*,*) '  prederiv1, semi-band width: rhs 4, lhs 2 (penta-CFD)'
  allocate ( pr1(-M:M), rh1(-N:N) ) !rh1(3),rh1(-3) is 0.d0
  allocate ( a1(my),b1(my),d1(my))
  allocate ( f1(my),g1(my),cf1(my),cg1(my) )
  !write(*,*) '  prederiv1: allocated'
  pr1 = 0d0; rh1 = 0d0
  a1  = 0d0; b1  = 0d0; d1 = 0d0; 
  !a1  = dcmplx(0d0,0.d0); b1  = dcmplx(0d0,0.d0); d1 = dcmplx(0d0,0.d0); 
  f1  = dcmplx(0d0,0.d0); g1  = dcmplx(0d0,0.d0)
  cf1 = dcmplx(0d0,0.d0); cg1 = dcmplx(0d0,0.d0)
  ! deallocate??
  N1  = N  ! N1???
  ! set B  ( penta-diagonal, B u' = A u )
  pr1(0) = 1.d0;
  pr1(1) = 8.d0/15.d0; 
  pr1(2) = 1.d0/15.d0;
  pr1(-2) = pr1(2);
  pr1(-1) = pr1(1);
  ! set A  ( B u' = A u )
  rh1(0) = 0.d0 
  rh1(1) = 154.d0/225.d0; 
  rh1(2) = 91.d0/450.d0; 
  rh1(3) = 2.d0/525.d0; 
  rh1(4) = -1.d0/12600.d0;
  rh1(-4) = -rh1(4);
  rh1(-3) = -rh1(3);
  rh1(-2) = -rh1(2);
  rh1(-1) = -rh1(1);

  !write(*,*) '  check pr1=',pr1 
  !write(*,*) '  check rh1=',rh1 

  !  ------ scale and prefactor 
  rh1 = rh1/h
  zrr = dcmplx(1.0,0.d0)
  call modChol5cyclic_cmplx(a1,b1,d1,f1,g1,cf1,cg1,pr1(-2:2),my,zrr)
  !call modChol7cyclic_cmplx(a1,b1,c1,d1,e1,f1,g1,ce1,cf1,cg1,pr1,my,zrr)
  ! shear will be included by calling modChol_addshear() in each time step. 

  !write(*,*) 'a1',a1 
  !write(*,*) 'b1',b1 
  !write(*,*) 'c1',c1 
  !write(*,*) 'd1',d1 
  !write(*,*) 'e1',e1 
  !write(*,*) 'f1',f1 
  !write(*,*) 'g1',g1 
  
end subroutine prederiv1

!**********************************************************************
! prederiv2:
!
! Computes the coefficients for the compact finite differences
!     for the second derivatives of bandwidth 2*Nb+1
!     with shear-periodic boundary condition. 
!
! References:  /shear/trunk/Doc/maple/test.ms   
!     and      /shear/trunk/Doc/diffyc.tex   
!      by A.S.
!
! Scheme is: B f' = A f
!            sum_{j=-Nb}^Na b_ij f''_{i+j} = 
!            sum_{j=-Na}^Nb a_ij f_{i+j}
! 
! Although this is optimized for Na=3, Nb=2, 
! we allocate pr2,rh2(-3:3) for computing lapmat in solving Helmholtz equation
! ( Laplacian solver for hepta-diagonal cyclic matrix ). 
!
!
!*********************************************************************
subroutine prederiv2
  use yderiv
  implicit none

  integer N, M
  integer,allocatable:: jj(:)
  real*8, allocatable:: mat(:,:),vec(:),kap(:)

  integer j,k,kk,kp
  real*8  pi
  complex*16 zrr,zii
  ! ----------- define overall parameters ----
  N=3    ! semi-band width for r.h.s
  M=2    ! semi-band width for l.h.s
  ! number of spectral points ! not implemented.
  pi     = 4d0*atan(1d0)
  !write(*,*) '  prederiv2, semi-band width: rhs 3, lhs 2 (penta-CFD)'
  allocate ( pr2(-N:N), rh2(-N:N), lapmat(-N:N) )
  allocate ( a2(my),b2(my),d2(my) )
  allocate ( f2(my),g2(my),cf2(my),cg2(my) )
  allocate ( ah(my),bh(my),ch(my),dh(my) )
  allocate ( eh(my),fh(my),gh(my),ceh(my),cfh(my),cgh(my) )
  !write(*,*) '  prederiv2: allocated'

  pr2 = 0d0; rh2 = 0d0; lapmat = 0d0;

  a2  = 0d0; b2  = 0d0; d2 = 0d0; 
  !a2  = dcmplx(0d0,0.d0); b2  = dcmplx(0d0,0.d0); d2 = dcmplx(0d0,0.d0); 
  f2  = dcmplx(0d0,0.d0); g2  = dcmplx(0d0,0.d0)
  cf2 = dcmplx(0d0,0.d0); cg2 = dcmplx(0d0,0.d0)

  ah  = 0d0; bh  = 0d0; ch  = 0d0; dh = 0d0; 
  eh  = dcmplx(0d0,0.d0); fh  = dcmplx(0d0,0.d0); gh = dcmplx(0d0,0.d0)
  ceh = dcmplx(0d0,0.d0); cfh = dcmplx(0d0,0.d0); cgh = dcmplx(0d0,0.d0)

  N2  = N
  
  ! B in B*D2 = A*Ru
  pr2(0) = 1.d0
  pr2(1) = 334.d0/899.d0; 
  pr2(2) = 43.d0/1798.d0;
  pr2(3) = 0.d0; ! dummy for computing lapmat
  pr2(-3) = pr2(3)
  pr2(-2) = pr2(2)
  pr2(-1) = pr2(1)
  
  rh2(0)= -14335.d0/8091.d0; 
  rh2(1)= 1065.d0/1798.d0; 
  rh2(2)= 519.d0/1798.d0; 
  rh2(3)= 79.d0/16182.d0;
  rh2(-3) = rh2(3)
  rh2(-2) = rh2(2)
  rh2(-1) = rh2(1)

  !write(*,*) '  check pr2=',pr2 
  !write(*,*) '  check rh2=',rh2 

  !  ------ scale and prefactor 
  rh2 = rh2/h**2
  zrr = dcmplx(1.d0,0.d0)
  call modChol5cyclic_cmplx(a2,b2,d2,f2,g2,cf2,cg2,pr2(-2:2),my,zrr)
  !call modChol7cyclic_cmplx(a2,b2,c2,d2,e2,f2,g2,ce2,cf2,cg2,pr2,my,zrr)
  ! shear will be included by calling modChol_addshear() in each time step.

end subroutine prederiv2

subroutine visc3d(u,visc,ddu,shp,shm,m,rk,scal,what,ias)
  use yderiv
  implicit none

  integer m,ias
  real*8 u(m,my),visc(m,my),ddu(m,my) ! ddu is tmp2d
  real*8 rk,scal
  complex*16 shp,shm
  character*1 what

! 2nd derivative in 'y':
  call derivyc(u,ddu,shp,shm,2,m,ias)

! 2nd derivative in 'x' and 'z': d^2(u)/dx^2 + d^2(u)/dz^2 -> visc 
  if (what.eq.'n') then
     visc = scal*(-rk*u + ddu)
  else
     visc = visc + scal*(-rk*u + ddu)
  endif

end subroutine visc3d

!************************************************************************
!                                                                       !
!    subroutine derivyc                                                 !
!                                                                       !
!    computes first or second derivatives in y of u, with size(u)=(my)  !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be complex or real, size my  !
!    Output:                                                            !
!        du: First derivative of u  (may be the same as u)              !
!                                                                       !
!************************************************************************

subroutine derivyc(u,du,shp,shm,iord,m,ias)
  use yderiv

  implicit none
  integer i,j,m,jm3,iord,na,ias
  real*8 u(m,my),du(m,my)
  complex*16 shp,shm,const

  ! wrapping to complex*16
  if (m.eq.1) then
     uc(:) = dcmplx(u(1,:),0.d0)
  else
     uc(:) = dcmplx(u(1,:),u(2,:))
  endif

  if (iord.eq.1) then    !!!!!  first derivative

     !write(*,*) 'wk1'
     na=4
     call matvecNcyclic(rh1(-na:na),uc,ur,wk1,shp,shm,my,m,na)

     if (m.eq.1) then
        call solve5cyclic(a1,b1,d1,f1,g1,cf1,cg1,wk1,duc,wk2,my)
        !call solve7cyclic(a1,b1,c1,d1,e1,f1,g1,ce1,cf1,cg1,wk1,duc,wk2,my)
        du(1,:) = dreal(duc)
     else
        !
        if (ias.eq.1) then
           dh=d1;fh=f1;gh=g1;
           call modChol5_addshear(a1,b1,dh,fh,gh,cfh,cgh,pr1(-2:2),my,shp)
        end if
        call solve5cyclic(a1,b1,dh,fh,gh,cfh,cgh,wk1,du,wk2,my)

     endif
  else                 !!!!!  second derivative 

     na=3
     call matvecNcyclic(rh2(-na:na),uc,ur,wk1,shp,shm,my,m,na)

     if (m.eq.1) then
        call solve5cyclic(a2,b2,d2,f2,g2,cf2,cg2,wk1,duc,wk2,my)
        !call solve7cyclic(a2,b2,c2,d2,e2,f2,g2,ce2,cf2,cg2,wk1,duc,wk2,my)
        du(1,:) = dreal(duc)
     else
        if (ias.eq.1) then
           dh=d2;fh=f2;gh=g2;
           call modChol5_addshear(a2,b2,dh,fh,gh,cfh,cgh,pr2(-2:2),my,shp)
        end if
        call solve5cyclic(a2,b2,dh,fh,gh,cfh,cgh,wk1,du,wk2,my)
     endif
  endif

end subroutine derivyc


!************************************************************************
!                                                                       !
!    subroutine lapcdy                                                  !
!                                                                       !
!    solves laplace's equation (d_yy-k2)p = u (helmholtz equation)      !
!                    shear-periodic  b.c.                               !
!    Input:                                                             !
!         u: rhs. Assumed to be complex(m=2) or real(m=1) of size my    !
!    Output:                                                            !
!         p: solution  (can be same as u)                               !
!                                                                       !
!************************************************************************

subroutine lapcdy(u,k2,p,shp,shm,m)
  use yderiv

  implicit none
  integer i,j,jm3,m,na
  real*8 k2,bnorm
  complex*16 shp,shm
  real*8 u(m,my),p(m,my)

  lapmat=rh2-k2*pr2
  !if ((k2.gt.-1.d0-tiny).and.(k2.lt.-1.d0+tiny)) then
  !   write(*,*) 'lapcdy: singular!'
  !   stop
  !endif
  !write(*,*) 'lapcdy: check shp=',shp

  ! shifting by shear is included

  !-------PREPARO EL TERMINO INDEPENDIENTE [pr2]*{u} 
  if (m.eq.1) then
     uc(:) = dcmplx(u(1,:),0.d0)
  else
     uc(:) = dcmplx(u(1,:),u(2,:))
  endif

  !bnorm = sqrt(sum(uc*dconjg(uc))) 
  !if (bnorm.lt.tiny) then
  !   write(*,*) 'diffy_shear: skip lapcdy, bnorm=', bnorm
  !   p = 0.d0
  !   return
  !endif

  call modChol7cyclic_cmplx(ah,bh,ch,dh,eh,fh,gh,ceh,cfh,cgh,lapmat,my,shp)
  na=2
  call matvecNcyclic(pr2(-na:na),uc,ur,wk1,shp,shm,my,m,na)

  if (m.eq.1) then
     ! no need to modify Cholesky elements  
     call solve7cyclic(ah,bh,ch,dh,eh,fh,gh,ceh,cfh,cgh,wk1,pc,wk2,my)
     p(1,:) = dreal(pc)
  else
     ! no need to modify Cholesky elements  
     call solve7cyclic(ah,bh,ch,dh,eh,fh,gh,ceh,cfh,cgh,wk1,p,wk2,my)
  end if  

end subroutine lapcdy


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    prepares commons and things for using cfdiff and laps later
!                             jj/may/05
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine derivadas(n,dy)  
  use yderiv
  implicit none
  integer n,my8
  real*8  dy,Length

  my = n
  h  = dy
!     ---  start doing something -------------
  my8 = my+4+4
  allocate (wk1(my), wk2(my), ur(my8))
  allocate (uc(my),duc(my),pc(my))
  wk1 = 0.d0; wk2 = 0.d0; ur = 0.d0
  uc = 0.d0; duc = 0.d0; pc = 0.d0
  ! ur is a working vector for matvec  
  !write(*,*) '  in derivadas: wk1 wk2 ur, allocated'
  call prederiv1
  !write(*,*) '  prederiv1 in derivadas: OK'
  call prederiv2
  !write(*,*) '  prederiv2 in derivadas: OK'

end subroutine derivadas

!***********************************************************************
!    subroutine for cyclicMatrix - vector product 
!    including shifting by mean shear
!***********************************************************************
subroutine matvecNcyclic(A,uc,ur,bb,shp,shm,my,m,na)
  !
  ! A is cyclic (2*na+1) band matrix
  !
  implicit none
  integer i,j,m,my,na,jmN 
  real*8 A(-na:na)
  complex*16 shp,shm,tmp
  complex*16 uc(my),ur(my+8) 
  complex*16 bb(my)

  if (m.eq.1) then
     ! for zeromode in x-dir 
     ur(1:na)=uc(my-na+1:my);
     ur(na+1:my+na)=uc(:);
     ur(my+na+1:my+2*na)=uc(1:na);

  else
     shm=dconjg(shp)
     ur(1:na)=uc(my-na+1:my)*shm;
     ur(na+1:my+na)=uc(:);
     ur(my+na+1:my+2*na)=uc(1:na)*shp;

  endif

  do j=1+na,my+na
     jmN = j-na
     tmp=dcmplx(0.d0,0.d0)
     do i=-na,na
        tmp = tmp+ dcmplx(A(i),0.d0)*ur(j+i)
     end do
     bb(jmN) = tmp
     !bb(jmN) = sum(dcmplx(A(:),0.d0)*ur(j-na:j+na))
     !write(*,*) bb(jmN)
  enddo

end subroutine matvecNcyclic

!***********************************************************************
!    subroutine for modified Cholesky decomposition
!***********************************************************************
subroutine modChol5cyclic_cmplx(La,Lb,Dd,Lf,Lg,cLf,cLg,a,N,shp)
  integer N
  real*8 La(N),Lb(N),Lc(N),Dd(N)
  !complex*16 La(N),Lb(N),Dd(N)
  complex*16 Lf(N),Lg(N),cLf(N),cLg(N)
  real*8 a(-2:2) ! penta-diagonal cyclic matrix
  complex*16 shp,tmp ! note that shp is complex*16
  ! do not include real*8 value in the last arguments, (,... ,1.d0) 
  integer i,j,k
  real*8 det

  j=1;
  Dd(j)=a(0);
  La(j+1)= 1.d0/Dd(j)*a(1);
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*a(2);
  !%
  j=2;
  Dd(j)=a(0) - (La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*(a(1) - (Lb(j+1)*Dd(j-1)*La(j)));
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*a(2);
  !%
  do j=3,N-4
     Dd(j) = a(0) - (Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
     La(j+1)= 1.d0/Dd(j)*(a(1) - Lb(j+1)*Dd(j-1)*La(j) );
     !%i=j+2
     Lb(j+2)= 1.d0/Dd(j)*a(2);
  end do
  !
  j=N-3;
  Dd(j) = a(0) - ( Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*(a(1) - Lb(j+1)*Dd(j-1)*La(j) );
  j=N-2;
  Dd(j) = a(0) - ( Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  !
  Lf(1)=1.d0/Dd(1)*a(2)*shp;
  i=2;
  Lf(i)=1.d0/Dd(i)*(-Lf(i-1)*Dd(i-1)*La(i));
  do i=3,N-4
     Lf(i)=1.d0/Dd(i)*(-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  end do
  i=N-3;
    Lf(i)=1.d0/Dd(i)*(a(2)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  i=N-2;
    Lf(i)=1.d0/Dd(i)*(a(1)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  !%%%
  j=N-1;
  tmp=0.d0;
  do k=1,(N-2)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Lf(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !
  Lg(1)=1.d0/Dd(1)*a(1)*shp;
  i=2;
    Lg(i)=1.d0/Dd(i)*(a(2)*shp-Lg(i-1)*Dd(i-1)*La(i));
  do i=3,(N-3)
     Lg(i)=1.d0/Dd(i)*(-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  end do
  i=N-2;
    Lg(i)=1.d0/Dd(i)*(a(2)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%
  i=N-1;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lf(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(1) - tmp); 
  !%%%
  j=N;
  tmp=0.d0;
  do k=1,N-1
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lg(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  
  cLf=dconjg(Lf);
  cLg=dconjg(Lg);
  ! check 
  !det = sum(Dd)
  !write(*,*) 'modChol5cyclic_cmplx: det = ', det 

end subroutine modChol5cyclic_cmplx

!***********************************************************************
subroutine modChol7cyclic_cmplx(La,Lb,Lc,Dd,Le,Lf,Lg,cLe,cLf,cLg,a,N,shp)
  integer N
  real*8 La(N),Lb(N),Lc(N),Dd(N)
  !complex*16 La(N),Lb(N),Lc(N),Dd(N)
  complex*16 Le(N),Lf(N),Lg(N),cLe(N),cLf(N),cLg(N)
  real*8 a(-3:3)
  complex*16 shp,tmp ! note that shp is complex*16
  ! do not include real*8 value in the last arguments, (,... ,1.d0) 
  integer i,j,k
  real*8 det

  j=1;
  Dd(j)=a(0);
  La(j+1)= 1.d0/Dd(j)*a(1);
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*a(2);
  !%i=j+3
  Lc(j+3)= 1.d0/Dd(j)*a(3);
  !%
  j=2;
  Dd(j)=a(0) - (La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*(a(1) - (Lb(j+1)*Dd(j-1)*La(j)));
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*(a(2) - Lc(j+2)*Dd(j-1)*La(j) );
  !%i=j+3
  Lc(j+3)= 1.d0/Dd(j)*a(3);
  !%
  j=3;
  Dd(j)=a(0) - (Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*( a(1) - (Lc(j+1)*Dd(j-2)*Lb(j) + Lb(j+1)*Dd(j-1)*La(j)));
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*( a(2) - Lc(j+2)*Dd(j-1)*La(j) );
  !%i=j+3
  Lc(j+3)= 1.d0/Dd(j)*a(3);
  !%
  do j=4,N-6
     Dd(j) = a(0) - (Lc(j)*Dd(j-3)*Lc(j)+ &
                     Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
     La(j+1)= 1.d0/Dd(j)*(a(1) - &
                    (Lc(j+1)*Dd(j-2)*Lb(j) + Lb(j+1)*Dd(j-1)*La(j)));
     !%i=j+2
     Lb(j+2)= 1.d0/Dd(j)*(a(2) - Lc(j+2)*Dd(j-1)*La(j) );
     !%i=j+3
     Lc(j+3)= 1.d0/Dd(j)*a(3);
  end do
  j=N-5;
  Dd(j) = a(0) - (Lc(j)*Dd(j-3)*Lc(j)+ &
                  Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*(a(1) - &
               (Lc(j+1)*Dd(j-2)*Lb(j) + Lb(j+1)*Dd(j-1)*La(j)));
  !%i=j+2
  Lb(j+2)= 1.d0/Dd(j)*(a(2) - Lc(j+2)*Dd(j-1)*La(j) );
  j=N-4;
  Dd(j) = a(0) - (Lc(j)*Dd(j-3)*Lc(j)+ &
                  Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  La(j+1)= 1.d0/Dd(j)*(a(1) - &
                 (Lc(j+1)*Dd(j-2)*Lb(j) + Lb(j+1)*Dd(j-1)*La(j)));
  j=N-3;
  Dd(j) = a(0) - (Lc(j)*Dd(j-3)*Lc(j)+ &
                  Lb(j)*Dd(j-2)*Lb(j)+La(j)*Dd(j-1)*La(j));
  !
  Le(1)=dcmplx(1.d0,0.d0)/Dd(1)*a(3)*shp; !%^\ast
  i=2;
  Le(2)=1.d0/Dd(i)*(-Le(i-1)*Dd(i-1)*La(i));
  i=3;
  Le(3)=1.d0/Dd(i)*(-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  do i=4,N-6
     Le(i)=1.d0/Dd(i)*(-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  end do
  i=N-5;
  Le(i)=1.d0/Dd(i)*(a(3)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  i=N-4;
  Le(i)=1.d0/Dd(i)*(a(2)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  i=N-3;
  Le(i)=1.d0/Dd(i)*(a(1)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  !%%%
  j=N-2;
  tmp=0.d0;
  do k=1,(N-3)
     tmp=tmp+ Le(k)*Dd(k)*dconjg(Le(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%%%
  Lf(1)=1.d0/Dd(1)*a(2)*shp;
  i=2;
  Lf(i)=1.d0/Dd(i)*(a(3)*shp-Lf(i-1)*Dd(i-1)*La(i));
  i=3;
  Lf(i)=1.d0/Dd(i)*(-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  do i=4,N-5
     Lf(i)=1.d0/Dd(i)*(-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  end do
  i=N-4;
  Lf(i)=1.d0/Dd(i)*(a(3)-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  i=N-3;
  Lf(i)=1.d0/Dd(i)*(a(2)-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  i=N-2;
  tmp=0.d0;
  do k=1,(N-3)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Le(k));
  end do
  Lf(i)= 1.d0/Dd(j)*(a(1) - tmp);
  !%%%
  j=N-1;
  tmp=0.d0;
  do k=1,(N-2)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Lf(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !
  Lg(1)=1.d0/Dd(1)*a(1)*shp;
  i=2;
  Lg(i)=1.d0/Dd(i)*(a(2)*shp-Lg(i-1)*Dd(i-1)*La(i));
  i=3;
  Lg(i)=1.d0/Dd(i)*(a(3)*shp-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  do i=4,(N-4)
     Lg(i)=1.d0/Dd(i)*(-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  end do
  i=N-3;
  Lg(i)=1.d0/Dd(i)*(a(3)-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%
  i=N-2;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Le(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(2) - tmp); 
  !%
  i=N-1;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lf(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(1) - tmp); 
  !%%%
  j=N;
  tmp=0.d0;
  do k=1,N-1
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lg(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  
  cLe=dconjg(Le);
  cLf=dconjg(Lf);
  cLg=dconjg(Lg);
  ! check 
  !det = sum(Dd)
  !write(*,*) 'modChol7cyclic_cmplx: det = ', det 

end subroutine modChol7cyclic_cmplx

!***********************************************************************
!    subroutine for addshear to modChol decomposition
!***********************************************************************
subroutine modChol5_addshear(La,Lb,Dd,Lf,Lg,cLf,cLg,a,N,shp)
  integer N
  real*8  La(N),Lb(N),Dd(N)
  !complex*16 La(N),Lb(N),Dd(N)
  complex*16 Lf(N),Lg(N),cLf(N),cLg(N)
  real*8 a(-2:2)
  complex*16 shp,tmp
  integer i,j,k
  real*8 det

  !%%%
  !%Lf(1)=1/Dd(1)*a(2)*shp;
  !%i=2;
  !%Lf(i)=1/Dd(i)*(a(3)*shp-Lf(i-1)*Dd(i-1)*La(i));
  !%i=3;
  !%Lf(i)=1/Dd(i)*(-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  !%for i=4:N-5
  !%  Lf(i)=1/Dd(i)*(-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*L%a(i));
  !%end
  Lf(1:N-4)=Lf(1:N-4)*shp;
  !i=N-4;
  !Lf(i)=1.d0/Dd(i)*(a(3)-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i) &
  !                      -Lf(i-1)*Dd(i-1)*La(i));
  i=N-3;
  Lf(i)=1.d0/Dd(i)*(a(2)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  i=N-2;
  Lf(i)=1.d0/Dd(i)*(a(1)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  !  tmp=0.d0;
  !  do k=1,(N-3)
  !     tmp=tmp+Lf(k)*Dd(k)*dconjg(Le(k));
  !  end do
  !  Lf(i)= 1.d0/Dd(j)*(a(1) - tmp);
  !  !%%%
  j=N-1;
  tmp=0.d0;
  do k=1,(N-2)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Lf(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%
  !%Lg(1)=1/Dd(1)*a(1)*shp;
  !%i=2;
  !%  Lg(i)=1/Dd(i)*(a(2)*shp-Lg(i-1)*Dd(i-1)*La(i));
  !%i=3;
  !%  Lg(i)=1/Dd(i)*(a(3)*shp-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%for i=4:(N-4)
  !%  Lg(i)=1/Dd(i)*(-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%end
  Lg(1:N-3)=Lg(1:N-3)*shp;
  !i=N-3;
  !Lg(i)=1.d0/Dd(i)*(a(3)-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i) &
  !                      -Lg(i-1)*Dd(i-1)*La(i));
  !%
  i=N-2;
  Lg(i)=1.d0/Dd(i)*(a(2)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
!  tmp=0.d0;
!  do k=1,(i-1)
!     tmp=tmp+Lg(k)*Dd(k)*dconjg(Le(k));
!  end do
!  Lg(i)= 1.d0/Dd(i)*(a(2) - tmp); !%^\ast
  !%
  i=N-1;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lf(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(1) - tmp); !%^\ast
  !%
  j=N;
  tmp=0.d0;
  do k=1,N-1
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lg(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%
  !cLe=dconjg(Le);
  cLf=dconjg(Lf);
  cLg=dconjg(Lg);
  !%%% check 
  !det = sum(Dd)
  !write(*,*) 'modChol_addshear: det = ', det 

end subroutine modChol5_addshear
!
!***********************************************************************
subroutine modChol7_addshear(La,Lb,Lc,Dd,Le,Lf,Lg,cLe,cLf,cLg,a,N,shp)
  integer N
  real*8  La(N),Lb(N),Lc(N),Dd(N)
  !complex*16 La(N),Lb(N),Lc(N),Dd(N)
  complex*16 Le(N),Lf(N),Lg(N),cLe(N),cLf(N),cLg(N)
  real*8 a(-3:3)
  complex*16 shp,tmp
  integer i,j,k
  real*8 det

  !%%%
  !%
  !%  Le(1)=1/Dd(1)*a(3)*shp; %^\ast
  !%  i=2;
  !%  Le(2)=1/Dd(i)*(-Le(i-1)*Dd(i-1)*La(i));
  !%  i=3;
  !%  Le(3)=1/Dd(i)*(-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  !%  for i=4:N-6
  !%    Le(i)=1/Dd(i)*(-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  !%  end
  Le(1:N-6)=Le(1:N-6)*shp;
  i=N-5;
  Le(i)=1.d0/Dd(i)*(a(3)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  i=N-4;
  Le(i)=1.d0/Dd(i)*(a(2)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  i=N-3;
  Le(i)=1.d0/Dd(i)*(a(1)-Le(i-3)*Dd(i-3)*Lc(i)-Le(i-2)*Dd(i-2)*Lb(i)-Le(i-1)*Dd(i-1)*La(i));
  !%%%
  j=N-2;
  tmp=0.d0;
  do k=1,(N-3)
     tmp=tmp+Le(k)*Dd(k)*dconjg(Le(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%%%
  !%Lf(1)=1/Dd(1)*a(2)*shp;
  !%i=2;
  !%Lf(i)=1/Dd(i)*(a(3)*shp-Lf(i-1)*Dd(i-1)*La(i));
  !%i=3;
  !%Lf(i)=1/Dd(i)*(-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*La(i));
  !%for i=4:N-5
  !%  Lf(i)=1/Dd(i)*(-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i)-Lf(i-1)*Dd(i-1)*L%a(i));
  !%end
  Lf(1:N-5)=Lf(1:N-5)*shp;
  i=N-4;
  Lf(i)=1.d0/Dd(i)*(a(3)-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i) &
                        -Lf(i-1)*Dd(i-1)*La(i));
  i=N-3;
  Lf(i)=1.d0/Dd(i)*(a(2)-Lf(i-3)*Dd(i-3)*Lc(i)-Lf(i-2)*Dd(i-2)*Lb(i) &
                        -Lf(i-1)*Dd(i-1)*La(i));
  i=N-2;
  tmp=0.d0;
  do k=1,(N-3)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Le(k));
  end do
  Lf(i)= 1.d0/Dd(j)*(a(1) - tmp);
  !%%%
  j=N-1;
  tmp=0.d0;
  do k=1,(N-2)
     tmp=tmp+Lf(k)*Dd(k)*dconjg(Lf(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%
  !%Lg(1)=1/Dd(1)*a(1)*shp;
  !%i=2;
  !%  Lg(i)=1/Dd(i)*(a(2)*shp-Lg(i-1)*Dd(i-1)*La(i));
  !%i=3;
  !%  Lg(i)=1/Dd(i)*(a(3)*shp-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%for i=4:(N-4)
  !%  Lg(i)=1/Dd(i)*(-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i)-Lg(i-1)*Dd(i-1)*La(i));
  !%end
  Lg(1:N-4)=Lg(1:N-4)*shp;
  i=N-3;
  Lg(i)=1.d0/Dd(i)*(a(3)-Lg(i-3)*Dd(i-3)*Lc(i)-Lg(i-2)*Dd(i-2)*Lb(i) &
                        -Lg(i-1)*Dd(i-1)*La(i));
  !%
  i=N-2;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Le(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(2) - tmp); !%^\ast
  !%
  i=N-1;
  tmp=0.d0;
  do k=1,(i-1)
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lf(k));
  end do
  Lg(i)= 1.d0/Dd(i)*(a(1) - tmp); !%^\ast
  !%
  j=N;
  tmp=0.d0;
  do k=1,N-1
     tmp=tmp+Lg(k)*Dd(k)*dconjg(Lg(k));
  end do
  Dd(j)=a(0) - dreal(tmp);
  !%
  cLe=dconjg(Le);
  cLf=dconjg(Lf);
  cLg=dconjg(Lg);
  !%%% check 
  !det = sum(Dd)
  !write(*,*) 'modChol_addshear: det = ', det 

end subroutine modChol7_addshear


!***********************************************************************
!    subroutine for forward- and backward-substitution 
!    solve A x = bb, by modified Cholesky decomposition A=(LDL*)
!
! La, ..., is the vectors which stores the elements of L as follows
! for i=4:nx-3
!   Lc(i)=L(i,i-3); 
! end
! for i=3:nx-3
!   Lb(i)=L(i,i-2);
! end
! for i=2:nx-3
!   La(i)=L(i,i-1);
! end
! for i=1:nx
!   Dd(i)=D(i,i);
! end
! for i=1:nx-3
!   Le(i)=L(N-2,i);
! end
! for i=1:nx-2
!   Lf(i)=L(N-1,i);
! end
! for i=1:nx-1
!   Lg(i)=L(N,i);
! end               
!***********************************************************************
subroutine solve5cyclic(La,Lb,Dd,Lf,Lg,cLf,cLg,bb,x,y,nx)
  implicit none

  integer nx
  real*8 La(nx),Lb(nx),Dd(nx)
  !complex*16 La(nx),Lb(nx),Dd(nx)
  complex*16 Lf(nx),Lg(nx),cLf(nx),cLg(nx)
  complex*16 tmp
  complex*16 bb(nx),x(nx),y(nx)
  integer i,j,k

  ! forwardsubstitution
  y(1)=bb(1);
  y(2)=bb(2)-La(2)*y(1);
  do i=3,nx-2
     y(i)=bb(i)-(Lb(i)*y(i-2)+La(i)*y(i-1));
  end do

  i=nx-1;
  tmp=dcmplx(0.d0,0.d0);
  do k=1,i-1
     tmp=tmp+Lf(k)*y(k);
  end do
  y(i)=bb(i)-tmp;
  
  i=nx;
  tmp=dcmplx(0.d0,0.d0);
  do k=1,i-1
     tmp=tmp+Lg(k)*y(k);
  end do
  y(i)=bb(i)-tmp;
  !
  ! backsubstitution
  x(nx)  =y(nx)/Dd(nx);
  x(nx-1)=y(nx-1)/Dd(nx-1) - cLg(nx-1)*x(nx);
  x(nx-2)=y(nx-2)/Dd(nx-2) - cLf(nx-2)*x(nx-1) - cLg(nx-2)*x(nx);
  x(nx-3)=y(nx-3)/Dd(nx-3) - La(nx-2)*x(nx-2) & 
                           - cLf(nx-3)*x(nx-1) - cLg(nx-3)*x(nx);
  do i=nx-4,1,-1
     x(i)=y(i)/Dd(i) - La(i+1)*x(i+1) - Lb(i+2)*x(i+2) &
                     - cLf(i)*x(nx-1) - cLg(i)*x(nx);
  end do

end subroutine solve5cyclic

!***********************************************************************

subroutine solve7cyclic(La,Lb,Lc,Dd,Le,Lf,Lg,cLe,cLf,cLg,bb,x,y,nx)
  implicit none

  integer nx
  real*8 La(nx),Lb(nx),Lc(nx),Dd(nx)
  !complex*16 La(nx),Lb(nx),Lc(nx),Dd(nx)
  complex*16 Le(nx),Lf(nx),Lg(nx),cLe(nx),cLf(nx),cLg(nx)
  complex*16 tmp
  complex*16 bb(nx),x(nx),y(nx)
  integer i,j,k

  ! forwardsubstitution
  y(1)=bb(1);
  y(2)=bb(2)-La(2)*y(1);
  y(3)=bb(3)-(Lb(3)*y(1)+La(3)*y(2));
  do i=4,nx-3
     y(i)=bb(i)-(Lc(i)*y(i-3)+Lb(i)*y(i-2)+La(i)*y(i-1));
  end do
  !extra three lines
  i=nx-2;
  tmp=dcmplx(0.d0,0.d0);
  do k=1,i-1
     tmp=tmp+Le(k)*y(k);
  end do
  y(i)=bb(i)-tmp;
  
  i=nx-1;
  tmp=dcmplx(0.d0,0.d0);
  do k=1,i-1
     tmp=tmp+Lf(k)*y(k);
  end do
  y(i)=bb(i)-tmp;
  
  i=nx;
  tmp=dcmplx(0.d0,0.d0);
  do k=1,i-1
     tmp=tmp+Lg(k)*y(k);
  end do
  y(i)=bb(i)-tmp;
  !
  ! backsubstitution
  x(nx)  =y(nx)/Dd(nx);
  x(nx-1)=y(nx-1)/Dd(nx-1) - cLg(nx-1)*x(nx);
  x(nx-2)=y(nx-2)/Dd(nx-2) - cLf(nx-2)*x(nx-1) - cLg(nx-2)*x(nx);
  x(nx-3)=y(nx-3)/Dd(nx-3) - cLe(nx-3)*x(nx-2) - cLf(nx-3)*x(nx-1) - cLg(nx-3)*x(nx);
  x(nx-4)=y(nx-4)/Dd(nx-4) - La(nx-3)*x(nx-3) &
       - cLe(nx-4)*x(nx-2) - cLf(nx-4)*x(nx-1) - cLg(nx-4)*x(nx);
  x(nx-5)=y(nx-5)/Dd(nx-5) - La(nx-4)*x(nx-4) - Lb(nx-3)*x(nx-3) &
       - cLe(nx-5)*x(nx-2) - cLf(nx-5)*x(nx-1) - cLg(nx-5)*x(nx);
  do i=nx-6,1,-1
     x(i)=y(i)/Dd(i) - La(i+1)*x(i+1) - Lb(i+2)*x(i+2) - Lc(i+3)*x(i+3) &
          - cLe(i)*x(nx-2) - cLf(i)*x(nx-1) - cLg(i)*x(nx);
  end do

end subroutine solve7cyclic

!***********************************************************************
!    SUBRUTINAS PARA TRABAJAR CON n-DIAGONALES
!    (Numerical Recipes in FORTRAN)
!      bandec, banbks
!      Solo para liso---
!***********************************************************************

SUBROUTINE banbks7(a,n,b)
  INTEGER n
  REAL*8 a(7,n),b(n)
  INTEGER i,k

  do k=1,n-3
     b(k+1) = b(k+1)-a(5,k)*b(k)
     b(k+2) = b(k+2)-a(6,k)*b(k)
     b(k+3) = b(k+3)-a(7,k)*b(k)
  enddo

  !     n-2

  b(n-1) = b(n-1)-a(5,n-2)*b(n-2)
  b(n)   = b(n)  -a(6,n-2)*b(n-2)

  !     n-1

  b(n) = b(n)    -a(5,n-1)*b(n-1)

  !     back substitution

  b(n) = b(n)*a(1,n)
  b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
  b(n-2) = (b(n-2)-a(2,n-2)*b(n-1)-a(3,n-2)*b(n))*a(1,n-2)

  do i=n-3,1,-1
     b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i)-a(4,i)*b(3+i))*a(1,i)
  enddo

END SUBROUTINE banbks7


!!!!!!!!!! -----------------------------------
subroutine bandec7(a,n)
  INTEGER n
  REAL*8 a(7,n)
  INTEGER j,k
  !
  do j=1,4
     a(j,1)=a(j+3,1)
  enddo

  do j=1,5
     a(j,2)=a(j+2,2)
  enddo

  do j=1,6
     a(j,3)=a(j+1,3)
  enddo


  !     LU

  do k=1,n-3
     a(1,k)   = 1d0/a(1,k)

     a(5,k)   = a(1,k+1)*a(1,k)
     a(1,k+1) = a(2,k+1)-a(5,k)*a(2,k)
     a(2,k+1) = a(3,k+1)-a(5,k)*a(3,k)
     a(3,k+1) = a(4,k+1)-a(5,k)*a(4,k)
     a(4,k+1) = a(5,k+1)

     a(6,k)   = a(1,k+2)*a(1,k)
     a(1,k+2) = a(2,k+2)-a(6,k)*a(2,k)
     a(2,k+2) = a(3,k+2)-a(6,k)*a(3,k)
     a(3,k+2) = a(4,k+2)-a(6,k)*a(4,k)
     a(4,k+2) = a(5,k+2)
     a(5,k+2) = a(6,k+2)

     a(7,k)   = a(1,k+3)*a(1,k)
     a(1,k+3) = a(2,k+3)-a(7,k)*a(2,k)
     a(2,k+3) = a(3,k+3)-a(7,k)*a(3,k)
     a(3,k+3) = a(4,k+3)-a(7,k)*a(4,k)
     a(4,k+3) = a(5,k+3)
     a(5,k+3) = a(6,k+3)
     a(6,k+3) = a(7,k+3)

  enddo

  !     k=n-2
  a(1,n-2) = 1d0/a(1,n-2)

  a(5,n-2) = a(1,n-1)*a(1,n-2)

  a(1,n-1) = a(2,n-1)-a(5,n-2)*a(2,n-2)
  a(2,n-1) = a(3,n-1)-a(5,n-2)*a(3,n-2)
  a(3,n-1) = a(4,n-1)-a(5,n-2)*a(4,n-2)
  a(4,n-1) = a(5,n-1)

  a(6,n-2) = a(1,n)*a(1,n-2)

  a(1,n)   = a(2,n)-a(6,n-2)*a(2,n-2)
  a(2,n)   = a(3,n)-a(6,n-2)*a(3,n-2)
  a(3,n)   = a(4,n)-a(6,n-2)*a(4,n-2)
  a(4,n)   = a(5,n)
  a(5,n)   = a(6,n)

  !     k=n-1

  a(1,n-1) = 1d0/a(1,n-1)
  a(5,n-1) = a(1,n)*a(1,n-1)

  a(1,n)   = a(2,n)-a(5,n-1)*a(2,n-1)
  a(2,n)   = a(3,n)-a(5,n-1)*a(3,n-1)
  a(3,n)   = a(4,n)-a(5,n-1)*a(4,n-1)
  a(4,n)   = a(5,n)

  !     the next loop will be used in banbk

  a(1,n)=1d0/a(1,n)

END SUBROUTINE bandec7


!**************************************************
!*    SUBRUTINA PARA RESOLVER SIST LINEALES
!*    (Numerical Recipes in FORTRAN)
!**************************************************
SUBROUTINE gaussj(a,n,np,b,m,mp)
  INTEGER m,mp,n,np,NMAX
  REAL*8 a(np,np),b(np,mp)
  PARAMETER (NMAX=50)
  INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  REAL*8 big,dum,pivinv
  do j=1,n
     ipiv(j)=0
  enddo
  do i=1,n
     big=0d0
     do j=1,n
        if(ipiv(j).ne.1)then
           do  k=1,n
              if (ipiv(k).eq.0) then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 pause 'singular matrix in gaussj'
              endif
           enddo
        endif
     enddo

     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do  l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do  l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
     pivinv=1d0/a(icol,icol)
     a(icol,icol)=1d0
     do  l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo
     do  l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo
     do  ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0d0
           do  l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo
           do  l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo
        endif
     enddo
  enddo

  do  l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do  k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  enddo

END subroutine gaussj


















