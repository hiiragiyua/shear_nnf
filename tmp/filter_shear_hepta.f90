module yfilter
  ! 
  real*8, allocatable:: pr1(:), rh1(:), rh1h(:) ! (constant) B u' = A u

  real*8, allocatable:: a1(:),b1(:),c1(:),d1(:),d1m1(:)
  complex*16, allocatable:: e1(:),f1(:),g1(:) 
  complex*16, allocatable:: ce1(:),cf1(:),cg1(:)

  real*8, allocatable:: ah(:),bh(:),ch(:),dh(:),dhm1(:)
  complex*16, allocatable:: eh(:),fh(:),gh(:) 
  complex*16, allocatable:: ceh(:),cfh(:),cgh(:)

  complex*16, allocatable:: wk1(:), wk2(:), ur(:)
  complex*16, allocatable:: uc(:), duc(:), pc(:)
  real*8  h
  integer my

end module yfilter


subroutine prefilter
  use yfilter
  implicit none

  integer N, M

  integer j,k
  real*8  pi
  complex*16 zrr,zii

  ! ----------- define overall parameters ----
  N=4    ! semi-band width for r.h.s (na)
  M=3    ! semi-band width for l.h.s (nb)
  ! number of spectral points ! not implemented.
  pi=4d0*atan(1d0)
  !write(*,*) '  prederiv1, semi-band width: rhs 4, lhs 2 (penta-CFD)'
  allocate ( pr1(-M:M), rh1(-N:N), rh1h(-N:N)  ) 
  allocate ( a1(my),b1(my),c1(my),d1(my),d1m1(my))
  allocate ( e1(my),f1(my),g1(my),ce1(my),cf1(my),cg1(my) )
  allocate ( ah(my),bh(my),ch(my),dh(my),dhm1(my) )
  allocate ( eh(my),fh(my),gh(my),ceh(my),cfh(my),cgh(my) )
  !write(*,*) '  prederiv1: allocated'
  pr1 = 0d0; rh1 = 0d0; rh1h = 0d0
  a1  = 0d0; b1  = 0d0; c1 = 0d0; d1 = 0d0; d1m1 = 0.d0;  
  e1  = dcmplx(0d0,0d0); f1  = dcmplx(0d0,0d0); g1  = dcmplx(0d0,0d0)
  ce1 = dcmplx(0d0,0d0); cf1 = dcmplx(0d0,0d0); cg1 = dcmplx(0d0,0d0)

  ah  = 0d0; bh  = 0d0; ch = 0d0; dh = 0d0; dhm1 = 0d0; 
  eh  = dcmplx(0d0); fh  = dcmplx(0d0); gh = dcmplx(0d0)
  ceh = dcmplx(0d0); cfh = dcmplx(0d0); cgh = dcmplx(0d0)
  ! deallocate??
  ! set B  ( hepta-diagonal, B u' = A u )
  pr1(0) = 1.d0;
  pr1(1) = 0.738981897704866d0; 
  pr1(2) = 0.315003368623046d0;
  pr1(3) = 0.053003478448814d0;
  pr1(-3) = pr1(3);
  pr1(-2) = pr1(2);
  pr1(-1) = pr1(1);
  ! set A  ( B u' = A u )
  rh1(0) = 0.987412035368317d0;
  rh1(1) = 0.749052269410218d0; 
  rh1(2) = 0.309968182770373d0; 
  rh1(3) = 0.054442102978149d0; 
  rh1(4) =-0.000179828066167d0;
  rh1(-4) = rh1(4); ! symmetric for filtering
  rh1(-3) = rh1(3);
  rh1(-2) = rh1(2);
  rh1(-1) = rh1(1);

  !write(*,*) '  check pr1=',pr1 
  !write(*,*) '  check rh1=',rh1 
  !
  !  ------ scale and prefactor, do not scale for filtering  
  !rh1 = rh1h 
  !rh1 = rh1/h
  zrr = dcmplx(1.d0,0.d0)
  !call modChol5cyclic_cmplx(a1,b1,d1,f1,g1,cf1,cg1,pr1(-2:2),my,zrr)
  call modChol7cyclic_cmplx(a1,b1,c1,d1,d1m1,e1,f1,g1,ce1,cf1,cg1,pr1,my,zrr)

  ! shear will be included by calling modChol_addshear() in each time step. 
  !write(*,*) 'a1',a1 
  !write(*,*) 'b1',b1 
  !write(*,*) 'c1',c1 
  !write(*,*) 'd1',d1 
  !write(*,*) 'e1',e1 
  !write(*,*) 'f1',f1 
  !write(*,*) 'g1',g1 

!  ==== THIS IS MATLAB SCRIPT to get coefficient === (see ./Doc/matlab/prederiv.m (icase==8) and transfer_function.m)
!  N=3;M=4; %nb na
!  kap(1)=1.0;    G(1)=0.0;
!  kap(2)=6/9; G(2)=0.72; 
!  kap(3)=8.5/9; G(3)=0.001; 
!  kap(4)=5.5/9; G(4)=0.91;
!  kap=kap*pi;
!  mat=zeros(N+M+1,N+M+1);
!  vec=zeros(N+M+1,1);
!  jj=(1:N); jjm=(1:M);
!  %
!  mat(1,1:N) = -1.0; % order 0
!  mat(1,N+2:N+M+1) = 1.0; % order 0
!  mat(1,N+1) =0.5;
!  for k=2:N+1 % order 2*k-2 % k=2,3,4
!    kk=2*k-2; % 2,4,6
!    mat(k,1:N) =  jj.^kk;
!    mat(k,N+2:N+M+1) = -jjm.^kk;
!  end
!  vec(1)=0.5;
!  for k=1:M
!    mat(N+1+k,1:N)= -G(k)*cos(kap(k)*jj);
!    mat(N+1+k,N+1) = 0.5; % a0
!    mat(N+1+k,N+2:N+M+1) = cos(kap(k)*jjm);
!    vec(N+1+k) = 0.5*G(k);
!  end
!  %
!  % --- check --
!  % mat
!  % vec
!  % ---
!  tmp=inv(mat)*vec;
!  b=tmp(1:N);a=tmp(N+2:N+M+1);
!  b0=1.0;a0=tmp(N+1);        
!  
end subroutine prefilter


subroutine filterc(u,du,shp,shm,m,ias)
  use yfilter

  implicit none
  integer i,j,m,na,ias
  real*8 u(m,my),du(m,my)
  complex*16 shp,shm,const

  ! wrapping to complex*16
  if (m.eq.1) then
     uc(:) = dcmplx(u(1,:),0.d0)
  else
     uc(:) = dcmplx(u(1,:),u(2,:))
  endif

  na=4
  call matvecNcyclic(rh1(-na:na),uc,ur,wk1,shp,shm,my,m,na)
  
  if (m.eq.1) then
     call solve7cyclic(a1,b1,c1,d1,d1m1,e1,f1,g1,ce1,cf1,cg1,wk1,duc,wk2,my)
     du(1,:) = dreal(duc)
  else
     !
     if (ias.eq.1) then
        dh=d1; dhm1=d1m1;

        eh(1:(my-6))=dreal(e1(1:(my-6)))*shp; eh(my-5:my)=e1(my-5:my) 
        fh(1:(my-5))=dreal(f1(1:(my-5)))*shp; fh(my-4:my)=f1(my-4:my) 
        gh(1:(my-4))=dreal(g1(1:(my-4)))*shp; gh(my-3:my)=g1(my-3:my) 
        call modChol_addshear(a1,b1,c1,dh,dhm1,eh,fh,gh,ceh,cfh,cgh,pr1(-3:3),my,shp)
     end if

     call solve7cyclic(a1,b1,c1,dh,dhm1,eh,fh,gh,ceh,cfh,cgh,wk1,du,wk2,my)
     
  endif
    
end subroutine filterc


subroutine filteras(n)  
  use yfilter
  implicit none
  integer n,my8
  real*8  dy,Length

  my = n
  h  = dy
  ! ---  start doing something ---
  my8 = my+4+4
  allocate (wk1(my), wk2(my), ur(my8))
  allocate (uc(my),duc(my),pc(my))
  wk1 = 0.d0; wk2 = 0.d0; ur = 0.d0
  uc = 0.d0; duc = 0.d0; pc = 0.d0

  call prefilter

end subroutine filteras

