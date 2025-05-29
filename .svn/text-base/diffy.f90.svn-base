!=======================================================================
!       PAQUETE DE SUBRUTINAS PARA DIFERENCIAS FINITAS COMPACTAS
!       O.F, nov 2002
!       S.H. may 2005
!       J.J. july 2007,  uniform grid for disk
!
!=======================================================================
module yderiv
  real*8, allocatable:: pr1(:,:), rh1(:,:), endder(:)
  real*8, allocatable:: pr2(:,:), pr2i(:,:), rh2(:,:),lapmat(:,:)
  real*8, allocatable:: wk1(:), wk2(:)
  real*8  h      
  integer N1,N2, N3, my

end module yderiv

!**********************************************************************
! prederiv1:
!
! Computes the coefficients for the compact finite differences
!     for the first derivatives of bandwidth 2*N+1
!     including the schemes for the end points. 
! References:  Disk notebook ctr 2007, pp. 9-14
!     and maple \javi\wall\sergio\pade7.ms   
!     and       \javi\wall\sergio\pade7end_sp.ms   
!
! Scheme is: sum_{j=-N}^N a_ij f'_{i+j} = 
!            sum_{j=-N}^N b_ij f_{i+j}
!
! At the ends, the inner band-width is maintained,
!    but outer limit is j=-i (for j=0:N-1)
! 
! Although this is a general routine it is optimized for N=3
!    It uses in general M spectral points to dispersion relation
!       at khat = pi*kap(1:M).  Here M=3.
!
!*********************************************************************
subroutine prederiv1
  use yderiv
  implicit none

  integer N, M
  integer,allocatable:: jj(:),nsp(:)
  real*8, allocatable:: mat(:,:),vec(:),kap(:)

  integer j,k,kk,kp
  real*8  pi

  ! ----------- define overall parameters ----
  N=3    ! semi-band width
  M=3    ! number of spectral points
  allocate (kap(M),nsp(N))
  kap(1) = 0.5d0
  kap(2) = 0.9d0
  kap(3) = 5d-1*(kap(1)+kap(2))
  pi=4d0*atan(1d0)
  kap = kap*pi

  allocate ( pr1(-N:N,my), rh1(-N:N,my) )
  pr1 = 0d0
  rh1 = 0d0
  N1  = N
  !-----------------------------------------------------------------------
  !     central part
  ! ----------------------------------------------------------------------
  allocate (mat(2*N, 2*N), vec(2*N), jj(N))
  vec = 0d0
  mat = 0d0
  do j=1,N
     jj(j) = j
  enddo

  do k=1,2*N-M            ! order 2*k-1 
     kk = 2*k-1
     mat(k,1:N)     = kk*jj**(kk-1)    
     mat(k,N+1:2*N) = -jj**kk      
  enddo
  vec(1)            = -0.5d0 ! the first-order derivative at j=0

  do k=1,M                   ! the spectral points
     mat(2*N-M+k,1:N)     =  kap(k)*cos(kap(k)*jj)    
     mat(2*N-M+k,N+1:2*N) = -sin(kap(k)*jj)
     vec(2*N-M+k)         = -0.5d0*kap(k)  
  enddo
     
  call gaussj(mat,2*N,2*N,vec,1,1)

  do j=1,N
     rh1( j,N+1:my-N) =  vec(N+j)
     rh1(-j,N+1:my-N) = -vec(N+j)
     pr1( j,N+1:my-N) =  vec(j)
     pr1(-j,N+1:my-N) =  vec(j)
  enddo
  pr1(0,:) = 1d0 

  deallocate (mat,vec,jj)

  !-----------------------------------------------------------------------
  !
  !       end points
  !       mixed order and spectral, look in notebook
  !       only here for M=3 spectral points
  !-----------------------------------------------------------------------

  if (M.ne.3 .or. N.ne.3) then
     write(*,*) 'the endpoints in prederiv1 only optimized for M=N=3; M=',M,'  N=',N
  endif

  nsp(1)    = 2   ! how many spectral relations to use at endpoint kp=1:N
  nsp(2)    = 2   
  nsp(3)    = 3 

  allocate (mat(0:4*N+1,-N:3*N+1), vec(0:4*N+1) , jj(-N:N))
  do j=-N,N
     jj(j) = j
  enddo

  !    ------------  kp is the point
  do kp=1,N 

     mat=0d0
     vec=0d0
     
     do k=0,N-kp                !  the coefficients that are zero      
        mat(2*k,k-N)     = 1d0
        mat(2*k+1,N+k+1) = 1d0
     enddo
     kk=2*(N-kp+1)  
     mat(kk,0) = 1d0            !  a(0)=1
     vec(kk)   = 1d0
     do k=1,nsp(kp)             !  -- the spectral relations
        mat(kk+1,-N:N)       =  kap(k)*cos(kap(k)*jj)
        mat(kk+1,N+1:3*N+1)  = -sin(kap(k)*jj)
        mat(kk+2,-N:N)       =  kap(k)*sin(kap(k)*jj)
        mat(kk+2,N+1:3*N+1)  =  cos(kap(k)*jj)
        kk=kk+2
     enddo

     kk=kk+1
     mat(kk,N+1:3*N+1)  = 1d0        !  -- the zeroth order   
     do k=1,4*N
        kk=kk+1
        if (kk>4*N+1) exit
        mat(kk,-N:N)      =  k*jj**(k-1)
        if (k.eq.1)  mat(kk,0)=1d0   !  -- fix 0**0 in first order
        mat(kk,N+1:3*N+1) = -jj**k 
     enddo

     call gaussj(mat,4*N+2,4*N+2,vec,1,1)

     pr1(:,kp) =  vec(0:2*N)
     rh1(:,kp) =  vec(2*N+1:4*N+1)
  enddo

  ! -----------------  the other end
  do j=0,N-1
     pr1(:,my-j) =  pr1(N:-N:-1,j+1) 
     rh1(:,my-j) = -rh1(N:-N:-1,j+1)
  enddo

  !  ------ scale and prefactor 
  rh1 = rh1/h
  call bandec7(pr1,my)

  deallocate (kap,jj,nsp,mat,vec) 

end subroutine prederiv1

!**********************************************************************
! prederiv2:
!
! Computes the coefficients for the compact finite differences
!     for the second derivatives of bandwidth 2*N+1
!     including the schemes for the end points. 
! References:  Disk notebook ctr 2007, pp. 9-14
!     and maple \javi\wall\sergio\pade7.ms   
!     and       \javi\wall\sergio\pade7end_sp.ms   
!
! Scheme is: sum_{j=-N}^N a_ij f'_{i+j} = 
!            sum_{j=-N}^N b_ij f_{i+j}
!
! At the ends, the inner band-width is maintained,
!    but outer limit is j=-i (for j=0:N-1)
! 
! Although this is a general routine it is optimized for N=3
!    It uses in general M spectral points to dispersion relation
!       at khat = pi*kap(1:M).  Here M=3.
!
!*********************************************************************
subroutine prederiv2
  use yderiv
  implicit none

  integer N, M
  integer,allocatable:: jj(:),nsp(:)
  real*8, allocatable:: mat(:,:),vec(:),kap(:)

  integer j,k,kk,kp
  real*8  pi

  ! ----------- define overall parameters ----
  N=3    ! semi-band width
  M=2    ! number of spectral points
  allocate (kap(M),nsp(N))
  kap(1) = 0.9d0
  kap(2) = 0.5d0
  pi     = 4d0*atan(1d0)
  kap    = kap*pi

  allocate ( pr2(-N:N,my), pr2i(-N:N,my), rh2(-N:N,my), lapmat(-N:N,my)  )
  pr2 = 0d0
  rh2 = 0d0
  N2  = N
  !-----------------------------------------------------------------------
  !     central part
  ! ----------------------------------------------------------------------
  allocate (mat(2*N+1, 2*N+1), vec(2*N+1), jj(N))
  vec = 0d0
  mat = 0d0
  do j=1,N
     jj(j) = j
  enddo

  mat(1,N+2:2*N+1)      = 1d0   ! order 0
  mat(1,N+1)          = 5d-1
  do k=2,2*N+1-M              ! order 2*k-2 
     kk = 2*k-2
     mat(k,1:N)       =  kk*(kk-1)*jj**(kk-2)    
     mat(k,N+2:2*N+1) = -jj**kk      
  enddo
  vec(2)            = -1d0 ! the second-order derivative of y^2 at j=0
  
  do k=1,M                       ! the spectral points
     mat(2*N+1-M+k,1:N)       =  kap(k)**2*cos(kap(k)*jj)    
     mat(2*N+1-M+k,N+2:2*N+1) =  cos(kap(k)*jj)
  enddo
  mat(2*N+2-M:2*N+1,N+1)  = 5d-1              ! the cos(0) on the rhs      
  vec(2*N+2-M:2*N+1)      = -0.5d0*kap**2     ! the derivative at 0 on the lhs      

  call gaussj(mat,2*N+1,2*N+1,vec,1,1)

  do j=1,N
     rh2( j,N+1:my-N) =  vec(N+1+j)
     rh2(-j,N+1:my-N) =  vec(N+1+j)
     pr2( j,N+1:my-N) =  vec(j)
     pr2(-j,N+1:my-N) =  vec(j)
  enddo
  rh2(0,N+1:my-N) = vec(N+1)
  pr2(0,:)        = 1d0 

  deallocate (mat,vec,jj)

  !-----------------------------------------------------------------------
  !
  !       end points
  !       mixed order and spectral, look in notebook
  !       only here for M=3 spectral points
  !-----------------------------------------------------------------------

  if (M.ne.2 .or. N.ne.3) then
     write(*,'(x,a,i4,a,i4)') &
        '% the endpoints in prederiv2 only optimized for M=2, N=3; M=',M,'  N=',N
  endif

  nsp(1)    = 0   ! how many spectral relations to use at endpoint kp=1:N
  nsp(2)    = 0   
  nsp(3)    = 1 

  allocate (mat(0:4*N+1,-N:3*N+1), vec(0:4*N+1) , jj(-N:N))
  do j=-N,N
     jj(j) = j
  enddo

  !    ------------  kp is the point
  do kp=1,N 

     mat=0d0
     vec=0d0
     
     do k=0,N-kp                !  the coefficients that are zero      
        mat(2*k,k-N)     = 1d0
        mat(2*k+1,N+k+1) = 1d0
     enddo
     kk=2*(N-kp+1)
     if (kp.eq.1) then
        mat(kk,N) =1d0  
        kk=kk+1
        mat(kk,3*N) =1d0  
        kk=kk+1
        mat(kk,3*N-1) =1d0  
        kk=kk+1
     endif
     if (kp.eq.2) then
        mat(kk,3*N) =1d0  
        kk=kk+1
     endif

     mat(kk,0) = 1d0            !  a(0)=1
     vec(kk)   = 1d0
     do k=1,nsp(kp)             !  -- the spectral relations
        mat(kk+1,-N:N)       =  kap(k)**2*cos(kap(k)*jj)
        mat(kk+1,N+1:3*N+1)  =  cos(kap(k)*jj)
        mat(kk+2,-N:N)       =  kap(k)**2*sin(kap(k)*jj)
        mat(kk+2,N+1:3*N+1)  =  sin(kap(k)*jj)
        kk=kk+2
     enddo

     kk=kk+1
     mat(kk,N+1:3*N+1)  = 1d0        !  -- the zeroth order   
     kk=kk+1
     mat(kk,N+1:3*N+1)  = jj         !  -- the first order   
     do k=2,4*N
        kk=kk+1
        if (kk>4*N+1) exit
        mat(kk,-N:N)      =  k*(k-1)*jj**(k-2)
        if (k.eq.2)  mat(kk,0)=2d0   !  -- fix 0**0 in second order
        mat(kk,N+1:3*N+1) = -jj**k 
     enddo

     call gaussj(mat,4*N+2,4*N+2,vec,1,1)

     pr2(:,kp) =  vec(0:2*N)
     rh2(:,kp) =  vec(2*N+1:4*N+1)
  enddo

  ! -----------------  the other end
  do j=0,N-1
     pr2(:,my-j) =  pr2(N:-N:-1,j+1) 
     rh2(:,my-j) =  rh2(N:-N:-1,j+1)
  enddo

  !  ------ scale and prefactor 
  rh2 = rh2/h**2
  pr2i= pr2
  call bandec7(pr2i,my)

  deallocate (kap,jj,nsp,mat,vec) 


  !!!!!!!!!!!!!!!!!!  compute explicit first derivative for end-point
  allocate(endder(0:N), mat(0:N,0:N))
  mat(0,:) = 1d0
  endder   = 0d0
  endder(1)= 1d0
  
  do j=0,N
     do k=1,N
        mat(k,j) = j**k
     enddo
  enddo
  call gaussj(mat,N+1,N+1,endder,1,1)

  deallocate(mat)

end subroutine prederiv2



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

subroutine derivyc(u,du,iord,m)
  use yderiv

  implicit none
  integer i,j,iord,m
  real*8 u(m,my),du(m,my)

  if (iord.eq.1) then    !!!!!  first derivative
     
     do i=1,m
        do j=1,N1
           wk1(j) = sum(rh1(1-j:N1,j)*u(i,1:j+N1))
        enddo

        do j=N1+1,my-N1
           wk1(j) = sum(rh1(:,j)*u(i,j-N1:j+N1))
        enddo

        do j=my-N1+1,my
           wk1(j) = sum(rh1(-N1:my-j,j)*u(i,j-N1:my))
        enddo

        call banbks7(pr1,my,wk1)
        du(i,:)=wk1
     enddo

  else                 !!!!!  second derivative 

     do i=1,m
        do j=1,N2
           wk1(j) = sum(rh2(1-j:N2,j)*u(i,1:j+N2))
        enddo

        do j=N2+1,my-N2
           wk1(j) = sum(rh2(:,j)*u(i,j-N2:j+N2))
        enddo

        do j=my-N2+1,my
           wk1(j) = sum(rh2(-N2:my-j,j)*u(i,j-N2:my))
        enddo

        call banbks7(pr2i,my,wk1)
        du(i,:)=wk1
     enddo

  endif

end subroutine derivyc


!************************************************************************
!                                                                       !
!    subroutine lapcdy                                                  !
!                                                                       !
!    solves laplace's equation (d_yy-k2)p = u                           !
!                    Dirichlet or Neumann b.c.                          !
!    Input:                                                             !
!         u: rhs. Assumed to be complex(m=2) or real(m=1) of size my    !
!    Output:                                                            !
!         p: solution  (can be same as u)                               !
!                                                                       !
!************************************************************************

subroutine lapcdy(u,k2,p,m,bc1,bc2,flag)
  use yderiv

  implicit none
  integer i,j,m
  real*8 u(m,my),p(m,my),k2,bc1(m),bc2(m)
  character*1 flag

  lapmat=rh2-k2*pr2
  if (flag.eq.'d') then   ! --- Dirichlet ----
     lapmat(:,1)  = 0d0
     lapmat(:,my) = 0d0
     lapmat(0,1)  = 1d0
     lapmat(0,my) = 1d0
  else                    ! --- Neumann -------
     if  (k2.eq.0) then      ! --- avoid singularity of Laplace 
        lapmat(:,1)  = 0d0
        lapmat(0,1)  = 1d0     
     else
        lapmat(0:N2,1)  =  endder
     endif
     lapmat(-N2:0,my)= -endder(N2:0:-1)
  endif

  call bandec7(lapmat,my)

  !-------PREPARO EL TERMINO INDEPENDIENTE [pr2]*{u} 
  
  do i=1,m

     do j=2,N2
        wk1(j) = sum(pr2(1-j:N2,j)*u(i,1:j+N2))
     enddo
     
     do j=N2+1,my-N2
        wk1(j) = sum(pr2(:,j)*u(i,j-N2:j+N2))
     enddo
     
     do j=my-N2+1,my-1
        wk1(j) = sum(pr2(-N2:my-j,j)*u(i,j-N2:my))
     enddo
     wk1(1)  = bc1(i)
     wk1(my) = bc2(i)
     
     call banbks7(lapmat,my,wk1)
     p(i,:) = wk1

  enddo

end subroutine lapcdy


!************************************************************************
!                                                                       !
!    subroutine lapc                                                    !
!                                                                       !
!    solves laplace's equation (d_yy-k2)p = u                           !
!                              p(0)=p(L) = 0                            !
!    Input:                                                             !
!         u: rhs. Assumed to be complex of size my                      !
!    Output:                                                            !
!         p: solution  (can be same as u)                               !
!                                                                       !
!************************************************************************
subroutine lapc(u,k2,p)
  use yderiv

  implicit none
  integer i,j
  real*8 u(2,my),p(2,my),k2


  lapmat=rh2-k2*pr2
  lapmat(:,1)  = 0d0
  lapmat(:,my) = 0d0
  lapmat(0,1)  = 1d0
  lapmat(0,my) = 1d0
  call bandec7(lapmat,my)

  !-------PREPARO EL TERMINO INDEPENDIENTE [pr2]*{u} 
  
  do i=1,2

     do j=2,N2
        wk1(j) = sum(pr2(1-j:N2,j)*u(i,1:j+N2))
     enddo
     
     do j=N2+1,my-N2
        wk1(j) = sum(pr2(:,j)*u(i,j-N2:j+N2))
     enddo
     
     do j=my-N2+1,my-1
        wk1(j) = sum(pr2(-N2:my-j,j)*u(i,j-N2:my))
     enddo
     wk1(1)  = 0d0
     wk1(my) = 0d0
     
     call banbks7(lapmat,my,wk1)
     p(i,:) = wk1

  enddo

end subroutine lapc




!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    prepares commons and things for using cfdiff and laps later
!                             jj/may/05
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine derivadas(n,dy)  
  use yderiv
  implicit none
  integer n
  real*8  dy

  my = n
  h  = dy
!     ---  start doing something -------------
  allocate (wk1(my), wk2(my))

  call prederiv1
  call prederiv2


end subroutine derivadas


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


















