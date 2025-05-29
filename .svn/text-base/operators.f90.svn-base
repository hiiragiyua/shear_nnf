subroutine add_force(rhv,iopt)

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: rhv  
  integer iopt

  uprim=force_roll ! take care that uprim is overwritten here

  if (iopt.eq.1) then
     ! 1pair-streamwise-roll force F=(Fy,Fz), div.F = 0
     !  then hv += (d^2/dz^2 + d^2/dy^2)*Fy
     !       hg += 0
     ! for k=1:mgalz
     !    hvF(k,:) = amp*cos(gam*z_k)*(-4+12*y.^2) 
     !            + amp*((1-y.^2).^2)*gam*gam*(-cos(gam*z_k));
     ! end do
     if ((kb.le.1).and.(ke.ge.1)) then
        !write(*,*) myid,': adding roll-force'
        rhv(:,0,1) = rhv(:,0,1) + &
             dcmplx(-0.5d0,0.d0)*uprim*( -4.d0*(1.d0-3.d0*y*y) &
             -gam*gam*(1.d0-y*y)*(1.d0-y*y) )
     endif
  elseif (iopt.eq.2) then
     ! 2pair-streamwise-roll
     !  for k=1:mgalz
     !     rhv(k,:) = amp*cc(k)*4*y.*(-3+5*y.^2) 
     !              + amp*y.*((1-y.^2).^2)*gam*gam*(-cc(k));
     !  end
     if ((kb.le.1).and.(ke.ge.1)) then
        rhv(:,0,1) = rhv(:,0,1) + dcmplx(0.5d0,0.d0)*uprim*( &
             4.d0*y*(-3.d0+5.d0*y*y) - y*(1.d0-y*y)*(1.d0-y*y)*gam*gam )
     end if
  end if

end subroutine add_force

subroutine add_damping(phi,vor,u00,w00,rhv,rhg,rf0u,rf0w,iopt)

  use ctes
  use running
  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: phi,vor,rhv,rhg
  real(8),dimension(0:my1) :: u00,w00,rf0u,rf0w
  integer iopt,i,k
  real(8) aa,up
  real(8),dimension(:),allocatable :: Hw ! masking function or filter function

  allocate(Hw(0:my1))

  !aa=Ly/8.d0;
  !aa=Lz/2.d0;

  ! set masked-damping parameter scaled by Lz
  aa=Lz*damp_aa
  up=damp_up  ! small permeability

  Hw = 1.d0 - 0.25*( 1.d0+tanh(6.d0*(aa-y)/Lz + 3.d0))*(1.d0+tanh(6.d0*(aa+y)/Lz+3.d0) );

  if (iopt.eq.1) then
    ! do k=kb,ke      
    !    do i=0,mx1
    !
    !       ! damping
    !       rhv(:,i,k) = rhv(:,i,k) - up*Hw(:)*phi(:,i,k) 
    !       rhg(:,i,k) = rhg(:,i,k) - up*Hw(:)*vor(:,i,k) 
    ! 
    !    enddo
    ! enddo
    ! 
     rf0u = rf0u - up*Hw(:)*u00(:)
    ! rf0w = rf0w - up*Hw(:)*w00(:)     

  end if

  deallocate(Hw)

end subroutine add_damping

subroutine moving_frame(vor,phi,xoff,iopt)
!subroutine moving_frame(phi,mtime,iopt)
  ! only works for (j,i,k)
  ! projection onto the moving frame of mean shear flow in fourier [exp(-i k_x ss y Dt)]

  use ctes
  use running
  ! do not add module 'use bcs'

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  integer i,j,k,iopt
  real(8) mtime,xoff

  complex(8) sydt(0:my1,0:mx1)

  !if (myid.eq.0) write(*,*) 'moving_frame', iopt 
  !if (iadd_mode.le.4) then
  mtime = xoff ! nondimensional time (including shear) 
  !else (iadd_mode.eq.5)
  !   mtime = xoff
  !end if


  if (iopt.eq.1) then 
     if (mtime*s.lt.0.d0) write(*,*) 'moving_frame error!', & 
          'mtime should be positive'
  elseif (iopt.eq.-1) then 
     if (mtime*s.gt.0.d0) write(*,*) 'moving_frame error!', &
          'mtime should be negative'
  endif

  do i=0,mx1
     do j=0,my1
        sydt(j,i) = cdexp(xalp(i)*y(j)*mtime)
     end do
  end do

  do k=kb,ke
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
     do i=0,mx1
        !do j=0,my1
           vor(:,i,k) = vor(:,i,k)*sydt(:,i)
        !end do
     end do
     !$OMP END PARALLEL DO
  end do
  do k=kb,ke
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
     do i=0,mx1
        !do j=0,my1
           phi(:,i,k) = phi(:,i,k)*sydt(:,i)
        !end do
     end do
     !$OMP END PARALLEL DO
  end do

end subroutine moving_frame

subroutine moving_frame_t(temp,xoff,iopt)
  ! only works for (j,i,k)
  ! projection onto the moving frame of mean shear flow in fourier [exp(-i k_x ss y Dt)]

  use ctes
  use running
  ! do not add module 'use bcs'

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: temp
  integer i,j,k,iopt
  real(8) mtime,xoff

  complex(8) sydt(0:my1,0:mx1)

  !if (myid.eq.0) write(*,*) 'moving_frame', iopt 
  !if (iadd_mode.le.4) then
  mtime = xoff ! nondimensional time (including shear) 
  !else (iadd_mode.eq.5)
  !   mtime = xoff
  !end if


  if (iopt.eq.1) then 
     if (mtime*s.lt.0.d0) write(*,*) 'moving_frame error!', & 
          'mtime should be positive'
  elseif (iopt.eq.-1) then 
     if (mtime*s.gt.0.d0) write(*,*) 'moving_frame error!', &
          'mtime should be negative'
  endif

  do i=0,mx1
     do j=0,my1
        sydt(j,i) = cdexp(xalp(i)*y(j)*mtime)
     end do
  end do

  do k=kb,ke
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
     do i=0,mx1
        !do j=0,my1
           temp(:,i,k) = temp(:,i,k)*sydt(:,i)
        !end do
     end do
     !$OMP END PARALLEL DO
  end do

end subroutine moving_frame_t

subroutine add_extshear(vor,phi,fac)

  ! for iadd_mode = 5

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  integer i,j,k,iopt
  real(8) fac

  do k=kb,ke
     do i=0,mx1
        do j=0,my1  
           phi(j,i,k) = phi(j,i,k) + fac*vor(j,i,k)*xgam(k)
        end do
     end do
  end do

end subroutine add_extshear

subroutine deriv_xfou(vor,phi,iopt)

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  integer i,j,k,iopt

  if (iopt.eq.1) then 
     ! the first derivatives in x
     do k=kb,ke
        do i=0,mx1
           do j=0,my1
              vor(j,i,k) = xalp(i)*vor(j,i,k)
           end do
        end do
     end do
     do k=kb,ke
        do i=0,mx1
           do j=0,my1
              phi(j,i,k) = xalp(i)*phi(j,i,k)
           end do
        end do
     end do
  else
     ! not implemented
  end if

end subroutine deriv_xfou

subroutine deriv_zfou(vor,phi,iopt)

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  integer i,j,k,iopt

  if (iopt.eq.1) then 
     ! the first derivatives in x
     do k=kb,ke
        do i=0,mx1
           do j=0,my1
              vor(j,i,k) = xgam(k)*vor(j,i,k)
           end do
        end do
     end do
     do k=kb,ke
        do i=0,mx1
           do j=0,my1
              phi(j,i,k) = xgam(k)*phi(j,i,k)
           end do
        end do
     end do
  else
     ! not implemented
  end if

end subroutine deriv_zfou

subroutine deriv_yfou(vorwk,phiwk,u00wk,w00wk,iopt)

  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vorwk, phiwk
  real(8),dimension(0:my1) :: u00wk, w00wk

  complex(8),dimension(:,:,:),allocatable :: dvordy, dphidy
  real(8),dimension(:),allocatable :: du00dy, dw00dy
  integer i,j,k,iopt
  complex(8) shp,shm
  ! set boundary condition outside, not here (rev1395)
  ! xoffb = xoffb_ini 
  ! xofft = xofft_ini
  do i=0,mx1
     shwkb(i)  = zexp(-xalp(i)*xoffb) ! negative shift
     shwkt(i)  = zexp(-xalp(i)*xofft) ! positive shift
  enddo  

  allocate(dvordy(0:my1,0:mx1,kb:ke),dphidy(0:my1,0:mx1,kb:ke))
  allocate(du00dy(my),dw00dy(my))
  dvordy=0.d0; dphidy=0.d0;     
  du00dy=0.d0; dw00dy=0.d0;     

  ! no need to chang to (j,i,k)

  if (iopt.eq.1) then 
     ! first derivative in y
     do k=kb,ke    
        do i=0,mx1
           shp = shwkt(i)
           shm = shwkb(i)
           ! --- d(ome_2)/ dy
           call derivyc(vorwk(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
           ! --- dv/dy 
           call derivyc(phiwk(0,i,k),dphidy(0,i,k),shp,shm,1,2,0) ! skip _addshear
        enddo
     enddo

     !--- Note: we need special treatments for the mean flow
     shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
     call derivyc(u00wk,du00dy,shp,shm,1,1,0)    ! ----- ome3_0 = -du0/dy ---
     !--- Note: derivative of fluctuation of u is pure-periodic
     call derivyc(w00wk,dw00dy,shp,shm,1,1,0)    ! ----- ome1_0 =  dw0/dy ---
  else
     ! not implemented
  end if
  
  vorwk=dvordy;
  phiwk=dphidy;

  u00wk=du00dy;
  w00wk=dw00dy;

  deallocate(du00dy,dw00dy)
  deallocate(dvordy,dphidy)

end subroutine deriv_yfou

subroutine dealiasing_y(vorwk,phiwk,u00wk,w00wk, & 
     &                  u00p,w00p,u00c,w00c,fac)

  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vorwk, phiwk
  real(8),dimension(0:my1) :: u00wk, w00wk

  ! 1d tmp arrays
  real(8),dimension(0:my1+2) :: u00p, w00p  ! extended buffer for rfty
  complex(8),dimension(0:my1) :: u00c, w00c

  integer i,j,k,iopt
  real(8) fac,djda,djdb
  complex(8) zero
  complex(8) shp,shm

  if ((iydealiasing.eq.1).or.(iydealiasing.eq.2)) then
     ! zero-zero modes 
     u00p=0.d0; w00p=0.d0;
     u00p(0:my1)=u00wk; w00p(0:my1)=w00wk;
     
     call rfty(u00p,my+2,1,-1) ! phy --> fou
     call rfty(w00p,my+2,1,-1) ! phy --> fou
     do j=my23,my1+2
        u00p(j)=0.d0
        w00p(j)=0.d0
     enddo
     call rfty(u00p,my+2,1,+1) ! phy <-- fou
     call rfty(w00p,my+2,1,+1) ! phy <-- fou
     u00wk=u00p(0:my1); w00wk=w00p(0:my1);
     
     zero=dcmplx(0.d0,0.d0)
     u00c=zero; w00c=zero;
  end if
  if (iydealiasing.eq.1) then
     do k=kb,ke
        ! dealiasing only for kx=0
        i=0
        u00c=vorwk(0:my1,i,k);
        w00c=phiwk(0:my1,i,k);
        call cfty(u00c,2,1,1,-1) ! phy --> fou
        call cfty(w00c,2,1,1,-1) ! phy --> fou
        !djd=fac*dfloat(ny1*i)*(Ly/Lx)/dfloat(mx1)
        !djd=fac*dfloat(i)*(Ly/Lx)
        ! from Kida Tanaka
        !djda=dfloat(ny1)/dfloat(mx1)*(-fac-Lx/Ly)*Ly/Lx*dfloat(i)
        !djdb=dfloat(ny1)/dfloat(mx1)*(-fac)*Ly/Lx*dfloat(i)
        
        djda=dfloat(ny1)/dfloat(mx1)*(-fac-Lx/Ly)*Ly/Lx*dfloat(i)
        djdb=dfloat(ny1)/dfloat(mx1)*(-fac)*Ly/Lx*dfloat(i)
        do j=0,my1
           if ( (j.gt.(ny1 + djdb)).and.(j.lt.(ny2+djdb)) ) then 
              ! if(myid.eq.0) write(*,*) 'djd cut',j,djdb,fac
              u00c(j)=zero
              w00c(j)=zero
              !elseif ( (j.gt.(ny1 + djda)).and.(j.lt.(ny2+djda)) ) then 
              !  ! if(myid.eq.0) write(*,*) 'djd cut',j,djdb,fac
              !   u00c(j)=0.d0
              !   w00c(j)=0.d0
           end if
        enddo
        call cfty(u00c,2,1,1,+1) ! phy <-- fou
        call cfty(w00c,2,1,1,+1) ! phy <-- fou
        vorwk(0:my1,i,k)=u00c;
        phiwk(0:my1,i,k)=w00c; 
     end do
  elseif (iydealiasing.eq.2) then
     ! dealiasing considering in moving frame ! does not work well
     !write(*,*) 'do not use iydealiasing ... stop'
     !stop
     ! 
     do k=kb,ke
        do i=0,mx1
           u00c=vorwk(0:my1,i,k);
           w00c=phiwk(0:my1,i,k);
           call cfty(u00c,2,1,1,-1) ! phy --> fou
           call cfty(w00c,2,1,1,-1) ! phy --> fou
           !
           ! from Kida Tanaka
           if (fac.le.0.d0) then
              djda=dfloat(ny1)/dfloat(mx1)*(-fac-Lx/Ly)*Ly/Lx*dfloat(i) ! -
              djdb=dfloat(ny1)/dfloat(mx1)*(-fac)*Ly/Lx*dfloat(i)       ! +
           else
              !djdb=dfloat(ny1)/dfloat(mx1)*(-fac+Lx/Ly)*Ly/Lx*dfloat(i) ! +
              !djda=dfloat(ny1)/dfloat(mx1)*(-fac)*Ly/Lx*dfloat(i)       ! -
              djda=dfloat(ny1)/dfloat(mx1)*(-fac+Lx/Ly)*Ly/Lx*dfloat(i) ! +
              djdb=dfloat(ny1)/dfloat(mx1)*(-fac)*Ly/Lx*dfloat(i)       ! -
           end if

           do j=0,my1
              if ( (j.ge.(ny1 + djdb)).and.(j.le.(ny2+djdb)) ) then ! just on moving frame
              !if ( (j.gt.(ny1 + djda)).and.(j.lt.(ny2 + djdb)) ) then ! remeshing-like
                 !if ((myid.eq.0).and.(i.eq.mx1).and.(kb.eq.0)) then 
                 !   write(*,*) 'djd cut',j,djdb,-fac/(Lx/Ly)
                 !end if
                 u00c(j)=zero
                 w00c(j)=zero
              elseif ( (j.gt.(ny1 + djda)).and.(j.lt.(ny2+djda)) ) then 
                 !  ! if(myid.eq.0) write(*,*) 'djd cut',j,djdb,fac
                 u00c(j)=0.d0
                 w00c(j)=0.d0
              end if
           enddo
           call cfty(u00c,2,1,1,+1) ! phy <-- fou
           call cfty(w00c,2,1,1,+1) ! phy <-- fou
           vorwk(0:my1,i,k)=u00c;
           phiwk(0:my1,i,k)=w00c;              
        end do        
     end do

  elseif (iydealiasing.eq.3) then
           
     shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
     call filterc(u00wk,u00wk,shp,shm,1,0)
     call filterc(w00wk,w00wk,shp,shm,1,0)
     ! apply a compact filter
     do k=kb,ke    
        do i=0,mx1
           shp = shwkt(i)
           shm = shwkb(i)
           call filterc(vorwk(0,i,k),vorwk(0,i,k),shp,shm,2,1)
           call filterc(phiwk(0,i,k),phiwk(0,i,k),shp,shm,2,0) ! skip _addshear
        enddo
     enddo


  end if

end subroutine dealiasing_y

subroutine phase_shift_fou(vor,phi,shiftx,shiftz)

  ! take care: the sign is defined to apply xf: phase_shift()*xf(0,Tp) ==> x0
  ! then positive shift of x0 ==> xf

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  integer i,j,k,iopt,iproc,ierr
  real(8) shiftx,shiftz
  real(8),dimension(0:numerop-1) :: checkx, checkz
  complex(8),dimension(:,:),allocatable :: shift

  allocate(shift(0:mx1,kb:ke))
  shift=0.d0

  ! check the consistency of shift paramters
  call MPI_GATHER(shiftx,1,MPI_REAL8,checkx,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(shiftz,1,MPI_REAL8,checkz,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  if (myid.eq.0) then
     do iproc=1,numerop-1
        if ((abs(shiftx-checkx(iproc)).gt.tiny).or. &
             (abs(shiftz-checkz(iproc)).gt.tiny)) then
           write(*,*) ' ERROR in phase shift: CHECK shift parameters for all proc.'
           stop
        endif
     end do
  endif

  do k=kb,ke
     do i=0,mx1
        shift(i,k)=exp(dcmplx(xalp(i)*(-shiftx)))*exp(dcmplx(xgam(k)*(-shiftz)))
     end do
  end do

  do k=kb,ke
     do i=0,mx1
        do j=0,my1
           vor(j,i,k) = vor(j,i,k)*shift(i,k)
        end do
     end do
  end do
  do k=kb,ke
     do i=0,mx1
        do j=0,my1
           phi(j,i,k) = phi(j,i,k)*shift(i,k)
        end do
     end do
  end do

  deallocate(shift)

end subroutine phase_shift_fou

subroutine phase_shifty(vorwk,phiwk,u00wk,w00wk, & 
     &                  shifty,bcfac)

  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vorwk, phiwk
  real(8),dimension(0:my1) :: u00wk, w00wk
  real(8),dimension(0:numerop-1) :: checky
  ! 1d tmp arrays
  !real(8),dimension(:),allocatable :: u00p,w00p ! extended buffer for rfty
  complex(8),dimension(:),allocatable :: u00c,w00c ! for cfty


  integer i,j,k,iopt,iproc,ierr
  real(8) bcfac,djda,djdb,bet,shifty
  complex(8) zero
  complex(8) shp,shm
  complex(8),dimension(:),allocatable :: shift

 ! check the consistency of shift paramters
  call MPI_GATHER(shifty,1,MPI_REAL8,checky,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  allocate(shift(0:my1))
  shift=0.d0
  allocate(u00c(0:my1),w00c(0:my1))
  u00c=0.d0; w00c=0.d0;
  bet  = 2.d0*pi/Ly ! note that the factor of 2
  do j=0,my/2-1
     shift(j)=exp(dcmplx(cii*bet*dfloat(j)*(-shifty)))
  end do

  do j=my/2,my1 ! bug fixed at rev.1565
     shift(j)=exp(dcmplx(cii*(-bet)*dfloat(my-j)*(-shifty)))
  end do

  ! rev.1564
  !do j=ny1+1,ny2-1
  !   shift(j) = 0.d0 ! 2/3 dealiasing in the shifted grid..
  !end do
  ! zero-zero modes 
  u00c=0.d0; w00c=0.d0;
  do j=0,my1
     u00c(j)=dcmplx(u00wk(j),0d0); 
     w00c(j)=dcmplx(w00wk(j),0d0);
  end do
  call cfty(u00c,2,1,1,-1) ! phy -> fou
  call cfty(w00c,2,1,1,-1) ! phy --> fou
  do j=1,my1
     u00c(j)=u00c(j)*shift(j)
     w00c(j)=w00c(j)*shift(j)
  end do
  call cfty(u00c,2,1,1,+1) ! phy <-- fou
  call cfty(w00c,2,1,1,+1) ! phy <-- fou
  u00wk=dreal(u00c(0:my1)); 
  w00wk=dreal(w00c(0:my1));
 

  !bcfac = (xwkb/Lx - int(xwkb/Lx) )*Lx/Ly

  call moving_frame(vorwk,phiwk,-bcfac,+1)    
  !if (myid.eq.0) write(*,*)  '      vor phi are distorted to satisfy B.C.', xwkb_ini,fac
  do k=kb,ke
     do i=0,mx1
        u00c(:)=vorwk(:,i,k)
        call cfty(u00c,2,1,1,-1) ! phy -> fou
        do j=1,my1
           u00c(j) = u00c(j)*shift(j)
        end do
        call cfty(u00c,2,1,1,+1) ! phy <-- fou
        vorwk(:,i,k) = u00c(:)
     end do
  end do
  do k=kb,ke
     do i=0,mx1
        w00c(:)=phiwk(:,i,k)
        call cfty(w00c,2,1,1,-1) ! phy -> fou
        do j=1,my1
           w00c(j) = w00c(j)*shift(j)
        end do
        call cfty(w00c,2,1,1,+1) ! phy <-- fou
        phiwk(:,i,k) = w00c(:)
     end do
  end do
  call moving_frame(vorwk,phiwk,bcfac,-1)    

  deallocate (u00c,w00c)
  deallocate (shift)

end subroutine phase_shifty


subroutine get_phase(vor,phi,u00,w00,phasex,phasey,phasez,ireport)
  !
  ! Note:vor, phi, u00, w00 should not be broken
  ! phasex = Im[oy(1,0)(y=0)]/Re[oy(1,0)(y=0)], later shifted by atan(phasex)/(2*pi)*Lx
  ! phasez = Im[oy(0,1)(y=0)]/Re[oy(0,1)(y=0)]
  !        = -Im[oy(0,mz)(y=0)]/Re[oy(0,mz)(y=0)] (complex conjugate)
  ! phasey
  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"

  integer i,j,k,ierr,ireport
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  real(8),dimension(0:my1) :: u00, w00
  real(8),dimension(:),allocatable :: u00tmp
  real*8 phasex,phasey,phasez

  allocate(u00tmp(0:my1+2))
  u00tmp=0.d0
  j=my/2
  phasex=0.d0; phasez=0.d0

  !do k=kb,ke
  !   write(*,*) myid,k,'vor:',vor(j,0,k)
  !enddo
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !call system("sleep 0.5")
  !do k=kb,ke
  !   write(*,*) myid,k,'phi:',phi(j,0,k)
  !enddo
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !call system("sleep 0.5")
  if (myid.eq.0) then
     if (abs(dreal(vor(j,1,0))).lt.tiny) then
        phasex=0.d0
        if (ireport.eq.1) write(*,*) ' iget_phase in x error: phasex will be infty',dreal(vor(j,1,0))
     else
        phasex=dimag(vor(j,1,0))/dreal(vor(j,1,0))
     end if
     if (abs(dreal(vor(j,0,1))).lt.tiny) then
        phasez=0.d0
        if (ireport.eq.1) write(*,*) ' iget_phase in z error: phasez will be infty',dreal(vor(j,0,1))
     else
        phasez=dimag(vor(j,0,1))/dreal(vor(j,0,1))
     end if
     if (ke.eq.0) then 
        write(*,*) 'error in get_phase, master has only k=0'
        stop
     endif
     if (ireport.eq.1) write(*,*) myid,' iget_phase in x,z:',phasex,phasez 
     !,atan2(real(vor(j,1,0),imag(vor(j,1,0)) )
  endif
  call  MPI_BCAST(phasex,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(phasez,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !write(*,*) myid,' iget_phase in x,z:',atan(phasex),atan(phasez)
     
  u00tmp=0.d0
  u00tmp(0:my1)=u00
 
  call rfty(u00tmp,my+2,1,-1) ! phy -> fou
  if (abs(u00tmp(2)).lt.tiny) then
     phasey=0.d0;
  else
     phasey=u00tmp(3)/u00tmp(2) ! imag(u00c(1))/real(u00c(1))
     !phasey=u00tmp(2) ! real(u00c(1))
  end if
  if ((myid.eq.0).and.(ireport.eq.1)) write(*,*)' iget_phase in y:',phasey
  if ((myid.eq.0).and.(ireport.eq.1)) write(*,*)' check u00f', u00tmp(0:3)
  call  MPI_BCAST(phasey,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! bug fixed (rev.1113)

  ! FFTy check for rfty and cfty 
  !if(myid.eq.0) then
  !   do j=0,my1+2,2
  !      write(*,*) j/2,'u00f:',u00tmp(j),u00tmp(j+1)
  !   end do
  !endif

  !u00c(0:my1)=dcmplx(u00,0.d0)
  !call cfty(u00c,2,1,1,-1) ! phy -> fou
  !if(myid.eq.0) then
  !   do j=0,my1
  !      write(*,*) j,'u00c:',u00c(j)
  !   end do
  !endif

  deallocate(u00tmp)

end subroutine get_phase


subroutine getv(phi,vv,iopt)

  use ctes
  use bcs
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vv, phi
  integer i,j,k,k1,iopt
  real(8) rk2,rk,fac
  complex(8) shp, shm
  real(8),dimension(:),allocatable:: wk1dc  

  ! BUG fixed at rev.1395, set b.c. before calling getv
  !xoffb = xoffb_ini
  !xofft = xofft_ini

  xwkt = xofft; xwkb = xoffb

  do i=0,mx1
     shb(i) = cdexp(-xalp(i)*xwkb) ! negative shift
     sht(i) = cdexp(-xalp(i)*xwkt) ! positive shift
  enddo
  if (iopt.eq.1) then ! phi-> v
     ! --- calcula la v a partir de phi ---
     do k=kb,ke
        k1 = icx(k)
        do i=0,mx1
           if ((i.ne.0).or.(k.ne.0)) then
              rk = gam2(k)+alp2(i) 
              shp = sht(i)
              shm = shb(i)
              call lapcdy(phi(0,i,k),rk,vv(0,i,k),shp,shm,2)   ! -- v
           end if
        enddo
     enddo
  elseif (iopt.eq.-1) then ! phi <- v

     allocate( wk1dc(2*my))
     wk1dc = 0.d0;

     fac = 1.d0
     do k=kb,ke          ! --- computes  lap.(v)
        do i=0,mx1
           if ((i.ne.0).or.(k.ne.0)) then
              shp = sht(i) ! use the previous one before updated by (7)
              shm = shb(i)
              rk = gam2(k)+alp2(i)
              call visc3d(vv(0,i,k),phi(0,i,k),wk1dc,shp,shm,2,rk,fac,'n',1)
           endif
        enddo
     enddo

     deallocate(wk1dc)
  end if
     
end subroutine getv


subroutine window_y(vor,phi,u00,w00)
 
  use ctes
  use bcs
  use running

  implicit none
  include "mpif.h"

  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  real(8),dimension(0:my1) :: u00, w00
  integer i,j,k,iopt,ierr
  real(8) aa
  real(8),dimension(:),allocatable :: filt
  complex(8),dimension(:,:,:),allocatable ::   vel


  !write(*,*) myid,' window_y'

  allocate(filt(0:my1))
  ! implemented (rev.1188)
  ! see: W(z), Gibson Brandt (2014, JFM), but in y-dir for HST
  !aa=Lz;
!  aa=Ly/4.d0;
  aa=Lz*0.25; ! for window norm. ! rev.1613 in diff_upo.f90
  filt=(1.d0+tanh(6.d0*(aa-y)/Lz + 3.d0))*(1.d0+tanh(6.d0*(aa+y)/Lz+3.d0) );
  filt=filt/4.d0

  allocate(vel(0:my1,0:mx1,kb:ke))
  vel=0.d0
  call getv(phi,vel,1) ! --- phi ==> v
  do k=kb,ke    
     do i=0,mx1
        vel(:,i,k) = vel(:,i,k)*filt
     enddo
  enddo
  call getv(phi,vel,-1) ! --- phi <== v

  do k=kb,ke    
     do i=0,mx1
        vor(:,i,k) = vor(:,i,k)*filt
     enddo

  enddo
  
  u00(:) = u00(:)*filt
  w00(:) = w00(:)*filt

  deallocate(vel)
  deallocate(filt)

end subroutine window_y

!function intplin(f,x,xp,n)
!! copied from rslogvis.f90 (which is removed)
!! for linear interpolation
!implicit none
!
!integer n
!real*8 intplin,f(n),x(n),xp
!
!integer i,il,ig
!
!i=1
!do while (xp>x(i))
!   il=i
!   i=i+1
!enddo
!ig=i
!
!intplin=(f(ig)-f(il))/(x(ig)-x(il))*(xp-x(il))+f(il)
!
!end

subroutine energy(vor,phi,u00,w00,energ)

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi
  real(8),dimension(0:my1) :: u00, w00
  integer i,j,k,iopt,ierr

  real(8),dimension(1:4) :: enerl, energ

  enerl(1:4)=0.d0
  energ(1:4)=0.d0
  do k=kb,ke
     i=0
     do j=0,my1  
        enerl(1) = enerl(1) + vor(j,i,k)*dconjg(vor(j,i,k))
     end do
  end do
  do k=kb,ke
     do i=1,mx1
        do j=0,my1  
           enerl(1) = enerl(1) + 2.d0*vor(j,i,k)*dconjg(vor(j,i,k))
        end do
     end do
  end do
  !
  do k=kb,ke
     i=0
     do j=0,my1  
        enerl(2) = enerl(2) + phi(j,i,k)*dconjg(phi(j,i,k))
     end do
  end do
  do k=kb,ke
     do i=1,mx1
        do j=0,my1  
           enerl(2) = enerl(2) + 2.d0*phi(j,i,k)*dconjg(phi(j,i,k))
        end do
     end do
  end do
  !
  if (myid.eq.0) then
     do j=0,my1  
        enerl(3) = enerl(3) + u00(j)*u00(j)
     end do
     do j=0,my1  
        enerl(4) = enerl(4) + w00(j)*w00(j)
     end do
  end if
  call MPI_ALLREDUCE(enerl,energ,4,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine energy

subroutine apply_filtering(vor,phi,u00,w00,iopt)

  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi

  real(8),dimension(0:my1) :: u00, w00
  integer iopt,i,j,k,ierr,iproc,ini

  complex*16 shp,shm

  shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
  call filterc(u00,u00,shp,shm,1,0)
  call filterc(w00,w00,shp,shm,1,0)
  ! apply a compact filter
  do i=0,mx1
     shb(i)  = zexp(-xalp(i)*xoffb) ! negative shift
     sht(i)  = zexp(-xalp(i)*xofft) ! positive shift
  enddo
  
  do k=kb,ke    
     do i=0,mx1
        shp = sht(i)
        shm = shb(i)
        call filterc(vor(0,i,k),vor(0,i,k),shp,shm,2,1)
        call filterc(phi(0,i,k),phi(0,i,k),shp,shm,2,0) ! skip _addshear
     enddo
  enddo
end subroutine apply_filtering

subroutine add_sym(vor,phi,u00,w00,vors,phis,u00s,w00s,chwk,iopt)

  use ctes
  use running
  use bcs
  use timer

  implicit none
  include "mpif.h"
  
  complex(8),dimension(0:my1,0:mx1,kb:ke) :: vor, phi, vors, phis

  real(8),dimension(buffsize) :: chwk
  real(8),dimension(0:my1) :: u00, w00, u00s, w00s
  integer iopt,i,j,k,km,ierr,iproc,ini
  real(8) phasex,phasey,phasez ! they are also defined in gmres_module...

  addsymtimer = addsymtimer - MPI_WTIME()

  cii=dcmplx(0.d0,1.d0)
  
 ! if (ini.eq.1) then
 !    !write(*,*) myid, 'allocate buffers for addsym'
 !    allocate( vors(0:my1,0:mx1,0:mz1),phis(0:my1,0:mx1,0:mz1) )
 !    allocate( u00s(0:my1),w00s(0:my1) )
 !    vors=0.d0; phis=0.d0
 !    u00s=0.d0; w00s=0.d0;
 ! 
 !    allocate(sgcount(0:numerop-1),rgcount(0:numerop-1),rgdispl(0:numerop-1))
 !    iproc=0
 !    sgcount(0) = mx*my*(kend(iproc)-kbeg(iproc)+1)
 !    do iproc=1,numerop-1
 !       sgcount(iproc) = mx*my*(kend(iproc)-kbeg(iproc)+1)
 !    enddo
 !    rgcount = sgcount
 !    rgdispl(0) = 0
 !    do iproc=1,numerop-1
 !       rgdispl(iproc) = rgdispl(iproc-1) + rgcount(iproc-1)
 !    enddo
 ! elseif (ini.eq.-1) then
 !    deallocate(u00s,w00s)
 !    deallocate(vors,phis)
 !    return
 ! end if

  vors=vor; phis=phi;
  u00s=u00; w00s=w00;
  !
  if (((iopt.eq.5).or.(iopt.eq.12)).or.(iopt.eq.6)) then
     ! ALLGATHERV may be slow ...
     !call MPI_ALLGATHERV(vor,my*(mx1+1)*mmz*2,MPI_DOUBLE_PRECISION, &
     !     vors,rgcount,rgdispl,MPI_DOUBLE_PRECISION, &
     !     MPI_COMM_WORLD,ierr)
     !call MPI_ALLGATHERV(phi,my*(mx1+1)*mmz*2,MPI_DOUBLE_PRECISION, &
     !     phis,rgcount,rgdispl,MPI_DOUBLE_PRECISION, &
     !     MPI_COMM_WORLD,ierr)  
     !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !
     chwk=0.d0
     !write(*,*) myid, 'change to Y-cut'
     call chjik2ikj(vors,vors,chwk,chwk) ! Y-cut
     call chjik2ikj(phis,phis,chwk,chwk) ! Y-cut
     !
  end if
  
  if ((iopt.eq.5).or.(iopt.eq.12)) then
     ! (sym5-1) averaging the upper-half domain and the lower-half domain.
     ! using shift-reflection symmetry, in Y-cut 
     !
     u00(0)=0.d0
     do j=1,my1
        u00(j)= 0.5d0*(u00s(j) - u00s(my1-j+1))
     end do
     w00=0.d0
     !
     call z_sym(vors,phis,iopt)
     !
     !write(*,*) myid, 'change back to Z-cut',iopt
     call chikj2jik(vors,vors,chwk,chwk)
     call chikj2jik(phis,phis,chwk,chwk)
     !
     !write(*,*) myid, 'add shift-rotation',iopt
     ! (sym5-2) averaging the upper- and lower-half domain using shift-rotation
     !         Z-cut
     !if (myid.eq.0) then
     !   ! kz0 modes
     !   do i=0,mx1
     !      !j=0 is special 
     !      do j=1,my1
     !         vor(j,i,0) = 0.5d0*(vors(j,i,0)-dconjg(vors(my1-j+1,i,0)));
     !         phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0)));
     !      end do
     !   end do
     !end if
     !!
     !do k=max(kb,1),ke
     !   km = mz1-k+1
     !   do i=0,mx1
     !      do j=1,my1
     !         vor(j,i,k) = 0.5d0*(vors(j,i,k)-dconjg(vors(my1-j+1,i,km))*cdexp(xgam(k)*Lz2));
     !         phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,km))*cdexp(xgam(k)*Lz2));
     !      end do
     !   end do
     !end do
     !write(*,*) myid, 'add shift-rotation',iopt
     if (myid.eq.0) then
        ! averaging the upper- and lower-half domain using shift-rotation
        ! reflectioin symmetry
        ! kz0 modes
        do i=0,mx1
           ! j=0 is a special case 
           do j=1,my1
              ! Note: this is Lx/2 shifted one
              ! of S_3 symmetry by Gibson, Halcrow & Cvitanovic (JFM, 2009), CORRIGENDUM?
              vor(j,i,0) = 0.5d0*(vors(j,i,0)+dconjg(vors(my1-j+1,i,0))*cdexp(xalp(i)*Lx2) );
              phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0))*cdexp(xalp(i)*Lx2) );
           end do
        end do
     end if
     ! 
     do k=max(kb,1),ke
        do i=0,mx1
           ! j=0 is a special case
           do j=1,my1
              vor(j,i,k) = 0.5d0*(vors(j,i,k)+dconjg(vors(my1-j+1,i,k)) & 
                   *cdexp(xgam(k)*Lz2)*cdexp(xalp(i)*Lx2) );
              phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,k)) & 
                   *cdexp(xgam(k)*Lz2)*cdexp(xalp(i)*Lx2) );
           end do
        end do
     end do
     
     
     if (iopt.eq.12) then ! sy5 + sym7 (shiftx0.25

        call phase_shift_fou(vor,phi,0.25d0*Lx,0.d0*Lz)         
        vors=vor; phis=phi;
        u00s=u00; w00s=w00;

        u00(0)=0.d0
        do j=1,my1
           u00(j)= 0.5d0*(u00s(j) - u00s(my1-j+1))
        end do
        w00(0)=0.d0
        do j=1,my1
           w00(j)= 0.5d0*(w00s(j) - w00s(my1-j+1))
        end do
        !write(*,*) myid, 'add shift-rotation',iopt
        if (myid.eq.0) then
           ! averaging the upper- and lower-half domain using shift-rotation
           ! reflectioin symmetry
           ! kz0 modes
           do i=0,mx1
              ! j=0 is a special case 
              do j=1,my1
                 vor(j,i,0) = 0.5d0*(vors(j,i,0)+dconjg(vors(my1-j+1,i,0)));
                 phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0)));
              end do
           end do
        end if
        ! 
        do k=max(kb,1),ke
           do i=0,mx1
              ! j=0 is a special case
              do j=1,my1
                 vor(j,i,k) = 0.5d0*(vors(j,i,k)+dconjg(vors(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
                 phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
              end do
           end do
        end do
        
        call phase_shift_fou(vor,phi,-0.25d0*Lx,0.d0*Lz) 
     end if        
    
  else if (iopt.eq.6) then
     ! mirror-symmetric and shift-rotation
     ! averaging the upper-half domain and the lower-half domain.
     ! using shift-reflection symmetry
     !
     u00(0)=0.d0
     do j=1,my1
        u00(j)= 0.5*(u00s(j) - u00s(my1-j+1))
     end do
     w00=0.d0
     !
     call z_sym(vors,phis,iopt)
     !
     !write(*,*) myid, 'change back to Z-cut',iopt
     call chikj2jik(vors,vors,chwk,chwk)
     call chikj2jik(phis,phis,chwk,chwk)
     !
     !write(*,*) myid, 'add shift-rotation',iopt
     ! averaging the upper- and lower-half domain using shift-rotation
     !if (myid.eq.0) then
     !   ! kz0 modes
     !   do i=0,mx1
     !      ! j=0 is a special case
     !      do j=1,my1
     !         vor(j,i,0) = 0.5d0*(vors(j,i,0)-dconjg(vors(my1-j+1,i,0)));
     !         phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0)));
     !      end do
     !   end do
     !end if
     ! 
     !do k=max(kb,1),ke
     !   km = mz1-k+1
     !   do i=0,mx1
     !      ! j=0 is a special case
     !      do j=1,my1
     !         vor(j,i,k) = 0.5d0*(vors(j,i,k)-dconjg(vors(my1-j+1,i,km))*cdexp(xgam(k)*Lz2));
     !         phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,km))*cdexp(xgam(k)*Lz2));
     !      end do
     !   end do
     !end do
     !
     !
     !write(*,*) myid, 'add shift-rotation',iopt
     if (myid.eq.0) then
        ! averaging the upper- and lower-half domain using shift-rotation
        ! reflectioin symmetry
        ! kz0 modes
        do i=0,mx1
           ! j=0 is a special case 
           do j=1,my1
              vor(j,i,0) = 0.5d0*(vors(j,i,0)+dconjg(vors(my1-j+1,i,0)));
              phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0)));
           end do
        end do
     end if
     ! 
     do k=max(kb,1),ke
        do i=0,mx1
           ! j=0 is a special case
           do j=1,my1
              vor(j,i,k) = 0.5d0*(vors(j,i,k)+dconjg(vors(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
              phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
           end do
        end do
     end do
     
   else if (iopt.eq.7) then
      ! shift-rotation-reflection, which does not conflict with mirror-symmetry
      ! averaging the upper-half domain and the lower-half domain.
      ! using shift-reflection symmetry
      !
      u00(0)=0.d0
      do j=1,my1
         u00(j)= 0.5d0*(u00s(j) - u00s(my1-j+1))
      end do
      w00(0)=0.d0
      do j=1,my1
         w00(j)= 0.5d0*(w00s(j) - w00s(my1-j+1))
      end do

      !
      !write(*,*) myid, 'add shift-rotation',iopt
      if (myid.eq.0) then
         ! averaging the upper- and lower-half domain using shift-rotation
         ! reflectioin symmetry
         ! kz0 modes
         do i=0,mx1
            ! j=0 is a special case 
            do j=1,my1
               vor(j,i,0) = 0.5d0*(vors(j,i,0)+dconjg(vors(my1-j+1,i,0)));
               phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0)));
            end do
         end do
      end if
      ! 
      do k=max(kb,1),ke
         do i=0,mx1
            ! j=0 is a special case
            do j=1,my1
               vor(j,i,k) = 0.5d0*(vors(j,i,k)+dconjg(vors(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
               phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,k))*cdexp(xgam(k)*Lz2));
            end do
         end do
      end do

   else if (iopt.eq.8) then
      ! shift-rotation-reflection, which does not conflict with shift-reflection-symmetry
      ! averaging the upper-half domain and the lower-half domain.
      ! using shift-reflection symmetry
      !
      u00(0)=0.d0
      do j=1,my1
         u00(j)= 0.5d0*(u00s(j) - u00s(my1-j+1))
      end do
      w00(0)=0.d0
      do j=1,my1
         w00(j)= 0.5d0*(w00s(j) - w00s(my1-j+1))
      end do
      if (myid.eq.0) then
        ! averaging the upper- and lower-half domain using shift-rotation
        ! reflectioin symmetry
        ! kz0 modes
        do i=0,mx1
           ! j=0 is a special case 
           do j=1,my1
              ! Note: this is Lx/2 shifted one
              ! of S_3 symmetry by Gibson, Halcrow & Cvitanovic (JFM, 2009), CORRIGENDUM?
              vor(j,i,0) = 0.5d0*(vors(j,i,0)+dconjg(vors(my1-j+1,i,0))*cdexp(xalp(i)*Lx2) );
              phi(j,i,0) = 0.5d0*(phis(j,i,0)-dconjg(phis(my1-j+1,i,0))*cdexp(xalp(i)*Lx2) );
           end do
        end do
     end if
     ! 
     do k=max(kb,1),ke
        do i=0,mx1
           ! j=0 is a special case
           do j=1,my1
              vor(j,i,k) = 0.5d0*(vors(j,i,k)+dconjg(vors(my1-j+1,i,k)) & 
                   *cdexp(xgam(k)*Lz2)*cdexp(xalp(i)*Lx2) );
              phi(j,i,k) = 0.5d0*(phis(j,i,k)-dconjg(phis(my1-j+1,i,k)) & 
                   *cdexp(xgam(k)*Lz2)*cdexp(xalp(i)*Lx2) );
           end do
        end do
     end do

   !elseif (iopt.eq.10) then 
   !  first-Fourier slice testing ... do not use this ... 
   !   !if (iremove_singularity.eq.1) then
   !  call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,0)
   !  call phase_shift_fou(vor,phi,atan(phasex)/(2.d0*pi)*Lx,atan(phasez)/(2.d0*pi)*Lz )
   !  !phase shift in y
   !  call phase_shifty(vor,phi,u00,w00,atan(phasey)/(2.d0*pi)*Ly)
   !  !call get_phase(vor,phi,u00,w00,phasex,phasey,phasez,0) ! for checking zero phase 
   !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !   !endif
   end if ! iopt

   addsymtimer = addsymtimer+MPI_WTIME()

end subroutine add_sym

subroutine z_sym(vors,phis,iopt)

  use ctes
  use running
  use bcs
  !use addsym

  ! done in Y-cut plane

  implicit none
  include "mpif.h"
  
  !complex*16  vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke) ! Z-cut
  complex(8),dimension(0:mx1,0:mz1,jb:je) :: vors, phis ! Y-cut
  integer iopt,i,j,k,km,ierr,iproc

  if ((iopt.eq.5).or.(iopt.eq.12)) then
     ! (sym5-1) averaging the upper-half domain and the lower-half domain.
     ! using shift-reflection symmetry, in Y-cut 
     ! kz0 modes
     
     do j=jb,je
        do i=0,mx1,2
           ! add shift-reflection symmetry
           vors(i,0,j) = 0.d0;
           !-- because of the following relationship 
           !%%vors(i,1,j)= -vors(i,1,j)*exp(cii*(i-1)*pi);
           !%%phis(i,1,j)= phis(i,1,j)*exp(cii*(i-1)*pi);        
           !%%vrs(i,1,j)= vrs(i,1,j)*exp(cii*(i-1)*pi); 
           !%%dvdyrs(i,1,j)= dvdyrs(i,1,j)*exp(cii*(i-1)*pi);        
        end do
     end do
     do j=jb,je
        do i=1,mx1,2
           ! add shift-reflection symmetry
           phis(i,0,j) = 0.d0;
        end do
     end do
     !% nonzero-kz
     do j=jb,je                 
        do k=1,mz1/2
           ! set_conj for kx0 mode, real(vor) = 0, imag(phi) = 0 (rev.1158)
           vors(0,k,j) = dcmplx(0.0, dimag(vors(0,k,j)))
           phis(0,k,j) = dcmplx(dreal(phis(0,k,j)), 0.d0)
           do i=0,mx1
              !% add shift-reflection symmetry
              vors(i,mz1-k+1,j) =-vors(i,k,j)*cdexp(xalp(i)*Lx2);
              phis(i,mz1-k+1,j) = phis(i,k,j)*cdexp(xalp(i)*Lx2);
           end do
        end do
     end do
  else if (iopt.eq.6) then
     ! mirror-symmetric and shift-rotation
     ! averaging the upper-half domain and the lower-half domain.
     ! using shift-reflection symmetry
     !
     ! kz0 modes
     do j=jb,je
        do i=0,mx1
           ! add mirror-symmetry
           vors(i,0,j) = 0.d0;
        end do
     end do
     !% nonzero-kz
     do j=jb,je                 
        do k=1,mz1/2
           !% add mirror-symmetry: set_conj for kx0 mode, real(vor) = 0, imag(phi) = 0 (rev.1158)
           vors(0,k,j) = dcmplx(0.0, dimag(vors(0,k,j)))
           phis(0,k,j) = dcmplx(dreal(phis(0,k,j)), 0.d0)
           do i=0,mx1
              !% add mirror-symmetry: omgy, cos(kz); phi, sin(kz) 
              vors(i,mz1-k+1,j) = -vors(i,k,j);
              phis(i,mz1-k+1,j) =  phis(i,k,j);
           end do
        end do
     end do
  else
     write(*,*) 'z-sym: this is not implemented, iopt =', iopt
     stop
  end if

  !call set_conj(vors,phis,0)  ! rev(1145)
  
end subroutine z_sym


subroutine set_conj(vor,phi,ireport)

  ! set conjugate relation for kx=0. Done in y-cut plane...

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  !complex*16 vor(0:my1,0:mx1,kb:ke),phi(0:my1,0:mx1,kb:ke)
  complex*16 vor(0:mx1,0:mz1,jb:je), phi(0:mx1,0:mz1,jb:je)
  integer i,j,k,ireport,ierr
  real*8 fac, maxre, maxim, maxre_l, maxim_l
  
  maxre_l=0.d0;
  maxim_l=0.d0;
  maxre=0.d0;
  maxim=0.d0;
  do j=jb,je 
     do k=nz+1,mz1
        maxre_l=max(maxre_l, abs(dreal(vor(0,mz-k,j) - dconjg(vor(0,k,j)))))
        maxim_l=max(maxim_l, abs(dimag(vor(0,mz-k,j) - dconjg(vor(0,k,j)))))
     enddo
  enddo

  call MPI_ALLREDUCE(maxre_l,maxre,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(maxim_l,maxim,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)  
  if (ireport.eq.1) then
     if(myid.eq.0) write(*,*) 'set_conj check:',maxre,maxim 
  elseif ((ireport.eq.0).and.((maxre.gt.1.d-6).or.(maxim.gt.1.d-6))) then
     if(myid.eq.0) write(*,*) 'set_conj: max non-conj is large:',maxre,maxim 
  end if

  do j=jb,je 
     do k=nz+1,mz1 ! bug fixed at rev.1133
        ! icx(k)=mz-k
        
        vor(0,k,j) = dconjg(vor(0,mz-k,j))
     end do
  end do
  do j=jb,je 
     do k=nz+1,mz1
        phi(0,k,j) = dconjg(phi(0,mz-k,j))
     end do
  end do

end subroutine set_conj


subroutine set_conj_t(vor,ireport)

  ! set conjugate relation for kx=0. Done in y-cut plane...

  use ctes
  use running

  implicit none
  include "mpif.h"
  
  complex*16 vor(0:mx1,0:mz1,jb:je)
  integer i,j,k,ireport,ierr
  real*8 fac, maxre, maxim, maxre_l, maxim_l
  
  maxre_l=0.d0;
  maxim_l=0.d0;
  maxre=0.d0;
  maxim=0.d0;
  do j=jb,je 
     do k=nz+1,mz1
        maxre_l=max(maxre_l, abs(dreal(vor(0,mz-k,j) - dconjg(vor(0,k,j)))))
        maxim_l=max(maxim_l, abs(dimag(vor(0,mz-k,j) - dconjg(vor(0,k,j)))))
     enddo
  enddo

  call MPI_ALLREDUCE(maxre_l,maxre,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(maxim_l,maxim,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)  
  if (ireport.eq.1) then
     if(myid.eq.0) write(*,*) 'set_conj check:',maxre,maxim 
  elseif ((ireport.eq.0).and.((maxre.gt.1.d-6).or.(maxim.gt.1.d-6))) then
     if(myid.eq.0) write(*,*) 'set_conj: max non-conj is large:',maxre,maxim 
  end if

  do j=jb,je 
     do k=nz+1,mz1 ! bug fixed at rev.1133
        ! icx(k)=mz-k        
        vor(0,k,j) = dconjg(vor(0,mz-k,j))
     end do
  end do

end subroutine set_conj_t


subroutine compute_div(u1c,u3c,dvdy,div)

  use ctes
  use running

  implicit none
  include "mpif.h"

  integer k
  complex(8),dimension(0:mx1,0:mz1) :: dvdy, div, u1c, u3c 
  
  div(0,1:mz1) = dvdy(0,1:mz1) + u3c(0,1:mz1)*xgam(1:mz1)
  do k=0,mz1
     div(1:mx1,k) = dvdy(1:mx1,k) + u1c(1:mx1,k)*xalp(1:mx1) &
          + u3c(1:mx1,k)*xgam(k)
  end do

end subroutine compute_div

subroutine check_parameter(pp)
 
  use ctes,only: numerop,myid,tiny
  
  implicit none
  include "mpif.h"

  real(8),dimension(0:numerop-1) :: chkpp
  real(8) pp
  integer iproc,ierr

  call MPI_GATHER(pp,1,MPI_REAL8,chkpp,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  if (myid.eq.0) then
     do iproc=1,numerop-1
        if (abs(pp-chkpp(iproc)).gt.tiny) then
           write(*,*) iproc,' ERR dns: CHECK for all proc. the parameter to be the same.',pp
           stop
        endif
     enddo
  end if
end subroutine check_parameter

subroutine rdtfunc(u1,u2,u3,u0,v0,w0,k0x,k0y,k0z,nu,t,ierr)

  implicit none
  ! input: u0, v0, w0, k0x, k0y, k0z, nu, t
  ! out: u1, u2, u3 ! (including viscosity damping)
  ! 
  integer ierr
  complex(8) u0,v0,w0,div,u1,u2,u3,Du
  real(8) k0x,k0y,k0z,kx,ky,kz,s,nu,t
  real(8) alp20,alp0,alp2,psir,psir0
  !
  s=1.d0 ! assuming s=1
  ! wave-numbers
  kx = k0x ;
  ky = k0y - s*k0x*t;
  kz = k0z ;
  ! check the continuity
  !
  div = k0x*u0 + k0y*v0 + k0z*w0;
  if ( abs(div).gt.1.e-12 ) then 
     write(*,*) 'Check the continuity of initial condition for RDT:', div
     ierr=1;
     return;
  end if
  !
  alp20= k0x**2 + k0z**2;
  alp0 = sqrt(alp20);
  alp2 = alp20 + ky**2;
  !
  psir = atan(ky/alp0); 
  psir0= atan(k0y/alp0);
  !
  ! check the solution without normalization (check it by maple RDT.mw)
  u2=(k0x**2+k0y**2+k0z**2)/(kx**2+ky**2+kz**2)*v0;
  u1=kx*(-ky)/(kx**2+ky**2+kz**2)/alp20*v0 - kx*(-k0y)/(k0x**2+k0y**2+k0z**2)/alp20*v0 + &
       &  (kz**2)/kx/alp20*(psir-psir0)/alp0*v0;
  u1=u1*(k0x**2+k0y**2+k0z**2) + u0;
  !
  u3=kz*(-ky)/(kx**2+ky**2+kz**2)/alp20*v0 - kz*(-k0y)/(k0x**2+k0y**2+k0z**2)/alp20*v0 - & 
       &  kz*(psir-psir0)/alp20/alp0*v0;
  u3=u3*(k0x**2+k0y**2+k0z**2) + w0;
  !
  Du=exp(-nu*t*( (k0x**2+k0y**2+k0z**2)-s*t*k0y*k0x+1d0/3d0*(s*t*k0x)**2 ));
  !
  u1=u1*Du; u2=u2*Du; u3=u3*Du;

end subroutine rdtfunc


