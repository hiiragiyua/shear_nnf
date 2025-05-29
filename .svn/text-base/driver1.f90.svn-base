program driver1 
  !/********************************************************************/
  !/*                                                                  */     
  !/* main routine for debugging fourier transformation                */
  !/*                                                                  */
  !/********************************************************************/
!  use mpi
  use ctes
  use running

  implicit none
  include "mpif.h"

  integer iproc,itags,newtag,imess,i,j,idum
  real*8  val
  real*8, allocatable:: vor(:),phi(:),wk(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: hv(:),hg(:),vorwk(:),phiwk(:),dvordy(:),chwk(:)
  real*8, allocatable:: wk1dr(:),wk1dc(:)
 
  real*8  randu,dummy
  integer istat(MPI_STATUS_SIZE),ierr

  ! /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)
  
  ! version (revision number of Subversion)
  if(myid.eq.0) write(*,*) 'HSF by A.S., Version: rev.388_optimizing_vcopy'

  ! (should be corresponding to svn revision) 
  !--------------- initializes commons and things
  call initcr
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  
  vor = 0.d0; phi = 0.d0; chwk = 0.d0; u00 = 0.d0; w00 = 0.d0;
  !  
  !if(myid.eq.0) write(*,*) 'go to read file ...'

  if (readflag.eq.0) then
     ! generate i.c.
     !
     !-- original wall version (slowest to be turbulence)
     ! call getini_original(vor,phi,u00,w00,time)
     !
     !-- periodic in y (not checked, but faster to be turbulence than original)
     ! call getini(vor,phi,u00,w00,time)
     !
     !-- Taylor-Green (not checked)     
     call getini_TG(vor,phi,u00,w00,time)       
     !-- sekimoto (TG + some modes, i=2,3)
     !call getini_seki(vor,phi,u00,w00,time)
     !-- streak profile: U=s*y+du*cos(gam*z)
     !call getini_streak(vor,phi,u00,w00,time) 
     !   
     ! add ramdum disturbance
     !idum=-123456
     !do i=1,buffsize 
     !   vor(i)=vor(i)+randu(idum)*(1.d-6)
     !end do
     !do i=1,buffsize 
     !   phi(i)=phi(i)+randu(idum)*(1.d-6)
     !end do
  else
     call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
  endif

  time_ini=time
  etime=1.d10
  noescru=0
  !if(myid.eq.0) write(*,*) '... file read'
  
  !  -------  allocate rest of buffers 
  allocate( hv(buffsize),hg(buffsize) )
  hv = 0.d0; hg = 0.d0; 
  allocate( phiwk(buffsize),vorwk(buffsize),dvordy(buffsize) )
  phiwk = 0.d0; vorwk = 0.d0; dvordy = 0.d0;
  allocate( wk1dr(my), wk1dc(2*my))
  wk1dr = 0.d0; wk1dc = 0.d0;
  ! --- the time loop
  if(myid.eq.0) write(*,*) 'just before cross ...'
  call chjik2ikj(phi,phi,chwk,chwk)
  call chjik2ikj(vor,vor,chwk,chwk)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) call check(vor,phi)

  !call cross(vor,phi,u00,w00,hv,hg,phiwk,vorwk,dvordy,chwk,wk1dr,wk1dc)  
  ! ------------------------------------------------------------------
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)   !! used to be statement 200

end program driver1

subroutine check(vorc,phic)

  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"
  complex*16 phic(0:mx1,0:mz1,jb:je),  & 
       vorc(0:mx1,0:mz1,jb:je)

  integer i
  real*8, allocatable:: tmpx(:),tmpx2(:)

  allocate (tmpx(mgalx+2),tmpx2(mgalx+2)) 

  do i=0,mx1
     write(10,*) dreal(vorc(i,0,jb)),dimag(vorc(i,0,jb))
  end do

  tmpx(:) = 0.d0
  call dcopy((mx1+1)*2,vorc(0,0,jb),1,tmpx,1)
  call rft(tmpx,mgalx+2,1,1)

  do i=1,mgalx+2
     write(11,*) tmpx(i)
  enddo

  tmpx2(:) = 0.d0
  call dcopy((mx1+1)*2,vorc(0,0,jb),1,tmpx2,1)

  call rftw(tmpx2,mgalx+2,1,1)

  do i=1,mgalx+2
     write(12,*) tmpx2(i)
  enddo

  deallocate(tmpx,tmpx2)

end subroutine check
