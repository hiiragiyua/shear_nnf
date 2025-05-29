module interp
  ! ----- dimensions ---------
  integer mgalx2,mgalz2,my2,mx1_2,mz1_2
  integer mgalzp,myp,mzp,mgalz1p,nxymax,buffsize
  real*8, allocatable:: wkr(:,:,:)
  real*8, allocatable:: vor2, phi2(:,:,:)

end module interp 

program interpolation

  use ctes
  use ctes2
  use running
  use statistics 

  implicit none
  include "mpif.h"

  integer iproc,itags,newtag,imess,i,j
  real*8  val
  real*8, allocatable:: vor(:),phi(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: hv(:),hg(:),dvordy(:),chwk(:)

  real*8, allocatable:: vor2(:),phi2(:)
  real*8, allocatable:: u002(:),w002(:)
  real*8, allocatable:: hv2(:),hg2(:),dvordy2(:),chwk2(:)

  integer istat(MPI_STATUS_SIZE),ierr
  character*3 ext, fname*80
  !                              /*   initializes everything    */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)
  !--------------- initializes commons and things
  call initcr(0)
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
  !  -------  allocate rest of buffers 
  allocate( hv(buffsize),hg(buffsize),dvordy(buffsize))
  hv = 0.d0; hg = 0.d0; dvordy = 0.d0;

  if (myid.eq.0) then
     write(*,*) 'Note: interp and scaling has not implemented'
     stop
     if (numerop.eq.1) then
        write(*,*) 'Please use 1 processor'
     end if
     write(*,*) 'base file name', &
          ' (including path, default is file output root = [fileinp]'
     !read(*,*) fileinp
     write(*,*) 'read file index (from var.***, to var.*** )'
     read(*,*) istart
     read(*,*) iend
     ishear = 0;
     !write(*,*) 'Which fields do you need ([u v w o1 o2 o3 lap.v]: 0 or 1)'
     !read(*,*) ifil(1),ifil(2),ifil(3),ifil(4),ifil(5),ifil(6),ifil(7)
     ifil= 0;
     ! Note: do not forget to set the number of elements in MPI_BCAST
     !write(*,*) 'Do you want spectram? (0 or 1)'
     !read(*,*) igetspec
     igetspec = 0
     nosvemos = 1
  end if

  call  MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(ishear,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(ifil,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(igetspec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(nosvemos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (nosvemos.eq.0) iend=istart-1 ! force to stop

  !if(myid.eq.0) write(*,*) 'go to read file ...'
  call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
  !if(myid.eq.0) write(*,*) '... file read'

  ! ====================================================================
  buffsize2 = mx2*max(maxval(jend2-jbeg2+1)*mz2,maxval(kend2-kbeg2+1)*my2)

  ! --------------- allocates big buffers
  allocate(vor2(buffsize2), phi2(buffsize2), chwk2(buffsize2), u002(my2), w002(my2) )
  !  -------  allocate rest of buffers 
  allocate( hv2(buffsize2),hg2(buffsize2),dvordy2(buffsize2))
  hv2 = 0.d0; hg2 = 0.d0; dvordy2 = 0.d0;

  call intvorphi(vor,phi,u00,w00,vor2,phi2,u002,w002)

  ! write interporated fields
  call escru2(vor2,phi2,u002,w002,hv2,hg2,chwk2) 
  !
  !                    /* finalize procedure      */
  if (myid.ne.0) then
     itags=myid
     call MPI_SEND(1.,1,MPI_REAL8,0,itags,MPI_COMM_WORLD,ierr)
  else
     do iproc=1,numerop-1
        call MPI_RECV(val,1,MPI_REAL8,MPI_ANY_SOURCE,MPI_ANY_TAG, & 
                                      MPI_COMM_WORLD,istat,ierr)
        imess=istat(MPI_TAG)
!        write(*,*) 'process',imess,'over'
        newtag=100
        call MPI_SEND(1.,1,MPI_REAL8,imess,newtag,MPI_COMM_WORLD,ierr)
     enddo
  endif

  if (myid.ne.0) then
     call MPI_RECV(val,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  else
!     write(*,*) 'process',0,'over'
  endif
  call MPI_FINALIZE(ierr)   !! used to be statement 200

end program interpolation

subroutine intvorphi(vor,phi,u00,w00,vor2,phi2,u002,w002)


end subroutine intvorphi


subroutine escru2(vor2,phi2,u002,w002,hv2,hg2,chwk2)

  use ctes2

end subroutine 
