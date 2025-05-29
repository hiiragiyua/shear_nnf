      subroutine comm_setup

      use ctes

      implicit none
      include "mpif.h"

      integer istat(MPI_STATUS_SIZE),ierr
      integer i,j,k,iproc

      allocate (scount(0:numerop-1)) ! for AllToAll communication
      allocate (sdispl(0:numerop-1)) ! for AllToAll communication
      allocate (rcount(0:numerop-1)) ! for AllToAll communication
      allocate (rdispl(0:numerop-1)) ! for AllToAll communication

      iproc=0
      scount(0) = mx*mmy*(kend(iproc)-kbeg(iproc)+1)
      sdispl(0) = 0
      do iproc=1,numerop-1
         scount(iproc) = mx*mmy*(kend(iproc)-kbeg(iproc)+1)
         sdispl(iproc) = sdispl(iproc-1) + scount(iproc-1)
      enddo
       
      iproc=0
      rcount(0) = mx*(jend(iproc)-jbeg(iproc)+1)*mmz
      rdispl(0) = 0             
      do iproc=1,numerop-1    
         rcount(iproc) = mx*(jend(iproc)-jbeg(iproc)+1)*mmz
         rdispl(iproc) = rdispl(iproc-1) + rcount(iproc-1)
      enddo
      
      end

      subroutine chikj2jik(xy,xz,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses ALLTOALLV  (xy -> xz)  c/
c sends a block (mx,jb:je,mz) to a (mx,my,kb:ke) ???? 
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes
      use timer

      implicit none 
      include "mpif.h"
      real*8 xy(0:mx-1,0:mz1,jb:*),xz(0:2*my-1,0:mx1,kb:*),
     &       wk1(*),wk2(0:mx-1,kb:ke,0:*),isend(mx,mz*mmy)  
        
      integer istat(MPI_STATUS_SIZE),ierr
      integer iproc,i,j,k

c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c NOTE: wk2 dimensioned at least (mx*mzp*myp)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
     
      transtimer = transtimer-MPI_WTIME()
      call pretrans(xy,isend)
      transtimer = transtimer+MPI_WTIME() 
      
      commtimer = commtimer-MPI_WTIME()
      call MPI_ALLTOALLV(isend, scount, sdispl, MPI_REAL8,
     &     wk1, rcount, rdispl, MPI_REAL8,
     &     MPI_COMM_WORLD, ierr)
      commtimer  = commtimer+MPI_WTIME()

      transtimer = transtimer-MPI_WTIME()
      do k=kb,ke 
         do j=0,my1
            do i=0,mx1       
               xz(2*j  ,i,k) = wk2(2*i  ,k,j)
               xz(2*j+1,i,k) = wk2(2*i+1,k,j)
            enddo
         enddo 
      enddo
      transtimer = transtimer+MPI_WTIME()

      end
     
      subroutine pretrans(xy,isend)
      use ctes

      implicit none
      include "mpif.h"
      integer istat(MPI_STATUS_SIZE),ierr
      real*8 xy(0:mx-1,0:mz1,mmy),isend(0:mx-1,0:mz*mmy-1)
      integer i,j,k,k1,k2,iproc

      do iproc=0,numerop-1
         k=kend(iproc)-kbeg(iproc)+1
         do j=1,mmy
            k1=kbeg(iproc)+(j-1)*k+(mmy-1)*kbeg(iproc)
            k2=k1+k-1
            isend(:,k1:k2)=xy(:,kbeg(iproc):kend(iproc),j)
         enddo
      enddo
   
      end

      subroutine chjik2ikj(xz,xy,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV  (xz -> xy)   c/
c sends a block (my,mx,kb:ke) to a (mx,mk,jb:je)
c sends a block (mx,mz,jb:je) to a (mx,kb:ke,my) ???? c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
!      use mpi
      use ctes
      use timer

      implicit none 
      include "mpif.h"
      real*8 xz(0:2*my-1,0:mx1,kb:*),xy(0:mx-1,0:mz1,0:*)
      real*8 wk1(*),wk2(0:mx-1,kb:ke,0:*),irecv(mx,mz*mmy)
      integer istat(MPI_STATUS_SIZE),ierr
      integer i,j,k,iproc
      
      transtimer = transtimer-MPI_WTIME()
      do k=kb,ke
         do j=0,my1
            do i=0,mx1
               wk2(2*i  ,k,j) = xz(2*j  ,i,k)
               wk2(2*i+1,k,j) = xz(2*j+1,i,k)
            enddo
         enddo
      enddo
      transtimer = transtimer+MPI_WTIME()

      commtimer = commtimer-MPI_WTIME()      
      call MPI_ALLTOALLV(wk2, rcount, rdispl, MPI_REAL8,
     &     irecv, scount, sdispl, MPI_REAL8,
     &     MPI_COMM_WORLD, ierr)
      commtimer = commtimer+MPI_WTIME()

      transtimer = transtimer-MPI_WTIME()
      call posttrans(irecv,xy)
      transtimer = transtimer+MPI_WTIME()

      end

      subroutine posttrans(irecv,xy)
      use ctes

      implicit none
      include "mpif.h"
      integer istat(MPI_STATUS_SIZE),ierr
      real*8 xy(0:mx-1,0:mz1,mmy),irecv(0:mx-1,0:mz*mmy-1)
      integer i,j,k,k1,k2,iproc

      do iproc=0,numerop-1
         k=kend(iproc)-kbeg(iproc)+1
         do j=1,mmy
            k1=kbeg(iproc)+(j-1)*k+(mmy-1)*kbeg(iproc)
            k2=k1+k-1
            xy(:,kbeg(iproc):kend(iproc),j)=irecv(:,k1:k2)
         enddo
      enddo

      end     

