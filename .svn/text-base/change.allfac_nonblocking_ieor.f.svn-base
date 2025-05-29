       subroutine chikj2jik(xy,xz,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV  (xy -> xz)   c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
!      use mpi
      use ctes

      implicit none 
      include "mpif.h"
      real*8 xy(0:mx-1,0:mz1,0:*),xz(0:2*my-1,0:mx1,kb:*),wk1(*),
     &        wk2(0:mx-1,kb:ke,0:*)
        
      integer istat((numerop-1)*2*MPI_STATUS_SIZE+1),ierr
      integer irequests
      integer irequestr((numerop-1)+1)
      integer ireqr,iflag
      integer iproc,iproc_s,iproc_r,nsetotr,ipoxz,pproc
      integer i,j,k,mmz2

      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c NOTE: wk2 dimensioned at least (mx*mzp*myp)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/

       mmz2=mmz*mx

       ireqr=0
       iflag=0

       if (myid.eq.0) then
         commtimer = commtimer-MPI_WTIME()
       endif

       !do iproc=0,numerop-1
       do pproc=0,numerop-1

          !iproc_s = mod(myid+iproc,numerop)
          !iproc_r = mod(numerop+myid-iproc,numerop)
          
          iproc_s=ieor(myid,pproc)
          iproc_r=ieor(myid,pproc)

          if(iflag.gt.0)then
             call MPI_Wait(irequests,istat,ierr)
          endif
          if(iproc_s.ne.myid)then

             nsetotr=mmz2*(jend(iproc_r)-jbeg(iproc_r)+1)
             ipoxz=1+mmz2*(jbeg(iproc_r))

             iflag=1
             call MPI_ISEND(xy(0,kbeg(iproc_s),0),1,myslice(iproc_s), 
     .            iproc_s,0,MPI_COMM_WORLD,irequests,ierr)
             ireqr=ireqr+1
             call MPI_IRECV(wk1(ipoxz),nsetotr,MPI_REAL8,
     .            iproc_r,0,MPI_COMM_WORLD,irequestr(ireqr),ierr)
          else

            ipoxz=1+mmz2*(jb)
            call pacy2z(wk1(ipoxz),xy,iproc_s)

          endif

       enddo

       call MPI_Waitall(ireqr,irequestr,istat,ierr)
       call MPI_Wait(irequests,istat,ierr)

       if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
         transtimer = transtimer-MPI_WTIME()
       endif


       do k=kb,ke 
          do j=0,my1
             do i=0,mx1       
                xz(2*j  ,i,k) = wk2(2*i  ,k,j)
                xz(2*j+1,i,k) = wk2(2*i+1,k,j)
             enddo
          enddo
       enddo

       if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
       endif

       end


      subroutine chjik2ikj(xz,xy,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV               c/
c sends a block (mx,mz,jb:je) to a (mx,kb:ke,my) c/
c xz -> xy                                       c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
!      use mpi
      use ctes

      implicit none 
      include "mpif.h"
      real*8 xz(0:2*my-1,0:mx1,kb:*),xy(0:mx-1,0:mz1,0:*)
      real*8 wk1(*),wk2(0:mx-1,kb:ke,0:*)
      integer istat((numerop-1)*2*MPI_STATUS_SIZE+1),ierr,ireqs,iflag
      integer irequests((numerop-1)+1)
      integer irequestr
      integer iproc,iproc_s,iproc_r,nsetots,ipoxz,pproc
      integer i,j,k,mmy2,mmz2
      
      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

      ireqs=0
      iflag=0

      mmz2=mmz*mx
      mmy2=mmy*mx

      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif

      do k=kb,ke
         do j=0,my1
            do i=0,mx1
               wk2(2*i  ,k,j) = xz(2*j  ,i,k)
               wk2(2*i+1,k,j) = xz(2*j+1,i,k)
            enddo
         enddo
      enddo

      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
         commtimer = commtimer-MPI_WTIME()
      endif

       do pproc=0,numerop-1

          if(iflag.gt.0)then
             call MPI_Wait(irequestr,istat,ierr)
          endif
!          iproc_s = mod(myid+iproc,numerop)
!          iproc_r = mod(numerop+myid-iproc,numerop)

          iproc_s=ieor(myid,pproc)
          iproc_r=ieor(myid,pproc)

          if(iproc_s.ne.myid)then
             
             nsetots=(jend(iproc_s)-jbeg(iproc_s)+1)*mmz2
             ipoxz=1+mmz2*(jbeg(iproc_s))

c             call MPI_SENDRECV(wk1(ipoxz),nsetots,
c     .            MPI_REAL8,iproc_s,0,xy(0,kbeg(iproc_r),0),
c     .            1,myslice(iproc_r),
c     .            iproc_r,0,MPI_COMM_WORLD,istat,ierr)
             
             iflag=1
             ireqs=ireqs+1
             call MPI_ISEND(wk1(ipoxz),nsetots,MPI_REAL8,
     .            iproc_s,0,MPI_COMM_WORLD,irequests(ireqs),ierr)
             call MPI_IRECV(xy(0,kbeg(iproc_r),0),1,myslice(iproc_r),
     .            iproc_r,0,MPI_COMM_WORLD,irequestr,ierr)

          else

             ipoxz=1+mmz2*(jb)
             call unpacz2y(wk1(ipoxz),xy,iproc_r)

          endif

       enddo

      call MPI_Waitall(ireqs,irequests,istat,ierr)
      call MPI_Wait(irequestr,istat,ierr)
      if (myid.eq.0) commtimer = commtimer+MPI_WTIME()

      end

      subroutine zcopy_one(n,a,b)
      implicit none 
      integer n,ia,ib
      complex*16 a(n),b(n)

      b=a

      end
      

      subroutine unpacz2y(xyi,xyo,ip)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use ctes

      implicit none 
      integer  ip
      real*8 xyi(mx,kend(ip)-kbeg(ip)+1,mmy),xyo(mx,0:mz1,mmy)

      xyo(1:mx,kbeg(ip):kend(ip),:)=xyi

      end


      subroutine pacy2z(xyo,xyi,ip)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use ctes

      implicit none
      integer ip
      real*8 xyi(mx,0:mz1,mmy),xyo(mx,kend(ip)-kbeg(ip)+1,mmy)

      xyo=xyi(1:mx,kbeg(ip):kend(ip),:)

      end








