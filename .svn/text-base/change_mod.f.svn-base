      subroutine comm_setup

      use ctes
      implicit none
      include "mpif.h"

      integer istat(MPI_STATUS_SIZE),ierr
      integer mybu,mzbu,i,j,k,iproc

      allocate (myslice(0:numerop-1))
      myslice = 0.d0
      
      do iproc=0,numerop-1
         if (iproc.ne.myid) then
            mzbu = kend(iproc)-kbeg(iproc)+1
            call MPI_TYPE_VECTOR(mmy,mx*mzbu,mx*mz,MPI_REAL8,
     &           myslice(iproc),ierr)
            call MPI_TYPE_COMMIT(myslice(iproc),ierr)
         endif
      enddo

      end subroutine

     
      subroutine chikj2jik(xy,xz,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV               c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes

      implicit none 
      include 'mpif.h'
       
      real*8 xy(0:mx-1,0:mz1,0:*),xz(0:2*my-1,0:mx1,kb:*),wk1(*),
     &     wk2(0:mx-1,kb:ke,0:*)
        
      integer istat(MPI_STATUS_SIZE),ierr
      integer iproc,pproc,nsetotr,ipoxz,mzbu
      integer i,j,k,mmz2
      integer pnodes
       
      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c NOTE:
cc  wk2 dimensioned at least (mx*mzp*myp)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/

      mmz2=mmz*mx
       
      !if (myid.eq.0) then
      commtimer = commtimer-MPI_WTIME()
      !endif
      
      if (mod(numerop,3)==0) then
         pnodes=numerop*4/3
      else
         pnodes=numerop
      endif
      
      do pproc=1,pnodes-1
         iproc=ieor(myid,pproc)
         
         if (iproc<numerop) then
            nsetotr=mmz2*(jend(iproc)-jbeg(iproc)+1)
            !ipoxz=1+mmz2*(jbeg(iproc)-1) ! for pointer starting from 1
            ipoxz=1+mmz2*(jbeg(iproc)) ! for pointer starting from 0
            
            call MPI_SENDRECV(xy(0,kbeg(iproc),0),1,myslice(iproc),
     .           iproc,0,wk1(ipoxz),nsetotr,
     .           MPI_REAL8,iproc,0,
     .           MPI_COMM_WORLD,istat,ierr)
            
         endif
      enddo
      
      !if (myid.eq.0) then
      commtimer = commtimer+MPI_WTIME()
      !endif
       
      mzbu=kend(myid)-kbeg(myid)+1
      !ipoxz=1+mmz2*(jb-1) ! for channel 
      ipoxz=1+mmz2*(jb)

      call pacy2z(wk1(ipoxz),xy)
      
      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif
      
      do j=0,my1
         do k=kb,ke 
            do i=0,mx1       
               xz(2*j  ,i,k) = wk2(2*i  ,k,j)
               xz(2*j+1,i,k) = wk2(2*i+1,k,j)
            enddo
         enddo
      enddo
      
      !if (myid.eq.0) then
      transtimer = transtimer+MPI_WTIME()
      !endif
      end
      
      subroutine chjik2ikj(xz,xy,wk1,wk2)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV               c/
c sends a block (mx,mz,jb:je) to a (mx,kb:ke,my) c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes

      implicit none 
      include "mpif.h"
      
      real*8 xz(0:2*my-1,0:mx1,kb:*),xy(0:mx-1,0:mz1,0:*) 
      real*8 wk1(*),wk2(0:mx-1,kb:ke,0:*)
        
      integer istat(MPI_STATUS_SIZE),ierr
      integer iproc,pproc,nsetots,ipoxz,mzbu
      integer i,j,k,mmy2,mmz2
      integer pnodes

      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

      mmz2=mmz*mx
      mmy2=mmy*mx

      !if (myid.eq.0) then
      transtimer = transtimer-MPI_WTIME()
      !endif
       
      do k=kb,ke
         do j=0,my1
            do i=0,mx1       
               wk2(2*i  ,k,j) = xz(2*j  ,i,k)
               wk2(2*i+1,k,j) = xz(2*j+1,i,k)
            enddo
         enddo
      enddo

      !if (myid.eq.0) then
      transtimer = transtimer+MPI_WTIME()
      commtimer = commtimer-MPI_WTIME()
      !endif

      if (mod(numerop,3)==0) then
         pnodes=numerop*4/3
      else
         pnodes=numerop
      endif


      do pproc=1,pnodes-1
         iproc=ieor(myid,pproc)
         
         if (iproc<numerop) then
            nsetots=(jend(iproc)-jbeg(iproc)+1)*mmz2
            !ipoxz=1+mmz2*(jbeg(iproc)-1) ! for channel
            ipoxz=1+mmz2*(jbeg(iproc))
            
            call MPI_SENDRECV(wk1(ipoxz),nsetots,
     .           MPI_REAL8,iproc,0,xy(0,kbeg(iproc),0),
     .           1,myslice(iproc),
     .           iproc,0,MPI_COMM_WORLD,istat,ierr)
         endif
      enddo

      !if (myid.eq.0) then
      commtimer = commtimer+MPI_WTIME()
      !endif

      mzbu=kend(myid)-kbeg(myid)+1
      !ipoxz=1+mmz2*(jb-1) ! for channel
      ipoxz=1+mmz2*(jb)
      call unpacz2y(wk1(ipoxz),xy)

      end



      subroutine unpacz2y(xyi,xyo)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use ctes
      
      implicit none 
      integer mzbu
      real*8 xyi(mx,kb:ke,*),xyo(mx,0:mz1,*)
      integer  i,j,k

      do j=1,mmy
         do k=kb,ke
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end


      subroutine pacy2z(xyo,xyi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use ctes

      implicit none
      integer mzbu
      real*8 xyi(mx,0:mz1,*),xyo(mx,kb:ke,*)
      integer i,j,k

      do j=1,mmy
         do k=kb,ke
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end
