#define ROTATE
#undef ROTATE
c/********************************************************************/
c/*   FAST TRANSFORMS INTERFACE PACKAGE                              */
c/*..................................................................*/
c/*  MODULE FOR MPI SP2                                              */
c/*..................................................................*/
c/*  CPP flags: ROTATE     transpose matrix before cft               */
c/*     this coule be effective for a large xz-plane, if             */
c/*     the tranpose function is used since it must be optimized     */
c/********************************************************************/
        subroutine fourxz(fou,phys,phys2,iopt,n)
        use ctes
        implicit none
        integer iopt,n
        real*8 fou(*),phys(*)
        real*8 phys2(*)
        integer j,ipy
c/********************************************************************/
c/*                                                                  */
c/*    makes a 3-dimensional real to complex fourier transform       */
c/*      from F-F-Physical   (array fou)                             */
c/*  to/ from Phys-Phys-Phys (array phys)                            */
c/*                                                                  */
c/*       iopt >=0    ===>  inverse transforms( fou ---> fis)        */
c/*       iopt < 0    ===>  direct  transforms( fis ---> fou)        */
c/*                                                                  */
c/*       does the expanding and compacting itself                   */
c/*                                                                  */
c/*  NOTE:                                                           */
c/*    fis ---> fou : supposes high frecuency modes are not          */
c/*                   needed (dealiasing) and throws them away       */
c/*                                                                  */
c/********************************************************************/

c
c       /* physical to fourier transforms  */
        if (iopt.lt.0) then

           call rft(phys,mgalx+2,n*mgalz,-1)
c
#ifdef ROTATE
           call cxz2zx(phys,phys2)
           call cft(phys2,2,2*mgalz,mx1+1,-1)
           call pack0rotate(phys2,fou,n)
#else
           call cft(phys,mgalx+2,2,mx1+1,-1) ! OK
           call dcopy((mgalx+2)*mgalz,phys,1,phys2,1)
           call pack0(phys2,fou,n)  ! affected by optimization between -O1 and -O2
#endif
c        /* fourier to physical transforms  */
         else

#ifdef ROTATE
           call fill00rotate(fou,phys2,n)
           call cft(phys2,2,2*mgalz,mx1+1,1)
           call czx2xz(phys,phys2)
           call rft(phys,mgalx+2,n*mgalz,1)
#else
           call fill00(fou,phys2,n) ! note: fou phys should not share the memory 
           call cft(phys2,mgalx+2,2,mx1+1,1)
           call rft(phys2,mgalx+2,n*mgalz,1)
           call dcopy((mgalx+2)*mgalz,phys2,1,phys,1)
#endif

         endif

         end


c/********************************************************************/
c/*                                                                  */
c/*         fills high armonics of a 3D fourier representation       */
c/*         (f-f-T) with  zeroes to get extra points enough to       */
c/*         perform dealiasing.                                      */
c/*                                                                  */
c/*     input:                                                       */
c/*       fou: variable at fourier space                             */
c/*                                                                  */
c/*    output:                                                       */
c/*      dfou: variable at fourier space ready to dealiasing         */
c/*                                                                  */
c/********************************************************************/
      subroutine fill00(fou,foud,n)
      use ctes
      implicit none
      integer n
      complex*16 fou(0:mx1,0:mz1,n),foud(0:mgx,0:mgalz1,n),zero

      integer i,j,k,kk

      zero = dcmplx(0.d0,0.d0)

      !foud = zero ! do not set zero, since fou and foud may share its memory
      !call dcopy(2*(mgx+1)*mgalz,0.d0,0,foud,1)

      do j=1,n

         do k=0,nz1
            do i=0,mx1
               foud(i,k,j)=fou(i,k,j)
            enddo
         enddo

         do k=nz2,mgalz1
            kk = k-mgalz1+mz1
            do i=0,mx1
               foud(i,k,j)=fou(i,kk,j)
            enddo
         enddo
 
         do k=nz1+1,nz2-1
            do i=0,mx1
               foud(i,k,j)=zero
           enddo
         enddo

         do k=0,mgalz1
            do i=mx1+1,mgx
               foud(i,k,j)=zero
            enddo
         enddo

      enddo

      end

      subroutine fill00_original_fixed(fou,foud,n)
      use ctes
      implicit none
      integer n
      complex*16 fou(0:mx1,0:mz1,n),foud(0:mgx,0:mgalz1,n),zero

      integer i,j,k,kk

      zero = dcmplx(0.d0,0.d0)
      !foud = zero ! do not set zero, since fou and foud may share its memory
      !DEC$ NOVECTOR 
      do 10 j=n,1,-1
         !DEC$ NOVECTOR 
         do 70 k=mgalz1,nz2,-1
            kk = k-mgalz1+mz1
            !DEC$ NOVECTOR 
            do 80 i=mx1,0,-1
               foud(i,k,j)=fou(i,kk,j)
 80         continue
 70      continue

         !DEC$ NOVECTOR 
         do 20 k=nz1,0,-1
            !DEC$ NOVECTOR 
            do 30 i=mx1,0,-1
               foud(i,k,j)=fou(i,k,j)
 30         continue
 20      continue

         !DEC$ NOVECTOR 
         do 50 k=nz2-1,nz1+1,-1
            !DEC$ NOVECTOR 
            do 60 i=mx1,0,-1
               foud(i,k,j)=zero
 60         continue
 50      continue

         !DEC$ NOVECTOR 
         do 25 k=mgalz1,0,-1
            !DEC$ NOVECTOR 
            do 35 i=mgx,mx1+1,-1
               foud(i,k,j)=zero
 35         continue
 25      continue

 10   continue

      end

      subroutine fill00rotate(fou,foud,n)
      use ctes
      implicit none
      integer n
      complex*16 fou(0:mx1,0:mz1,n),foud(0:mgalz1,0:mgx,n),zero

      integer i,j,k,kk

      zero = dcmplx(0.d0,0.d0)
      do j=n,1,-1

         do i=0,mx1
            do k=0,nz1
               foud(k,i,j)=fou(i,k,j)
            enddo
         enddo

         do i=0,mx1
            do k=nz1+1,nz2-1
               foud(k,i,j)=zero
            enddo
         enddo

         do i=0,mx1
            do k=nz2,mgalz1
               kk = k-mgalz1+mz1
               foud(k,i,j)=fou(i,kk,j)
            enddo
         enddo

         do i=mx1+1,mgx
            do k=0,mgalz1
               foud(k,i,j)=zero
            enddo
         enddo
      enddo

      end



c/********************************************************************/
c/*                                                                  */
c/*         fills high armonics of a 3D fourier representation       */
c/*         (f-f-T) with  zeroes to get extra points enough to       */
c/*         perform dealiasing.                                      */
c/*                                                                  */
c/*     input:                                                       */
c/*      dfou: variable at fourier space after dealiasing with       */
c/*            zeroes at high frecuencies                            */
c/*                                                                  */
c/*    output:                                                       */
c/*       fou: packed variable at fourier space (if fou=fdou         */
c/*            fou overwrites fdou)                                  */
c/*                                                                  */
c/********************************************************************/
      subroutine pack0(foud,fou,n)
      use ctes
      implicit none
      integer n
      complex*16 fou(0:mx1,0:mz1,n),foud(0:mgx,0:mgalz1,n)
      integer i,j,k,kk

      do j=1,n

         do k=0,nz1

            do i=0,mx1
               fou(i,k,j)=foud(i,k,j)
            enddo
         enddo

         do k=nz2,mgalz1
            kk = k-mgalz1+mz1

            do i=0,mx1
               fou(i,kk,j)=foud(i,k,j)
            enddo
         enddo

      enddo


      end

      subroutine pack0rotate(foud,fou,n)
      use ctes
      implicit none
      integer n
      complex*16 fou(0:mx1,0:mz1,n),foud(0:mgalz1,0:mgx,n)
      integer i,j,k,kk
      do j=1,n

         do i=0,mx1
            do k=0,nz1
               fou(i,k,j)=foud(k,i,j)
            enddo
         enddo

         do i=0,mx1
            do k=nz2,mgalz1
            kk = k-mgalz1+mz1
               fou(i,kk,j)=foud(k,i,j)
            enddo
         enddo

      enddo

      end

c/********************************************************************/
c/*                                                                  */
c/*      computes the derivative of u respect ch and gives it        */
c/*      back en du (for an array of the form array(mx,mz,my)        */
c/*      in f-f-t space)                                             */
c/*                                                                  */
c/*     input:                                                       */
c/*         u: array with fourier coefficients                       */
c/*        ch: indicates the direction of the derivative             */
c/*            (options ch='x', or 'z')                              */
c/*    output:                                                       */
c/*        du: derivative coefficients (if du=u then du over-        */
c/*            write u)                                              */
c/*                                                                  */
c/********************************************************************/
      subroutine deriv(u,du,ch,n)
      use ctes
      implicit none
      integer n
      complex*16 u(0:mx1,0:mz1,n),du(0:mx1,0:mz1,n)
      character ch
      complex*16    dk

      complex*16 ci,zero
      integer i,j,k

      ci=dcmplx(0d0,1d0)
      zero = dcmplx (0d0,0d0)

      if(ch.eq.'x') then

         do 10 j=1,n
            do 20 k=0,mz1
               do 30 i=0,mx1
                  du(i,k,j)=xalp(i)*u(i,k,j)
 30            continue
 20         continue
 10      continue

      elseif(ch.eq.'z') then

         do 40 j=1,n
            do 50 k=0,mz1
               dk=xgam(k)
               do 60 i=0,mx1
                  du(i,k,j)=dk*u(i,k,j)
 60            continue
 50         continue
 40      continue

      else
          write(*,*) 'error en deriv'
      endif
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine cxz2zx(phys,phys2)
      use ctes
      implicit none

      integer i,k,kk
      complex*16 phys(mgx+1,mgalz), phys2(mgalz,mgx+1)

c         do k=1,mgalz
c      do i=1,mgx+1
c            phys2(k,i) = phys(i,k) 
c         enddo
c      enddo

      phys2 = transpose(phys)

c      do kk=1,mgalz,blockingik2ki
c         do i=1,mgx+1
c            do k=kk,kk+blockingik2ki-1
c            phys2(k,i) = phys(i,k) 
c            enddo
c         enddo
c      enddo

      end
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine czx2xz(phys,phys2)
      use ctes
      implicit none

      integer i,k,kk
      complex*16 phys(mgx+1,mgalz), phys2(mgalz,mgx+1)

c      do k=1,mgalz
c         do i=1,mgx+1
c            phys(i,k) = phys2(k,i) 
c         enddo
c      enddo

      phys = transpose(phys2)

c      do kk=1,mgalz,blockingki2ik
c         do i=1,mgx+1
c            do k=kk,kk+blockingki2ik-1
c            phys(i,k) = phys2(k,i) 
c            enddo
c         enddo
c      enddo

      end
    
