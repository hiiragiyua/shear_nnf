c ==============================================================
c              initialize cft_w
c     uses fftw_3.2,2 supposed to work as the old rftdouble
c     not tested 
c     compile with flags: 
c     
c ==============================================================
      subroutine cfti(n)
      use cfttmp
      implicit none
      include "fftw3.f"

      integer  n

      !if (nmax.lt.n) write(*,*) 'cftw error: set nmax bigger than',n
      
      if (allocated(fdum)) then
      
         call dfftw_destroy_plan(planf)
         call dfftw_destroy_plan(planb)
         !     write(*,*) 'destroy plans in cfti',planf,planb
         !     write(*,*) 'cleanup fftw in cfti',planf,planb
         call dfftw_cleanup()         
         deallocate(fdum,bdum) 
         !write(*,*) 'cleanup fftw in cfti',planf,planb
      
      end if

      if (.not. allocated(fdum)) allocate(fdum(n),bdum(n))
      
      nb=n
      nf=n
      dnf=1d0/dfloat(n)
      dnb=1d0

      fdum = 0.d0
      
      call dfftw_plan_dft_1d(planf,nf,fdum,fdum,
     .                                FFTW_FORWARD,FFTW_MEASURE )

      bdum = 0.d0
      call dfftw_plan_dft_1d(planb,nb,bdum,bdum,
     .                                FFTW_BACKWARD,FFTW_MEASURE )
      
      
      end


c ================================================================ c
c             the complex cft                                      c
c ================================================================ c

      subroutine cft(c,ise,isa,m,iopt)
      use cfttmp
      implicit none
      integer isa,m,iopt,i,j,ise,ise1,isa1,j0,j00
      complex*16 c(*)
      
      isa1 = isa/2
      ise1 = ise/2
 
      j00  = 1
            
      if (iopt<0) then

         do j=1,m
            j0 = j00
            do i=1,nf 
               fdum(i)=c(j0)
               j0 = j0 + ise1
            enddo
            
            call dfftw_execute(planf)
            
            j0 = j00
            do i=1,nf
               c(j0)=dnf*fdum(i)
               j0 = j0 + ise1
            enddo
            j00 = j00 + isa1
         enddo
         
      else   
         
         do j=1,m
            j0 = j00
            do i=1,nb 
               bdum(i)=c(j0)
               j0 = j0 + ise1
            enddo
            
            call dfftw_execute(planb)
            j0 = j00
            do i=1,nb 
               c(j0)=bdum(i)
               j0 = j0 + ise1
            enddo   
            j00 = j00 + isa1
         enddo
         
      endif

      end

! fftw in y-direction

      subroutine cfti_y(n)
      use cftytmp
      implicit none
      include "fftw3.f"

      integer  n

      
      if (allocated(fdum)) then
      
         call dfftw_destroy_plan(planf)
         call dfftw_destroy_plan(planb)
         !     write(*,*) 'destroy plans in cfti_y',planf,planb
         !     write(*,*) 'cleanup fftw in cfti_y',planf,planb
         call dfftw_cleanup()         
         deallocate(fdum,bdum) 
         !write(*,*) 'cleanup fftw in cfti_y',planf,planb
      
      end if

      if (.not. allocated(fdum)) allocate(fdum(n),bdum(n))
      
      nb=n
      nf=n
      dnf=1d0/dfloat(n)
      dnb=1d0

      fdum = 0.d0
      
      call dfftw_plan_dft_1d(planf,nf,fdum,fdum,
     .                                FFTW_FORWARD,FFTW_MEASURE )

      bdum = 0.d0
      call dfftw_plan_dft_1d(planb,nb,bdum,bdum,
     .                                FFTW_BACKWARD,FFTW_MEASURE )
      
      
      end


c ================================================================ c
c             the complex cft                                      c
c ================================================================ c

      subroutine cfty(c,ise,isa,m,iopt)
      use cftytmp
      implicit none
      integer isa,m,iopt,i,j,ise,ise1,isa1,j0,j00
      complex*16 c(*)
      
      isa1 = isa/2
      ise1 = ise/2
 
      j00  = 1
            
      if (iopt<0) then

         do j=1,m
            j0 = j00
            do i=1,nf 
               fdum(i)=c(j0)
               j0 = j0 + ise1
            enddo
            
            call dfftw_execute(planf)
            
            j0 = j00
            do i=1,nf
               c(j0)=dnf*fdum(i)
               j0 = j0 + ise1
            enddo
            j00 = j00 + isa1
         enddo
         
      else   
         
         do j=1,m
            j0 = j00
            do i=1,nb 
               bdum(i)=c(j0)
               j0 = j0 + ise1
            enddo
            
            call dfftw_execute(planb)
            j0 = j00
            do i=1,nb 
               c(j0)=bdum(i)
               j0 = j0 + ise1
            enddo   
            j00 = j00 + isa1
         enddo
         
      endif

      end
