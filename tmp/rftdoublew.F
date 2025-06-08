c ==============================================================
c              initialize rft_w
c     uses fftw_3.3.2 supposed to work as the old rftdouble 
c   
c ==============================================================
      subroutine rfti(n)
      use rfttmp
      implicit none
      include "fftw3.f"

      integer  n

      !if (nmax.lt.n) write(*,*) 'rftw error: set nmax bigger than',n

      if (allocated(fdum)) then
      
         call dfftw_destroy_plan(planf)
         call dfftw_destroy_plan(planb)
         !     write(*,*) 'destroy plans in rfti',planf,planb
         !     write(*,*) 'cleanup fftw in rfti',planf,planb
         call dfftw_cleanup()         
         deallocate(fdum,bdum) 
         !write(*,*) 'cleanup fftw in rfti',planf,planb
      
      end if

      !if ( allocated(fdum)) deallocate(fdum,bdum) ! for the case of mgalx=my

      allocate(fdum(n+2),bdum(n+2)) 

      nb=n; nf=n
      dnf=1.d0/dfloat(n)
      dnb=1.d0
      
      fdum = 0.d0
      bdum = 0.d0

      call dfftw_plan_dft_r2c_1d(planf,nf,fdum,fdum,FFTW_MEASURE )
      call dfftw_plan_dft_c2r_1d(planb,nb,bdum,bdum,FFTW_MEASURE )
      
      deallocate (fdum,bdum)

      end


c ================================================================
c             the real  rft
c ================================================================
      subroutine rft(c,isa,m,iopt)
      use rfttmp
      implicit none
      integer isa,m,iopt,i,j
      real*8 c(isa,*)

      if (iopt.lt.0) then

         do j=1,m
            !do i=1,nf
            !   fdum(i)=c(i,j)
            !enddo
            !call dfftw_execute_dft_r2c(planf,fdum,fdum) 
            ! skip copying (rev.397)
            call dfftw_execute_dft_r2c(planf,c(1,j),c(1,j))
            do i=1,nf+2
               !c(i,j)=dnf*fdum(i)
               c(i,j)=dnf*c(i,j)
            enddo
         enddo

      else   

         do j=1,m
            !do i=1,nb+2 
            !   bdum(i)=c(i,j)
            !enddo
            !call dfftw_execute_dft_c2r(planb,bdum,bdum)
            call dfftw_execute_dft_c2r(planb,c(1,j),c(1,j))
            !do i=1,nb
            !   c(i,j)=bdum(i)
            !enddo            
            c(nb+1,j)=0.d0
            c(nb+2,j)=0.d0
         enddo

      endif

      end
!
! y direction
!

      subroutine rfti_y(n)
      use rftytmp
      implicit none
      include "fftw3.f"

      integer  n

      if (allocated(fdum)) then
      
         call dfftw_destroy_plan(planf)
         call dfftw_destroy_plan(planb)
         !     write(*,*) 'destroy plans in rfti_y',planf,planb
         !     write(*,*) 'cleanup fftw in rfti_y',planf,planb
         call dfftw_cleanup()         
         deallocate(fdum,bdum) 
         !write(*,*) 'cleanup fftw in rfti_y',planf,planb
      
      end if

      !if ( allocated(fdum)) deallocate(fdum,bdum) ! for the case of mgalx=my

      allocate(fdum(n+2),bdum(n+2)) 

      nb=n; nf=n
      dnf=1.d0/dfloat(n)
      dnb=1.d0
      
      fdum = 0.d0
      bdum = 0.d0

      call dfftw_plan_dft_r2c_1d(planf,nf,fdum,fdum,FFTW_MEASURE )
      call dfftw_plan_dft_c2r_1d(planb,nb,bdum,bdum,FFTW_MEASURE )
      
      deallocate (fdum,bdum)

      end


c ================================================================
c             the real rft in y-direction
c ================================================================
      subroutine rfty(c,isa,m,iopt)
      use rftytmp
      implicit none
      integer isa,m,iopt,i,j
      real*8 c(isa,*)

      if (iopt.lt.0) then

         do j=1,m
            call dfftw_execute_dft_r2c(planf,c(1,j),c(1,j))
            do i=1,nf+2
               c(i,j)=dnf*c(i,j)
            enddo
         enddo

      else   

         do j=1,m
            call dfftw_execute_dft_c2r(planb,c(1,j),c(1,j))
            c(nb+1,j)=0.d0
            c(nb+2,j)=0.d0
         enddo

      endif

      end
