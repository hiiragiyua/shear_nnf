      subroutine ffti(m,n)
      use ffttmp
      use omp_lib

      implicit none
      include "fftw3.f"

      integer m,n,err

      allocate(fdum(m+2,n),bdum(m+2,n))
      fdum=0.d0; bdum=0.d0 
      norm=1.d0/dfloat(m*n)

      call dfftw_init_threads(err)
      call dfftw_plan_with_nthreads(OMP_GET_MAX_THREADS())
      call dfftw_plan_dft_r2c_2d(planf,m,n,fdum,fdum,FFTW_MEASURE)
      call dfftw_plan_dft_c2r_2d(planb,m,n,bdum,bdum,FFTW_MEASURE)
   
      deallocate(fdum,bdum)
      end subroutine ffti

      subroutine fourxz(fou,phys,tmpxzr,iopt,n1)
      use ctes,only: mgalx,mgalz,mx1,mz1
      use ffttmp
      use omp_lib

      implicit none
      include  "mpif.h"

      integer    iopt,n1
      complex*16 fou(0:mx1,0:mz1)
      real*8     phys(mgalx+2,mgalz)
      real*8     tmpxzr(mgalx+2,mgalz)
      
      real*8  pack_t,fill_t,rft_t,cft_t
      common/timers_fft/ pack_t,fill_t,rft_t,cft_t
      save/timers_fft/      
  
      real*8  t1,t2

      if (iopt.lt.0) then
         !$OMP PARALLEL WORKSHARE
         tmpxzr=phys
         !$OMP END PARALLEL WORKSHARE
         t1=MPI_WTIME()         
         call dfftw_execute_dft_r2c(planf,tmpxzr,tmpxzr)
         t2=MPI_WTIME()
         rft_t=rft_t+t2-t1
         !$OMP PARALLEL WORKSHARE
         tmpxzr=norm*tmpxzr
         !$OMP END PARALLEL WORKSHARE
         t1=MPI_WTIME()
         call pack00(tmpxzr,fou)
         t2=MPI_WTIME()
         pack_t=pack_t+t2-t1
      else
         t1=MPI_WTIME()
         call fill00(fou,tmpxzr)
         t2=MPI_WTIME()
         fill_t=fill_t+t2-t1
         t1=t2
         call dfftw_execute_dft_c2r(planb,tmpxzr,tmpxzr)
         t2=MPI_WTIME()
         cft_t=cft_t+t2-t1
         !$OMP PARALLEL WORKSHARE
         tmpxzr(mgalx+1:mgalx+2,:)=0.d0
         phys=tmpxzr
         !$OMP END PARALLEL WORKSHARE
      endif

      end subroutine fourxz

      subroutine fill00(fou,foud)
      use ctes
      use omp_lib

      implicit none

      complex*16 fou(0:mx1,0:mz1),foud(0:mgx,0:mgalz1),zero
      integer i,k,kk

      zero = dcmplx(0.d0,0.d0)
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,kk)

      !$OMP DO SCHEDULE(STATIC)
      do k=0,nz1
         do i=0,mx1
            foud(i,k)=fou(i,k)
         enddo
      enddo
      !$OMP END DO NOWAIT
      
      !!$OMP PARALLEL WORKSHARE
      !foud(0:mx1,0:nz1)      = fou(0:mx1,0:nz1)
      !foud(0:mx1,nz1+1:nz2-1)= zero
      !foud(0:mx1,nz2:mgalz1) = fou(0:mx1,nz2-mgalz1+mz1
      !&                              :mgalz1-mgalz1+mz1)
      !foud(mx1:mgx,:)        = zero
      !!$OMP END PARALLEL WORKSHARE 
        
      !$OMP DO SCHEDULE(STATIC)
      do k=nz1+1,nz2-1
         do i=0,mx1
            foud(i,k)=zero
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO SCHEDULE(STATIC)
      do k=nz2,mgalz1
         kk = k-mgalz1+mz1
         do i=0,mx1
            foud(i,k)=fou(i,kk)
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO SCHEDULE(STATIC)
      do k=0,mgalz1
         do i=mx1+1,mgx
            foud(i,k)=zero
         enddo
      enddo

      !$OMP END PARALLEL
      end

      subroutine pack00(foud,fou)
      use ctes
      use omp_lib

      implicit none

      complex*16 fou(0:mx1,0:mz1),foud(0:mgx,0:mgalz1)
      integer i,k,kk

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,kk)

      !$OMP DO SCHEDULE(STATIC)
      do k=0,nz1
         do i=0,mx1
            fou(i,k)=foud(i,k)
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO SCHEDULE(STATIC)
      do k=nz2,mgalz1
         kk = k-mgalz1+mz1
         do i=0,mx1
            fou(i,kk)=foud(i,k)
         enddo
      enddo    

      !$OMP END PARALLEL

      !!$OMP PARALLEL WORKSHARE
      !fou(:,0:nz1)=foud(0:mx1,0:nz1)
      !fou(:,nz2-mgalz1+mz1:mgalz1-mgalz1+mz1)=foud(0:mx1,nz2:mgalz1)
      !!$OMP END PARALLEL WORKSHARE
      end

