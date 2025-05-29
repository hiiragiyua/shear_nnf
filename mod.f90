module ctes
  ! ****************************************************************
  !       things which are fixed for the whole run (in cross)
  ! ****************************************************************

  ! ----- dimensions ---------
  integer mgalx,mgalz,my,mx,mz,mgalx1,mgalz1,mx1,my1,mz1,mgx,nz,nz1,nz2
  integer mgy,my23,ny,ny1,ny2

  ! ---- coordinates and parameters ---------
  real(8) Re,alp,gam,Ly,s,chi,dx,dz,dy,Lx,Lz,Lx2,Lz2 
  real(8) CFL,CFLv,cflx,cfly,cflz,cflvv,pmesp,uprim,umaxl(0:3),umax(0:3)
  real(8) cflkk,CFLt,CFLk,CFLb,tbmaxl,tbmax
  real(8) vbulk, xforce, zforce, force_roll
  real(8) damp_up,damp_aa

  ! for checking parameters (for newton-continuation parameters)
  real(8) Ree,alpe,game,Lye,se,chie
  real(8) CFLe, uprime, vbulke, xforcee, zforcee, force_rolle
  real(8) damp_upe,damp_aae


  integer  readflag
  integer,dimension(:),allocatable :: icx
  real(8),dimension(:),allocatable :: alp2, gam2, y
  complex(8),dimension(:),allocatable :: xalp, xgam

  !    --------  parallelization and blocking -------
  integer  numerop,myid,jb,je,kb,ke,mmy,mmz,blockingik2ki,blockingki2ik
  integer  mgalzp,myp,mzp,mgalz1p,nxymax,buffsize
  integer,dimension(:),allocatable :: jbeg,jend,kbeg,kend
  !integer, allocatable:: scount(:),sdispl(:),rcount(:),rdispl(:) 
  ! Note: these array affects myslice (rev.400) 

  integer,dimension(:),allocatable :: myslice

  !  --------  junk -----
  real(8) pi, pi2
  complex(8) cii
  real(8) tiny,big
  parameter(tiny = 1.d-13)
  parameter(big = 1.d10) 
  
  integer irev,nopt,nparams
  parameter(irev=1742)  ! revision information. 
  integer ireve

  parameter(maxParams=100)  ! the max. number of options
  !parameter(nopt=9) ! the number of switch options 
  !parameter(nparams=13) ! the number of input parameters

end module ctes

module running
! ****************************************************************
!       things that used to be running commons
! ****************************************************************

  real(8)  Deltat,time, etime, time_ini, dtimag, dumptime, timep_sh
  real(8)  dtr, dump2timep, dump2tint
  real(8)  dtrv, dtrk ! itemperature
  real(8) sym_shiftx,sym_shiftz

  integer nimag,nstep,nhist,ihist,icfl,ifix_dt,iwrote,idump_mode,iydealiasing
  integer iget_cfy ! to dump cfy 
  integer istart,istep,iend,ifatal,istop,idump,noescru,nohist,nocf,iskip_screenout
  ! some option for post processing
  integer explicit, iadd_force, iadd_mode, iadd_sym, iadd_damping, itemperature, init_temp
  integer explicite, iadd_forcee, iadd_modee, iadd_syme, iadd_dampinge, itemperaturee
  integer iuse_newton,iread_footer, inorm
  integer iread_hdf5, iwrite_hdf5
  !  ------- file things ----
  integer iinp,iout,id22,isn,ispf,ifile,iocf,iocfb,iocfy,iohre !31 -- 40,19
  integer iotfb
  character*128 filinp,filinm1,filout,filstt
  character*4 filext  ! from rev.1153

  character*128,dimension(:),allocatable :: hrefile
  integer hrelines,ihre_error

end module running

module timer
  ! timer
  real(8) commtimer, transtimer, totaltimer
  real(8) iter_time, comm_time, wtime
  real(8) addsymtimer
  ! 
  real(8) t_derivyc, t_add_shear, t_uvw, t_copy, t_fourxz, &
        t_uome, t_nl, t_hvhg, t1, t2, t_hvhg0, t_otra, t_waiting  
 
end module timer

module statistics
! ****************************************************************
!                statistics and diagnostics
! ****************************************************************
  !  ---- spectra -------
  integer  nacumsp,nspec,nspecy,nspect
  !parameter (nspec=8, nspect=4)
  real(8),dimension(:,:,:),allocatable :: sp, spl
  real(8),dimension(:,:,:),allocatable :: sp2, spl2 ! for sta4
  real(8),dimension(:,:,:),allocatable :: spt, splt ! for sta4
  !  ----- one-point statistics ---
  real(8),dimension(:,:),allocatable :: stats, stats2
  real(8),dimension(:,:),allocatable :: tstats
  integer nstat,ntstat,nacum,nner,nner2,ntner
  !parameter (nstat = 23, ntstat = 5)

  parameter (nner=10,nner2=4)
  real(8) ener(nner),ener2(nner2),stime0, total_dt
  parameter (ntner=4)
  real(8) enert(ntner)

end module statistics

module bcs
! ****************************************************************
!                boundary condistions in y direcition
!                (shear-periodic information)
! ****************************************************************
  real(8) xoffb,xofft, xwkb,xwkt
  real(8) zoffb,zofft, zwkb,zwkt

  real(8) u00b,u00t, w00b,w00t
  real*8 xoffb_ini, xofft_ini  ! move from the module 'gmres'

  complex(8),dimension(:),allocatable :: shb, sht, shwkb, shwkt
!  complex(8),dimension(:,:),allocatable :: shb, sht, shwkb, shwkt
end module bcs

module rfttmp
! ****************************************************************
!   fftw3.3.2 tmp arrays and plans for rft 
! ****************************************************************  

  real*8,allocatable:: fdum(:),bdum(:)
  real*8  dnf,dnb
  integer*4  nf,nb
  integer*8  planf,planb
  
end module rfttmp

module cfttmp
! ****************************************************************
!   fftw3.3.2 tmp arrays and plans for cft 
! 
! ****************************************************************  

  complex*16,allocatable:: fdum(:),bdum(:)
  real*8     dnf,dnb 
  integer  nf,nb
  integer*8  planf,planb

end module cfttmp

module ffttmp
  ! this is for threaded 2D fft (fou3D_thread.f) 
  complex*16,allocatable:: fdum(:,:),bdum(:,:)
  real*8     norm
  integer*8  planf,planb

end module ffttmp



! ****************************************************************
!   fftw3.3.2 tmp arrays and plans for cft in y-dir
! 
! ****************************************************************  
module rftytmp

  real*8,allocatable:: fdum(:),bdum(:)
  real*8  dnf,dnb
  integer*4  nf,nb
  integer*8  planf,planb
  
end module rftytmp

module cftytmp

  complex*16,allocatable:: fdum(:),bdum(:)
  real*8     dnf,dnb 
  integer  nf,nb
  integer*8  planf,planb

end module cftytmp

module addsym
  ! ***
  ! for adding symmetry
  ! ***
  complex*16,dimension(:,:,:),allocatable:: vors,phis
  real*8,dimension(:),allocatable :: u00s,w00s
  integer,dimension(:),allocatable :: sgcount,rgcount,rgdispl 
  
end module addsym

module LES
  integer     iuse_LES,idynamic,iadd_visc,ifix_CsDeltag
  real*8      Cles,Csy,Deltag,max_nut,nutmax,cutoffx,cutoffz,LESflt(3)
  real*8      Clese,Deltage,cutoffxe,cutoffze,LESflte(3)
  real*8      CsDeltag_fix
  ! tmp 2d plane array
  real*8,     allocatable:: ux(:,:),uy(:,:),uz(:,:),vx(:,:),vy(:,:), &
                            vz(:,:),wx(:,:),wy(:,:),wz(:,:)
  complex*16, allocatable:: uxc(:,:),uyc(:,:),uzc(:,:),vxc(:,:),     &
                            vyc(:,:),vzc(:,:),wxc(:,:),wyc(:,:),wzc(:,:)
  real*8,     allocatable:: u11r(:,:),u22r(:,:),u33r(:,:),u12r(:,:), &
                            u13r(:,:),u23r(:,:)
  complex*16, allocatable:: u11c(:,:),u22c(:,:),u33c(:,:),u12c(:,:), &
                            u13c(:,:),u23c(:,:)
  real*8,     allocatable:: nuSGS(:,:),Sij2(:,:),Tll(:,:),Cast(:,:)
  integer,    allocatable:: cutoff(:,:)

  ! SGS is mut.
  real*8,dimension(:,:,:),allocatable:: nut
  real*8,dimension(:),allocatable:: txy00,tyz00
  ! buffsize arrays for change...
  !real*8,dimension(:),allocatable:: rhvc1, rhvc2, rhgc1, rhgc2  

end module LES

module temp
  !integer itemperature
  real*8      fkappa,  bym,  gbeta
  !real*8,dimension(:,:,:)
  real*8,     allocatable:: uTr(:,:),vTr(:,:),wTr(:,:)
  complex*16, allocatable:: uTc(:,:),vTc(:,:),wTc(:,:)
end module temp


! modules for error-free-transformation...

module eft
  ! double-double precision
  integer, parameter :: nk=8 
  real(nk), parameter :: bfactor=2**27+1 ! 64 bit
  !real(nk), parameter :: bfactor=2.d0**32.d0+1.d0 ! extende 80 bit
contains

  subroutine TwoSum(a,b,x,y)
    implicit none
    real(nk) a,b,c,d,x,y
    c=0d0; d=a; ! in case that a and x share the memory
    x=a+b; c=x-a; 
    y=(d-(x-c)) + (b-c);
  end subroutine TwoSum

  subroutine Split(a,aH,aL)
    implicit none
    real(nk) a,aH,aL,c
    aH=0d0; aL=0d0; c=0.d0;
    !c=bfactor*a;
    c=(2.d0**27.d0+1.d0)*a;
    aH=c-(c-a); aL=a-aH;
  end subroutine Split

  subroutine TwoProduct(a,b,x,y)
    implicit none
    real(nk) a,b,aH,aL,bH,bL,x,y
    x=a*b;
    call Split(a,aH,aL); call Split(b,bH,bL);
    y=aL*bL-(((x-aH*bH)-aL*bH)-aH*bL);
  end subroutine TwoProduct

  subroutine VecSum(n,p) 
    ! p(n) is the sum, not self-preservative
    implicit none 
    integer i,n
    real(nk) p(n)
    do i=2,n
       call TwoSum(p(i),p(i-1),p(i),p(i-1))
    end do
  end subroutine VecSum

  subroutine SumK(n,p,K,res)
    ! note; p(n) is the numerical sum
    ! res(ult): is the error-free transformation
    implicit none
    integer i,n,K,l
    real(nk) p(n)
    real(nk) res

    do l=1,K-1
       call VecSum(n,p)
    end do
    
    res=0.d0
    do i=1,n-1
       res = res + p(i)
    end do    
    res = p(n)+res
    
  end subroutine SumK

  subroutine DotK(n,x,y,res)
    implicit none
    integer i,n,K,nn
    real(nk) x(n),y(n)
    real(nk) res,p,h
    real(nk), dimension(:), allocatable :: r
    K=3
    nn=2*n
    allocate(r(nn))
    r(:)=0d0
    p=0.d0; h=0.d0
    call TwoProduct(x(1),y(1),p,r(1))
    !write(*,*) 'TwoProd',p,r(1)
    do i=2,n
       call TwoProduct(x(i),y(i),h,r(i))
       call TwoSum(p,h,p,r(n+i-1))
    end do
    r(nn)=p;
    call SumK(nn,r,K-1,res)
    deallocate(r)

  end subroutine DotK

    ! for complex, see Graillat, et al. (2012)
  subroutine Sum2(n,p,res)
    implicit none
    integer i,n,K,l
    real(nk) p(n)
    real(nk), allocatable :: q(:),s(:),tmp(:)
    real(nk) res
    
    !allocate(q(n),s(n),tmp(n))
    !tmp(1)=p(1); s(1)=0.d0;
    !do i=2,n
    !   call TwoSum(tmp(i-1),p(i),tmp(i),q(i))
    !   s(i)=s(i-1)+q(i)
    !end do
    !res=tmp(n)+s(n); 
    !deallocate(q,s,tmp)
    ! this is VecSum
    do i=2,n
       call TwoSum(p(i),p(i-1),p(i),p(i-1))
    end do
    res=0.d0
    do i=1,n-1
       res = res + p(i)
    end do
    res = p(n)+res

  end subroutine Sum2

  subroutine Sum2cmplx(n,p,cres)
    implicit none
    integer n
    complex(nk) p(n)
    real(nk),allocatable :: pr(:),pm(:)
    real(nk) res(2)
    complex(nk) cres
    allocate(pr(n),pm(n))
    res=0.d0;
    pr=dreal(p(:)); pm=dimag(p(:))
    call Sum2(n,pr,res(1)) ! real part
    call Sum2(n,pm,res(2)) ! imaginary part
    cres=dcmplx(res(1),res(2))
    deallocate(pr,pm)
  end subroutine Sum2cmplx

  subroutine Dot2(n,x,y,res)
    implicit none
    integer i,n
    real(nk) x(n),y(n)
    real(nk) res,p,q,s,h,r
    call TwoProduct(x(1),y(1),p,s)
    do i=2,n
       call TwoProduct(x(i),y(i),h,r)
       call TwoSum(p,h,p,q)
       s=(s+(q+r))
    end do
    res=(p+s)
  end subroutine Dot2
  
  subroutine Dot2cmplx(n,x,y,cres)
    ! complex dot product: cres=conj(x)*y
    implicit none
    integer i,n
    !real(nk) x(2,n),y(2,n)
    complex(nk) x(n),y(n)
    complex(nk) cres
    real(nk) res(2),res1,res2,res3,res4
    real(nk), dimension(:), allocatable :: a,b,c,d

    allocate(a(n),b(n),c(n),d(n))
    a=dreal(x(:)); b=dimag(x(:)); 
    c=dreal(y(:)); d=dimag(y(:))
    
    call Dot2(n,a,c,res1)
    call Dot2(n,b,d,res2)

    call Dot2(n,a,d,res3)
    call Dot2(n,b,-c,res4)

    deallocate(a,b,c,d)

    cres=dcmplx(res1+res2, res3+res4)

  end subroutine Dot2cmplx

end module eft
