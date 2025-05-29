! ====================================================================
!
!                   get data field from file
!                           
! ====================================================================
subroutine get_header()

  use ctes
  use running
  use bcs
  !use LES, only: iuse_LES

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,irec,recl
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  integer ihalf,iskip,ihalfy,ncopy

  real*8, allocatable:: ye(:)

  integer  mgalxe,mye,mgalze,numerope
  real*8  tiempo,ubulk,wbulk


  integer nopte,nparamse ! rev.1433
  integer ioptions(maxParams) ! 100 is maximum ... see in mod.f90
  real(8) params(maxParams)

  ! reading header
  if (myid.eq.0) then     ! ----- this is done by the master

     open(iinp,file=filinp,status='old',form='unformatted',action='read')

     read(iinp) mgalxe,mye,mgalze,numerope
     read(iinp) tiempo,Ree,alpe,Lye,game,se,chie,me ! me=mx,my,mz

     Lye=Lye*pi
     !write(*,*) 'read Ree',Ree ! do not forget to BCAST
     !allocate(ye(1:me(2))) ! myid = 0, it starts from ye(1), but y(0)
     !read(iinp) ye(1:me(2))
     !read(iinp) xoffb,xofft,zoffb,zofft
     
     close(iinp)
     !deallocate(ye)
     ! write(*,*) 'headers:',mgalxe, Ree, me, xoffb
  end if

  mgalx=mgalxe; my=mye; mgalz=mgalze

  call MPI_BCAST(mgalx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(my,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mgalz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(numerope,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     
  call MPI_BCAST(Ree,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Lye,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(game,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(se,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(chie,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  write(*,*) 'CFL is set as 0.6, in get_header ...'
  cfl=0.6 ! 

end subroutine get_header


subroutine getfil(vor,phi,u00,w00,wk,tiempo)

  use ctes
  use running
  use bcs
  use LES, only: iuse_LES

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  integer ihalf,iskip,ihalfy,ncopy, icheck_write

  real(8),dimension(2,mx1+1,mz,jb:je) :: vor, phi
  real(8),dimension(buffsize) :: wk
  real(8),dimension(my) :: u00, w00

  real*8, allocatable:: wk00(:,:), wk2(:,:,:,:), ye(:)

  integer  mgalxe,mye,mgalze,numerope
  real*8  tiempo,ubulk,wbulk

  integer nopte,nparamse ! rev.1433
  integer ioptions(maxParams) ! 100 is maximum ... see in mod.f90
  real(8) params(maxParams)

  !      ---------------  zero everything
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  wk(buffsize)=0.d0

  icheck_write=0 ! to skip writing on screen

  if (myid.eq.0) then     ! ----- this is done by the master

     open(iinp,file=filinp,status='unknown',form='unformatted',action='read')

     read(iinp) mgalxe,mye,mgalze,numerope

     if((mgalx.ne.mgalxe).or.(my.ne.mye).or.(mgalz.ne.mgalze)) then
        
        write(*,*) ' different physical grid points!'
        write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,mye,mgalze, &
             ' =>', mgalx, my, mgalz 
        
        if ((my.ne.mye).and.(mgalx.ne.mgalxe)) then
           write(*,*) ' Do not interpolate in my at the same time...'
           stop
        end if
        !stop
     endif

     read(iinp) tiempo,Ree,alpe,Lye,game,se,chie,me ! me=mx,my,mz

     Lye=Lye*pi
     !write(*,*) 'read Ree',Ree ! do not forget to BCAST
     if (abs(Lye-Ly).gt.tiny) write(*,*)' getfil: Ly is different !!',Lye, &
          ' ==>', Ly
     if (abs(se-s).gt.tiny) then 
        if (abs(se).lt.tiny) then
           write(*,*)' Now the shear is started'
        elseif (abs(s).lt.tiny) then
           write(*,*)' Now the shear is removed'
        else
           write(*,*)' getfil: s is different, STOP !!',se
           !stop ! rev.1426 for reverse run.
        end if
     endif
     allocate(ye(1:me(2))) ! myid = 0, it starts from ye(1), but y(0)
     read(iinp) ye(1:me(2))
     read(iinp) xoffb,xofft,zoffb,zofft
     if (icheck_write.eq.1) then
        write(*,*) 'reading input file ..., xoffb =',xoffb 
        write(*,*) '                   ..., xofft =',xofft
     endif
     if (abs(Lye + ye(1)*2.d0).gt.tiny) then
        write(*,*)' getfil: Lye is different from ye(1), reading data is broken, STOP!!',Lye,ye(1)
        stop
     end if
     !rescaling offset by shear
     ihalf=0; ihalfy=0;
     if (alpe.ne.alp) then
   
        if ((int((alp+0.00001)/alpe).eq.2).and.(mgalxe/mgalx.eq.2)) then
           xoffb=xoffb*alpe/alp*2.d0
           xofft=xofft*alpe/alp*2.d0
           write(*,*) ' '
           write(*,*) ' make the box to be half using only even modes...'
           write(*,*) '   keep the offset xoffb,xofft = ',xoffb,xofft
           ihalf=1
        else
           xoffb=xoffb*alpe/alp
           xofft=xofft*alpe/alp
           write(*,*) ' rescaling shear-periodic b.c. offset: factor =', &
                alpe/alp, '  xoffb,xofft = ',xoffb,xofft
        endif

        if (me(2).ne.my) then
           write(*,*) ' Do not interpolate in my at the same time...'
           stop
        end if

     end if
 
     if (icheck_write.eq.1) then
        write(*,*)
        write(*,*) 'reading input file ..., time=',tiempo
        write(*,'(a,3(1X,I5))') '                    ..., (mx,my,mz)=',me
        write(*,*)
     endif
 
     do iproc=1,numerop-1
        call MPI_SEND(me,3,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ye   ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(tiempo,1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr) 

        call MPI_SEND(ihalf,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ihalfy,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
     enddo

  else       ! -------- this is done by all slaves

     call MPI_RECV(me,3,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     allocate(ye(1:me(2))) ! slaves allocate, 
     call MPI_RECV(ye  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(tiempo,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

     call MPI_RECV(ihalf,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(ihalfy,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  endif


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !write(*,*) 'allocate wk00', me(2), my
  allocate( wk00(me(2),2) )
  !     ------------- 00 modes ------------------
  if (myid.eq.0) then     ! ----- this is done by the master

     if ((me(2).lt.my).and.(abs(Ly-Lye).lt.tiny)) then
        if(myid.eq.0) then
           write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' => my=',my
        end if
     end if
     if (my.eq.me(2)/2) then
        ihalfy=1
        if (myid.eq.0) write(*,*) 'reducing my to be half '
     end if

     !write(*,*) 'allocated wk00', me(2), my
     wk00 = 0.d0
     !write(*,*) 'reading 00 modes'
     read(iinp) wk00
     mym = min(my,me(2))
     if (ihalfy.eq.1) then
        u00(1:my) = wk00(1:me(2):2,1)
        w00(1:my) = wk00(1:me(2):2,2)
        !write(*,*) 'read 00 mode'
     else
        u00(1:mym) = wk00(1:mym,1)
        w00(1:mym) = wk00(1:mym,2)
     end if
     do iproc=1,numerop-1
        call MPI_SEND(u00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
     end do
  else       ! -------- this is done by all slaves
     call MPI_RECV(u00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(w00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  endif
  deallocate (wk00)


  ! ---- the other modes ----
  ntotr=me(1)*me(3)*2
  allocate ( wk2(2,2,me(1)/2,me(3)) )
  wk2 = 0.d0
  if (ihalf.eq.1) then
     iskip=2
  else
     iskip=1
  endif

  mxm=min(mx,me(1))
  mxm=mxm/2 ! ??
  mym=min(my,me(2))
  mzm=min(mz,me(3))

  klen=(mzm+1)/2
  kini1=me(3) - klen + 1
  kini2=mz    - klen + 1
  
  ! first expand the fourier mode in x-z and then use cubic-spilne to interpolate in y
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then     ! ----- this is done by the master, first vor

     !write(*,*) 'reading other modes'
    
     do iproc=0,numerop-1   
        do j=jbeg(iproc),jend(iproc)
           if (j.le.(me(2)-1)) then ! the case for my is larger than mye
              read(iinp) wk2
           else
              wk2=0.d0
           endif

           if (iproc.eq.0) then  ! master stores its thing
              vor(:,1:mxm,1:klen,j)    =wk2(1,:,1:mxm*iskip:iskip,1:klen)
              vor(:,1:mxm,1+kini2:mz,j)=wk2(1,:,1:mxm*iskip:iskip,1+kini1:me(3))
              phi(:,1:mxm,1:klen,j)    =wk2(2,:,1:mxm*iskip:iskip,1:klen)
              phi(:,1:mxm,1+kini2:mz,j)=wk2(2,:,1:mxm*iskip:iskip,1+kini1:me(3))
           else      ! sends everything else to nodes
              call MPI_SEND(wk2,ntotr,MPI_REAL8,iproc,iproc, & 
                   MPI_COMM_WORLD,ierr)
           endif

           if (ihalfy.eq.1) then
              ! the case for my being half
              read(iinp) wk2
              write(*,*) 'skip one plane by just reading '              
           endif

         enddo
     enddo
     
     if (iread_footer.eq.1) then ! default is iread_footer=1

        ! !!! read some extra parameters (testing)
        ! write(iout) CFL, uprim, vbulk, xforce ! rev934
        ! !!!write(iout) irev,explicit,iadd_force,iadd_mode,iadd_sym ! rev1153  
        ! write(iout) irev,explicit,iadd_force,iadd_mode,iadd_sym,iadd_damping ! rev1200      
        if (iread_hdf5.eq.1) then
           
        else
           read(iinp,err=980) CFLe, uprime, vbulke, xforcee ! real*8
           read(iinp,err=980) ireve,nopte,nparamse
           read(iinp,err=980) ioptions(1:nopte)
           read(iinp,err=980) params(1:nparamse)
        end if

        if (ireve.lt.1200) write(*,*) 'check the revision, the file has ireve=',ireve
        
        if (mye.ne.my) then 
           write(*,*) 'check reading extra options for my.ne.mye, STOP!, ireve=',ireve
           write(*,*) 'check reading , CFLe, etc. =',CFLe, uprime,vbulke,xforcee
           !stop
        end if

     else if (iread_footer.eq.0) then
        
        if (icheck_write.eq.1) write(*,*) 'iread_footer =',iread_footer 
        nopte=0; nparamse=0; ireve=0;

     end if

!980  write(*,*) myid,'end of file read'     
980     close(iinp)

  else          ! -------- the slaves receive higher modes ----------

     do j=jb,je
        if ((ihalfy.eq.1).and.(mod(j+1,2).eq.0)) then
           ! Skip reading
           continue
        endif
        call MPI_RECV(wk2,ntotr,MPI_REAL8,0,MPI_ANY_TAG, & 
             MPI_COMM_WORLD,istat,ierr)
        vor(:,1:mxm,1:klen,j)    =wk2(1,:,1:mxm*iskip:iskip,1:klen)
        vor(:,1:mxm,1+kini2:mz,j)=wk2(1,:,1:mxm*iskip:iskip,1+kini1:me(3))
        phi(:,1:mxm,1:klen,j)    =wk2(2,:,1:mxm*iskip:iskip,1:klen)
        phi(:,1:mxm,1+kini2:mz,j)=wk2(2,:,1:mxm*iskip:iskip,1+kini1:me(3)) 
     enddo

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! reading vb vt is removed (rev260, 2012/Feb/02)
  deallocate (wk2)
  call  MPI_BCAST(Ree,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(Lye,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! BCAST are forgotten (2013/10/14) alpe,Lye,game,se,chie,me
  call  MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(game,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(se,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(chie,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  !call  MPI_BCAST(CFLe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(uprime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(vbulke,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(xforcee,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(ireve,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nopte,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nparamse,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ioptions,nopte,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(params,nparamse,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  call check_parameters(nopte,ioptions,nparamse,params)
  ! set conjugate relation for kx=0
  !if(myid.eq.0) write(*,*) myid,'set_conj' 
  call set_conj(vor,phi,-1)
  ! ----- before going any further exchange cuts in y with cuts in z!
  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)

  ! ------ you should not interpolate from a small number of grids to ------
  ! ------- a much larger, interpolation itself has error --------
  ncopy=0
  if ((me(2).eq.my/3)) then
      ncopy=3
  elseif ((me(2).eq.my/5)) then
      ncopy=5
  endif

  if ((me(2).lt.my).and.(ncopy.eq.0)) then
     ! !note not to interpolated with different Ly.
     if (abs(Lye-Ly).lt.tiny) then
        write(*,*) myid,': interpolating in y-dir'
        call interp_y(vor,phi,u00,w00,me,ye)
     else
        if (myid.eq.0) write(*,*) ' Do not interpolate in y-dir, with different Ly from Lye in hre2.dat'
        stop
     end if

  elseif ((me(2).lt.my).and.(ncopy.gt.0)) then
     write(*,*) myid,'getfil, ncopy=',ncopy 
     call ncopy_y(vor,phi,u00,w00,me,ncopy)

  endif 

  deallocate(ye) ! deallocated by all proc.
  ! cleanup work array for change (to avoid an un-expected bug, rev(635))
  wk(buffsize)=0.d0

  ! set bulk velocity to be zero
  ubulk=sum(u00(1:my))/dfloat(my)
  wbulk=sum(w00(1:my))/dfloat(my)
  if ((myid.eq.0).and.(icheck_write.eq.1)) then 
     write(*,*) 'check ubulk, wbulk=', ubulk, wbulk
  end if
  !!u00(:)=u00(:)-ubulk
  !!w00(:)=w00(:)-wbulk

  xwkb=xoffb; xwkt=xofft; ! bug fixed, 2014/June/24 sekimoto
  zwkb=zoffb; zwkt=zofft;


  if (abs(game-gam).gt.tiny) then ! rev1140
     write(*,*) 'rescale velocity ... and lap.v: game/gam=',game/gam
     u00=u00*(game/gam)
     w00=w00*(game/gam)
     phi=phi*(gam/game) ! (game/gam*(gam/game)^2)
     if (iuse_LES.eq.0) then
        vor=vor*(game/gam)*sqrt(Re/Ree) ! rev1157
     end if
  elseif (abs(Ree-Re).gt.tiny) then ! rev1157
     if ((iuse_LES.eq.0).and.(iuse_newton.eq.0)) then
        write(*,*) 'rescale vorticity  sqrt(Re/Ree)',sqrt(Re/Ree)
        vor=vor*(game/gam)*sqrt(Re/Ree) ! rev1157
     end if
  end if

end subroutine getfil

subroutine ncopy_y(vor,phi,u00,w00,me,ncopy)
  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr  
  integer i,j,k,ncopy,inc,jst,jen
  integer me(3)
  complex*16 shp
  real*8 ye(my)
  !real*8 phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  complex*16 phi(1:my,0:mx1,kb:ke),vor(1:my,0:mx1,kb:ke)
  real*8 u00(1:my),w00(1:my)

  ! ---------------- note: my>me(2) ---------------------
  do k=kb,ke
     do i=0,mx1

        do inc=1,ncopy-1
           shp = zexp(-xalp(i)*xofft*inc) ! positive shift
           jst=1+me(2)*inc
           jen=me(2)+me(2)*inc
           vor(jst:jen,i,k)  = vor(1:me(2),i,k)*shp
           phi(jst:jen,i,k)  = phi(1:me(2),i,k)*shp
        
        enddo
     enddo
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! update boundary condition
  xofft=xofft*dfloat(ncopy)
  xoffb=xoffb*dfloat(ncopy)

  ! --------------------- 00 mode -----------------------
  do inc=1,ncopy-1
     jst=1+me(2)*inc
     jen=me(2)+me(2)*inc
     u00(jst:jen) = u00(1:me(2))
     w00(jst:jen) = w00(1:me(2))
  enddo

end subroutine ncopy_y


subroutine ncopy_y_t(tb,tb00,me,ncopy)
  use ctes
  use running
  use bcs
  
  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr  
  integer i,j,k,ncopy,inc,jst,jen
  integer me(3)
  complex*16 shp
  real*8 ye(my)
  complex*16 tb(1:my,0:mx1,kb:ke)
  real*8 tb00(1:my)

  ! ---------------- note: my>me(2) ---------------------
  do k=kb,ke
     do i=0,mx1

        do inc=1,ncopy-1
           shp = zexp(-xalp(i)*xofft*inc) ! positive shift
           jst=1+me(2)*inc
           jen=me(2)+me(2)*inc
           tb(jst:jen,i,k)  = tb(1:me(2),i,k)*shp
        
        enddo
     enddo
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! update boundary condition
  xofft=xofft*dfloat(ncopy)
  xoffb=xoffb*dfloat(ncopy)

  ! --------------------- 00 mode -----------------------
  do inc=1,ncopy-1
     jst=1+me(2)*inc
     jen=me(2)+me(2)*inc
     tb00(jst:jen) = tb00(1:me(2))
  enddo

end subroutine ncopy_y_t

!===================real*4==>real*8=======================!
subroutine getfil_convert(vor,phi,u00,w00,wk,tiempo)

  use ctes
  use running
  use bcs
  use LES, only: iuse_LES

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2

  real*8 vor(mx,mz,jb:je),phi(mx,mz,jb:je)
  real*8 wk(*)
  real*8 u00(my),w00(my)
  real*8, allocatable:: wk1(:,:)
  real*4, allocatable:: wk2(:,:,:)

  integer  mgalxe,mye,mgalze,numerope
  !real*8  Ree,alpe,Lye,game,se,chie
  real*8  ye(my)
  real*8  tiempo

  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0

  if (myid.eq.0) then    

     open(iinp,file=filinp,status='unknown',form='unformatted',action='read')

     read(iinp) mgalxe,mye,mgalze,numerope

     if((mgalx.ne.mgalxe).or.(my.ne.mye).or.(mgalz.ne.mgalze)) then
        if(myid.eq.0) then
           write(*,*) ' different physical grid points!'
           write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,mye,mgalze, &
                                                     ' =>', mgalx, my, mgalz 
        endif
     endif

     read(iinp) tiempo,Ree,alpe,Lye,game,se,chie,me 
     read(iinp) ye(1:me(2))
     read(iinp) xoffb,xofft,zoffb,zofft
     
     write(*,*) 'reading input file ..., xoffb =',xoffb 
     write(*,*) '                   ..., xofft =',xofft
     
    !------rescaling offset by shear------
     if (alpe.ne.alp) then
        xoffb=xoffb*alpe/alp
        xofft=xofft*alpe/alp
        write(*,*) ' rescaling shear-periodic b.c. offset: factor =', &
                     alpe/alp, '  xoffb,xofft = ',xoffb,xofft
     endif
     write(*,*)
     write(*,*) 'reading input file ..., time=',tiempo
     write(*,'(a,3(1X,I5))') '                    ..., (mx,my,mz)=',me
     write(*,*)
   
     !------------- 00 modes ---------------
     if (me(2).ne.my) then
        if(myid.eq.0) then
           write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' => my=',my
        end if
     end if
     allocate(wk1(me(2),2))
     wk1 = 0.d0
     read(iinp) wk1
     mym = min(my,me(2))
     u00(1:mym) = wk1(1:mym,1)
     w00(1:mym) = wk1(1:mym,2)
     deallocate (wk1)
     do iproc=1,numerop-1
        call MPI_SEND(me,3,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(y    ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ye   ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(u00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
     enddo

  else      

     call MPI_RECV(me,3,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(y    ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(ye   ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(u00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(w00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

  endif

  ntotr=me(1)*me(3)*2
  allocate ( wk2(2,me(1),me(3)) )
  wk2=0.0
  mxm=min(mx,me(1))
  mym=min(my,me(2))
  mzm=min(mz,me(3))

  klen=(mzm+1)/2
  kini1=me(3) - klen + 1
  kini2=mz    - klen + 1
  ! first expand the fourier mode in x-z and then use cubic-spilne to
  ! interpolate in y
  if (myid.eq.0) then 

     do iproc=0,numerop-1   
        do j=jbeg(iproc),jend(iproc)
           if (j.le.(me(2)-1)) then
              read(iinp) wk2
           else
              wk2=0.0
           endif
          
           if (iproc.eq.0) then
              vor(1:mxm,1:klen,j)    =real(wk2(1,1:mxm,1:klen),kind=8)
              vor(1:mxm,1+kini2:mz,j)=real(wk2(1,1:mxm,1+kini1:me(3)),kind=8)
              phi(1:mxm,1:klen,j)    =real(wk2(2,1:mxm,1:klen),kind=8)
              phi(1:mxm,1+kini2:mz,j)=real(wk2(2,1:mxm,1+kini1:me(3)),kind=8)
           else      
              call MPI_SEND(wk2,ntotr,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,ierr)
           endif
         enddo
     enddo

     close(iinp)

  else         
     do j=jb,je
        call MPI_RECV(wk2,ntotr,MPI_REAL4,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
        vor(1:mxm,1:klen,j)    =real(wk2(1,1:mxm,1:klen),kind=8)
        vor(1:mxm,1+kini2:mz,j)=real(wk2(1,1:mxm,1+kini1:me(3)),kind=8)     
        phi(1:mxm,1:klen,j)    =real(wk2(2,1:mxm,1:klen),kind=8)
        phi(1:mxm,1+kini2:mz,j)=real(wk2(2,1:mxm,1+kini1:me(3)),kind=8)     
     enddo

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  deallocate (wk2)
  
  ! ----- before going any further exchange cuts in y with cuts in z!
  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)

  if (me(2).lt.my) then
     call interp_y(vor,phi,u00,w00,me,ye)
  elseif (me(2).gt.my) then
     if (myid.eq.0) write(*,*) 'Sorry, it does not support my less than me(2),:('
     stop
  endif

end subroutine getfil_convert



subroutine interp_y(vor,phi,u00,w00,me,ye)
  use ctes
  use running
  use bcs

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr  
  integer i,j,k
  integer me(3)
  complex*16 shp
  real*8 ye(my)
  real*8 phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  real*8 u00(0:my1),w00(0:my1)

  ! --------- wk3 and wk4 for second derivatives ---------- !
  ! -------- wk1 and wk2 for original phi and vor --------- !
  real*8, allocatable:: wk1(:),wk2(:)
  real*8, allocatable:: wk3(:),wk4(:)
  allocate(wk1(0:2*me(2)+1),wk2(0:2*me(2)+1))
  allocate(wk3(0:2*me(2)+1),wk4(0:2*me(2)+1))
  wk1=0.d0;wk2=0.d0;wk3=0.d0;wk4=0.d0

  ye(me(2)+1)=Ly/2.d0

  ! ---------------- note: my>me(2) ---------------------
  do k=kb,ke
     do i=0,mx1
        shp = zexp(-xalp(i)*xofft) ! positive shift

        wk1(0:2*me(2)-1)=vor(0:2*me(2)-1,i,k)
        
        wk1(2*me(2))  =dreal(dcmplx(vor(0,i,k),vor(1,i,k))*shp)
        wk1(2*me(2)+1)=dimag(dcmplx(vor(0,i,k),vor(1,i,k))*shp)

        wk2(0:2*me(2)-1)=phi(0:2*me(2)-1,i,k)
        wk2(2*me(2))  =dreal(dcmplx(phi(0,i,k),phi(1,i,k))*shp)
        wk2(2*me(2)+1)=dimag(dcmplx(phi(0,i,k),phi(1,i,k))*shp)   
        call cubic_spline(ye(1:me(2)+1),wk1,me(2)+1,wk3,2)  
        call cubic_spline(ye(1:me(2)+1),wk2,me(2)+1,wk4,2)
        do j=0,my1
           call interpolate(y(j),ye(1:me(2)+1),wk1,wk3,me(2)+1,vor(2*j:2*j+1,i,k),2)
           call interpolate(y(j),ye(1:me(2)+1),wk2,wk4,me(2)+1,phi(2*j:2*j+1,i,k),2)
        enddo
     enddo
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! ---------------------00 mode -----------------------
  wk1(0:me(2)-1)=u00(0:me(2)-1); wk1(me(2))=u00(0)
  wk2(0:me(2)-1)=w00(0:me(2)-1); wk2(me(2))=w00(0)
  call cubic_spline(ye(1:me(2)+1),wk1(0:me(2)),me(2)+1,wk3(0:me(2)),1)  
  call cubic_spline(ye(1:me(2)+1),wk2(0:me(2)),me(2)+1,wk4(0:me(2)),1)
  do j=0,my1
     call interpolate(y(j),ye(1:me(2)+1),wk1(0:me(2)),wk3(0:me(2)),me(2)+1,u00(j),1)
     call interpolate(y(j),ye(1:me(2)+1),wk2(0:me(2)),wk4(0:me(2)),me(2)+1,w00(j),1)
     ! BUG fixed: 2012/July/13, (wk3 -> wk4) by Siwei
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! -after this we get the vor and phi interpolated in y-
  deallocate(wk1,wk2,wk3,wk4)

end subroutine interp_y


! -------------- numerical recipes ---------------
subroutine cubic_spline(x,y,N,dy2,m)
  integer N,i,j,k,m
  real *8 x(N),y(m,N),dy2(m,N),u(m,N)
  real *8 p,sig
  if (m.eq.1) then
     dy2(1,1)=0.d0
     dy2(1,N)=0.d0
     u(1,1)=0.d0
     u(1,N)=0.d0
     do i=2,N-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*dy2(1,i-1)+2.d0
        dy2(1,i)=(sig-1.d0)/p
        u(1,i)=(6.d0*((y(1,i+1)-y(1,i))/(x(i+1)-x(i))-(y(1,i)-y(1,i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(1,i-1))/p
     enddo
     
     do k=N-1,2,-1
        dy2(1,k)=dy2(1,k)*dy2(1,k+1)+u(1,k)
     enddo
  else 
     dy2(:,1)=0.d0
     dy2(:,N)=0.d0
     u(:,1)=0.d0
     u(:,N)=0.d0
     do i=2,N-1
        do j=1,2
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*dy2(j,i-1)+2.d0
           dy2(j,i)=(sig-1.d0)/p
           u(j,i)=(6.d0*((y(j,i+1)-y(j,i))/(x(i+1)-x(i))-(y(j,i)-y(j,i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(j,i-1))/p
        enddo
    enddo

    do k=N-1,2,-1   
        do j=1,2
           dy2(j,k)=dy2(j,k)*dy2(j,k+1)+u(j,k)
        enddo
    enddo

  endif
end subroutine cubic_spline

subroutine getfil_yshift(vor,phi,u00,w00,wk,tiempo)

  use ctes
  use running
  use bcs
  use LES, only: iuse_LES

  implicit none
  include "mpif.h"

  integer istat(MPI_STATUS_SIZE),ierr

  integer i,j,k,jj
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  integer ihalf,iskip,ihalfy

  real*8 vor(2,mx1+1,mz,jb:je),phi(2,mx1+1,mz,jb:je),wk(*)
  real*8 u00(my),w00(my)

  real*8, allocatable:: wk00(:,:), wk2(:,:,:,:), ye(:)
  real*8, allocatable:: u00tmp(:), w00tmp(:)
  real*8, allocatable:: wksh(:,:,:,:)

  integer  mgalxe,mye,mgalze,numerope
  !real*8  alpe,game,se,chie ! Ree, Lye is now in mod.f90 to share
  !real*8  ye(my)
  real*8  tiempo

  ! only for 1 proc.
  if (numerop.ne.1) then
     write(*,*) 'getfil_yshift ERROR; use only 1 proc...STOP '
     stop
  endif

  !      ---------------  zero everything, fitf
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  wk(buffsize)=0.d0

  if (myid.eq.0) then     ! ----- this is done by the master

     open(iinp,file=filinp,status='unknown',form='unformatted',action='read')

     read(iinp) mgalxe,mye,mgalze,numerope

     if((mgalx.ne.mgalxe).or.(my.ne.mye).or.(mgalz.ne.mgalze)) then
        
        write(*,*) ' different physical grid points!'
        write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,mye,mgalze, &
             ' =>', mgalx, my, mgalz 
        
        if ((my.ne.mye).and.(mgalx.ne.mgalxe)) then
           write(*,*) ' Do not interpolate in my at the same time...'
           stop
        end if
        !stop
     endif

     read(iinp) tiempo,Ree,alpe,Lye,game,se,chie,me ! me=mx,my,mz

     !write(*,*) 'read Ree',Ree ! do not forget to BCAST
     Lye=Lye*pi
     if (abs(Lye-Ly).gt.tiny) write(*,*)' getfil: Ly is different !!',Lye
     if (abs(se-s).gt.tiny) then 
        if (abs(se).lt.tiny) then
           write(*,*)' Now the shear is started'
        elseif (abs(s).lt.tiny) then
           write(*,*)' Now the shear is removed'
        else
           write(*,*)' getfil: s is different, STOP !!',se
           stop
        end if
     endif
     allocate(ye(1:me(2))) ! myid = 0
     read(iinp) ye(1:me(2))
     read(iinp) xoffb,xofft,zoffb,zofft
     write(*,*) 'reading input file ..., xoffb =',xoffb 
     write(*,*) '                   ..., xofft =',xofft
     
     !rescaling offset by shear
     ihalf=0; ihalfy=0;
     if (alpe.ne.alp) then
   
        if ((int((alp+0.00001)/alpe).eq.2).and.(mgalxe/mgalx.eq.2)) then
           xoffb=xoffb*alpe/alp*2.d0
           xofft=xofft*alpe/alp*2.d0
           write(*,*) ' '
           write(*,*) ' make the box to be half using only even modes...'
           write(*,*) '   keep the offset xoffb,xofft = ',xoffb,xofft
           ihalf=1
        else
           xoffb=xoffb*alpe/alp
           xofft=xofft*alpe/alp
           write(*,*) ' rescaling shear-periodic b.c. offset: factor =', &
                alpe/alp, '  xoffb,xofft = ',xoffb,xofft
        endif

        if (me(2).ne.my) then
           write(*,*) ' Do not interpolate in my at the same time...'
           stop
        end if

     end if
 
     write(*,*)
     write(*,*) 'reading input file ..., time=',tiempo
     write(*,'(a,3(1X,I5))') '                    ..., (mx,my,mz)=',me
     write(*,*)

     do iproc=1,numerop-1
        call MPI_SEND(me,3,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
        !call MPI_SEND(y    ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ye   ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zoffb, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(zofft, 1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(tiempo,1,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr) 
        ! fixed (rev.544)

        call MPI_SEND(ihalf,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ihalfy,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
     enddo

  else       ! -------- this is done by all slaves

     call MPI_RECV(me,3,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     !call MPI_RECV(y   ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     allocate(ye(1:me(2))) ! slaves allocate, 
     call MPI_RECV(ye  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(xofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zoffb, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(zofft, 1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(tiempo,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

     call MPI_RECV(ihalf,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(ihalfy,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  endif


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !write(*,*) 'allocate wk00', me(2), my
  allocate( wk00(me(2),2) )
  !     ------------- 00 modes ------------------
  if (myid.eq.0) then     ! ----- this is done by the master

     if (me(2).lt.my) then
        if(myid.eq.0) then
           write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' => my=',my
        end if
     end if
     if (my.eq.me(2)/2) then
        ihalfy=1
        if (myid.eq.0) write(*,*) 'reducing my to be half '
     end if

     write(*,*) 'allocated wk00', me(2), my
     wk00 = 0.d0
     !write(*,*) 'reading 00 modes'
     read(iinp) wk00
     mym = min(my,me(2))
     if (ihalfy.eq.1) then
        u00(1:my) = wk00(1:me(2):2,1)
        w00(1:my) = wk00(1:me(2):2,2)
        !write(*,*) 'read 00 mode'
     else
        u00(1:mym) = wk00(1:mym,1)
        w00(1:mym) = wk00(1:mym,2)
     end if

     ! yshift 
     allocate(u00tmp(my),w00tmp(my))
     do j=1,my
        jj= mod(j+my/2-1,my)+1
        u00tmp(j) = u00(jj) 
        w00tmp(j) = w00(jj) 
     end do
     u00=u00tmp
     w00=w00tmp
     deallocate(u00tmp,w00tmp)

     do iproc=1,numerop-1
        call MPI_SEND(u00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w00  ,my,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,ierr)
     end do
  else       ! -------- this is done by all slaves
     call MPI_RECV(u00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
     call MPI_RECV(w00  ,my,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  endif
  deallocate (wk00)


  ! ---- the other modes ----
  ntotr=me(1)*me(3)*2
  allocate ( wk2(2,2,me(1)/2,me(3)) )
  wk2 = 0.d0
  if (ihalf.eq.1) then
     iskip=2
  else
     iskip=1
  endif

  mxm=min(mx,me(1))
  mxm=mxm/2
  mym=min(my,me(2))
  mzm=min(mz,me(3))

  klen=(mzm+1)/2
  kini1=me(3) - klen + 1
  kini2=mz    - klen + 1
  
  ! first expand the fourier mode in x-z and then use cubic-spilne to interpolate in y
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then     ! ----- this is done by the master, first vor

     !write(*,*) 'reading other modes'
    
     do iproc=0,numerop-1   
        do j=jbeg(iproc),jend(iproc)
           if (j.le.(me(2)-1)) then ! the case for my being increased
              read(iinp) wk2
           else
              wk2=0.d0
           endif

           if (iproc.eq.0) then  ! master stores its thing
              vor(:,1:mxm,1:klen,j)    =wk2(1,:,1:mxm*iskip:iskip,1:klen)
              vor(:,1:mxm,1+kini2:mz,j)=wk2(1,:,1:mxm*iskip:iskip,1+kini1:me(3))
              phi(:,1:mxm,1:klen,j)    =wk2(2,:,1:mxm*iskip:iskip,1:klen)
              phi(:,1:mxm,1+kini2:mz,j)=wk2(2,:,1:mxm*iskip:iskip,1+kini1:me(3))
           else      ! sends everything else to nodes
              call MPI_SEND(wk2,ntotr,MPI_REAL8,iproc,iproc, & 
                   MPI_COMM_WORLD,ierr)
           endif

           if (ihalfy.eq.1) then
              ! the case for my being half
              read(iinp) wk2
              write(*,*) 'skip one plane by just reading '              
           endif

         enddo

     enddo

     close(iinp)

     ! yshift 
     allocate(wksh(2,mx1+1,mz,jb:je)) ! j starts from zero
     do j=0,my1
        jj= mod(j+my/2,my)
        wksh(:,:,:,j) = vor(:,:,:,jj) 
     end do
     vor=wksh;
     do j=0,my1
        jj= mod(j+my/2,my)
        wksh(:,:,:,j) = phi(:,:,:,jj) 
     end do
     phi=wksh;


     deallocate(wksh)
     
  else          ! -------- the slaves receive higher modes ----------

     do j=jb,je
        if ((ihalfy.eq.1).and.(mod(j+1,2).eq.0)) then
           ! Skip reading
           continue
        endif
        call MPI_RECV(wk2,ntotr,MPI_REAL8,0,MPI_ANY_TAG, & 
             MPI_COMM_WORLD,istat,ierr)
        vor(:,1:mxm,1:klen,j)    =wk2(1,:,1:mxm*iskip:iskip,1:klen)
        vor(:,1:mxm,1+kini2:mz,j)=wk2(1,:,1:mxm*iskip:iskip,1+kini1:me(3))
        phi(:,1:mxm,1:klen,j)    =wk2(2,:,1:mxm*iskip:iskip,1:klen)
        phi(:,1:mxm,1+kini2:mz,j)=wk2(2,:,1:mxm*iskip:iskip,1+kini1:me(3)) 
     enddo

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! reading vb vt is removed (rev260, 2012/Feb/02)
  deallocate (wk2)
  call  MPI_BCAST(Ree,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(Lye,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 ! BCAST are forgotten (2013/10/14) alpe,Lye,game,se,chie,me
  call  MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(game,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(se,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(chie,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! ----- before going any further exchange cuts in y with cuts in z!
  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)

  ! ------ you should not interpolate from a small number of grids to ------
  ! ------- a much larger, interpolation itself has error --------
  if (me(2).lt.my) then
     call interp_y(vor,phi,u00,w00,me,ye)
  endif 

  deallocate(ye) ! deallocated by all proc.
  ! cleanup work array for change (to avoid an un-expected bug, rev(635))
  wk(buffsize)=0.d0

  xwkb=xoffb; xwkt=xofft; ! bug fixed, 2014/June/24 sekimoto
  zwkb=zoffb; zwkt=zofft;

end subroutine getfil_yshift


subroutine interpolate(c,x,y,dy2,N,y2,m)
  integer N,klo,khi,k,m
  real*8 x(N),y(m,N),dy2(m,N),h,a,b,y2(m),c
  klo=1
  khi=N
  do while (khi-klo.gt.1) 
     k=(khi+klo)/2
     if(x(k).gt.c) then
        khi=k
     else
        klo=k
     endif
  enddo
   
  h=x(khi)-x(klo)
  a=(x(khi)-c)/h
  b=(c-x(klo))/h
  if (m.eq.1) then
     y2(1)=a*y(1,klo)+b*y(1,khi)+((a**3-a)*dy2(1,klo)+(b**3-b)*dy2(1,khi))*h**2/6.d0
  else 
     y2(1)=a*y(1,klo)+b*y(1,khi)+((a**3-a)*dy2(1,klo)+(b**3-b)*dy2(1,khi))*h**2/6.d0
     y2(2)=a*y(2,klo)+b*y(2,khi)+((a**3-a)*dy2(2,klo)+(b**3-b)*dy2(2,khi))*h**2/6.d0
  endif

end subroutine interpolate

!/********************************************************************/
!/*                                                                  */
!/*        writes statistics and an intermediate solution            */
!/*                                                                  */
!/*    single  jjs  Aug./07                                          */
!/*    double  A.S.                                                  */
!/********************************************************************/
subroutine escru(vor,phi,u00,w00,wk1,wk2,wk)

  use ctes
  use running
  use statistics
  use bcs
  use LES

  implicit none
  include "mpif.h"

  integer i,j,k, iproc, ipo, leng ! , nopt, nparams
  real(8),dimension(mx*mz,jb:je) :: wk1, wk2 
  real(8),dimension(buffsize) :: vor, phi, wk
  real(8),dimension(0:my1) :: u00, w00

  real(8) write_time 
  character*4 ext1
  character*128 filnam
  integer istat(MPI_STATUS_SIZE),ierr
  !integer ista3, ista4, ivar,

  write_time = -MPI_WTIME()
  if(myid.eq.0) then         !!!!! coming from node 0

     if(id22.gt.9999.or.id22.lt.0) then
        write(*,*) 'number of images out of range',id22
        write(*,*) 'check ifile in hre.dat'
        ifatal=1 
        !call MPI_FINALIZE(ierr)
        stop
     else
        write(ext1,'(i4.4)') id22 ! from rev.1095
        id22 = id22+1
     endif

  endif

  ! /*       get and write statistics first       */
  if (nacum.ne.0) then

     if (myid.eq.0) write(*,*) 'stat esc', nacum
     if (myid.eq.0) then     
        do iproc=1,numerop-1
           ipo=jbeg(iproc)
           leng=(jend(iproc)-jbeg(iproc)+1)*nstat
           call MPI_RECV(stats(1,ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        enddo

        stats(17:19,:) = stats(17:19,:)/dble(mgalx*mgalz)  
        !!! triple products  ?????

        filnam = filout(1:index(filout,' ')-1)//'.'//ext1//'.sta3' 
        open(isn,file=filnam,status='unknown',form='unformatted')
        rewind(isn)
        write(isn) nacum
        write(isn) mgalx,my,mgalz
        write(isn) time,stime0,Re,alp,Ly/pi,gam,s,chi
        write(isn) y
        write(isn) nstat
        write(isn) stats
        !close(isn) ! for adding sp(:,:,:) later

     else

        call MPI_SEND(stats(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)

     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta3
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl,sp,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (myid.eq.0) then
        !write(*,*) 'spec esc', nacumsp ! = nacum
        write(*,*) 'writing spec in .sta3'
        write(isn) nspec
        write(isn) sp(:,:,:)

        write(isn) CFL, uprim, vbulk, xforce ! rev934 ! xforce is redundunt, but keep the format
        write(isn) irev,nopt,nparams
        write(isn) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
             &      iuse_LES,iadd_visc,idynamic ! rev1200 
        write(isn) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(1:3)

        close(isn)
     endif


     ! for testing dt-multiplied statistics in sta4 


     if (myid.eq.0) write(*,*) 'stat sta4 ', total_dt
     if (myid.eq.0) then     
        do iproc=1,numerop-1
           ipo=jbeg(iproc)
           leng=(jend(iproc)-jbeg(iproc)+1)*nstat
           call MPI_RECV(stats2(1,ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        enddo

        stats2(17:19,:) = stats2(17:19,:)/dble(mgalx*mgalz)  
        !!! triple products  ?????

        filnam = filout(1:index(filout,' ')-1)//'.'//ext1//'.sta4_fixed' 
        open(isn,file=filnam,status='unknown',form='unformatted')
        rewind(isn)
        write(isn) nacum,total_dt
        write(isn) mgalx,my,mgalz
        write(isn) time,stime0,Re,alp,Ly/pi,gam,s,chi
        write(isn) y
        write(isn) nstat
        write(isn) stats2
        !close(isn) ! for adding sp(:,:,:) later

     else

        call MPI_SEND(stats2(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)

     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta4
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl2,sp2,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (myid.eq.0) then
        !write(*,*) 'spec esc', nacumsp ! = nacum
        write(*,*) 'writing spec in .sta4'
        write(isn) nspec
        write(isn) sp2(:,:,:)

        write(isn) CFL, uprim, vbulk, xforce ! rev934 ! xforce is redundunt, but keep the format
        write(isn) irev,nopt,nparams
        write(isn) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
             &      iuse_LES,iadd_visc,idynamic ! rev1200 
        write(isn) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(1:3)

        close(isn)
     endif

  endif   ! ----- if nacum==0 ---


  ! /*       get and write spectra       */
  !if (nacumsp.ne.0) then
  !
  !   if (myid.eq.0) then 
  !      !write(*,*) 'spec esc', nacumsp
  !      filnam = filout(1:index(filout,' ')-1)//'.'//ext1//'.spe'
  !   endif
  !
  !   call write_spec(filnam)
  !
  !endif   ! ----- if nacumsp==0 ---

  ! reset the buffers for statistics
  if(myid.eq.0) write(*,*) '    note: reset stats = 0.0   '
  nacum = 0
  stats = 0.0
  ! reset the spectra 
  if(myid.eq.0) write(*,*) '    note: reset sp = 0.0, spl = 0.0   '
  nacumsp = 0
  sp = 0.d0; spl = 0.d0

  if(myid.eq.0) write(*,*) '    note: reset sp2 = 0.0, spl = 0.0   '
  stats2 = 0.0; sp2 = 0.d0; spl2 = 0.d0
  total_dt = 0;

  stime0 = time
  ! /*       write restart file       */

  call chjik2ikj(vor,wk1,wk,wk)       ! --- first  vor
  call chjik2ikj(phi,wk2,wk,wk)       ! --- phi 

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then

     !u00 = u00 + u00cut
     if(myid.eq.0) write(*,*) ' field dumped at St = ', xofft !time*s
     filnam = filout(1:index(filout,' ')-1)//'.'//ext1//trim(filext)
     open (iout,file=filnam,status='unknown',form='unformatted')
     rewind(iout)
     write(*,*) 'in escru restart: ',trim(filnam)  
     
     !write(iout) irev,iuseQL,iCouette,iadd_sym  
     write(iout) mgalx,my,mgalz,numerop
     write(iout) time,Re,alp,Ly/pi,gam,s,chi,mx,my,mz  
     write(iout) y
     write(iout) xoffb,xofft,zoffb,zofft
     !!!write(iout) xwkb,xwkt,zwkb,zwkt ! for the RK substep ! for debug
     write(iout) u00,w00       

     !------the master first writes its stuff,  
     !      then receives from everybody, and writes it
     !
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
           call MPI_RECV(wk1,leng,MPI_REAL8,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk2,leng,MPI_REAL8,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        endif
        do j=0,jend(iproc)-jbeg(iproc)
        !do j=jbeg(iproc),jend(iproc) !this is wrong, 
                                      !since wk1(mx,jb,je) of master
           write(iout) (wk1(i,j),wk2(i,j), i=1,mx*mz)
        enddo
     enddo
     write(iout) CFL, uprim, vbulk, xforce ! rev934 ! xforce is redundunt, but keep the format
     !write(iout) irev,explicit,iadd_force,iadd_mode,iadd_sym ! rev1153  
     !nopt = 9 ! the number of switch options
     !nparams = 13 ! the number of input parameters
     write(iout) irev,nopt,nparams
     write(iout) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
          &      iuse_LES,iadd_visc,idynamic ! rev1200 
     write(iout) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(1:3)
     !u00 = u00 - u00cut

  else             ! --- everybody else sends things to master 
     call MPI_SEND(wk1,mx*mz*mmy,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk2,mx*mz*mmy,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write_time = MPI_WTIME()+write_time
  if (myid.eq.0) then

     write(*,*) 'closing files: time to write:',write_time
     close(iout)
  endif

end subroutine escru


subroutine writefield_ppp(vor,wkr,ext)

  ! write instantaneous field 'vor' (real*4)

  use ctes
  use bcs
  use running ! iout

  implicit none
  include "mpif.h"

  integer i,j,k, iproc, ipo, leng
  real*4 wkr(mgalx+2,mgalz,jb:je)
  real*4 vor(mgalx+2,mgalz,jb:je)
  real*8 write_time
  character*90 filename
  integer istat(MPI_STATUS_SIZE),ierr
  character*4 ext  ! from rev.1095
  character*4 ext1

  wkr = vor
  iout = 99
  if (myid.eq.0) then

     !write(ext1,'(i3.3)') id22
     !
     filename = filinp(1:index(filinp,' ')-1)//ext(1:4)
     write(*,*) 'writing :',trim(filename)
     open (iout,file=filename,status='unknown',form='unformatted')
     rewind(iout)
     
     !write(iout) irev,iuseQL,iCouette,iadd_sym
     write(iout) mgalx,my,mgalz,numerop
     write(iout) time,Re,alp,Ly/pi,gam,s,chi,mx,my,mz  
     write(iout) y
     write(iout) xoffb,xofft,zoffb,zofft
     !!!write(iout) xwkb,xwkt,zwkb,zwkt ! for the RK substep ! debug

     !------the master first writes its stuff,  
     !      then receives from everybody, and writes it
     !
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*(mgalx+2)*mgalz
           call MPI_RECV(wkr,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)  
        endif
        do j=0,jend(iproc)-jbeg(iproc)
           write(iout) ((wkr(i,k,j), i=1,mgalx),k=1,mgalz)
        enddo
     enddo

  else             ! --- everybody else sends things to master

    call MPI_SEND(vor,(mgalx+2)*mgalz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  close(iout)

end subroutine writefield_ppp


subroutine write_spec(specname)
  ! modified for homogeneous shear flow 
  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"

  integer i,j,k, iproc, ipo, leng, idum

  character*4 ext1, specname*80
  integer istat(MPI_STATUS_SIZE),ierr

  !leng=nspec*(mx1+1)*(nz1+1)
  leng=nspec*(mx1+1)*(mz1+1)
  call MPI_ALLREDUCE(spl,sp,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then
     write(*,*) 'spec esc', nacumsp
     write(*,*) 'writing :',trim(specname)
     open(isn,file=specname,status='unknown',form='unformatted')
     !rewind(isn)
     write(isn) nacumsp
     idum = 0 ! this used to be 'nspecy' 
     write(isn) mx1,mz1,my,time,Re,alp,Ly/pi,gam,s,chi,idum
     !
     write(isn) nspec
     write(isn) sp(:,:,:)
     write(isn) stime0 ! added to know when start counting
     close(isn)
     
  endif

end subroutine write_spec

subroutine write_stat

  use ctes
  use running
  use statistics
  use bcs

  implicit none
  include "mpif.h"

  integer i,j,k, iproc, ipo, leng

  real*8 write_time
  character*4 ext1
  character*128 filnam
  integer istat(MPI_STATUS_SIZE),ierr
  !integer ista3, ista4, ivar,


  ! /*       get and write statistics first       */
  write(ext1,'(i4.4)') id22 ! from rev.1095

  if (nacum.ne.0) then

     if (myid.eq.0) write(*,*) 'stat esc', nacum
     if (myid.eq.0) then     
        do iproc=1,numerop-1
           ipo=jbeg(iproc)
           leng=(jend(iproc)-jbeg(iproc)+1)*nstat
           call MPI_RECV(stats(1,ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        enddo

        stats(17:19,:) = stats(17:19,:)/dble(mgalx*mgalz)  
        !!! triple products  ?????

        filnam = filout(1:index(filout,' ')-1)//'.'//ext1//'.sta3.check' 
        open(isn,file=filnam,status='unknown',form='unformatted')
        rewind(isn)
        write(isn) nacum
        write(isn) mgalx,my,mgalz
        write(isn) time,stime0,Re,alp,Ly/pi,gam,s,chi
        write(isn) y
        write(isn) nstat
        write(isn) stats
        !close(isn) ! for adding sp(:,:,:) later

     else

        call MPI_SEND(stats(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)

     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta3
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl,sp,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (myid.eq.0) then
        !write(*,*) 'spec esc', nacumsp ! = nacum
        write(*,*) 'writing spec in .sta3'
        write(isn) nspec
        write(isn) sp(:,:,:)
        close(isn)
     endif


     ! for testing dt-multiplied statistics in sta4 


     if (myid.eq.0) write(*,*) 'stat sta4 ', total_dt
     if (myid.eq.0) then     
        do iproc=1,numerop-1
           ipo=jbeg(iproc)
           leng=(jend(iproc)-jbeg(iproc)+1)*nstat
           call MPI_RECV(stats2(1,ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc, &
                MPI_COMM_WORLD,istat,ierr)
        enddo

        stats2(17:19,:) = stats2(17:19,:)/dble(mgalx*mgalz)  
        !!! triple products  ?????

        filnam = filout(1:index(filout,' ')-1)//'.'//ext1//'.sta4_fixed' 
        open(isn,file=filnam,status='unknown',form='unformatted')
        rewind(isn)
        write(isn) nacum,total_dt
        write(isn) mgalx,my,mgalz
        write(isn) time,stime0,Re,alp,Ly/pi,gam,s,chi
        write(isn) y
        write(isn) nstat
        write(isn) stats2
        !close(isn) ! for adding sp(:,:,:) later

     else

        call MPI_SEND(stats2(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)

     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta4
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl2,sp2,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (myid.eq.0) then
        !write(*,*) 'spec esc', nacumsp ! = nacum
        write(*,*) 'writing spec in .sta4'
        write(isn) nspec
        write(isn) sp2(:,:,:)
        close(isn)
     endif

  endif   ! ----- if nacum==0 ---


  ! reset the buffers for statistics
  if(myid.eq.0) write(*,*) '    note: reset stats = 0.0   '
  nacum = 0
  stats = 0.0
  ! reset the spectra 
  if(myid.eq.0) write(*,*) '    note: reset sp = 0.0, spl = 0.0   '
  nacumsp = 0
  sp = 0.d0; spl = 0.d0

  if(myid.eq.0) write(*,*) '    note: reset sp2 = 0.0, spl = 0.0   '
  stats2 = 0.0; sp2 = 0.d0; spl2 = 0.d0
  total_dt = 0;

  stime0 = time
 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine write_stat


subroutine check_parameters(nopte,ioptions,nparamse,params)

  use ctes
  use running
  use LES
  
  implicit none
  include "mpif.h"
  
  integer nparamse, nopte, ierr
  integer ioptions(1:nopte)
  real*8 params(1:nparamse)
  real*8 Cstmp

  if (ireve.ge.1235) then
     ! (rev.1235)
     !write(iout) irev,nopt,nparams
     !write(iout) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
     !     &      iuse_LES,iadd_visc,idynamic ! rev1200 
     !write(iout) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(3)
    
 
     if (myid.eq.0) write(*,*) myid,'checking footers (from rev.1235) nopte:',nopte,ioptions
     if (myid.eq.0) write(*,*) myid,'checking footers (from rev.1235) nparamse',nparamse,params
     if (iadd_force.ne.1) then
        force_rolle=params(2)
        xforcee = params(3)
        zforcee = params(4)
        if (abs(force_roll-force_rolle).gt.tiny) then
           if (myid.eq.0) write(*,*) ' force_rolle =',force_rolle, ' force_roll=',force_roll
        end if
        if (abs(xforcee-xforce).gt.tiny) then
           if (myid.eq.0) write(*,*) ' xforcee =',xforcee, ' xforce=',xforce
        end if
        if (abs(zforcee-zforce).gt.tiny) then
           if (myid.eq.0) write(*,*) ' zforcee =',zforcee, ' zforce=',zforce
        end if
     end if

     if (iadd_damping.eq.1) then
        damp_aae=params(5)
        damp_upe=params(6)
        if (abs(damp_aae-damp_aa).gt.tiny) then
           if (myid.eq.0) write(*,*) ' damp_aae =',damp_aae, ' damp_aa=',damp_aa
        end if
        if (abs(damp_upe-damp_up).gt.tiny) then
           if (myid.eq.0) write(*,*) ' damp_upe =',damp_upe, ' damp_up=',damp_up
        end if 
     end if

     if (iuse_LES.eq.1) then
        Clese=params(7)
        if (abs(Clese-Cles).gt.tiny) then ! This is Cs (Smagorinsky constant) 
           if (myid.eq.0) write(*,*) ' Clese =',Clese, ' Cles=',Cles
        end if 
        Deltage=params(8) ! from rev.1299

        if (ifix_CsDeltag.eq.1) then
           Cstmp=Clese*Deltage ! ref.1348
           ! error check
           if (CsDeltag_fix.gt.0.d0) then
              if (abs(CsDeltag_fix-Cstmp).gt.tiny) then
                 write(*,*) 'something wrong for Cs*Deltag fixed'
                 stop
              end if
           end if
           CsDeltag_fix=Cstmp

           Cles=CsDeltag_fix/Deltag
           if(myid.eq.0) write(*,*) ' ifix_CsDeltag: Cs is adjusted', Cles
        end if

        if (abs(Deltage-Deltag).gt.tiny) then
           if (myid.eq.0) then
              write(*,*) ' Deltage =',Deltage, ' Deltag=',Deltag
              write(*,*) ' Deltag is not a parameter, this should not happen!'
              if (abs(Deltag*Cles - Deltage*Clese).gt.tiny) then
                 write(*,*) ' stop if follows are differnt Delta*Cs: ',Deltag*Cles,Deltage*Clese
                 !stop
              else
                 write(*,*) '  Cs is correctly adjusted ==> run '
              end if
           end if
        end if
     end if     

  else
!     if (myid.eq.0) write(*,*) 'skip checking footers in the previous field (from rev.1235)'

  end if

end subroutine check_parameters
