module ppp
  ! 3D arrays for output 
  real*4, allocatable:: ur(:,:,:), vr(:,:,:), wr(:,:,:)
  real*4, allocatable:: oxr(:,:,:),oyr(:,:,:),ozr(:,:,:)
  real*4, allocatable:: wkr(:,:,:)

  !real*4, allocatable(:,:,:):: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  real*4, allocatable:: Qa(:,:,:), Ra(:,:,:) 
  real*4, allocatable:: Qs(:,:,:), Qw(:,:,:)
  real*4, allocatable:: Rs(:,:,:)

  real*4, allocatable:: tempr(:,:,:)

  integer ifil(20),ifils(20),ishear,igetspec,nosvemos,isave_hdf5 ! for conv.f90

  ! 2D tmp plane for derivatives
  real*8,     allocatable:: ux(:,:),uy(:,:),uz(:,:),vx(:,:),vy(:,:), &
       vz(:,:),wx(:,:),wy(:,:),wz(:,:)
  complex*16, allocatable:: uxc(:,:),uyc(:,:),uzc(:,:),vxc(:,:),     &
       vyc(:,:),vzc(:,:),wxc(:,:),wyc(:,:),wzc(:,:)
  real*8,     allocatable:: ttt(:,:),ddd(:,:)
  real*8,     allocatable:: s12(:,:),s23(:,:),s13(:,:)

end module ppp 

program conv

  ! convert omega2, lap.v ==> velocities and vorticities (real*4)

  use ctes
  use running
  use statistics 
  use hdf5 ! for f90
  use ppp

  implicit none
  include "mpif.h"
  !include "hdf5.h" ! do not use this header file 

  integer iproc,itags,newtag,imess,i,j,iskip,idigit
  real*8  val, hard, t1, t2
  real*8, allocatable:: vor(:),phi(:),chwk(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: hv(:),hg(:),dvordy(:)

  ! -- set buffers for itemperature --
  real(8),dimension(:),allocatable :: tb
  real(8),dimension(:),allocatable :: tb00

  integer istat(MPI_STATUS_SIZE),ierr
  character*3 ext3, fname*256, filin_base*256, filout_base*256
  character*4 ext4 

  ! --- hdf5 ---
  integer h5err
  !  
  
  call MPI_INIT(ierr)

  call h5open_f(h5err) ! 

  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)
  !--------------- initializes commons and things
  call set_options()
  call initcr(0)
  write(*,*) myid, 'allocates buffers for convh5'
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), u00(my), w00(my) )
  vor = 0.d0; phi = 0.d0
  u00 = 0.d0; w00 = 0.d0
  !  -------  allocate rest of buffers 
  allocate( hv(buffsize),hg(buffsize),dvordy(buffsize))
  hv = 0.d0; hg = 0.d0; dvordy = 0.d0;

  if (itemperature.eq.1) then
     allocate(tb(buffsize),tb00(my))
     tb=0.d0; tb00=0.d0
  end if

  write(*,*) myid, 'buffers are allocated'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then
     ! initialize input option
     ! Note: do not forget to set the number of elements in MPI_BCAST     
     write(*,*) 'conv (fftw version, change_ieor, hre2, fixed v, error check, rev.824)'
     write(*,*) 'base file name to read and write', &
          ' (including relative path, like in hre.dat [filout]'
     read(*,'(a)') filin_base
     read(*,'(a)') filout_base
     write(*,*) 'read file index (from var.***, to var.*** ) [istart,iend,iskip]'
     read(*,*) istart, iend, iskip, idigit
     write(*,*) 'Which fields do you need ([(u [and t] v w) (ox oy oz) (Qa Ra) (Qs Rs) Qw]: (0 or 1)'
     read(*,*) ifil(1),ifil(2),ifil(3),ifil(4),ifil(5),ifil(6),ifil(7),ifil(8),ifil(9),ifil(10),ifil(11)
     write(*,*) 'Do you want spectram? (0 or 1)'
     read(*,*) igetspec
     !write(*,*) 'read in hdf5? (iread_hdf5 = 0 or 1)'
     !read(*,*) iread_hdf5
     write(*,*) 'going to write in hdf5 (isave_hdf5 = 1)'
     isave_hdf5 = 1
     write(*,*) 'Do you want them shifted in moving frame? (0 or 1)' 
     ! Now we also dump the shifted fields by S*y (ishear=1), (rev.554)
     read(*,*) ishear
     ! estimate data size (MB)
     hard = dfloat((sum(ifil)+1))*dfloat( (mgalx+2)*mgalz*my*(iend-istart+1)/(1000000/4) )
     write(*,*) 'need: ',hard, ' MByte for phisical fields (+alpha for spectra)'
     !write(*,*) 'We need ',(mx1+1)*(nz1+1)*7*nspec*8.d0/1000000.0, ' MByte for spectra'
     write(*,*) '(0: No  1: OK)'
     read(*,*) nosvemos
  end if

  call MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ishear,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(idigit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ifil,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  ifils=ifil; ! for a big HDF5 
  call MPI_BCAST(igetspec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(isave_hdf5,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_BCAST(iread_hdf5,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nosvemos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (nosvemos.eq.0) iend=istart-1 ! stop
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do istep = istart,iend,iskip
     !
     if (myid.eq.0) write(*,*) myid,'istep, itemperature=',istep,itemperature
     ! reinitialize the buffers
     vor = 0.d0; phi = 0.d0; 
     u00 = 0.d0; w00 = 0.d0;
     if (itemperature.eq.1) tb=0.d0
     !
     ! the shuffled buffers by change are dangerous to reuse
     hv = 0.d0; hg = 0.d0; dvordy = 0.d0;
     !
     !
     if (idigit.eq.3) then
        write(ext3,'(i3.3)') istep
     elseif (idigit.eq.4) then
        write(ext4,'(i4.4)') istep
     else 
        write(*,*) 'set idigit= 3 or 4', idigit
        stop
     endif
     nacumsp = nacumsp +1
     nacum   = nacum+1

     allocate(chwk(buffsize)); chwk=0.d0

     if (idigit.eq.3) then 
        filout =  filout_base(1:index(filout_base,' ')-1)//'.'//ext3(1:3)
     elseif (idigit.eq.4) then
        filout =  filout_base(1:index(filout_base,' ')-1)//'.'//ext4(1:4)
     endif
     if (iread_hdf5.eq.0) then
        if (idigit.eq.3) then 
           filinp =  filin_base(1:index(filin_base,' ')-1)//'.'//ext3(1:3)
        elseif (idigit.eq.4) then
           filinp =  filin_base(1:index(filin_base,' ')-1)//'.'//ext4(1:4)
        endif
        call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
     elseif(iread_hdf5.eq.1) then
        if (idigit.eq.3) then 
           filinp =  filin_base(1:index(filin_base,' ')-1)//'.'//ext3(1:3)//'.h5'
        elseif (idigit.eq.4) then
           filinp =  filin_base(1:index(filin_base,' ')-1)//'.'//ext4(1:4)//'.h5'
        endif
        if (itemperature.eq.1) then
           call getfil_t_hdf5(vor,phi,tb,u00,w00,tb00,chwk,time)
        else
           call getfil_hdf5(vor,phi,u00,w00,chwk,time) 
           ! only for real*4 data of Madrid HST, note that Ly is devided by pi in the original data format
        end if
     else
        if (myid.eq.0) write(*,*) 'set iread_hdf5 = 1 or 0'
        stop
     end if

     if (abs(Lye-Ly).gt.tiny) then 
        write(*,*) myid,' getfil: Ly is different, STOP!!',Lye,Ly
        stop
     end if
     if (abs(Ree-re).gt.tiny) then
        write(*,*) myid,' getfil: re is different, STOP!!',ree,re
        stop
     end if
     if (abs(alpe-alp).gt.tiny) then
        write(*,*) myid,' getfil: alp is different, STOP!!',alpe,alp
        stop
     end if
      if (abs(game-gam).gt.tiny) then
        write(*,*) myid,' getfil: gam is different, STOP!!',game,gam
        stop
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     if(myid.eq.0) write(*,*) '... file read:', filinp,istep,'/',iend

     t1= MPI_WTIME()
     call pre_uvwp(vor,phi,u00,w00,hv,hg,dvordy,chwk) ! precompute for uvwp...
     if (itemperature.eq.1) then
        !if(myid.eq.0) write(*,*) 'chjik2ikj, 1 / 6'
        call chjik2ikj(tb,tb,chwk,chwk)
     end if
     t2= MPI_WTIME()
     if (myid.eq.0) write(*,*) '  --- time: pre_uvwp =',(t2-t1)
     deallocate(chwk)
     
     if (iread_hdf5.eq.1) then
        ! separate into some processes for a big HDF5 data (H32).
        ! maximum 3 buffers for writing_ppp_hdf5 
        ifil=0; 
        if (((ifils(1).eq.1).or.(ifils(2).eq.1).or.(ifils(3).eq.1)).or.(ifils(11).eq.1)) then
           ifil(1)=ifils(1); ifil(2)=ifils(2); ifil(3)=ifils(3); ! uvw     
           if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
           t1= MPI_WTIME()
           call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
           t2= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: uvw =',(t2-t1)
           if (itemperature.eq.1) call get_t(tb,tb00)
           t1= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: t =',(t1-t2)
        end if

        ifil=0;
        if ((ifils(4).eq.1).or.(ifils(11).eq.1)) then
           ifil(4)=1; ifil(5)=1; ifil(6)=1; ifil(11)=1;! oxoyoz     
           if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
           t1= MPI_WTIME()
           call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
           t2= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: oxoyoz =',(t2-t1)          
        end if 

        ifil=0; 
        if ((ifils(7).eq.1).or.(ifils(8).eq.1)) then
           ifil(7)=1; ifil(8)=1; ! Qa Ra     
           if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
           t1= MPI_WTIME()
           call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
           t2= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: QR =',(t2-t1)      
        end if 

        ifil=0; 
        if (ifils(11).eq.1) then
           ! Qw     
           ifil(11)=1
           if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
           t1= MPI_WTIME()
           call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
           t2= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: Qw =',(t2-t1)      
        end if 

        ifil=0; 
        if ((ifils(9).eq.1).or.(ifils(10).eq.1)) then
           ifil(9)=1; ifil(10)=1; ifil(11)=1 ! QsRs     
           if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
           t1= MPI_WTIME()
           call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
           t2= MPI_WTIME()
           if (myid.eq.0) write(*,*) '  --- time: QsRs =',(t2-t1)      
        end if 

     else
        if (myid.eq.0) write(*,'(a,20(1X,I1))') 'ifil=',ifil
        call uvwp(vor,phi,u00,w00,hv,hg,dvordy) 
        ! write uvwp...
     end if
  end do

  deallocate( hv, hg, dvordy )
  deallocate( vor, phi, u00, w00 )
  if (itemperature.eq.1) then
     deallocate(tb,tb00)
  end if
  ! write spectra
  if (igetspec.eq.1) then
     if (myid.eq.0) then 
        fname = filout(1:index(filout,' ')-1)//'.spe3.conv' ! buggy 
     endif
     call write_spec(fname) ! buggy 
  end if
  !
  ! 
  call h5close_f(h5err)
  call MPI_FINALIZE(ierr)

end program conv

subroutine pre_uvwp(vor,phi,u00,w00,hv,hg,dvordy,chwk)

  use ctes
  use running
  use statistics 
  use bcs
  use ppp

  implicit none
  include "mpif.h"

  integer i,j,k
  real*8  rk,rk2

  complex*16 shp,shm
  real*8  u00(0:my1),w00(0:my1)
  real*8  phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  real*8  hv(0:2*my-1,0:mx1,kb:ke) 
  real*8  hg(0:2*my-1,0:mx1,kb:ke)
  real*8  dvordy(0:2*my-1,0:mx1,kb:ke)
  real*8  chwk(buffsize)
  
  xwkt = xofft; xwkb = xoffb 
  zwkt = zofft; zwkb = zoffb 
  
  if (s.gt.0.5d0) then ! s shoud be 0 or 1     
     do i=0,mx1
        shb(i)  = cdexp(-xalp(i)*xwkb) ! negative shift
        sht(i)  = cdexp(-xalp(i)*xwkt) ! positive shift
     enddo
  end if

  !   ---   calcula la v a partir de phi */
  if(myid.eq.0) write(*,*) 'calcula v'
  do k=kb,ke
     do i=0,mx1
        if ((i.ne.0).or.(k.ne.0)) then
           rk = gam2(k)+alp2(i)
           shp = sht(i)
           shm = shb(i)
           call lapcdy(phi(0,i,k),rk,hg(0,i,k),shp,shm,2)   ! -- v
        end if
     enddo
  enddo
  
  u00b = s*(-Ly) ! negative
  u00t = s*Ly ! positive
  w00b = 0.d0
  w00t = 0.d0

  if(myid.eq.0) write(*,*) 'computes ome2/dy, dvdy'
  do k=kb,ke          ! --- computes 
     do i=0,mx1
        shp = sht(i)
        shm = shb(i)
        ! --- d(ome_2)/ dy
        call derivyc(vor(0,i,k),dvordy(0,i,k),shp,shm,1,2,1)
        ! --- dv/dy 
        call derivyc(hg(0,i,k),hv(0,i,k),shp,shm,1,2,0) ! skip add_shear
        if (itemperature.eq.1) then
           
        end if
     enddo
  enddo

  !if(myid.eq.0) write(*,*) 'chjik2ikj, 1 / 5'
  call chjik2ikj(phi,phi,chwk,chwk)
  !if(myid.eq.0) write(*,*) 'chjik2ikj, 2 / 5'
  call chjik2ikj(vor,vor,chwk,chwk)
  !if(myid.eq.0) write(*,*) 'chjik2ikj, 3 / 5'
  call chjik2ikj(hv,hv,chwk,chwk)
  !if(myid.eq.0) write(*,*) 'chjik2ikj, 4 / 5'
  call chjik2ikj(hg,hg,chwk,chwk)
  !if(myid.eq.0) write(*,*) 'chjik2ikj, 5 / 5'
  call chjik2ikj(dvordy,dvordy,chwk,chwk)

end subroutine pre_uvwp

subroutine get_t(tempc,tb00)
  use ctes
  use running
  use statistics 
  use bcs
  use ppp

  implicit none
  include "mpif.h"

  integer i,j,k
  real*8  rk,rk2

  complex*16 shp,shm

  complex(8),dimension(0:mx1,0:mz1,jb:je) :: tempc
  real*8  tb00(0:my1)

  integer iproc,istat(MPI_STATUS_SIZE),ierr
  real*8, allocatable:: tbr(:,:) !, tbc(:,:) ! bug
  complex*16, allocatable:: tbc(:,:)

  ! tmp 2d arrays
  real*8, allocatable:: tmpxzr(:,:) ! for fouxz, phys2 

  allocate(tmpxzr(mgalx+2,mgalz))
  tmpxzr = 0.d0; 

  allocate(tbr(mgalx+2,mgalz),tbc(0:mx1,0:mz1) ) ! real*8 and complex*16
  tbr=0.d0; tbc=dcmplx(0.d0,0.d0);

  allocate(tempr(mgalx+2,mgalz,jb:je)); 
  tempr = 0.d0; ! real*4

  DO J = JB,JE 
     call dcopy((mx1+1)*(mz1+1)*2,tempc(0,0,j),1,tbc(0,0),1)
     tbc(0,0) = tb00(j)
     call fourxz(tbc,tbr,tmpxzr,1,1)     
     tempr(:,:,j) = tbr
  enddo

  deallocate(tbc,tbr)
  deallocate(tmpxzr)

  if (myid.eq.0) write(*,*) 'writing temperature'
  ! write physical fields
  allocate(wkr(mgalx+2,mgalz,jb:je)) ! this should be jsizemaxl
  wkr = 0.d0

  if (isave_hdf5.eq.1) then     
     call writefield_ppp_hdf5(tempr,wkr,'.tr'); 
  else
     call writefield_ppp(tempr,wkr,'.tr');
  end if

  deallocate (tempr)
  deallocate (wkr) 

end subroutine get_t

subroutine uvwp(vor,phi,u00,w00,hv,hg,dvordy)

  use ctes
  use running
  use statistics 
  use bcs
  use ppp

  implicit none
  include "mpif.h"

  integer i,j,k
  real*8  rk,rk2

  complex*16 shp,shm
  real*8  u00(0:my1),w00(0:my1)

  ! real*8  phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
  ! real*8  hv(0:2*my-1,0:mx1,kb:ke) 
  ! real*8  hg(0:2*my-1,0:mx1,kb:ke)
  ! real*8  dvordy(0:2*my-1,0:mx1,kb:ke)

  real*8 phi(buffsize), vor(buffsize), hv(buffsize), hg(buffsize), dvordy(buffsize)
  integer iproc,istat(MPI_STATUS_SIZE),ierr

  real*8, allocatable:: rf0u(:),rf0w(:)
  real*8, allocatable:: u1r(:,:),u2r(:,:),u3r(:,:),o1r(:,:),o2r(:,:),o3r(:,:)

  ! tmp 2d arrays
  real*8, allocatable:: tmpxzr(:,:) ! for fouxz, phys2 
  real*8, allocatable:: tmpx(:)
  real*8, allocatable:: wk1dr(:),wk1dc(:)
  character ext1*3

  ! ---- temporal for writing every field ----
  !character fil*90
  ! ---- -------------------------------------

  real*8 commtimer,transtimer,totaltimer
  common/timers/ commtimer, transtimer, totaltimer
  save/timers/
  real*8 iter_time,comm_time,t1,t2

  comm_time = 0d0
  commtimer = 0d0
  transtimer= 0d0
  totaltimer= 0d0

  ! ------------------- allocate everything ------------  
  allocate(u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz) )
  allocate(o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz) )
  allocate(rf0u(0:my1), rf0w(0:my1))
  u1r = 0.d0; u2r = 0.d0; u3r = 0.d0;
  o1r = 0.d0; o2r = 0.d0; o3r = 0.d0;
  rf0u= 0.d0; rf0w= 0.d0

  ! ---- allocate tmp 2D arrays ----
  allocate(tmpxzr(mgalx+2,mgalz))
  tmpxzr = 0.d0; 
  allocate(tmpx(mgalx+2))
  tmpx = 0.d0;
  allocate( wk1dr(my), wk1dc(2*my))
  wk1dr = 0.d0; wk1dc = 0.d0;

  if (myid.eq.0) write(*,*) myid,' 2D arrays are allocated in uvwp' 

  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !
  if (myid.eq.0) then
     totaltimer = totaltimer-MPI_WTIME()
     iter_time  = -MPI_WTIME()
  endif

  ! allocate tmp 3D arrays
  if(myid.eq.0) write(*,*) ' allocating tmp 3D arrays'
  if (ifil(1).eq.1) then
     allocate(ur(mgalx+2,mgalz,jb:je)); ur = 0.d0; 
  end if
  if (ifil(2).eq.1) then
     allocate(vr(mgalx+2,mgalz,jb:je)); vr = 0.d0; 
  end if
  if (ifil(3).eq.1) then 
     allocate(wr(mgalx+2,mgalz,jb:je)); wr = 0.d0;
  end if

  if (ifil(4).eq.1) then 
     allocate(oxr(mgalx+2,mgalz,jb:je)); oxr = 0.d0;
  end if
  if (ifil(5).eq.1) then 
     allocate(oyr(mgalx+2,mgalz,jb:je)); oyr = 0.d0; 
  end if
  if (ifil(6).eq.1) then
     allocate(ozr(mgalx+2,mgalz,jb:je)); ozr = 0.d0;
  end if

  if (ifil(7).eq.1) then 
     allocate(Qa(mgalx+2,mgalz,jb:je)); Qa = 0.d0;
  end if
  if (ifil(8).eq.1) then
     allocate(Ra(mgalx+2,mgalz,jb:je)); Ra = 0.d0;  
  end if

  if (ifil(9).eq.1) then
     allocate(Qs(mgalx+2,mgalz,jb:je)); Qs = 0.d0; 
  end if
  if (ifil(10).eq.1) then
     allocate(Rs(mgalx+2,mgalz,jb:je)); Rs = 0.d0;
  end if

  if (ifil(11).eq.1) then 
     allocate(Qw(mgalx+2,mgalz,jb:je)); Qw = 0.d0; 
  end if
  !Rw = 0.d0;
  t1= MPI_WTIME()
  if (myid.eq.0) write(*,*) 'hvhg_uvwp'
  call hvhg_uvwp(phi,vor,u00,w00,hv,hg,rf0u,rf0w,dvordy, & 
       u1r,u2r,u3r,o1r,o2r,o3r,u1r,u2r,u3r,o1r,o2r,o3r, &
       tmpxzr,tmpx,wk1dr)
  t2= MPI_WTIME()
  if (myid.eq.0) write(*,*) 'hvhg_uvwp, done, time=',t2-t1
  deallocate (tmpxzr,tmpx,wk1dr,wk1dc)
  deallocate (rf0u,rf0w)
  deallocate (u1r,u2r,u3r,o1r,o2r,o3r)
  !
  ! write physical fields
  allocate(wkr(mgalx+2,mgalz,jb:je)) ! this should be jsizemaxl
  wkr = 0.d0

  if (isave_hdf5.eq.1) then
     
     if (ifil(1).eq.1) then 
        call writefield_ppp_hdf5(ur,wkr,'.ur');deallocate(ur);
     end if
     if (ifil(2).eq.1) then
        call writefield_ppp_hdf5(vr,wkr,'.vr');deallocate(vr);
     end if
     if (ifil(3).eq.1) then 
        call writefield_ppp_hdf5(wr,wkr,'.wr');deallocate(wr);
     end if

     if (ifil(4).eq.1) then 
        call writefield_ppp_hdf5(oxr,wkr,'.ox');deallocate(oxr);
     end if
     if (ifil(5).eq.1) then 
        call writefield_ppp_hdf5(oyr,wkr,'.oy');deallocate(oyr);
     end if
     if (ifil(6).eq.1) then 
        call writefield_ppp_hdf5(ozr,wkr,'.oz');deallocate(ozr);
     end if

     if (ifil(7).eq.1) then
        call writefield_ppp_hdf5(Qa,wkr,'.Qa'); deallocate(Qa);
     end if
     if (ifil(8).eq.1) then 
        call writefield_ppp_hdf5(Ra,wkr,'.Ra'); deallocate(Ra);
     end if

     if (ifil(9).eq.1) then 
        call writefield_ppp_hdf5(Qs,wkr,'.Qs');deallocate(Qs);
     end if
     if (ifil(10).eq.1) then 
        call writefield_ppp_hdf5(Rs,wkr,'.Rs');deallocate(Rs);
     end if
     if (ifil(11).eq.1) then 
        call writefield_ppp_hdf5(Qw,wkr,'.Qw');deallocate(Qw);
     end if
  else

     if (ifil(1).eq.1)  call writefield_ppp(ur,wkr,'.ur');
     if (ifil(2).eq.1)  call writefield_ppp(vr,wkr,'.vr');
     if (ifil(3).eq.1)  call writefield_ppp(wr,wkr,'.wr');

     if (ifil(4).eq.1)  call writefield_ppp(oxr,wkr,'.ox');
     if (ifil(5).eq.1)  call writefield_ppp(oyr,wkr,'.oy');
     if (ifil(6).eq.1)  call writefield_ppp(ozr,wkr,'.oz');

     if (ifil(7).eq.1)  call writefield_ppp(Qa,wkr,'.Qa');
     if (ifil(8).eq.1)  call writefield_ppp(Ra,wkr,'.Ra');

     if (ifil(9).eq.1)  call writefield_ppp(Qs,wkr,'.Qs');
     if (ifil(10).eq.1)  call writefield_ppp(Rs,wkr,'.Rs');
     if (ifil(11).eq.1)  call writefield_ppp(Qw,wkr,'.Qw');

  end if

  deallocate (wkr)  
  
end subroutine uvwp

subroutine hvhg_uvwp(phic,ome2c,u00,w00,rhvc,rhgc,rf0u,rf0w,ome1c, &
     u1r,u2r,u3r,o1r,o2r,o3r, &
     u1c,u2c,u3c,o1c,o2c,o3c,tmpxzr,tmpx,tmpy)
  
  use ctes
  use running
  use bcs
  use ppp

  implicit none
  include "mpif.h"
  
  complex*16 phic(0:mx1,0:mz1,jb:je), & 
             ome1c(0:mx1,0:mz1,jb:je),ome2c(0:mx1,0:mz1,jb:je), &
             rhgc (0:mx1,0:mz1,jb:je),rhvc (0:mx1,0:mz1,jb:je)

  real*8 rf0u(0:my1),rf0w(0:my1),u00(0:my1),w00(0:my1)

  !------ 6 * (mgalx+2) * mgalz buffer planes, THESE PLANES SHARE STORAGE !!
  real*8 u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz),  &
         o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz)
  complex*16 u1c(0:mx1,0:mz1),u2c(0:mx1,0:mz1),u3c(0:mx1,0:mz1),  &
             o1c(0:mx1,0:mz1),o2c(0:mx1,0:mz1),o3c(0:mx1,0:mz1)
  ! -------------------------------------------------------------------------
  real*8 divr(mgalx+2,mgalz)
  complex*16 divc(0:mx1,0:mz1)
  real*8 tmpxzr(mgalx+2,mgalz)
  real*8 tmpx(mgalx+2)
  real*8 tmpy(0:my1)
  ! -------------------------------------------------------------------------
  integer i,j,k,iproc
  integer istat(MPI_STATUS_SIZE),ierr

  real*8 aa, v00, div_max, sytime, diff_xoff
  complex*16  shp,shm,cfac
  complex(8) kkwm1(0:mx1) 
  !complex*16  sydt(0:mx1,0:my1),sydb(0:mx1,0:my1)

  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
  !  at this point:
  !  rhv is dv / dy -- F-F-P
  !  rhg is v -- F-F-P
  !  phi is nabla^2(v) -- F-F-P
  !  ome1 is d(omega_2)/dy -- F-F-P
  !  ome2 is omega_2 --F-F-P
  ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

  ! ---------------- Initialize variables outside of the y loop
  !
  !--- Note: we need special treatments for the mean flow
  shp=dcmplx(1.d0,0.d0); shm=dcmplx(1.d0,0.d0)
  call derivyc(u00,rf0u,shp,shm,1,1,0)    ! ----- ome3_0 = -du0/dy ---
  !--- Note: derivative of fluctuation of u is pure-periodic
  call derivyc(w00,rf0w,shp,shm,1,1,0)    ! ----- ome1_0 =  dw0/dy ---

  umaxl= 0d0
  umax = 0d0

  !sytime = xoffb/(s*Ly) ! negative time
  ! Take care ... check Ly ...
  !sytime = (xoffb/Lx - int(xoffb/Lx) )*Lx/Ly
  !sytime = -10*Lx/(s*Ly)/20.d0
  !
  !if (myid.eq.0) then
  !   write(*,*) 'sytime =', sytime, xoffb/Ly
  !   do i=0,mx1
  !      do j=0,my1
  !         write(*,*) i,j,'diff time', cdexp(xalp(i)*y(j)*sytime) - cdexp(xalp(i)*y(j)*xoffb/Ly)
  !      end do
  !   enddo
  !end if
  !
  !do j=0,my1
  !   do i=0,mx1
  !      sydt(i,j) = cdexp(xalp(i)*y(j)*sytime) ! for negative shift
  !      sydb(i,j) = cdexp(xalp(i)*y(j)*xoffb/Ly) ! for positive shift
  !      diff_xoff = diff_xoff + (abs(sydt(i,j)-sydb(i,j)))**2
  !   end do
  !end do
  !stop
  !write(*,*) 'diff_xoff', diff_xoff/(mx1+1)/my
 
  allocate(ux(mgalx+2,mgalz),uy(mgalx+2,mgalz),uz(mgalx+2,mgalz))
  allocate(vx(mgalx+2,mgalz),vy(mgalx+2,mgalz),vz(mgalx+2,mgalz))
  allocate(wx(mgalx+2,mgalz),wy(mgalx+2,mgalz),wz(mgalx+2,mgalz))

  allocate(ttt(mgalx+2,mgalz),ddd(mgalx+2,mgalz))
  allocate(s12(mgalx+2,mgalz),s23(mgalx+2,mgalz),s13(mgalx+2,mgalz))  

  allocate(uxc(0:mx1,0:mz1),uyc(0:mx1,0:mz1),uzc(0:mx1,0:mz1))
  allocate(vxc(0:mx1,0:mz1),vyc(0:mx1,0:mz1),vzc(0:mx1,0:mz1))
  allocate(wxc(0:mx1,0:mz1),wyc(0:mx1,0:mz1),wzc(0:mx1,0:mz1))

  ux=0.d0; uy=0.d0; uz=0.d0;
  vx=0.d0; vy=0.d0; vz=0.d0;
  wx=0.d0; wy=0.d0; wz=0.d0;

  uxc=dcmplx(0.d0,0.d0); uyc=dcmplx(0.d0,0.d0); uzc=dcmplx(0.d0,0.d0)
  vxc=dcmplx(0.d0,0.d0); vyc=dcmplx(0.d0,0.d0); vzc=dcmplx(0.d0,0.d0)
  wxc=dcmplx(0.d0,0.d0); wyc=dcmplx(0.d0,0.d0); wzc=dcmplx(0.d0,0.d0)

  DO J = JB,JE       
     !---------------------- start operating by planes
     !---------------------- 00 modes, u,w,ome1,ome3 ! total velocity and vorticity with Sy
     u1c(0,0) =  u00(j) + s*y(j) ! add_shear, only for post-processr
     u3c(0,0) =  w00(j)

     o3c(0,0) = -rf0u(j) - s    ! add_shear, only for post-process
     o1c(0,0) =  rf0w(j)
     
     ! --------------------- computes non 0 modes of ome1, ome3, u1, u3 
     o3c(0,1:mz1) = -ome1c(0,1:mz1,j)/xgam(1:mz1)
     o1c(0,1:mz1) = - phic(0,1:mz1,j)/xgam(1:mz1)

     kkwm1(0:mx1)=0.d0
     do k=0,mz1
        kkwm1(1:mx1)=1.d0/(alp2(1:mx1)+gam2(k))
        o3c(1:mx1,k) = (ome1c(1:mx1,k,j)*xgam(k)-phic(1:mx1,k,j)*xalp(1:mx1)) &
                       *kkwm1(1:mx1)
        o1c(1:mx1,k) = (ome1c(1:mx1,k,j)*xalp(1:mx1)+phic(1:mx1,k,j)*xgam(k)) &
                       *kkwm1(1:mx1)
     end do

     u3c(0,1:mz1) = - rhvc(0,1:mz1,j)/xgam(1:mz1)
     u1c(0,1:mz1) =  ome2c(0,1:mz1,j)/xgam(1:mz1)
     do k=0,mz1
        kkwm1(1:mx1)=1.d0/(alp2(1:mx1)+gam2(k))
        u3c(1:mx1,k) = (rhvc(1:mx1,k,j)*xgam(k)+ome2c(1:mx1,k,j)*xalp(1:mx1)) &
                       *kkwm1(1:mx1)
        u1c(1:mx1,k) = (rhvc(1:mx1,k,j)*xalp(1:mx1)-ome2c(1:mx1,k,j)*xgam(k)) &
                       *kkwm1(1:mx1)
     enddo

     ! -------------- copy v and omega_2 into their planes
     u2c = rhgc(:,:,j)
     o2c = ome2c(:,:,j)  

     ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
     !  at this point:
     !  u1 is u,  u2 is v, u3 is w
     !  o1 is omega_1,  o2 is omega_2,  o3 is omega_3
     !  all variables in Fourierx -- Fourier z -- Physical y
     ! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

     ! ---------------- do spectra some day
     ! ---------------- only if my plane contains spectra information
     ! -- fourier statistics
     if (igetspec.eq.1) call hist(u1r,u2r,u3r,o1r,o2r,o3r, & 
                                  u1c,u2c,u3c,o1c,o2c,o3c, & 
                                  rhvc,j,'f')  

     ! COMPUTE VGT HERE for QR(7,8) QsRs(9,10) Qw(11)
     if ( (ifil(7).eq.1).or.(ifil(9).eq.1).or.(ifil(11).eq.1) ) then
        do k=0,mz1
           cfac=xgam(k)
           do i=0,mx1
              uzc(i,k) = u1c(i,k)*cfac     ! du/dz
              uxc(i,k) = u1c(i,k)*xalp(i)  ! du/dx
              vzc(i,k) = u2c(i,k)*cfac     ! dv/dz
              vxc(i,k) = u2c(i,k)*xalp(i)  ! dv/dx
              uyc(i,k) = vxc(i,k)-o3c(i,k) ! du/dy=dv/dx-o3c ! total o3c
              wyc(i,k) = vzc(i,k)+o1c(i,k) ! dw/dy=dv/dz+o1c
              wzc(i,k) = u3c(i,k)*cfac     ! dw/dz
              wxc(i,k) = u3c(i,k)*xalp(i)  ! dw/dx
           end do
        end do
        call fourxz(uxc,ux,tmpxzr,1,1)
        call fourxz(uzc,uz,tmpxzr,1,1)
     
        call fourxz(vzc,vz,tmpxzr,1,1)
        call fourxz(vxc,vx,tmpxzr,1,1)
        call fourxz(uyc,uy,tmpxzr,1,1)
        call fourxz(wyc,wy,tmpxzr,1,1)
       
        call fourxz(wxc,wx,tmpxzr,1,1)
        call fourxz(wzc,wz,tmpxzr,1,1)
        do k=0,mz1
           do i=0,mx1
              vyc(i,k) = rhvc(i,k,j)        ! dv/dy
           end do
        end do
        call fourxz(vyc,vy,tmpxzr,1,1)

     end if

     if (ifil(7).eq.1) then
        do k=1,mgalz
           do i=1,mgalx
              ttt(i,k) = ux(i,k)*ux(i,k) + vy(i,k)*vy(i,k) + wz(i,k)*wz(i,k) + & 
                   2.d0*(vz(i,k)*wy(i,k) + uy(i,k)*vx(i,k) + uz(i,k)*wx(i,k))
           enddo
        enddo
        Qa(:,:,j) = -0.5d0*ttt
     end if

     if ((ifil(8).eq.1)) then
        do k=1,mgalz
           do i=1,mgalx
              ddd(i,k) = ux(i,k)*(vy(i,k)*wz(i,k)-wy(i,k)*vz(i,k)) - &
                   &     uy(i,k)*(vx(i,k)*wz(i,k)-wx(i,k)*vz(i,k)) + &
                   &     uz(i,k)*(vx(i,k)*wy(i,k)-wx(i,k)*vy(i,k))
           enddo
        enddo
        Ra(:,:,j) = -ddd
     end if


     if (((ifil(9).eq.1).or.(ifil(10).eq.1)).or.(ifil(12).eq.1)) then
        do k=1,mgalz
           do i=1,mgalx
              s12(i,k) = 0.5d0*(uy(i,k) + vx(i,k)) 
              s23(i,k) = 0.5d0*(vz(i,k) + wy(i,k)) 
              s13(i,k) = 0.5d0*(uz(i,k) + wx(i,k))
           enddo
        enddo
        if (((ifil(9).eq.1).or.(ifil(10).eq.1)).or.(ifil(12).eq.1)) then
           do k=1,mgalz
              do i=1,mgalx
                 ttt(i,k) = -0.5d0*(ux(i,k)*ux(i,k) + vy(i,k)*vy(i,k) + wz(i,k)*wz(i,k)) - & ! bug fixed at rev1679 
                      ( s12(i,k)*s12(i,k) + s13(i,k)*s13(i,k) + s23(i,k)*s23(i,k) ) ! bug fixed at rev1680
                 ddd(i,k) = s12(i,k)*s23(i,k)*s13(i,k) + & 
                      ux(i,k)*(s12(i,k)*s12(i,k) + s13(i,k)*s13(i,k)) + &
                      vy(i,k)*(s12(i,k)*s12(i,k) + s23(i,k)*s23(i,k)) + & 
                      wz(i,k)*(s13(i,k)*s13(i,k) + s23(i,k)*s23(i,k)) + &
                      ux(i,k)*vy(i,k)*wz(i,k)
              enddo
           enddo
           Qs(:,:,j) = ttt
           Rs(:,:,j) = -ddd ! bug fixed at rev1681
        end if
     end if


     if ((ifil(1).eq.1).and.(ifil(2).eq.1).and.(ifil(3).eq.1)) then
        ! 
        call fourxz(u1c,u1r,tmpxzr,1,1)         ! u
        if (ifil(1).eq.1) ur(:,:,j) = u1r 
        call fourxz(u2c,u2r,tmpxzr,1,1)         ! v
        if (ifil(2).eq.1) vr(:,:,j) = u2r
        call fourxz(u3c,u3r,tmpxzr,1,1)         ! w
        if (ifil(3).eq.1) wr(:,:,j) = u3r
     end if
     if (((ifil(4).eq.1).and.(ifil(5).eq.1).and.(ifil(6).eq.1)).or.(ifil(11).eq.1)) then
        call fourxz(o1c,o1r,tmpxzr,1,1)         !omega_1
        if (ifil(4).eq.1) oxr(:,:,j) = o1r
        call fourxz(o2c,o2r,tmpxzr,1,1)         !omega_2
        if (ifil(5).eq.1) oyr(:,:,j) = o2r
        call fourxz(o3c,o3r,tmpxzr,1,1)         !omega_3 
        if (ifil(6).eq.1) ozr(:,:,j) = o3r        
        !
        if (ifil(11).eq.1) then
           do k=1,mgalz
              do i=1,mgalx
                 ttt(i,k) = 0.25d0*(o1r(i,k)*o1r(i,k) + o2r(i,k)*o2r(i,k) + o3r(i,k)*o3r(i,k))
              end do
           end do
           Qw(:,:,j) = ttt
        end if
     end if

     if (ifil(12).eq.1) then
        ! Rw = Ra - Rs
     end if

     ! -- physical statistics
     if (igetspec.eq.1) call hist(u1r,u2r,u3r,o1r,o2r,o3r, &
                                  u1c,u2c,u3c,o1c,o2c,o3c, & 
                                  rhvc,j,'p')
 
     umaxl(0) = max(umaxl(0),maxval(dabs(u1r(1:mgalx,1:mgalz)-s*y(j))))
     umaxl(1) = max(umaxl(1),maxval(dabs(u1r(1:mgalx,1:mgalz)))) 
     umaxl(2) = max(umaxl(2),maxval(dabs(u2r(1:mgalx,1:mgalz))))
     umaxl(3) = max(umaxl(3),maxval(dabs(u3r(1:mgalx,1:mgalz))))

  enddo

  deallocate(s12,s23,s13)
  deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz,ttt,ddd)
  deallocate(uxc,uyc,uzc,vxc,vyc,vzc,wxc,wyc,wzc)

  call MPI_ALLREDUCE(umaxl,umax,4,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  if (myid.eq.0) write(*,*) ' up,u,v,wmax check: ', umax(0:3)

end subroutine hvhg_uvwp


! --------------- read and write by HDF5 -------------------
subroutine writefield_ppp_hdf5(vor,wkr,ext)

  ! write instantaneous field 'vor' (real*4)

  use ctes
  use bcs
  use running 
  use hdf5
  use ppp,only:ishear

  implicit none
  include "mpif.h"
  !include "hdf5.h"

  integer i,j,k, iproc, ipo, leng, jpl
  real*4 wkr(mgalx+2,mgalz,jb:je)
  real*4 vor(mgalx+2,mgalz,jb:je)
  real*8 write_time
  character*256 filename
  integer istat(MPI_STATUS_SIZE),ierr
  integer, dimension(4):: mpe
  integer, dimension(3):: me
  real*8, dimension(7):: dum
  real*8, dimension(4):: buffsh

  real*8  xgrid(mgalx),ygrid(my),zgrid(mgalz) ! scaled by Lz

  character*3 ext
  character*3 ext1
  integer h5err

  ! --- hdf5 ---
  integer(hid_t):: fid, info_id
  integer(hsize_t), dimension(3):: dims, cdims, totaldims, start, offset
  integer(hsize_t), dimension(3):: dim0

  integer(hid_t):: dset,attr
  integer(hid_t):: dspace,mspace,aspace

  ! -----------

  jpl=je-jb+1;

  dim0 = (/1,0,0/)
  dims = (/ mgalx+2,mgalz,jpl /)   ! local data size on memory
  cdims = (/ mgalx,mgalz,jpl /)    ! local cropped date size (mgalx+2 => mgalx)
  totaldims = (/ mgalx,mgalz,my /) ! the size of output dataset (u,v,w ...)
  start = (/ 0, 0, 0/)
  offset = (/ 0, 0, 0/)

  wkr = vor
  iout = 99
  fid = iout

  write_time = MPI_WTIME()

  if (myid.eq.0) then
     !
     if (ishear.eq.1) then
        ext1='_sh'
        filename = filout(1:index(filout,' ')-1)//ext(1:3)//ext1(1:3)//'.h5'
     else
        filename = filout(1:index(filout,' ')-1)//ext(1:3)//'.h5'
     end if
     write(*,*) 'writing : ',trim(filename)
     !open (iout,file=filename,status='unknown',form='unformatted')
     !rewind(iout)
     call h5fcreate_f(filename,H5F_ACC_TRUNC_F,fid,h5err)
     !  header 
     !write(iout) mgalx,my,mgalz,numerop               ! mpe(1:4) 
     mpe=(/mgalx,my,mgalz,numerop/)
     !write(iout) time,Re,alp,Ly/pi,gam,s,chi,mx,my,mz ! dum(1:7)  
     dum=(/time,Re,alp,Ly/pi,gam,s,chi/)
     me=(/mx,my,mz/)
     !write(iout) y                                    ! ye(1:my)
     !!!write(iout) xoffb,xofft,zoffb,zofft
     !write(iout) xwkb,xwkt,zwkb,zwkt                  ! buff(1:4)
     buffsh=(/xwkb,xwkt,zwkb,zwkt/)
     !
     dim0(1)=4;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'mpe',H5T_STD_I32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_INTEGER,mpe,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     dim0(1)=3;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'me',H5T_STD_I32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_INTEGER,me,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     dim0(1)=7;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'dum',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,dum,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     dim0(1)=4;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'buff',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,buffsh,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     dim0(1)=1;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'Axz',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,Lx/Lz,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     dim0(1)=1;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'Ayz',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,Ly/Lz,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     dim0(1)=1;
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'Rez',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,Re*Lz*Lz,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     dim0(1)=1; ! add info. of the revision number of the code
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'shear_revision',H5T_STD_I32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_INTEGER,irev,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     ! write the scaled grid x,y,z 
     dim0(1)=mgalx
     do i=1,mgalx
        xgrid(i)= dfloat(i-1)*dx/Lz
     end do
!!$     !xp=2.0*pi/alp*dfloat(0:mgalx-1)/dfloat(mgalx);
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'xgrid',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,xgrid,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     dim0(1)=my
     ygrid=y/Lz
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'ygrid',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ygrid,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     dim0(1)=mgalz
     do k=1,mgalz
        zgrid(k)= dfloat(k-1)*dz/Lz
     end do
!!$     !zp=2.0*pi/gam*dfloat(0:mgalz-1)/dfloat(mgalz);
     call h5screate_simple_f(1,dim0,dspace,h5err)
     call h5dcreate_f(fid,'zgrid',H5T_IEEE_F32LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,zgrid,dim0,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
     !
     ! --- main data ----------------------------
     call h5screate_simple_f(3,totaldims,dspace,h5err)
     call h5dcreate_f(fid,ext(2:3),H5T_IEEE_F32LE,dspace,dset,h5err)
     !------the master first writes its stuff,  
     !      then receives from everybody, and writes it
     !
     start(3)=0
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*(mgalx+2)*mgalz
           call MPI_RECV(wkr,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)  
        endif
        !  do j=0,jend(iproc)-jbeg(iproc)
        !!!write(iout) ((wkr(i,k,j), i=1,mgalx),k=1,mgalz)
        !  enddo
        cdims(3) = jend(iproc)-jbeg(iproc)+1 ! fixed dims -> cdims
        ! dims should be the maxjsize
        ! Create the local data set
        call h5screate_simple_f(3,dims,mspace,h5err)
        !write(*,*) myid,'dcreate_f:',h5err
        call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F, & 
             offset,cdims,h5err)
        !write(*,*) myid,'sselect mspace:',h5err
        start(3)=jbeg(iproc)
        !start(3)=jbeg(iproc)-1 ! check it
        call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F, & 
             offset+start,cdims,h5err)
        !write(*,*) myid,'sselect dspace:',h5err
        !
        ! Commit the memspace to the disk
        call h5dwrite_f(dset,H5T_NATIVE_REAL,wkr,dims,h5err, & 
             mspace,dspace,H5P_DEFAULT_F)
        call h5fflush_f(fid,H5F_SCOPE_LOCAL_F,ierr)
        call h5sclose_f(mspace,h5err)
        
     enddo

     !Close datasets and global dataspace
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)        
     call h5fclose_f(fid,h5err)

  else             ! --- everybody else sends things to master

    call MPI_SEND(vor,(mgalx+2)*mgalz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)

  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write_time =  MPI_WTIME() - write_time
  if (myid.eq.0) write(*,*) '  ---  --- write_time = ',write_time

end subroutine writefield_ppp_hdf5

