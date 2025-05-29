
module xp_check
  real*8, allocatable:: xp(:,:)  
  real*8, allocatable:: xref(:,:)
  real*8, allocatable:: diff_xp(:,:,:),time_xp(:),diff_ph(:,:)
  real*8, allocatable:: phase_xp(:,:),phase_ref(:,:)
end module xp_check

program diff_upo_phase

  ! convert omega2, lap.v ==> velocities and vorticities (real*4)

  use ctes
  use running
  use statistics
  use xp_check
  use LES
  !use hdf5
  use gmres

  implicit none
  include "mpif.h"

  integer iproc,itags,newtag,imess,i,j,ibase,ii,ind,ixp,ixrefp,iTp,mTp
  integer nbinx, nbinz, iphx,iphz, ix,iz
  real*8  val, hard, xrnorm,xpnorm,sca(4), sca_l(4), fac
  real*8, allocatable:: vor(:),phi(:)
  real*8, allocatable:: u00(:),w00(:),u00wk(:),w00wk(:),dump(:),dumc(:)
  real*8, allocatable:: hv(:),hg(:),dvordy(:),chwk(:)
  integer istat(MPI_STATUS_SIZE),ierr,kmax,iph,ist_u,ien_u,idigit_u,idigit

  character*3 ext
  character*4 ext4
  character*128 fname_base,fname1,fname2,upofile

  !  /*   initializes everything    */
  
  call MPI_INIT(ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)

  if(myid.eq.0) write(*,*) ' %--- Diffp for finding initial upo by A.S., Version: rev.', irev

  !--------------- initializes commons and things
  call initcr(0)
  call set_options()
  iread_footer=0;

  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my), dump(my+2), dumc(my))
  vor = 0.d0; phi = 0.d0; chwk = 0.d0
  u00 = 0.d0; w00 = 0.d0
  !  -------  allocate rest of buffers 
  allocate( hv(buffsize),hg(buffsize),dvordy(buffsize))
  hv = 0.d0; hg = 0.d0; dvordy = 0.d0;  
  !
  iuse_v=1
  iuse_us=1
  iuse_ffty=0
  iset_bulk0=0;
  ifilter=0; ! not implemented...
  iydealiasing=0; ! this does not change the norm...??

  if(myid.eq.0) write(*,*) 'iuse_v,iuse_us=',iuse_v,iuse_us
  if(myid.eq.0) write(*,*) 'iuse_ffty=',iuse_ffty
  if(myid.eq.0) write(*,*) 'iset_bulk0=',iset_bulk0
  if(myid.eq.0) write(*,*) 'ifilter,iydealiasing',ifilter,iydealiasing
  iuse_fou_mode=2;
  if (myid.eq.0) then
     !nmax=mx*mz*my*2 - 2*my*2 + my*2 ! this is buggy, fixed (rev.1529) at 2015/Sep/16 
     nmax=mx*mz*my*2 - 2*my*2 + my*2 - (mz1-nz)*2*my*2      
     nmax1=nmax
  else
     nmax=1
     nmax1=1;
  end if
  nf=nmax1; nl=nf;

  if (myid.eq.0) then
     ! initialize input option
     ! Note: do not forget to set the number of elements in MPI_BCAST     
     write(*,*) 'diff for revisiting to the reference UPO'
     write(*,*) 'base file name of turbulence to read'
     read(*,'(a)') filout
     write(*,*) 'read file index (from var.***, to var.*** )'
     read(*,*) istart
     read(*,*) iend
     read(*,*) idigit
     write(*,*) 'refrence file name of UPO (assume that the phase of B.C. is the same with turb.'
     read(*,'(a)') upofile
     write(*,*) 'read file index (from upo.***, to upo.*** )'
     read(*,*) ist_u
     read(*,*) ien_u
     read(*,*) idigit_u
     write(*,*) ' output filename including path, automatically added .diffp'
     read(*,'(a)') fname_base
     write(*,*) '(0: No  1: OK)'
     read(*,*) nosvemos
  end if
  nbinx=30; nbinz=10;

  call  MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(idigit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call  MPI_BCAST(ist_u,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(ien_u,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(idigit_u,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call  MPI_BCAST(nosvemos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (nosvemos.ne.1) then
     iend=istart-1 ! stop
  end if

  kmax=ien_u-ist_u+1
  mTp=kmax/10;

  if ((mod(kmax,10).ne.0).and.(iuse_LES.eq.0)) then
     if (myid.eq.0) write(*,*) 'kmax error ... kmax/10 should be interger... stop'
     stop
  elseif (iuse_LES.eq.1) then
     if (myid.eq.0) write(*,*) 'LES-EQ mode.'
     mTp=1;
  end if


  if (myid.eq.0) then
     ! for checking a file
     fname1 = fname_base(1:index(fname_base,' ')-1)//'.diffuph'
     open(iout,file=fname1,status='unknown',form='unformatted')
     write(iout) mTp,istart,iend
  
!     fname2 = fname1(1:index(fname1,' ')-1)//'_phasex'
!     open(34,file=fname2,status='unknown',form='unformatted')
!     write(34) mTp,istart,iend

!     fname2 = fname1(1:index(fname1,' ')-1)//'_phasey'
!     open(35,file=fname2,status='unknown',form='unformatted')
!     write(35) mTp,istart,iend

!     fname2 = fname1(1:index(fname1,' ')-1)//'_phasez'
!     open(36,file=fname2,status='unknown',form='unformatted')
!     write(36) mTp,istart,iend

  end if

  sca_u00=s*Lz; sca_w00=s*Lz
  if (iuse_v.eq.1) then
     sca_vor=s; sca_phi=s*Lz ! here iuse_v =1 is assumed
     if (iuse_compensatory_norm.eq.1) then
        if (myid.eq.0) write(*,*) 'compensatory norm for vor and v',comp_norm
        sca_vor=s*sqrt(s*re*Lz*Lz)/comp_norm; sca_phi=s*Lz ! here iuse_v =1 is assumed
     end if
  else
     sca_vor=s; sca_phi=1.d0  ! the scaling for phi is not implemented yet.
  endif
  
  allocate(xp(nmax,kmax),xref(nmax,kmax))
  xp=0.d0; xref=0.d0; 
  allocate(phase_xp(3,kmax),phase_ref(3,kmax));
  phase_xp=0.d0;
  phase_ref=0.d0;
  allocate(xtmp(nmax))
  xtmp=0.d0
  allocate(diff_xp(mTp,nbinx,nbinz),time_xp(kmax))
  diff_xp=0.d0;
  allocate(diff_ph(3,mTp))
  diff_ph=0.d0;
  ind=0

  do iphz=1,nbinz;! consider phase-shift in x, the increments is fixed: 0.1*Lz (nbinz=10)
     do iphx=1,nbinx; ! consider phase-shift in x, the increments is fixed: 0.1*Lx (nbinx=10*Axz)
        phasex=dfloat(iphx-1)/dfloat(nbinx);
        phasez=dfloat(iphz-1)/dfloat(nbinz);
        ! ---- read the reference UPO time series ----
        do iph=ist_u,ien_u
           ! reinitialize the buffers, (rev565)
           vor = 0.d0; phi = 0.d0; chwk = 0.d0
           u00 = 0.d0; w00 = 0.d0
           ! the shuffled buffers by change are dangerous to reuse
           hv = 0.d0;      
           hg = 0.d0; 
           dvordy = 0.d0;
           if (idigit_u.eq.3) then
              write(ext,'(i3.3)') iph
              filinp =  upofile(1:index(upofile,' ')-1)//'.'//ext(1:3)
           else
              write(ext4,'(i4.4)') iph
              filinp =  upofile(1:index(upofile,' ')-1)//'.'//ext4(1:4)
           end if

           call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
           !call dealiasing_y(vor,phi,u00,w00,dump,dump,dumc,dumc,fac)
           !call get_phase(vor,phi,u00,w00,phasex,phasey,phasez)
           !call phase_shift_fou(vor,phi,atan(phasex)/(2*pi)*Lx,atan(phasez)/(2*pi)*Lz)
           !phase_ref(1,iph) = phasex
           !phase_ref(2,iph) = phasey
           !phase_ref(3,iph) = phasez

           !call phase_shift_fou(vor,phi,atan(phase_ref(1,1))/(2*pi)*Lx,atan(phase_ref(3,1))/(2*pi)*Lz)
           call phase_shift_fou(vor,phi,phasex*Lx,phasez*Lz)
           call throwing(vor,phi,u00,w00,chwk,xref(1,iph),1)

        end do
  
        if (myid.eq.0) write(*,*) '=== read turb fields',iphz,'/',nbinz,' ===='
        ! ---- read turbulence time series ---- 
        ind = 0
        do istep = istart,iend
           ind = ind +1
           ! reinitialize the buffers, (rev565)
           vor = 0.d0; phi = 0.d0; chwk = 0.d0
           u00 = 0.d0; w00 = 0.d0
           !
           ! the shuffled buffers by change are dangerous to reuse
           hv = 0.d0;      
           hg = 0.d0;   ! only hg should be reset (rev565)
           dvordy = 0.d0;
           if (idigit.eq.3) then
              write(ext,'(i3.3)') istep
              filinp =  filout(1:index(filout,' ')-1)//'.'//ext(1:3)
           else
              write(ext4,'(i4.4)') istep
              filinp =  filout(1:index(filout,' ')-1)//'.'//ext4(1:4)
           end if
           
           !if(myid.eq.0) write(*,*) 'go to read file ...'
           call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
           !call dealiasing_y(vor,phi,u00,w00,dump,dump,dumc,dumc,fac)
           !if(myid.eq.0) write(*,*) '... file read'
           ! 
           !write(*,*) myid,'throwing',nmax
           !if (mod(istep-istart,kmax).eq.0) then
           !   call get_phase(vor,phi,u00,w00,phasex,phasey,phasez)
           !   !call phase_shift_fou(vor,phi,atan(phasex)/(2*pi)*Lx,atan(phasez)/(2*pi)*Lz )
           !end if
           !if (myid.eq.0) write(34,*) phasex,phasey,phasez
           !
           !write(*,*) myid,'set diff,nmax,time,ibase,istep',nmax,time,ibase,istep
           if ((istep-istart+1).le.kmax) then 
              !phase_xp(1,istep)=phasex
              !phase_xp(2,istep)=phasey
              !phase_xp(3,istep)=phasez
              call throwing(vor,phi,u00,w00,chwk,xp(1,ind),1)
              time_xp(ind)=time
           else
              ibase=mod((istep-istart),kmax)
              ! istep=41, istart=1 ==> ibase=0
              !write(*,*) myid,'get norm'
              if (ifilter.eq.1) then
                 call scaling(xp(1,ibase+1),xtmp,1)
                 call norm(nmax,xtmp,xpnorm)
              else
                 call norm2(nmax,xp(1,ibase+1),xpnorm)
              end if
              diff_xp=0.d0
              diff_ph=0.d0
              do iTp=0,mTp-1  ! consider phase in time (only for the multiple of STs) 
                 ! /* the time integration of L2 norm of distance (X_turb - X_upo) */
                 do ii=0,kmax-1 
                    ixp=mod(ibase+ii,kmax)+1
                    ixrefp=mod(ibase+ii+iTp*10,kmax)+1
                    xtmp = xp(:,ixp) - xref(:,ixrefp)
                    if (ifilter.eq.1) then
                       call scaling(xtmp,xtmp,1)
                       call norm(nmax,xtmp,xpnorm)
                    else
                       call norm2(nmax,xtmp,xrnorm)
                    end if
                    diff_xp(iTp+1,iphx,iphz) = diff_xp(iTp+1,iphx,iphz) + xrnorm 
                    !atan(phasex)/(2*pi)
                    !diff_ph(1,iTp+1) = diff_ph(1,iTp+1) + (phase_xp(1,ixp) - phase_ref(1,ixrefp))**2
                    !diff_ph(2,iTp+1) = diff_ph(2,iTp+1) + (phase_xp(2,ixp) - phase_ref(2,ixrefp))**2
                    !diff_ph(3,iTp+1) = diff_ph(3,iTp+1) + (phase_xp(3,ixp) - phase_ref(3,ixrefp))**2
                 end do
                 diff_xp(iTp+1,iphx,iphz) = diff_xp(iTp+1,iphx,iphz)/dfloat(kmax)
                 !diff_ph(1,iTp+1) = diff_ph(1,iTp+1)/dfloat(kmax)
                 !diff_ph(2,iTp+1) = diff_ph(2,iTp+1)/dfloat(kmax)
                 !diff_ph(3,iTp+1) = diff_ph(3,iTp+1)/dfloat(kmax)
              end do
              if(myid.eq.0) then
                 write(iout) dfloat(ind),time_xp(ibase+1),xpnorm, & 
                      (diff_xp(ii,iphx,iphz),ii=1,mTp)
              end if
              !phase_xp(1,ibase+1)=phasex
              !phase_xp(2,ibase+1)=phasey
              !phase_xp(3,ibase+1)=phasez
              call throwing(vor,phi,u00,w00,chwk,xp(1,ibase+1),1)
              time_xp(ibase+1)=time
           end if
           
        end do
 
     end do; 
  end do; ! phase xbin,zbin loop
  
  if (myid.eq.0) then
     close(iout)
     !close(34)
     !close(35)
     !close(36)
  end if
     
  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)   !! used to be statement 200

end program diff_upo_phase

