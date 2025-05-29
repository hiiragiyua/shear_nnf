
module xp_check
  real*8, allocatable:: xp(:,:)  
  real*8, allocatable:: xref(:)
  real*8, allocatable:: diff_xp(:),time_xp(:)

end module xp_check

program diff_period

  ! convert omega2, lap.v ==> velocities and vorticities (real*4)

  use ctes
  use running
  use statistics
  use xp_check
  !use hdf5
  use gmres

  implicit none
  include "mpif.h"

  integer iproc,itags,newtag,imess,i,j,ibase,ii,ind,ixp
  real*8  val, hard, xrnorm,xpnorm,sca(4), sca_l(4)
  real*8, allocatable:: vor(:),phi(:)
  real*8, allocatable:: u00(:),w00(:)
  real*8, allocatable:: hv(:),hg(:),dvordy(:),chwk(:)
  integer istat(MPI_STATUS_SIZE),ierr,kmax
  !character*3 ext
  character*4 ext
  character*80 fname

  !  /*   initializes everything    */
  
  call MPI_INIT(ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numerop,ierr)

  if(myid.eq.0) write(*,*) ' %--- Diffp for finding initial upo by A.S., Version: rev.', irev

  !--------------- initializes commons and things
  call initcr(0)
  call set_options()
  ! --------------- allocates buffers
  allocate(vor(buffsize), phi(buffsize), chwk(buffsize), u00(my), w00(my) )
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
  if(myid.eq.0) write(*,*) 'iuse_v,iuse_us=',iuse_v,iuse_us
  if(myid.eq.0) write(*,*) 'iuse_ffty=',iuse_ffty
  if(myid.eq.0) write(*,*) 'iset_bulk0=',iset_bulk0

  iuse_fou_mode=2;
  if (myid.eq.0) then
     nmax=mx*mz*my*2 - 2*my*2 + my*2 
  else
     nmax=1
  end if

  if (myid.eq.0) then
     ! initialize input option
     ! Note: do not forget to set the number of elements in MPI_BCAST     
     write(*,*) 'diff to find initial condition'
     write(*,*) 'base file name to read', &
          ' is from hre2.dat [filinp,filout]'
     read(*,'(a)') filout
     write(*,*) 'read file index (from var.***, to var.*** )'
     read(*,*) istart
     read(*,*) iend
     write(*,*) 'maximum file-index intervals'
     read(*,*) kmax
     write(*,*) ' output filename including path, automatically added .diffp'
     read(*,'(a)') fname
     write(*,*) '(0: No  1: OK)'
     read(*,*) nosvemos
  end if

  call  MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(kmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call  MPI_BCAST(nosvemos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (nosvemos.ne.1) then
     iend=istart-1 ! stop
  end if

  if (myid.eq.0) then
     ! for checking a file
     fname = fname(1:index(fname,' ')-1)//'.diffp'
     open(iout,file=fname,status='unknown',form='unformatted')
     write(iout) kmax,istart,iend
  end if

  allocate(xp(nmax,kmax))
  xp=0.d0
  allocate(xref(nmax),xtmp(nmax))
  xref=0.d0; xtmp=0.d0
  allocate(diff_xp(kmax),time_xp(kmax))
  diff_xp=0.d0;
  ind=0
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

     !write(ext,'(i3.3)') istep
     write(ext,'(i4.4)') istep

     filinp =  filout(1:index(filout,' ')-1)//'.'//ext(1:4)

     !if(myid.eq.0) write(*,*) 'go to read file ...'
     call getfil(vor,phi,u00,w00,chwk,time)  ! read i.c.
     !if(myid.eq.0) write(*,*) '... file read'
     !
     ! 
     sca_l(3)=sum(vor*vor)
     
     if (iuse_v.eq.1) then 
       
        call getv(phi,hv,1) ! phi -> v
        sca_l(4)=sum(hv*hv)

     else
        sca_l(4)=sum(phi*phi)
     end if

     call MPI_ALLREDUCE(sca_l(3:4),sca(3:4),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     !
     sca_u00 = sqrt(sum(u00*u00))
     sca_w00 = sca_u00
     sca_vor = sqrt(sca(3))
     sca_phi = sqrt(sca(4))
     
     !sca_u00=1.d0; sca_w00=1.d0
     !sca_vor=1.d0; sca_phi=1.d0

     sca_u00=s*Lz; sca_w00=s*Lz
     sca_vor=s*sqrt(s*re*Lz*Lz)/2.5d0; ! set effective Reynolds number for LES 
     sca_phi=s*Lz    
  
     if (myid.eq.0) then
        write(*,*) ' sca_u00,w00,vor,phi:', sca_u00,sca_w00,sca_vor,sca_phi
        write(*,*) ' sca_time,sx,sz,re:', sca_time,sca_sx,sca_sy,sca_sz,sca_Ly,sca_re 
     endif

     !write(*,*) myid,'throwing',nmax
     call throwing(vor,phi,u00,w00,chwk,xref,1)
     !
     !write(*,*) myid,'set 1'     
     ibase=mod((istep-istart),kmax)
     !
     !write(*,*) myid,'set diff,nmax,time,ibase,istep',nmax,time,ibase,istep
     if ((istep-istart+1).le.kmax) then        

        xp(:,ibase+1)=xref
        time_xp(ibase+1)=time

     else

        !write(*,*) myid,'get norm'
        call norm(nmax,xp(1,ibase+1),xpnorm)
        diff_xp=0.d0
        do ii=0,kmax-1
           ixp=mod(ibase+ii,kmax)+1
           xtmp = xp(:,ixp) - xp(:,ibase+1)
           call norm(nmax,xtmp,diff_xp(ii+1))
        end do
        
        if(myid.eq.0) then
           !write(10,'(I6,1X,F14.6,32(1X,E14.6) )') &
           write(iout) &
                dfloat(istep-1),time_xp(ibase+1),xpnorm,(diff_xp(ii),ii=1,kmax)
        end if
        xp(:,ibase+1)=xref; time_xp(ibase+1)=time
     end if
     
  end do
  
  if (myid.eq.0) then
     close(iout)
  end if

  !----------------------------------------------------------------------- 
  ! /* finalize procedure      */
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)   !! used to be statement 200

end program diff_period

