subroutine getfil_hdf5(vor,phi,u00,w00,wk,tiempo)
  use ctes
  use running
  use bcs
  use hdf5
 
  implicit none
  include "mpif.h"
  
  integer i,j,k
  integer istat(MPI_STATUS_SIZE)
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  real*8  vor(mx,mz,jb:je),phi(mx,mz,jb:je),wk(*)
  real*8  u00(my),w00(my)
  
  real*4, allocatable:: vore(:,:),phie(:,:)
 
  integer  mgalxe,mye,mgalze,numerope
  !real*8   Ree,alpe,Lye,game,se,chie
  real*8   ye(my)
  real*8   tiempo,t1,t2
  integer  mpee(6)
  real*8   dume(7),buffe(4)
  real*8,  allocatable:: u00e(:),w00e(:),yee(:)
 
  integer(hsize_t) dimf(3),dimh1(1),dimh2(1),dimh3(1),dimh4(1),dimh5(1),dimh6(1)
  
  character(3), parameter :: dsetname1 = "vor"
  character(3), parameter :: dsetname2 = "phi"
  character(3), parameter :: dsetname3 = "mpe"
  character(3), parameter :: dsetname4 = "dum"
  character(4), parameter :: dsetname5 = "buff"
  character(1), parameter :: dsetname6 = "y"
  character(3), parameter :: dsetname7 = "u00"
  character(3), parameter :: dsetname8 = "w00"
  
  integer(hid_t) plist_id, file_id
  integer(hid_t) fspace1, fspace2, fspace3, fspace4, fspace5, fspace6, fspace7, fspace8
  integer(hid_t) mspace1, mspace2, mspace3, mspace4, mspace5, mspace6, mspace7, mspace8
  integer(hid_t) dset_id1,dset_id2,dset_id3,dset_id4,dset_id5,dset_id6,dset_id7,dset_id8
 
  integer(hsize_t)  countf(3), counth(1)
  integer(hssize_t) offset(3)
 
  integer, parameter :: rankf = 3
  integer, parameter :: rankh = 1
  integer  error,ierr
  
  dimf  = (/mx,mz,my/)
  dimh1 = (/6/)
  dimh2 = (/7/)
  dimh3 = (/4/)
  dimh4 = (/my/)
  dimh5 = (/my/)
  dimh6 = (/my/)

  t1= MPI_WTIME()

  if (myid.eq.0) then
     write(*,*) 'reading:',trim(filinp)
     if ((filinp(index(filinp,' ')-3:index(filinp,' ')-1)).ne.'.h5') then
        write(*,*) 'reading HDF5, but this is not hdf5 data, set reading option correctly'
        stop
     end if
     call h5fopen_f(filinp,H5F_ACC_RDONLY_F,file_id,error,H5P_DEFAULT_F)
     if (error.ge.1) stop
 
     call h5dopen_f(file_id, dsetname1, dset_id1, error)
     call h5dopen_f(file_id, dsetname2, dset_id2, error)
     call h5dopen_f(file_id, dsetname3, dset_id3, error)
     call h5dopen_f(file_id, dsetname4, dset_id4, error)
     call h5dopen_f(file_id, dsetname5, dset_id5, error)
     call h5dopen_f(file_id, dsetname6, dset_id6, error)
     call h5dopen_f(file_id, dsetname7, dset_id7, error)
     call h5dopen_f(file_id, dsetname8, dset_id8, error)
  
     call h5dread_f(dset_id3,H5T_NATIVE_INTEGER,mpee,dimh1,error)
     mgalxe=mpee(1);mgalze=mpee(3)
     me(1) =mpee(5);me(2)=mpee(2);me(3)=mpee(6)
     if((mgalx.ne.mgalxe).or.(my.ne.me(2)).or.(mgalz.ne.mgalze)) then
        write(*,*) ' different physical grid points!'
        write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,me(2),mgalze,&
                                                    ' =>', mgalx, my, mgalz
     endif
     call h5dread_f(dset_id4,H5T_NATIVE_DOUBLE,dume,dimh2,error)
     tiempo=dume(1);Ree=dume(2);alpe=dume(3);Lye=dume(4)*pi
     game  =dume(5);se =dume(6);chie=dume(7)
     !write(*,*) 'dume',dume, Ree, alpe, Lye, game, se, chie
     
     call h5dread_f(dset_id5, H5T_NATIVE_DOUBLE, buffe,dimh3, error)
     xoffb =buffe(1); xofft=buffe(2); zoffb=buffe(3); zofft=buffe(4)
     if (alpe.ne.alp) then
        xoffb=xoffb*alpe/alp
        xofft=xofft*alpe/alp
        write(*,*) ' rescaling shear-periodic b.c. offset: factor =', &
                    alpe/alp, '  xoffb,xofft = ',xoffb,xofft
     endif
     write(*,*)
     write(*,*) 'reading input file ..., time=',tiempo
     write(*,'(a,3(1X,I5))') '                    ...,(mx,my,mz)=',me
     write(*,*)
 
     allocate(yee(me(2)))
     call h5dread_f(dset_id6, H5T_NATIVE_DOUBLE, yee,  dimh4, error)
     ye(1:me(2))=yee

     if (me(2).ne.my) then
        write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' =>my=',my
     endif
     allocate(u00e(me(2)),w00e(me(2)))
     call h5dread_f(dset_id7, H5T_NATIVE_DOUBLE, u00e, dimh5, error)
     call h5dread_f(dset_id8, H5T_NATIVE_DOUBLE, w00e, dimh6, error)

     mym = min(my,me(2))
     u00(1:mym) = u00e(1:mym)
     w00(1:mym) = w00e(1:mym)
 
  endif
 
  call MPI_BCAST(me,3,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(ye,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(u00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(w00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  
  ntotr = me(1)*me(3)
  mxm   = min(mx,me(1))
  mym   = min(my,me(2))
  mzm   = min(mz,me(3))
 
  klen  = (mzm+1)/2
  kini1 = me(3) - klen + 1
  kini2 = mz    - klen + 1
 
  allocate(vore(me(1),me(3)),phie(me(1),me(3)))
  vore=0.0;phie=0.0
 
  if (myid.eq.0) then
 
     do iproc=0,numerop-1
        do j=jbeg(iproc),jend(iproc)
           countf(1) = me(1)
           countf(2) = me(3)
           countf(3) = 1
           offset(1) = 0
           offset(2) = 0
           offset(3) = j
           if (j.le.(me(2)-1)) then
              call h5dget_space_f(dset_id1,fspace1,error)
              call h5sselect_hyperslab_f(fspace1,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace1, error)

              call h5dget_space_f(dset_id2,fspace2,error)
              call h5sselect_hyperslab_f(fspace2,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace2, error) 
     
              call h5dread_f(dset_id1,H5T_NATIVE_REAL,vore,countf,error,mspace1,fspace1)
              call h5dread_f(dset_id2,H5T_NATIVE_REAL,phie,countf,error,mspace2,fspace2)
           else
              vore=0.0
              phie=0.0
           endif

           if (iproc.eq.0) then
              vor(1:mxm,1:klen,j)    =real(vore(1:mxm,1:klen),kind=8)
              vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm,1+kini1:me(3)),kind=8)
              phi(1:mxm,1:klen,j)    =real(phie(1:mxm,1:klen),kind=8)
              phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm,1+kini1:me(3)),kind=8)
           else
              call MPI_SEND(vore,ntotr,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,error)
              call MPI_SEND(phie,ntotr,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,error)
           endif
 
        enddo
     enddo
     call h5sclose_f(fspace1,error)
     call h5sclose_f(mspace1,error)
  
     call h5sclose_f(fspace2,error)
     call h5sclose_f(mspace2,error)

     call h5dclose_f(dset_id1,error)
     call h5dclose_f(dset_id2,error)
     call h5fclose_f(file_id,error)
  else
     do j=jb,je
        call MPI_RECV(vore,ntotr,MPI_REAL4,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,error)
        call MPI_RECV(phie,ntotr,MPI_REAL4,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,error)
        vor(1:mxm,1:klen,j)    =real(vore(1:mxm,1:klen),kind=8)
        vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm,1+kini1:me(3)),kind=8)
        phi(1:mxm,1:klen,j)    =real(phie(1:mxm,1:klen),kind=8)
        phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm,1+kini1:me(3)),kind=8)
     enddo
 
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,error)
 
  call  MPI_BCAST(Ree,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(Lye,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! BCAST are forgotten (2013/10/14) alpe,Lye,game,se,chie,me
  call  MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(game,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(se,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(chie,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  t2= MPI_WTIME()
  if (myid.eq.0) write(*,*) '  --- time: read HDF5 =',(t2-t1)

  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)
 
  if (me(2).lt.my) then
     call interp_y(vor,phi,u00,w00,me,ye)
  elseif (me(2).gt.my) then
     if (myid.eq.0)  write(*,*) 'Sorry, it does not support my less than me(2),:('
     stop
  endif
  if(myid.eq.0) deallocate(yee,u00e,w00e)
  deallocate(vore,phie)

  t2= MPI_WTIME()
  if (myid.eq.0) write(*,*) '  --- time: read HDF5, including comm.=',(t2-t1)

end subroutine getfil_hdf5

subroutine getfil_hdf5_convert(vor,phi,u00,w00,wk,tiempo)
  use ctes
  use running
  use bcs
  use hdf5
 
  implicit none
  include "mpif.h"
  
  integer i,j,k
  integer istat(MPI_STATUS_SIZE)
  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  real*8  vor(mx,mz,jb:je),phi(mx,mz,jb:je),wk(*)
  real*8  u00(my),w00(my)
  real*4, allocatable:: vore(:,:),phie(:,:)
 
  integer  mgalxe,mye,mgalze,numerope
  !real*8   Ree,alpe,Lye,game,se,chie
  real*8   ye(my)
  real*8   tiempo
  integer  mpee(6)
  real*8   dume(7),buffe(4)
  real*8,  allocatable:: u00e(:),w00e(:),yee(:)
 
  integer(hsize_t) dimf(3),dimh1(1),dimh2(1),dimh3(1),dimh4(1),dimh5(1),dimh6(1)
  
  character(3), parameter :: dsetname1 = "vor"
  character(3), parameter :: dsetname2 = "phi"
  character(3), parameter :: dsetname3 = "mpe"
  character(3), parameter :: dsetname4 = "dum"
  character(4), parameter :: dsetname5 = "buff"
  character(1), parameter :: dsetname6 = "y"
  character(3), parameter :: dsetname7 = "u00"
  character(3), parameter :: dsetname8 = "w00"
  
  integer(hid_t) plist_id, file_id
  integer(hid_t) fspace1, fspace2, fspace3, fspace4, fspace5, fspace6, fspace7, fspace8
  integer(hid_t) mspace1, mspace2, mspace3, mspace4, mspace5, mspace6, mspace7, mspace8
  integer(hid_t) dset_id1,dset_id2,dset_id3,dset_id4,dset_id5,dset_id6,dset_id7,dset_id8
 
  integer(hsize_t)  countf(3), counth(1)
  integer(hssize_t) offset(3)
 
  integer, parameter :: rankf = 3
  integer, parameter :: rankh = 1
  integer  error,ierr
  
  dimf  = (/mx,mz,my/)
  dimh1 = (/6/)
  dimh2 = (/7/)
  dimh3 = (/4/)
  dimh4 = (/my/)
  dimh5 = (/my/)
  dimh6 = (/my/)

  if (myid.eq.0) then
     call h5fopen_f(filinp,H5F_ACC_RDWR_F,file_id,error)
 
     call h5dopen_f(file_id, dsetname1, dset_id1, error)
     call h5dopen_f(file_id, dsetname2, dset_id2, error)
     call h5dopen_f(file_id, dsetname3, dset_id3, error)
     call h5dopen_f(file_id, dsetname4, dset_id4, error)
     call h5dopen_f(file_id, dsetname5, dset_id5, error)
     call h5dopen_f(file_id, dsetname6, dset_id6, error)
     call h5dopen_f(file_id, dsetname7, dset_id7, error)
     call h5dopen_f(file_id, dsetname8, dset_id8, error)
  
     call h5dread_f(dset_id3,H5T_NATIVE_INTEGER,mpee,dimh1,error)
     mgalxe=mpee(1);mgalze=mpee(3)
     me(1) =mpee(5);me(2)=mpee(2);me(3)=mpee(6)
     if((mgalx.ne.mgalxe).or.(my.ne.me(2)).or.(mgalz.ne.mgalze)) then
        write(*,*) ' different physical grid points!'
        write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,me(2),mgalze,&
                                                    ' =>', mgalx, my, mgalz
     endif
     call h5dread_f(dset_id4,H5T_NATIVE_REAL,dume,dimh2,error)
     tiempo=real(dume(1),kind=8);Ree=real(dume(2),kind=8)
     alpe  =real(dume(3),kind=8);Lye=real(dume(4),kind=8)
     game  =real(dume(5),kind=8);se =real(dume(6),kind=8)
     chie  =real(dume(7),kind=8)
  
     call h5dread_f(dset_id5, H5T_NATIVE_REAL, buffe,dimh3, error)
     xoffb =real(buffe(1),kind=8); xofft=real(buffe(2),kind=8)
     zoffb =real(buffe(3),kind=8); zofft=real(buffe(4),kind=8)
     
     if (alpe.ne.alp) then
        xoffb=xoffb*alpe/alp
        xofft=xofft*alpe/alp
        write(*,*) ' rescaling shear-periodic b.c. offset: factor =', &
                    alpe/alp, '  xoffb,xofft = ',xoffb,xofft
     endif
     write(*,*)
     write(*,*) 'reading input file ..., time=',tiempo
     write(*,'(a,3(1X,I5))') '                    ...,(mx,my,mz)=',me
     write(*,*)
 
     allocate(yee(me(2)))
     call h5dread_f(dset_id6, H5T_NATIVE_REAL, yee,  dimh4, error)
     ye(1:me(2))=real(yee,kind=8)

     if (me(2).ne.my) then
        write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' =>my=',my
     endif
     allocate(u00e(me(2)),w00e(me(2)))
     call h5dread_f(dset_id7, H5T_NATIVE_REAL, u00e, dimh5, error)
     call h5dread_f(dset_id8, H5T_NATIVE_REAL, w00e, dimh6, error)

     mym = min(my,me(2))
     u00(1:mym) = real(u00e(1:mym),kind=8)
     w00(1:mym) = real(w00e(1:mym),kind=8)
 
  endif
 
  call MPI_BCAST(me,3,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(ye,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(u00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  !call MPI_BCAST(u00,my,MPI_REAL8,0,MPI_COMM_WORLD,error) ! bug 
  call MPI_BCAST(w00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  
  ntotr = me(1)*me(3)
  mxm   = min(mx,me(1))
  mym   = min(my,me(2))
  mzm   = min(mz,me(3))
 
  klen  = (mzm+1)/2
  kini1 = me(3) - klen + 1
  kini2 = mz    - klen + 1
 
  allocate(vore(me(1),me(3)),phie(me(1),me(3)))
  vore=0.0;phie=0.0
 
  if (myid.eq.0) then
 
     do iproc=0,numerop-1
        do j=jbeg(iproc),jend(iproc)
           countf(1) = me(1)
           countf(2) = me(3)
           countf(3) = 1
           offset(1) = 0
           offset(2) = 0
           offset(3) = j
           if (j.le.(me(2)-1)) then
              call h5dget_space_f(dset_id1,fspace1,error)
              call h5sselect_hyperslab_f(fspace1,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace1, error)

              call h5dget_space_f(dset_id2,fspace2,error)
              call h5sselect_hyperslab_f(fspace2,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace2, error) 
     
              call h5dread_f(dset_id1,H5T_NATIVE_REAL,vore,countf,error,mspace1,fspace1)
              call h5dread_f(dset_id2,H5T_NATIVE_REAL,phie,countf,error,mspace2,fspace2)
           else
              vore=0.0
              phie=0.0
           endif

           if (iproc.eq.0) then
              vor(1:mxm,1:klen,j)    =real(vore(1:mxm,1:klen),kind=8)
              vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm,1+kini1:me(3)),kind=8)
              phi(1:mxm,1:klen,j)    =real(phie(1:mxm,1:klen),kind=8)
              phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm,1+kini1:me(3)),kind=8)
           else
              call MPI_SEND(vore,ntotr,MPI_REAL4,iproc,iproc, &
                            MPI_COMM_WORLD,error)
              call MPI_SEND(phie,ntotr,MPI_REAL4,iproc,iproc, &
                            MPI_COMM_WORLD,error)
           endif

        enddo
     enddo
     call h5sclose_f(fspace1,error)
     call h5sclose_f(mspace1,error)
  
     call h5sclose_f(fspace2,error)
     call h5sclose_f(mspace2,error)

     call h5dclose_f(dset_id1,error)
     call h5dclose_f(dset_id2,error)

     call h5fclose_f(file_id,error)
  else
     do j=jb,je
        call MPI_RECV(vore,ntotr,MPI_REAL4,0,MPI_ANY_TAG, &
                      MPI_COMM_WORLD,istat,error)
        call MPI_RECV(phie,ntotr,MPI_REAL4,0,MPI_ANY_TAG, &
                      MPI_COMM_WORLD,istat,error)
        vor(1:mxm,1:klen,j)    =real(vore(1:mxm,1:klen),kind=8)
        vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm,1+kini1:me(3)),kind=8)
        phi(1:mxm,1:klen,j)    =real(phie(1:mxm,1:klen),kind=8)
        phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm,1+kini1:me(3)),kind=8)
     enddo
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,error)
 
  call  MPI_BCAST(Ree,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(Lye,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! BCAST are forgotten (2013/10/14) alpe,Lye,game,se,chie,me
  call  MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(game,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(se,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call  MPI_BCAST(chie,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)
  
  if (me(2).lt.my) then
     call interp_y(vor,phi,u00,w00,me,ye)
  elseif (me(2).gt.my) then
     if (myid.eq.0)  write(*,*) 'Sorry, it does not support my less than me(2),:('
     stop
  endif
  if(myid.eq.0) deallocate(yee,u00e,w00e)
  deallocate(vore,phie)

end subroutine getfil_hdf5_convert

subroutine writehdf5(vor,phi,u00,w00,wk1,wk2,wk)
  use ctes
  use running !,only:time,id22,filout
  use statistics
  use bcs
  use hdf5
  use LES
        
  implicit none  
  include 'mpif.h'

  real*8  phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,mx,kb:ke)
  real*8  wk1(mx,mz,jb:je), wk2(mx,mz,jb:je)
  !real*4, allocatable::  wk11(:,:,:),wk22(:,:,:)
  real*8  wk(*)
  real*8  u00(my),w00(my),dum(7),buff(4)
  integer mpe(6)
  real*8  w_time
  integer iproc,leng,ipo

  integer nopte,nparamse ! rev.1433
  integer ioptions(maxParams) ! 100 is maximum ... see in mod.f90
  real(8) params(maxParams)
  

  integer istat(MPI_STATUS_SIZE),ierr
  character ext*4, filnam*256, filst*256
  
  integer(hsize_t) dimf(3),dimh1(1),dimh2(1),dimh3(1),dimh4(1),dimh5(1),dimh6(1)

  character(3), parameter :: dsetname1 = "vor"  
  character(3), parameter :: dsetname2 = "phi" 
  character(3), parameter :: dsetname3 = "mpe"
  character(3), parameter :: dsetname4 = "dum"
  character(4), parameter :: dsetname5 = "buff"
  character(1), parameter :: dsetname6 = "y"
  character(3), parameter :: dsetname7 = "u00"
  character(3), parameter :: dsetname8 = "w00"

  integer(hid_t) plist_id, file_id, fid    
  integer(hid_t) fspace1, fspace2, fspace3, fspace4, fspace5, fspace6, fspace7, fspace8
  integer(hid_t) mspace1, mspace2, mspace3, mspace4, mspace5, mspace6, mspace7, mspace8
  integer(hid_t) dset_id1,dset_id2,dset_id3,dset_id4,dset_id5,dset_id6,dset_id7,dset_id8
    
  integer(hsize_t)  countf(3), counth(1)  
  integer(hssize_t) offset(3) 

  integer, parameter :: rankf = 3
  integer, parameter :: rankh = 1
  integer  error, h5err

  ! for hdf5
  integer(hid_t) dspace, mspace, dset
  integer(hsize_t), dimension(2):: dimstat 
  integer(hsize_t), dimension(3):: dimspec
  

  !allocate(wk11(mx,mz,jb:je),wk22(mx,mz,jb:je))  ! real*4 version by siwei
 
  call chjik2ikj(vor,wk1,wk,wk)  
  call chjik2ikj(phi,wk2,wk,wk)  
  !wk11=real(wk1,kind=4)
  !wk22=real(wk2,kind=4)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  mpe =(/mgalx,my,mgalz,numerop,mx,mz/)
  dum =(/time,Re,alp,Ly/pi,gam,s,chi/)
  buff=(/xwkb,xwkt,zwkb,zwkt/)

  dimf  = (/mx,mz,my/)
  dimh1 = (/6/)
  dimh2 = (/7/)
  dimh3 = (/4/)
  dimh4 = (/my/)
  dimh5 = (/my/)
  dimh6 = (/my/)

  if(id22.gt.9999.or.id22.lt.0) then
     write(*,*) 'number of images out of range',id22
     write(*,*) 'check ifile in hre.dat'
     ifatal=1 
     !call MPI_FINALIZE(ierr)
     stop
  end if


  !write(ext,'(i3.3)') id22
  write(ext,'(i4.4)') id22
  id22   = id22+1

  ! /*       get and write statistics first       */
  if (nacum.ne.0) then
     if (myid.eq.0) write(*,*) 'stat esc', nacum

     ! dt-multiplied statistics  
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

        filst = filout(1:index(filout,' ')-1)//'.'//ext//'.sta.h5' 
        call h5fcreate_f(filst,H5F_ACC_TRUNC_F,fid,h5err)
        if(h5err.ne.0) then
           write(*,*) h5err, "ERROR: Problem openning ", trim(filst)
           write(*,*) "overwrite the file"
           call h5eprint_f(h5err,trim(filst))
           !stop
        endif
        call h5write_simple_int(fid,'nacum',nacum,1,h5err)
        call h5write_simple_int(fid,'nhist',nhist,1,h5err)
        call h5write_simple_double(fid,'total_dt',total_dt,1,h5err)
        call h5write_simple_int(fid,'mgalx',mgalx,1,h5err)
        call h5write_simple_int(fid,'my',my,1,h5err)
        call h5write_simple_int(fid,'mgalz',mgalz,1,h5err)
        call h5write_simple_double(fid,'time',time,1,h5err)
        call h5write_simple_double(fid,'stime0',stime0,1,h5err)
        call h5write_simple_double(fid,'Re',Re,1,h5err)
        call h5write_simple_double(fid,'alp',alp,1,h5err)
        call h5write_simple_double(fid,'Ly',Ly,1,h5err) ! Ly/pi is saved in the flow fields
        call h5write_simple_double(fid,'gam',gam,1,h5err)
        call h5write_simple_double(fid,'s',s,1,h5err)
        call h5write_simple_double(fid,'chi',chi,1,h5err)
        call h5write_simple_double(fid,'y',y,my,h5err)

        !call h5write_simple_int(fid,'nstat',nstat,1,h5err)
        dimstat=(/nstat,my/)
        call h5screate_simple_f(2,dimstat,dspace,h5err)
        call h5dcreate_f(fid,'stats2',H5T_IEEE_F64LE,dspace,dset,h5err)
        call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,stats2,dimstat,h5err)
        call h5dclose_f(dset,h5err)
        call h5sclose_f(dspace,h5err)

     else

        call MPI_SEND(stats2(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)

     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta.h5 together
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl2,sp2,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

     if (myid.eq.0) then
        dimspec=(/nspec,mx1+1,mz1+1/)
        call h5screate_simple_f(3,dimspec,dspace,h5err)
        call h5dcreate_f(fid,'sp2',H5T_IEEE_F64LE,dspace,dset,h5err)
        call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,sp2,dimspec,h5err)
        call h5dclose_f(dset,h5err)
        call h5sclose_f(dspace,h5err)        

        !write(isn) CFL, uprim, vbulk, xforce ! rev934 ! xforce is redundunt, but keep the format
        call h5write_simple_double(fid,"CFL",CFLe,1,h5err)
        call h5write_simple_double(fid,"uprim",uprim,1,h5err)
        call h5write_simple_double(fid,"vbulk",vbulk,1,h5err)
        call h5write_simple_double(fid,"xforce",xforce,1,h5err)

        !write(isn) irev,nopt,nparams
        call h5write_simple_int(fid,"irev",irev,1,h5err)
        call h5write_simple_int(fid,"nopt",nopt,1,h5err)        
        call h5write_simple_int(fid,"nparams",nparams,1,h5err)
           
        ioptions(1)=explicit; ioptions(2)=ifix_dt; ioptions(3)=iadd_force; 
        ioptions(4)=iadd_mode; ioptions(5)=iadd_sym; ioptions(6)=iadd_damping;
        ioptions(7)=iuse_LES; ioptions(8)=iadd_visc; ioptions(9)=idynamic; 
        call h5write_simple_int(fid,"ioptions",ioptions,nopt,h5err) 

        params(1)=Deltat; params(2)=force_roll; params(3)=xforce; params(4)=zforce;
        params(5)=damp_aa; params(6)=damp_up; params(7)=Cles; params(8)=Deltag; 
        params(9)=cutoffx; params(9)=cutoffz; params(10:12)=LESflt(1:3) 
        call h5write_simple_double(fid,"params",params,nparams,h5err) 

        !write(isn) explicit,ifix_dt,iadd_force,iadd_mode,iadd_sym,iadd_damping, & 
        !     &      iuse_LES,iadd_visc,idynamic ! rev1200 
        !write(isn) Deltat,force_roll,xforce,zforce,damp_aa,damp_up,Cles,Deltag,cutoffx,cutoffz,LESflt(1:3)

        !close(isn)

        call h5fclose_f(fid,h5err)
     endif

  endif   ! ----- if nacum==0 ---

  ! reset the buffers for statistics
  if(myid.eq.0) write(*,*) '    note: reset stats = 0.0   '
  nacum = 0
  stats = 0.0
  ! reset the spectra 
  if(myid.eq.0) write(*,*) '    note: reset sp = 0.0, spl = 0.0   '
  nacumsp = 0
  sp = 0.d0; spl = 0.d0; 
  
  stats2 = 0.0; sp2 = 0.d0; spl2 = 0.d0; 
  total_dt = 0;

  stime0 = time

  ! ---------------------------------------------------------------------
  filnam = filout(1:index(filout,' ')-1)//'.'//ext//'.h5'

  w_time =-MPI_WTIME()
  if (myid.eq.0) then

     call h5fcreate_f(filnam,H5F_ACC_TRUNC_F,file_id,error)

     counth = 6
     call h5screate_simple_f(rankh,dimh1,fspace3,error)
     call h5dcreate_f(file_id,dsetname3,H5T_STD_I32LE,fspace3,dset_id3,error)
     call h5screate_simple_f(rankh,counth,mspace3,error)
     call h5dwrite_f(dset_id3,H5T_NATIVE_INTEGER,mpe,dimh1,error)
     call h5sclose_f(fspace3,error)
     call h5sclose_f(mspace3,error)
     call h5dclose_f(dset_id3,error)

     counth = 7
     call h5screate_simple_f(rankh,dimh2,fspace4,error)
     call h5dcreate_f(file_id,dsetname4,H5T_IEEE_F64LE,fspace4,dset_id4,error)
     call h5screate_simple_f(rankh,counth,mspace4,error)
     call h5dwrite_f(dset_id4,H5T_NATIVE_DOUBLE,dum,dimh2,error)
     call h5sclose_f(fspace4,error)
     call h5sclose_f(mspace4,error)
     call h5dclose_f(dset_id4,error)

     counth = 4
     call h5screate_simple_f(rankh,dimh3,fspace5,error)
     call h5dcreate_f(file_id,dsetname5,H5T_IEEE_F64LE,fspace5,dset_id5,error)
     call h5screate_simple_f(rankh,counth,mspace5,error)
     call h5dwrite_f(dset_id5,H5T_NATIVE_DOUBLE,buff,dimh3,error)
     call h5sclose_f(fspace5,error)
     call h5sclose_f(mspace5,error)
     call h5dclose_f(dset_id5,error)
  
     counth = my
     call h5screate_simple_f(rankh,dimh4,fspace6,error)
     call h5dcreate_f(file_id,dsetname6,H5T_IEEE_F64LE,fspace6,dset_id6,error)
     call h5screate_simple_f(rankh,counth,mspace6,error)
     call h5dwrite_f(dset_id6,H5T_NATIVE_DOUBLE,y,dimh4,error)
     call h5sclose_f(fspace6,error)
     call h5sclose_f(mspace6,error)
     call h5dclose_f(dset_id6,error)
 
     counth = my
     call h5screate_simple_f(rankh,dimh5,fspace7,error)
     call h5dcreate_f(file_id,dsetname7,H5T_IEEE_F64LE,fspace7,dset_id7,error)
     call h5screate_simple_f(rankh,counth,mspace7,error)
     call h5dwrite_f(dset_id7,H5T_NATIVE_DOUBLE,u00,dimh5,error)
     call h5sclose_f(fspace7,error)
     call h5sclose_f(mspace7,error)
     call h5dclose_f(dset_id7,error)

     counth = my
     call h5screate_simple_f(rankh,dimh6,fspace8,error)
     call h5dcreate_f(file_id,dsetname8,H5T_IEEE_F64LE,fspace8,dset_id8,error)
     call h5screate_simple_f(rankh,counth,mspace8,error)
     call h5dwrite_f(dset_id8,H5T_NATIVE_DOUBLE,w00,dimh6,error)
     call h5sclose_f(fspace8,error)
     call h5sclose_f(mspace8,error)
     call h5dclose_f(dset_id8,error)

     call h5screate_simple_f(rankf,dimf,fspace1,error)
     call h5dcreate_f(file_id,dsetname1,H5T_IEEE_F32LE,fspace1,dset_id1,error)
     call h5sclose_f(fspace1,error)

     call h5screate_simple_f(rankf,dimf,fspace2,error)
     call h5dcreate_f(file_id,dsetname2,H5T_IEEE_F32LE,fspace2,dset_id2,error)
     call h5sclose_f(fspace2,error)
  
     offset(3) = 0
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
           !call MPI_RECV(wk11,leng,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           !call MPI_RECV(wk22,leng,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk1,leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk2,leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        endif
  
        countf(1) = dimf(1)
        countf(2) = dimf(2)
        countf(3) = jend(iproc)-jbeg(iproc)+1
        offset(1) = 0
        offset(2) = 0
        offset(3) = jbeg(iproc)  !  iproc*(my/numerop)

        call h5screate_simple_f(rankf,countf,mspace1,error)
        call h5screate_simple_f(rankf,countf,mspace2,error)

        call h5dget_space_f(dset_id1,fspace1,error)
        call h5sselect_hyperslab_f(fspace1,H5S_SELECT_SET_F,offset,countf,error)

        call h5dget_space_f(dset_id2,fspace2,error)
        call h5sselect_hyperslab_f(fspace2,H5S_SELECT_SET_F,offset,countf,error)
 
        !call h5dwrite_f(dset_id1,H5T_NATIVE_REAL,wk11,countf,error,file_space_id=fspace1,mem_space_id=mspace1)
        !call h5dwrite_f(dset_id2,H5T_NATIVE_REAL,wk22,countf,error,file_space_id=fspace2,mem_space_id=mspace2)

        call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,wk1,countf,error,file_space_id=fspace1,mem_space_id=mspace1)
        call h5dwrite_f(dset_id2,H5T_NATIVE_DOUBLE,wk2,countf,error,file_space_id=fspace2,mem_space_id=mspace2)
     enddo

     call h5sclose_f(fspace1,error)
     call h5sclose_f(mspace1,error)
     call h5sclose_f(fspace2,error)
     call h5sclose_f(mspace2,error) 
     call h5dclose_f(dset_id1,error)  
     call h5dclose_f(dset_id2,error)

     call h5fclose_f(file_id,error)

  else
     !call MPI_SEND(wk11,mx*mz*mmy,MPI_REAL4,0,myid,MPI_COMM_WORLD,ierr)
     !call MPI_SEND(wk22,mx*mz*mmy,MPI_REAL4,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk1,mx*mz*mmy,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk2,mx*mz*mmy,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr) 
  endif
  
  w_time = MPI_WTIME()+w_time
  if(myid.eq.0) write(*,*) 'time to write:',w_time
  
  !deallocate(wk11,wk22)

end subroutine writehdf5
