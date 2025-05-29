subroutine getfil_t_hdf5(vor,phi,tb,u00,w00,tb00,wk,tiempo)
  use ctes
  use running
  use bcs
  use hdf5
  use LES, only: iuse_LES

  implicit none
  include "mpif.h"
  
  integer i,j,k
  integer istat(MPI_STATUS_SIZE),ierr

  integer me(3),ntotr,iproc,mxm,mym,mzm,klen,kini1,kini2
  integer ihalf,iskip,kskip,ihalfy,ncopy,icheck_write

  real*8, allocatable:: ye(:)
  
  real*8  vor(mx,mz,jb:je),phi(mx,mz,jb:je)
  real*8  wk(buffsize)
  real*8  u00(my),w00(my)
  
  integer  mgalxe,mye,mgalze,numerope
  real*8   tiempo,ubulk,wbulk,t1,t2
  integer  mpee(6)
  real*8   dume(7),buffe(4)
  real*8,  allocatable:: u00e(:),w00e(:),yee(:)
  !
  ! itemperature
  real*8  tb(mx,mz,jb:je)
  real*8  tb00(my)
  real*8  tbulk
  !real*4, allocatable:: vore(:,:),phie(:,:) ! wk2 in the original getfil
  !real*4, allocatable:: tbe(:,:)
  real*8, allocatable:: vore(:,:),phie(:,:) ! wk2 in the original getfil
  real*8, allocatable:: tbe(:,:)
  real*8, allocatable:: tb00e(:) 

  integer nopte,nparamse ! rev.1433
  integer ioptions(maxParams) ! 100 is maximum ... see in mod.f90
  real(8) params(maxParams)

  integer(hsize_t) dimf(3),dimh1(1),dimh2(1),dimh3(1),dimh4(1),dimh5(1),dimh6(1)
  
  character(3), parameter :: dsetname1 = "vor"
  character(3), parameter :: dsetname2 = "phi"
  character(3), parameter :: dsetname3 = "mpe"
  character(3), parameter :: dsetname4 = "dum"
  character(4), parameter :: dsetname5 = "buff"
  character(1), parameter :: dsetname6 = "y"
  character(3), parameter :: dsetname7 = "u00"
  character(3), parameter :: dsetname8 = "w00"
  
  character(2), parameter :: dsetname9 = "tb"
  character(4), parameter :: dsetname10 = "tb00"

  integer(hid_t) plist_id, fid
  integer(hid_t) fspace1, fspace2, fspace3, fspace4, fspace5, fspace6, fspace7, fspace8, fspace9, fspace10
  integer(hid_t) mspace1, mspace2, mspace3, mspace4, mspace5, mspace6, mspace7, mspace8, mspace9, mspace10
  integer(hid_t) dset_id1,dset_id2,dset_id3,dset_id4,dset_id5,dset_id6,dset_id7,dset_id8,dset_id9,dset_id10
 
  integer(hsize_t)  countf(3), counth(1)
  integer(hssize_t) offset(3)
 
  integer, parameter :: rankf = 3
  integer, parameter :: rankh = 1
  integer  error,h5err

  ihalf=0; ihalfy=0; iskip=1;
  !      ---------------  zero everything
  vor = 0.d0
  phi = 0.d0 
  u00 = 0.d0
  w00 = 0.d0
  if (itemperature.eq.1) then
     tb=0.d0
     tb00=0.d0
  end if
  allocate(ye(my))
  wk(buffsize)=0.d0  

  icheck_write=0 ! to skip writing on screen

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
        write(*,*) 'reading HDF5, but this is not hdf5 data, set iread_hdf5 properly'
        stop
     end if
     !call h5fopen_f(filinp,H5F_ACC_RDWR_F,fid,error)
     ! Release    	Change
     ! 1.8.10	Removed H5F_ACC_RDWR_F and H5F_ACC_RDONLY_F from comments for access_flag field in Fortran subroutine, and changed “Possible values” to “Valid values”.
     call h5fopen_f(filinp,H5F_ACC_RDONLY_F,fid,error,H5P_DEFAULT_F)
     if (error.ge.1) then 
        write(*,*) 'h5fopen_f: error =', error 
        stop
     end if

     call h5dopen_f(fid, dsetname1, dset_id1, error)
     call h5dopen_f(fid, dsetname2, dset_id2, error)
     call h5dopen_f(fid, dsetname3, dset_id3, error)
     call h5dopen_f(fid, dsetname4, dset_id4, error)
     call h5dopen_f(fid, dsetname5, dset_id5, error)
     call h5dopen_f(fid, dsetname6, dset_id6, error)
     call h5dopen_f(fid, dsetname7, dset_id7, error)
     call h5dopen_f(fid, dsetname8, dset_id8, error)
     call h5dopen_f(fid, dsetname9, dset_id9, error)
     call h5dopen_f(fid, dsetname10, dset_id10, error)
 
     call h5dread_f(dset_id3,H5T_NATIVE_INTEGER,mpee,dimh1,error)
     mgalxe=mpee(1); mye=mpee(2);   mgalze=mpee(3)
     me(1) =mpee(5); me(2)=mpee(2); me(3)=mpee(6)
     if((mgalx.ne.mgalxe).or.(my.ne.mye).or.(mgalz.ne.mgalze)) then
        write(*,*) ' different physical grid points!'
        write(*,'(a,3(1X,I5),a,3(1X,I5))') '  reading ', mgalxe,me(2),mgalze,&
                                                    ' =>', mgalx, my, mgalz
     endif
     call h5dread_f(dset_id4,H5T_NATIVE_DOUBLE,dume,dimh2,error)
     tiempo=dume(1); Ree=dume(2); alpe=dume(3); Lye=dume(4)*pi
     game  =dume(5); se =dume(6); chie=dume(7)
     if (icheck_write.eq.1) write(*,*) 'dume',dume, Ree, alpe, Lye, game, se, chie
     
     call h5dread_f(dset_id5,H5T_NATIVE_DOUBLE,buffe,dimh3,error)
     xoffb =buffe(1); xofft=buffe(2); zoffb=buffe(3); zofft=buffe(4)
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

        if ((my.ne.mye).and.(mgalx.ne.mgalxe)) then
           write(*,*) ' Do not interpolate in my at the same time...'
           stop
        end if

     endif

     if (icheck_write.eq.1) then
        write(*,*)
        write(*,*) 'reading input file ..., time=',tiempo
        write(*,'(a,3(1X,I5))') '                    ...,(mx,my,mz)=',me
        write(*,*)
     end if

     allocate(yee(me(2)))
     call h5dread_f(dset_id6, H5T_NATIVE_DOUBLE, yee,  dimh4, error)
     ye(1:me(2))=yee ! copy from local yee to global ye

     if (me(2).ne.my) then
        write(*,*) 'Spline interpolation in y, me(2)= ',me(2),' =>my=',my
     endif
     allocate(u00e(me(2)),w00e(me(2)))
     call h5dread_f(dset_id7, H5T_NATIVE_DOUBLE, u00e, dimh5, error)
     call h5dread_f(dset_id8, H5T_NATIVE_DOUBLE, w00e, dimh6, error)

     mym = min(my,me(2))
     u00(1:mym) = u00e(1:mym)
     w00(1:mym) = w00e(1:mym)
     if (itemperature.eq.1) then
        allocate(tb00e(me(2)))
        call h5dread_f(dset_id10, H5T_NATIVE_DOUBLE, tb00e, dimh6, error)
        tb00(1:mym) = tb00e(1:mym)
     end if

  endif
  call MPI_BCAST(me,3,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(ye,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(xofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zoffb,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(zofft,1,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(u00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(w00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  call MPI_BCAST(tiempo,1,MPI_REAL8,0,MPI_COMM_WORLD,error) 
  call MPI_BCAST(ihalf,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)   
  call MPI_BCAST(ihalfy,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)   
  if (itemperature.eq.1) then
     call MPI_BCAST(tb00,my,MPI_REAL8,0,MPI_COMM_WORLD,error)
  end if

  ntotr = me(1)*me(3)

  mxm   = min(mx,me(1))
  mym   = min(my,me(2))
  mzm   = min(mz,me(3))
 
  klen  = (mzm+1)/2
  kini1 = me(3) - klen + 1
  kini2 = mz    - klen + 1
 
  allocate(vore(me(1),me(3)),phie(me(1),me(3)))
  vore=0.d0;phie=0.d0
  if (itemperature.eq.1) then
     allocate(tbe(me(1),me(3)))
     tbe=0.d0
  end if

  if (myid.eq.0) then
 
     do iproc=0,numerop-1
        do j=jbeg(iproc),jend(iproc)

           countf(1) = me(1)
           countf(2) = me(3)
           countf(3) = 1
           offset(1) = 0
           offset(2) = 0
           if (ihalfy.eq.1) then
              offset(3) = j*2 ! skip on plane to get half-coarse grid data in y
           else
              offset(3) = j
           end if
           if (j.le.(me(2)-1)) then
              call h5dget_space_f(dset_id1,fspace1,error)
              call h5sselect_hyperslab_f(fspace1,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace1, error)
              call h5dread_f(dset_id1,H5T_NATIVE_DOUBLE,vore,countf,error,mspace1,fspace1)
              !call h5dread_f(dset_id1,H5T_NATIVE_REAL,vore,countf,error,mspace1,fspace1)

              call h5dget_space_f(dset_id2,fspace2,error)
              call h5sselect_hyperslab_f(fspace2,H5S_SELECT_SET_F,offset,countf,error)
              call h5screate_simple_f(rankf,countf, mspace2, error)      
              call h5dread_f(dset_id2,H5T_NATIVE_DOUBLE,phie,countf,error,mspace2,fspace2)
              !call h5dread_f(dset_id2,H5T_NATIVE_REAL,phie,countf,error,mspace2,fspace2)
              
              if (itemperature.eq.1) then
                 call h5dget_space_f(dset_id9,fspace9,error)
                 call h5sselect_hyperslab_f(fspace9,H5S_SELECT_SET_F,offset,countf,error)
                 call h5screate_simple_f(rankf,countf, mspace9, error)
                 call h5dread_f(dset_id9,H5T_NATIVE_DOUBLE,tbe,countf,error,mspace9,fspace9)
                 !call h5dread_f(dset_id9,H5T_NATIVE_REAL,tbe,countf,error,mspace9,fspace9)
              end if
           else
              vore=0.d0
              phie=0.d0
              if (itemperature.eq.1) then
                 tbe=0.d0
              end if
           endif

           if (iproc.eq.0) then
              vor(1:mxm,1:klen,j)    =real(vore(1:mxm*iskip:iskip,1:klen),kind=8)
              vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm*iskip:iskip,1+kini1:me(3)),kind=8)
              phi(1:mxm,1:klen,j)    =real(phie(1:mxm*iskip:iskip,1:klen),kind=8)
              phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm*iskip:iskip,1+kini1:me(3)),kind=8)
              if (itemperature.eq.1) then
                 tb(1:mxm,1:klen,j)    =real(tbe(1:mxm*iskip:iskip,1:klen),kind=8)
                 tb(1:mxm,1+kini2:mz,j)=real(tbe(1:mxm*iskip:iskip,1+kini1:me(3)),kind=8)                 
              end if
           else
              call MPI_SEND(vore,ntotr,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,error)
              call MPI_SEND(phie,ntotr,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,error)
              if (itemperature.eq.1) then
                 call MPI_SEND(tbe,ntotr,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,error)
              endif
           endif
 
           !if (ihalfy.eq.1) then
           !   ! the case for my being half
           !   read(iinp) wk2
           !   write(*,*) 'skip one plane by just reading '              
           !endif           

        enddo
     enddo
     call h5sclose_f(fspace1,error)
     call h5sclose_f(mspace1,error)
     call h5dclose_f(dset_id1,error)
  
     call h5sclose_f(fspace2,error)
     call h5sclose_f(mspace2,error)
     call h5dclose_f(dset_id2,error)

     if (itemperature.eq.1) then
        call h5sclose_f(fspace9,error)
        call h5sclose_f(mspace9,error)
        call h5dclose_f(dset_id9,error)
     end if    

     if (iread_footer.eq.1) then ! default is iread_footer=1

        !write(*,*) 'iread_footer is not implemented for hdf5'
        !stop
        if (iread_hdf5.eq.1) then

           call h5read_simple_double(fid,"CFL",CFLe,1,h5err)
           call h5read_simple_double(fid,"uprim",uprime,1,h5err)
           call h5read_simple_double(fid,"vbulk",vbulke,1,h5err)
           call h5read_simple_double(fid,"xforce",xforcee,1,h5err)

           call h5read_simple_int(fid,"irev",ireve,1,h5err)
           call h5read_simple_int(fid,"nopt",nopte,1,h5err)
           call h5read_simple_int(fid,"nparams",nparamse,1,h5err)

           call h5read_simple_int(fid,"ioptions",ioptions,nopte,h5err)
           call h5read_simple_double(fid,"params",params,nparamse,h5err)

        else

           ! !!! read some extra parameters (testing)
           ! write(iout) CFL, uprim, vbulk, xforce ! rev934
           ! !!!write(iout) irev,explicit,iadd_force,iadd_mode,iadd_sym ! rev1153  
           ! write(iout) irev,explicit,iadd_force,iadd_mode,iadd_sym,iadd_damping ! rev1200      
           !read(iinp,err=980) CFLe, uprime, vbulke, xforcee ! real*8
           !read(iinp,err=980) ireve,nopte,nparamse
           !read(iinp,err=980) ioptions(1:nopte)
           !read(iinp,err=980) params(1:nparamse)
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

     call h5fclose_f(fid,error) ! fixed from rev.1720
     !call h5fflush_f(fid, H5F_SCOPE_GLOBAL_F, error)
     !write(*,*)  'getfil_t_hdf5 file closed and flushed. ', error

  else
     do j=jb,je
        if ((ihalfy.eq.1).and.(mod(j+1,2).eq.0)) then
           ! Skip reading
           continue
        endif
        call MPI_RECV(vore,ntotr,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,error)
        call MPI_RECV(phie,ntotr,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,error)
        
        vor(1:mxm,1:klen,j)    =real(vore(1:mxm,1:klen),kind=8)
        vor(1:mxm,1+kini2:mz,j)=real(vore(1:mxm,1+kini1:me(3)),kind=8)
        phi(1:mxm,1:klen,j)    =real(phie(1:mxm,1:klen),kind=8)
        phi(1:mxm,1+kini2:mz,j)=real(phie(1:mxm,1+kini1:me(3)),kind=8)
        if (itemperature.eq.1) then
           call MPI_RECV(tbe,ntotr,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,error)
           tb(1:mxm,1:klen,j)    =real(tbe(1:mxm,1:klen),kind=8)
           tb(1:mxm,1+kini2:mz,j)=real(tbe(1:mxm,1+kini1:me(3)),kind=8)
        end if
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
  !call set_conj(vor,phi,-1)

  t2= MPI_WTIME()
  if (myid.eq.0) write(*,*) '  --- time: read HDF5 =',(t2-t1)

  call chikj2jik(vor,vor,wk,wk)
  call chikj2jik(phi,phi,wk,wk)
  if (itemperature.eq.1) then
     call chikj2jik(tb,tb,wk,wk)
  end if

  ! ------ you should not interpolate from a small number of grids to ------
  ! ------- a much larger, interpolation itself has error --------
  ncopy=0
  if ((me(2).eq.my/2)) then
      ncopy=2
  elseif ((me(2).eq.my/3)) then
      ncopy=3
  elseif ((me(2).eq.my/5)) then
      ncopy=5
  elseif ((me(2).eq.my/6)) then
      ncopy=6 
  endif
 
  if ((me(2).lt.my).and.(ncopy.eq.0)) then
     ! !note not to interpolated with different Ly.
     if (abs(Lye-Ly).lt.tiny) then
        write(*,*) myid,': interpolating in y-dir, not supported in temperature'
        stop
        !call interp_y(vor,phi,u00,w00,me,ye)
     else
        if (myid.eq.0) write(*,*) ' Do not interpolate in y-dir, with different Ly from Lye in hre2.dat'
        stop
     end if

  elseif ((me(2).lt.my).and.(ncopy.gt.0)) then
     write(*,*) myid,'getfil, ncopy=',ncopy 
     call ncopy_y(vor,phi,u00,w00,me,ncopy)
     if (itemperature.eq.1) then
        call ncopy_y_t(tb,tb00,me,ncopy)
     end if
  endif 
  
  deallocate(ye) ! deallocated by all proc.
  if(myid.eq.0) then 
     deallocate(yee,u00e,w00e)
     if (itemperature.eq.1) then
        deallocate(tb00e)
     end if
  end if
  deallocate(vore,phie)
  if (itemperature.eq.1) then
     deallocate(tbe)
  end if
  ! set bulk velocity to be zero
  ubulk=sum(u00(1:my))/dfloat(my)
  wbulk=sum(w00(1:my))/dfloat(my)
  if ((myid.eq.0).and.(icheck_write.eq.1)) then 
     write(*,*) 'check ubulk, wbulk=', ubulk, wbulk
  end if
  if (itemperature.eq.1) then
     tbulk=sum(tb00(1:my))/dfloat(my)
     if ((myid.eq.0).and.(icheck_write.eq.1)) then 
        write(*,*) 'check tbulk =', tbulk
     end if
  end if


  xwkb=xoffb; xwkt=xofft; ! bug fixed, 2014/June/24 sekimoto
  zwkb=zoffb; zwkt=zofft;

  if (abs(game-gam).gt.tiny) then ! rev1140
     if (myid.eq.0) write(*,*) 'rescale velocity ... and lap.v: game/gam=',game/gam
     u00=u00*(game/gam)
     w00=w00*(game/gam)
     phi=phi*(gam/game) ! (game/gam*(gam/game)^2)
     if (iuse_LES.eq.0) then
        vor=vor*(game/gam)*sqrt(Re/Ree) ! rev1157
     end if
  elseif (abs(Ree-Re).gt.tiny) then ! rev1157
     if ((iuse_LES.eq.0).and.(iuse_newton.eq.0)) then
        if (myid.eq.0) write(*,*) 'rescale vorticity  sqrt(Re/Ree)',sqrt(Re/Ree)
        vor=vor*(game/gam)*sqrt(Re/Ree) ! rev1157
     end if
  end if

  t2= MPI_WTIME()
  if (myid.eq.0) write(*,*) '  --- time: read HDF5, including comm.=',(t2-t1)

end subroutine getfil_t_hdf5

! see readwrite.f90 for the routine ncopy()

subroutine write_t_hdf5(vor,phi,tb,u00,w00,tb00,wk1,wk2,wk3,wk)

  use hdf5
  use ctes
  use running !,only:time,id22,filout,itemperature
  use statistics
  use bcs
  use LES
  use temp, only: fkappa, bym, gbeta
        
  implicit none  
  include 'mpif.h'

  real*8  vor(0:2*my-1,0:mx1,kb:ke),phi(0:2*my-1,0:mx1,kb:ke) ! buffsize
  real*8  wk(buffsize)
  real*8  wk1(mx,mz,jb:je), wk2(mx,mz,jb:je), wk3(mx,mz,jb:je)
  !real*8, allocatable::  wk11(:,:,:),wk22(:,:,:),wk33(:,:,:)
  real*8  u00(my),w00(my),dum(7),buff(4)

  real*8 tb(0:2*my-1,0:mx1,kb:ke)
  !real*8 tb(buffsize)
  real*8 tb00(my)

  integer mpe(6)
  real*8  w_time
  integer iproc,leng,ipo

  integer nopte,nparamse ! rev.1433
  integer ioptions(maxParams) ! 100 is maximum ... see in mod.f90
  real(8) params(maxParams)

  integer istat(MPI_STATUS_SIZE),ierr
  character ext*4, fil*256, filst*256
  
  integer(hsize_t) dimf(3),dimh1(1),dimh2(1),dimh3(1),dimh4(1),dimh5(1),dimh6(1)

  character(3), parameter :: dsetname1 = "vor"  
  character(3), parameter :: dsetname2 = "phi" 
  character(3), parameter :: dsetname3 = "mpe"
  character(3), parameter :: dsetname4 = "dum"
  character(4), parameter :: dsetname5 = "buff"
  character(1), parameter :: dsetname6 = "y"
  character(3), parameter :: dsetname7 = "u00"
  character(3), parameter :: dsetname8 = "w00"
 
  character(2), parameter :: dsetname9 = "tb"
  character(4), parameter :: dsetname10 = "tb00"

  integer(hid_t) plist_id, fid

  ! hdf5 
  integer(hid_t) fspace1, fspace2, fspace3, fspace4, fspace5, fspace6, fspace7, fspace8, fspace9, fspace10
  integer(hid_t) mspace1, mspace2, mspace3, mspace4, mspace5, mspace6, mspace7, mspace8, mspace9, mspace10
  integer(hid_t) dset_id1,dset_id2,dset_id3,dset_id4,dset_id5,dset_id6,dset_id7,dset_id8,dset_id9,dset_id10
    
  integer(hsize_t) countf(3), counth(1)  
  integer(hsize_t) offset(3) 

  integer, parameter :: rankf = 3
  integer, parameter :: rankh = 1
  integer  h5err  

  ! for hdf5
  integer(hid_t) dspace, mspace, dset
  integer(hsize_t), dimension(2):: dimstat 
  integer(hsize_t), dimension(3):: dimspec

  !allocate(wk11(mx,mz,jb:je),wk22(mx,mz,jb:je),wk33(mx,mz,jb:je)) 
 
  call chjik2ikj(vor,wk1,wk,wk)  
  call chjik2ikj(phi,wk2,wk,wk)  
  call chjik2ikj(tb,wk3,wk,wk)  
  !wk11=real(wk1,kind=8)
  !wk22=real(wk2,kind=8)
  !wk33=real(wk3,kind=8)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  mpe =(/mgalx,my,mgalz,numerop,mx,mz/)
  dum =(/time,Re,alp,Ly/pi,gam,s,chi/)
  buff=(/xwkb,xwkt,zwkb,zwkt/)

  dimf  = (/mx,mz,my/) ! totaldims
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
           if (itemperature.eq.1) then
              leng=(jend(iproc)-jbeg(iproc)+1)*ntstat
              call MPI_RECV(tstats(1,ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc, &
                   MPI_COMM_WORLD,istat,ierr)
           end if
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
        if (itemperature.eq.1) then
           call h5write_simple_double(fid,'gbeta',gbeta,1,h5err)
           call h5write_simple_double(fid,'fkappa',fkappa,1,h5err)
           call h5write_simple_double(fid,'bym',bym,1,h5err)
        end if
        call h5write_simple_double(fid,'y',y,my,h5err)

        !call h5write_simple_int(fid,'nstat',nstat,1,h5err)
        dimstat=(/nstat,my/)
        call h5screate_simple_f(2,dimstat,dspace,h5err)
        call h5dcreate_f(fid,'stats2',H5T_IEEE_F64LE,dspace,dset,h5err)
        call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,stats2,dimstat,h5err)
        call h5dclose_f(dset,h5err)
        call h5sclose_f(dspace,h5err)
        if (itemperature.eq.1) then
           !call h5write_simple_int(fid,'ntstat',nstat,1,h5err)
           dimstat=(/ntstat,my/)
           call h5screate_simple_f(2,dimstat,dspace,h5err)
           call h5dcreate_f(fid,'tstats',H5T_IEEE_F64LE,dspace,dset,h5err)
           call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,tstats,dimstat,h5err)
           call h5dclose_f(dset,h5err)
           call h5sclose_f(dspace,h5err)
        end if

     else

        call MPI_SEND(stats2(1,jb),mmy*nstat,MPI_DOUBLE_PRECISION,0,myid, & 
             MPI_COMM_WORLD,ierr)
        if (itemperature.eq.1) then
           call MPI_SEND(tstats(1,jb),mmy*ntstat,MPI_DOUBLE_PRECISION,0,myid, & 
                MPI_COMM_WORLD,ierr)
        end if
     endif

     ! leng=nspec*(mx1+1)*(nz1+1)
     ! write spectra in .sta.h5 together
     leng=nspec*(mx1+1)*(mz1+1)
     call MPI_ALLREDUCE(spl2,sp2,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (itemperature.eq.1) then
        leng=nspect*(mx1+1)*(mz1+1)
        call MPI_ALLREDUCE(splt,spt,leng,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     if (myid.eq.0) then
        dimspec=(/nspec,mx1+1,mz1+1/)
        call h5screate_simple_f(3,dimspec,dspace,h5err)
        call h5dcreate_f(fid,'sp2',H5T_IEEE_F64LE,dspace,dset,h5err)
        call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,sp2,dimspec,h5err)
        call h5dclose_f(dset,h5err)
        call h5sclose_f(dspace,h5err)        

        if (itemperature.eq.1) then
           dimspec=(/nspect,mx1+1,mz1+1/)
           call h5screate_simple_f(3,dimspec,dspace,h5err)
           call h5dcreate_f(fid,'spt',H5T_IEEE_F64LE,dspace,dset,h5err)
           !call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,sp2,dimspec,h5err) ! BUG 2019/0913
           call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,spt,dimspec,h5err)
           call h5dclose_f(dset,h5err)
           call h5sclose_f(dspace,h5err)
        end if
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
  if (itemperature.eq.1) then
     tstats = 0.d0
     spt = 0.d0; splt = 0.d0;
  end if
  total_dt = 0;

  stime0 = time

  ! ---------------------------------------------------------------
  fil = filout(1:index(filout,' ')-1)//'.'//ext//'.h5'

  w_time =-MPI_WTIME()
  if (myid.eq.0) then

     call h5fcreate_f(fil,H5F_ACC_TRUNC_F,fid,h5err)
     if(h5err.ne.0) then
        write(*,*) h5err, "ERROR: Problem openning ", trim(fil)
        call h5eprint_f(h5err,trim(fil))
        write(*,*) "overwrite the file"
        !stop
     endif

     call h5write_simple_int(fid,dsetname3,mpe,dimh1,h5err)

     !counth = 6 ! mpe
     !call h5screate_simple_f(rankh,dimh1,fspace3,h5err)
     !call h5dcreate_f(file_id,dsetname3,H5T_STD_I32LE,fspace3,dset_id3,h5err)
     !call h5screate_simple_f(rankh,counth,mspace3,h5err)
     !call h5dwrite_f(dset_id3,H5T_NATIVE_INTEGER,mpe,dimh1,h5err)
     !call h5sclose_f(fspace3,h5err)
     !call h5sclose_f(mspace3,h5err)
     !call h5dclose_f(dset_id3,h5err)

     call h5write_simple_double(fid,dsetname4,dum,dimh2,h5err)
     !counth = 7 ! dum

     call h5write_simple_double(fid,dsetname5,buff,dimh3,h5err)
     !counth = 4 ! buff
  
     call h5write_simple_double(fid,dsetname6,y,dimh4,h5err)
     !counth = my ! y
 
     call h5write_simple_double(fid,dsetname7,u00,dimh5,h5err)
     !counth = my ! u00

     call h5write_simple_double(fid,dsetname8,w00,dimh6,h5err)
     !counth = my  ! w00

     if (itemperature.eq.1) then
        call h5write_simple_double(fid,dsetname10,tb00,dimh6,h5err)
        !counth = my  ! tb00
     end if

     call h5write_simple_double(fid,"CFL",CFLe,1,h5err)
     call h5write_simple_double(fid,"uprim",uprim,1,h5err)
     call h5write_simple_double(fid,"vbulk",vbulk,1,h5err)
     call h5write_simple_double(fid,"xforce",xforce,1,h5err)

     call h5write_simple_int(fid,"irev",irev,1,h5err)
     call h5write_simple_int(fid,"nopt",nopt,1,h5err)
     call h5write_simple_int(fid,"nparams",nparams,1,h5err)
           
     call h5write_simple_int(fid,"ioptions",ioptions,nopt,h5err) 
     call h5write_simple_double(fid,"params",params,nparams,h5err) 
     

     call h5screate_simple_f(rankf,dimf,fspace1,h5err)
     !call h5dcreate_f(file_id,dsetname1,H5T_IEEE_F32LE,fspace1,dset_id1,h5err)
     call h5dcreate_f(fid,dsetname1,H5T_IEEE_F64LE,fspace1,dset_id1,h5err)
     call h5sclose_f(fspace1,h5err)

     call h5screate_simple_f(rankf,dimf,fspace2,h5err)
     call h5dcreate_f(fid,dsetname2,H5T_IEEE_F64LE,fspace2,dset_id2,h5err)
     call h5sclose_f(fspace2,h5err)

     if (itemperature.eq.1) then 
        call h5screate_simple_f(rankf,dimf,fspace9,h5err)
        call h5dcreate_f(fid,dsetname9,H5T_IEEE_F64LE,fspace9,dset_id9,h5err)
        call h5sclose_f(fspace9,h5err)
     end if

     offset(3) = 0
     do iproc=0,numerop-1
        if (iproc.ne.0) then
           leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
           !call MPI_RECV(wk11,leng,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           !call MPI_RECV(wk22,leng,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,istat,ierr)           
           !if (itemperature.eq.1) call MPI_RECV(wk33,leng,MPI_REAL4,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk1,leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
           call MPI_RECV(wk2,leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)           
           if (itemperature.eq.1) call MPI_RECV(wk3,leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        endif
  
        countf(1) = dimf(1)
        countf(2) = dimf(2)
        countf(3) = jend(iproc)-jbeg(iproc)+1
        offset(1) = 0
        offset(2) = 0
        offset(3) = jbeg(iproc)  ! iproc*(my/numerop)

        call h5screate_simple_f(rankf,countf,mspace1,h5err)
        call h5screate_simple_f(rankf,countf,mspace2,h5err)
        if (itemperature.eq.1) call h5screate_simple_f(rankf,countf,mspace9,h5err)

        call h5dget_space_f(dset_id1,fspace1,h5err)
        call h5sselect_hyperslab_f(fspace1,H5S_SELECT_SET_F,offset,countf,h5err)

        call h5dget_space_f(dset_id2,fspace2,h5err)
        call h5sselect_hyperslab_f(fspace2,H5S_SELECT_SET_F,offset,countf,h5err)
        if (itemperature.eq.1) then
           call h5dget_space_f(dset_id9,fspace9,h5err)
           call h5sselect_hyperslab_f(fspace9,H5S_SELECT_SET_F,offset,countf,h5err)
        end if
        !call h5dwrite_f(dset_id1,H5T_NATIVE_REAL,wk11,countf,h5err,file_space_id=fspace1,mem_space_id=mspace1)
        !call h5dwrite_f(dset_id2,H5T_NATIVE_REAL,wk22,countf,h5err,file_space_id=fspace2,mem_space_id=mspace2)
        call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,wk1,countf,h5err,file_space_id=fspace1,mem_space_id=mspace1)
        call h5dwrite_f(dset_id2,H5T_NATIVE_DOUBLE,wk2,countf,h5err,file_space_id=fspace2,mem_space_id=mspace2)
        if (itemperature.eq.1) then
           !call h5dwrite_f(dset_id9,H5T_NATIVE_REAL,wk33,countf,h5err,file_space_id=fspace9,mem_space_id=mspace9)
           call h5dwrite_f(dset_id9,H5T_NATIVE_DOUBLE,wk3,countf,h5err,file_space_id=fspace9,mem_space_id=mspace9)
        end if

     enddo

     call h5sclose_f(fspace1,h5err)
     call h5sclose_f(mspace1,h5err)
     call h5sclose_f(fspace2,h5err)
     call h5sclose_f(mspace2,h5err) 
     call h5dclose_f(dset_id1,h5err)  
     call h5dclose_f(dset_id2,h5err)

     if (itemperature.eq.1) then
        call h5sclose_f(fspace9,h5err)
        call h5sclose_f(mspace9,h5err) 
        call h5dclose_f(dset_id9,h5err)
     end if
     call h5fclose_f(fid,h5err)

     !call writeheader_hdf5(fil,'vor_phi_temp') ! fopen and fclose are inside

  else
     !call MPI_SEND(wk11,mx*mz*mmy,MPI_REAL4,0,myid,MPI_COMM_WORLD,ierr)
     !call MPI_SEND(wk22,mx*mz*mmy,MPI_REAL4,0,myid,MPI_COMM_WORLD,ierr)
     !if (itemperature.eq.1) call MPI_SEND(wk33,mx*mz*mmy,MPI_REAL4,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk1,mx*mz*mmy,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(wk2,mx*mz*mmy,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
     if (itemperature.eq.1) call MPI_SEND(wk3,mx*mz*mmy,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
  w_time = MPI_WTIME()+w_time
  if(myid.eq.0) write(*,*) 'time to write:',w_time
  
  !deallocate(wk11,wk22)
  !if (itemperature.eq.1) then
  !   deallocate(wk33)
  !end if

end subroutine write_t_hdf5



subroutine h5write_simple_double(fid,var,array,dims,h5err)
  
  use hdf5
  implicit none
  
  character(len=*) :: var
  integer dims,h5err
  real(8) array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims
  call h5screate_simple_f(1,hdims,dspace,h5err)
  call h5dcreate_f(fid,trim(var),H5T_IEEE_F64LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,array,hdims,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

end subroutine h5write_simple_double

subroutine h5write_simple_int(fid,var,iarray,dims,h5err)
  
  use hdf5
  implicit none
  
  character(len=*) :: var
  integer dims,h5err
  integer iarray(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims
  call h5screate_simple_f(1,hdims,dspace,h5err)
  call h5dcreate_f(fid,trim(var),H5T_STD_I32LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_INTEGER,iarray,hdims,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

end subroutine h5write_simple_int


subroutine h5read_simple_double(fid,var,array,dims,h5err)
  
  use hdf5
  implicit none
  
  character(len=*) :: var
  integer dims,h5err
  real(8) array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims

  call h5dopen_f(fid,var,dset,h5err)
  call h5dread_f(dset,H5T_NATIVE_DOUBLE,array,hdims,h5err)
  call h5dclose_f(dset,h5err)

end subroutine h5read_simple_double

subroutine h5read_simple_int(fid,var,array,dims,h5err)
  
  use hdf5
  implicit none
  
  character(len=*) :: var
  integer dims,h5err
  integer array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims

  call h5dopen_f(fid,var,dset,h5err)
  call h5dread_f(dset,H5T_NATIVE_INTEGER,array,hdims,h5err)
  call h5dclose_f(dset,h5err)

end subroutine h5read_simple_int


subroutine interp_y_t(tb,tb00,me,ye)
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
  real*8 tb(0:2*my-1,0:mx1,kb:ke)
  real*8 tb00(0:my1)

  ! --------- wk3 and wk4 for second derivatives ---------- !
  ! -------- wk1 and wk2 for original phi and vor --------- !
  real*8, allocatable:: wk1(:),wk3(:)
  allocate(wk1(0:2*me(2)+1),wk3(0:2*me(2)+1))
  wk1=0.d0; wk3=0.d0;

  ye(me(2)+1)=Ly/2.d0

  ! ---------------- note: my>me(2) ---------------------
  do k=kb,ke
     do i=0,mx1
        shp = zexp(-xalp(i)*xofft) ! positive shift
        wk1(0:2*me(2)-1)=tb(0:2*me(2)-1,i,k)        
        wk1(2*me(2))  =dreal(dcmplx(tb(0,i,k),tb(1,i,k))*shp)
        wk1(2*me(2)+1)=dimag(dcmplx(tb(0,i,k),tb(1,i,k))*shp)
        call cubic_spline(ye(1:me(2)+1),wk1,me(2)+1,wk3,2)  
        do j=0,my1
           call interpolate(y(j),ye(1:me(2)+1),wk1,wk3,me(2)+1,tb(2*j:2*j+1,i,k),2)
        enddo
     enddo
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! ---------------------00 mode -----------------------
  wk1(0:me(2)-1)=tb00(0:me(2)-1); wk1(me(2))=tb00(0)
  call cubic_spline(ye(1:me(2)+1),wk1(0:me(2)),me(2)+1,wk3(0:me(2)),1)  
  do j=0,my1
     call interpolate(y(j),ye(1:me(2)+1),wk1(0:me(2)),wk3(0:me(2)),me(2)+1,tb00(j),1)
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! -after this we get the tb interpolated in y-
  deallocate(wk1,wk3)

end subroutine interp_y_t
