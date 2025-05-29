subroutine timing(t,iflag)
  implicit real*8(a-h,o-z)
  !     
  !     
  integer COUNT,COUNT_RATE,COUNT_MAX,icount
  save icount
  !-----|---------------------------------------------------------------|
  call SYSTEM_CLOCK(COUNT,COUNT_RATE,COUNT_MAX)
  
  if(iflag.ne.0)then        !no init
     if(count.lt.icount)then !count wrapped around COUNT_MAX
        count=count_max-icount+count
     endif
  endif
  t=dfloat(count)/dfloat(COUNT_RATE)
  
  icount=count
  !$$$      write(*,*)'timing: ',COUNT,COUNT_RATE,COUNT_MAX,icount,t
  !
  return
end subroutine timing
!-----|---------------------------------------------------------------|

