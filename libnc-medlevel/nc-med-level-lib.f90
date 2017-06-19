!** nc-med-level-lib.f90
!** SUPPORT LIBRARY for fortran NetCDF applications

!** Author: GB, 20.10.2010

SUBROUTINE readNetCDF_4dvar(fileNetCDF,varname,im,jm,km,tm,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer im,jm,km,tm
real MATRIX(im,jm,km,tm)
integer mycount

mycount = 0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_get_var (ncid,VARid,MATRIX)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_close(ncid)                         ; call handle_err1(stat,mycount,fileNetCDF )


end SUBROUTINE readNetCDF_4dvar

!****************************************************************************
SUBROUTINE readNetCDF_3dvar(fileNetCDF,varname,im,jm,km,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer im,jm,km
real MATRIX(im,jm,km)
integer counter

counter = 0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_inq_varid (ncid, varname, VARid)     ;
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_get_var (ncid,VARid,MATRIX)
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_close(ncid)                          ; call handle_err1(stat,counter,fileNetCDF )

end SUBROUTINE readNetCDF_3dvar





!****************************************************************************
!
SUBROUTINE readNetCDF_2dvar(fileNetCDF,varname,im,jm,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid, mycount
integer im,jm
real MATRIX(im,jm)

mycount = 0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_get_var (ncid,VARid,MATRIX)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )

stat = nf90_close(ncid)                         ; call handle_err1(stat,mycount,fileNetCDF )


end SUBROUTINE readNetCDF_2dvar


!****************************************************************************
SUBROUTINE readNetCDF_1dvar(fileNetCDF,varname,im,ARRAY)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer im
real ARRAY(im)
integer counter
counter=0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_get_var (ncid,VARid,ARRAY)
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_close(ncid)                          ; call handle_err1(stat,counter,fileNetCDF)


end SUBROUTINE readNetCDF_1dvar

SUBROUTINE readNetCDF_0dvar(fileNetCDF,varname,VALUE)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid, counter
real VALUE

counter= 0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,counter,fileNetCDF)
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_get_var (ncid,VARid,VALUE)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_close(ncid)                         ; call handle_err1(stat,counter,fileNetCDF)


end SUBROUTINE readNetCDF_0dvar


!***************************************************************************
SUBROUTINE getDIMENSION(fileNetCDF,stringname,n)
use netcdf
implicit none
character  fileNetCDF*(*)
integer n
integer DIMid,ncid,stat
character  dim_name*30, stringname*(*)
integer counter

counter = 0
stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)        ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_inq_dimid (ncid, stringname, DIMid)         ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_Inquire_Dimension (ncid, DIMid, dim_name, n); call handle_err1(stat,counter,fileNetCDF)
stat = nf90_close(ncid)                                 ; call handle_err1(stat,counter,fileNetCDF)
END SUBROUTINE getDIMENSION


SUBROUTINE MODIFY_NC_4D(fileNetCDF,Varname,im,jm,km,tm,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer ncid, stat, VARid
integer im,jm,km,tm
real MATRIX(im,jm,km,tm)

integer counter
counter = 0

stat = nf90_open( fileNetCDF, NF90_WRITE, ncid)   ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)       ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  MATRIX  )
call handle_err2(stat, fileNetCDF, Varname)       ; call handle_err1(stat,counter,fileNetCDF )

stat=nf90_close(ncid)                             ; call handle_err1(stat,counter,fileNetCDF )

END SUBROUTINE MODIFY_NC_4D


! ***********************************************************  

SUBROUTINE MODIFY_NC_3D(fileNetCDF,Varname,im,jm,km,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer ncid, stat, VARid
integer im,jm,km
real MATRIX(im,jm,km)
integer mycount

mycount = 0;
stat = nf90_open( fileNetCDF, NF90_WRITE, ncid) ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  MATRIX  )
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )

stat=nf90_close(ncid)                           ; call handle_err1(stat,mycount,fileNetCDF )

END SUBROUTINE MODIFY_NC_3D


! ***********************************************************  

SUBROUTINE MODIFY_NC_2D(fileNetCDF,Varname,im,jm,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer mycount
integer ncid, stat, VARid
integer im,jm
real MATRIX(im,jm)

mycount = 0
stat = nf90_open( fileNetCDF, NF90_WRITE, ncid) ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)        ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  MATRIX  )
call handle_err2(stat, fileNetCDF, Varname)        ; call handle_err1(stat,mycount,fileNetCDF )

stat=nf90_close(ncid)                           ;call handle_err1(stat,mycount,fileNetCDF )

END SUBROUTINE MODIFY_NC_2D



SUBROUTINE MODIFY_NC_1D(fileNetCDF,Varname,im, ARRAY)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer ncid, stat, VARid, mycount
integer im
real ARRAY(im)

mycount = 0
stat = nf90_open( fileNetCDF, NF90_WRITE, ncid) ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  ARRAY  );
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )

stat=nf90_close(ncid)                           ;call handle_err1(stat,mycount,fileNetCDF )


END SUBROUTINE MODIFY_NC_1D

! ***********************************************************************
SUBROUTINE MODIFY_NC_0D(fileNetCDF,Varname,VALUE)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer ncid, stat, VARid, mycount
real VALUE

stat = nf90_open( fileNetCDF, NF90_WRITE, ncid); call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  VALUE  );
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )

stat=nf90_close(ncid)                           ;call handle_err1(stat,mycount,fileNetCDF )


END SUBROUTINE MODIFY_NC_0D

! *************************************************************************
SUBROUTINE MODIFY_NC_3D_byte(fileNetCDF,Varname,im,jm,km,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer im,jm,km
integer(1) MATRIX(im,jm,km)
! local
integer ncid, stat, VARid
integer mycount


mycount = 0;
stat = nf90_open( fileNetCDF, NF90_WRITE, ncid) ; call handle_err1(stat,mycount,fileNetCDF)
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat, mycount,fileNetCDF)
stat = nf90_put_var(ncid, VARid,  MATRIX  )
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat, mycount,fileNetCDF)

stat=nf90_close(ncid)                           ; call handle_err1(stat, mycount,fileNetCDF)


END SUBROUTINE MODIFY_NC_3D_byte

!****************************************************************************

SUBROUTINE readNetCDF_WaterColumn(fileNetCDF,varname,im,jm,km,ARRAY)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer im,jm,km
real ARRAY(km)
integer START(3), COUNTER(3), mycount

mycount = 0
START(1)=im
START(2)=jm
START(3)=1

COUNTER(1)=1
COUNTER(2)=1
COUNTER(3)=km

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,mycount,fileNetCDF)
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat, mycount,fileNetCDF)

stat = nf90_get_var (ncid,VARid,ARRAY,START,COUNTER)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat, mycount,fileNetCDF)

stat = nf90_close(ncid)                         ; call handle_err1(stat, mycount,fileNetCDF)


end SUBROUTINE readNetCDF_WaterColumn


!****************************************************************************

SUBROUTINE readNetCDF_VALUE_in3dvar(fileNetCDF,varname,im,jm,km,NUMBER)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer im,jm,km
real NUMBER
integer START(3)

START(1)=im
START(2)=jm
START(3)=km

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
call handle_err(stat)
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err(stat)
stat = nf90_get_var (ncid,VARid,NUMBER,START)
call handle_err(stat)
stat = nf90_close(ncid)


end SUBROUTINE readNetCDF_VALUE_in3dvar





! *************************************************************
! * setta i fillValue di una matrice al valore fillValueToChange
! * per sapere chi è fillvalue usa un'altra matrice, detta BluePrint in cui
! * è noto il valore del fillvalue (es. 1.E+20)
SUBROUTINE setFillValue_3D(im,jm,km,Mblueprint,fillValueBP,MtoChange,fillValueTC)
implicit none
integer im,jm,km
integer i,j,k
real Mblueprint(im,jm,km)
real MtoChange( im,jm,km)
real fillValueBP
real fillValueTC




do i=1,im
do j=1,jm
do k=1,km

 if ( Mblueprint(i,j,k).eq.fillValueBP ) then
     Mtochange(i,j,k) = fillValueTC
 end if

end do
end do
end do

END SUBROUTINE setFillValue_3D

!****************************************************************************
SUBROUTINE setFillValue_4D(im,jm,km,tm,Mblueprint,fillValueBP,MtoChange,fillValueTC)
implicit none
integer im,jm,km,tm
integer i,j,k,t
real Mblueprint(im,jm,km,tm)
real MtoChange( im,jm,km,tm)
real fillValueBP
real fillValueTC


do i=1,im
do j=1,jm
do k=1,km
do t=1,tm
 if ( Mblueprint(i,j,k,t).eq.fillValueBP ) then
     Mtochange(i,j,k,t) = fillValueTC
 end if
end do
end do
end do
end do

END SUBROUTINE setFillValue_4D

!****************************************************************************
SUBROUTINE setFillValue_2D(im,jm,Mblueprint,fillValueBP,MtoChange,fillValueTC)
implicit none
integer im,jm
integer i,j
real Mblueprint(im,jm)
real MtoChange( im,jm)
real fillValueBP
real fillValueTC


do i=1,im
do j=1,jm


 if ( Mblueprint(i,j).eq.fillValueBP ) then
     Mtochange(i,j) = fillValueTC
 end if


end do
end do

END SUBROUTINE setFillValue_2D

!****************************************************************************

SUBROUTINE setFillValue_3D_byte(im,jm,km,Mblueprint,fillValueBP,MtoChange,fillValueTC)
implicit none
integer im,jm,km
integer i,j,k
real Mblueprint(im,jm,km)
integer(1) MtoChange( im,jm,km)
real fillValueBP
integer(1) fillValueTC


do i=1,im
do j=1,jm
do k=1,km

 if ( Mblueprint(i,j,k).eq.fillValueBP ) then
     Mtochange(i,j,k) = fillValueTC
 end if

end do
end do
end do

END SUBROUTINE setFillValue_3D_byte

!****************************************************************************
SUBROUTINE SET_ATT(fileNetCDF, varname, attname, attvalue)
use netcdf
implicit none

character fileNetCDF*(*), varname*(*), attname*(*), attvalue*(*)
integer s, nc, varID, mycount

mycount = 0
s = nf90_open(fileNetCDF, NF90_WRITE, nc)     ; call handle_err1(s, mycount,fileNetCDF)
s = nf90_inq_varid (nc, varname, varid)       ; call handle_err1(s, mycount,fileNetCDF)
s = nf90_redef(nc)                            ; call handle_err1(s, mycount,fileNetCDF)
s = nf90_put_att(nc,varid,attname, attvalue ) ; call handle_err1(s, mycount,fileNetCDF)
s = nf90_close(nc)                            ; call handle_err1(s, mycount,fileNetCDF)


END SUBROUTINE SET_ATT
!****************************************************************************
subroutine handle_err1(status,mycount, FileNetCDF)
USE netcdf
integer status,mycount
character fileNetCDF*(*)

mycount =mycount+1
if(status .ne. nf90_NoErr)  then
   write(*,*) 'netcdf call',mycount,'with status = ',status
   write(*,*)  'file :', fileNetCDF
   write(*,*) nf90_strerror(status)
   write(*,*) 'Stopped'
   call MPI_Abort(MPI_COMM_WORLD, -1, status)
endif
end


        subroutine handle_err2(status,fileNetCDF,varname)
        USE netcdf
        use mpi
        integer status
        character fileNetCDF*(*) ,varname*(*)
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2

!****************************************************************************
SUBROUTINE handle_err(stat)
! include 'netcdf.inc'
 use netcdf
 INTEGER stat
 IF (stat .NE. nf90_NoErr) THEN
    PRINT *, NF_STRERROR(stat)
    write(*,*) 'Stopped'
    STOP 1
 ENDIF
END SUBROUTINE handle_err
!****************************************************************************
