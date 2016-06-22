subroutine rdgrd

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!--------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!                                                                      !
! READ parameters of the MFS_16_72 grid                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use set_knd
  use drv_str
  use grd_str
  use netcdf
  use filenames

  implicit none

  INTEGER(i4)                    :: stat, ncid, idvar
  real(r4), ALLOCATABLE          :: x3(:,:,:), x2(:,:), x1(:)

  character                     :: cgrd

    write(cgrd,'(i1)') grd%grd_mod

! Vertical levels ------------------------------------------------
    stat = nf90_open(GRID_FILE, NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)


! Get dimensions 
      stat = nf90_inq_dimid (ncid, 'im', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%im)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'jm', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%jm)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'km', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%km)
    if (stat /= nf90_noerr) call netcdf_err(stat)

    write(drv%dia,*)'Grid dimensions are: ',grd%im,grd%jm,grd%km

!  Allocate grid arrays
     ALLOCATE ( grd%reg(grd%im,grd%jm))        ; grd%reg = huge(grd%reg(1,1))
     ALLOCATE ( grd%msk(grd%im,grd%jm,grd%km)) ; grd%msk = huge(grd%msk(1,1,1))
     ALLOCATE ( grd%ums(grd%im,grd%jm,grd%km), grd%vms(grd%im,grd%jm,grd%km))
     grd%ums = huge(grd%ums(1,1,1))  ; grd%vms = huge(grd%vms(1,1,1))
     ALLOCATE ( grd%f(grd%im,grd%jm)) ; grd%f = huge(grd%f(1,1));
      
     ALLOCATE ( grd%mdt(grd%im,grd%jm)) ; grd%mdt = huge(grd%mdt(1,1));
     ALLOCATE ( grd%hgt(grd%im,grd%jm)) ; grd%hgt = huge(grd%hgt(1,1))
     ALLOCATE ( grd%bx(grd%im,grd%jm))  ; grd%bx = huge(grd%bx(1,1))
     ALLOCATE ( grd%by(grd%im,grd%jm))  ; grd%by = huge(grd%by(1,1))
     ALLOCATE ( grd%lon(grd%im,grd%jm)) ; grd%lon = huge(grd%lon(1,1))
     ALLOCATE ( grd%lat(grd%im,grd%jm)) ; grd%lat = huge(grd%lat(1,1))
     ALLOCATE ( grd%dep(grd%km))        ; grd%dep = huge(grd%dep(1))
     ALLOCATE ( grd%dx(grd%im,grd%jm))  ; grd%dx  = huge(grd%dx(1,1))
     ALLOCATE ( grd%dy(grd%im,grd%jm))  ; grd%dy  = huge(grd%dy(1,1))
     ALLOCATE ( grd%dz(grd%km))         ; grd%dz  = huge(grd%dz(1))

     ALLOCATE ( grd%dxdy(grd%im,grd%jm))         ; grd%dxdy = huge(grd%dxdy(1,1))
     ALLOCATE ( grd%alx(grd%im,grd%jm) )         ; grd%alx  = huge(grd%alx(1,1))
     ALLOCATE ( grd%aly(grd%im,grd%jm) )         ; grd%aly  = huge(grd%aly(1,1))
     ALLOCATE ( grd%btx(grd%im,grd%jm) )         ; grd%btx  = huge(grd%btx(1,1))
     ALLOCATE ( grd%bty(grd%im,grd%jm) )         ; grd%bty  = huge(grd%bty(1,1))
     ALLOCATE ( grd%scx(grd%im,grd%jm) )         ; grd%scx  = huge(grd%scx(1,1))
     ALLOCATE ( grd%scy(grd%im,grd%jm) )         ; grd%scy  = huge(grd%scy(1,1))
     ALLOCATE ( grd%msr(grd%im,grd%jm,grd%km) )  ; grd%msr  = huge(grd%msr(1,1,1))
     ALLOCATE ( grd%imx(grd%km))                 ; grd%imx  = huge(grd%imx(1))
     ALLOCATE (  grd%jmx(grd%km))                ; grd%jmx  = huge(grd%jmx(1))
     ALLOCATE ( grd%istp(grd%im,grd%jm))         ; grd%istp = huge(grd%istp(1,1))
     ALLOCATE ( grd%jstp(grd%im,grd%jm))         ; grd%jstp = huge(grd%jstp(1,1))
     ALLOCATE ( grd%inx(grd%im,grd%jm,grd%km))   ; grd%inx  = huge(grd%inx(1,1,1))
     ALLOCATE ( grd%jnx(grd%im,grd%jm,grd%km))   ; grd%jnx  = huge(grd%jnx(1,1,1))
     ALLOCATE ( grd%fct(grd%im,grd%jm,grd%km) )  ; grd%fct  = huge(grd%fct(1,1,1))

     ALLOCATE ( Dump_chl(grd%im,grd%jm,grd%km) ) ; Dump_chl  = 0.0
     ALLOCATE ( Dump_msk(grd%im,grd%jm) )        ; Dump_msk  = 0.0
     ALLOCATE ( grd%chl(grd%im,grd%jm,grd%km,grd%nchl) )    ; grd%chl    = huge(grd%chl(1,1,1,1))
     ALLOCATE ( grd%chl_ad(grd%im,grd%jm,grd%km,grd%nchl) ) ; grd%chl_ad = huge(grd%chl_ad(1,1,1,1))

     ALLOCATE ( x3(grd%im,grd%jm,grd%km)) ;  x3 = huge(x3(1,1,1))
     ALLOCATE ( x2(grd%im,grd%jm))        ; x2 = huge(x2(1,1))
     ALLOCATE ( x1(grd%km) )              ;  x1 = huge(x1(1))




      stat = nf90_inq_varid (ncid, 'lon', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,grd%lon)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'lat', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,grd%lat)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'dx', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,grd%dx)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'dy', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,grd%dy)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'dz', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,grd%dz)
    if (stat /= nf90_noerr) call netcdf_err(stat)

    stat = nf90_inq_varid (ncid, 'dep', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    stat = nf90_get_var (ncid, idvar, x1)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    grd%dep(:) = x1(:)

      stat = nf90_inq_varid (ncid, 'topo', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x2)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    grd%hgt(:,:) = x2(:,:)

      stat = nf90_inq_varid (ncid, 'tmsk', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x3)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    grd%msk(:,:,:) = x3(:,:,:)

    stat = nf90_inq_varid (ncid, 'regs', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    stat = nf90_get_var (ncid, idvar, x2)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    grd%reg(:,:) = int(x2(:,:))

    stat = nf90_inq_varid (ncid, 'mdt', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    stat = nf90_get_var (ncid, idvar, x2)
    grd%mdt(:,:) = x2(:,:)
    if (stat /= nf90_noerr) call netcdf_err(stat)

    stat = nf90_close(ncid)

! ----------------------------------------------------------------

    DEALLOCATE ( x3, x2, x1 )

end subroutine rdgrd


