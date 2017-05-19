subroutine rdquotas

!---------------------------------------------------------------------------
!anna                                                                          !
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
! READ quotas for biological variables                                 !
!                                                                      !
! Version 1: A.Teruzzi 2012                                            !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use bio_str
  use grd_str
  use drv_str
!  use netcdf

  implicit none

  INTEGER(i4)      :: stat, nchl, ncmps, nlevs, jjm, iim
  INTEGER(i4)      :: ncid, idvar

!     stat = nf90_open('quotas.nc', NF90_NOWRITE, ncid)
!     if (stat /= nf90_noerr) call netcdf_err(stat)
!
!
!    ! Get dimensions
!     stat = nf90_inq_dimid (ncid, 'nphy', idvar)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inquire_dimension (ncid, idvar, len = nchl)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inq_dimid (ncid, 'ncom', idvar)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inquire_dimension (ncid, idvar, len = ncmps)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inq_dimid (ncid, 'nlev', idvar)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inquire_dimension (ncid, idvar, len = nlevs)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inq_dimid (ncid, 'jm', idvar)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inquire_dimension (ncid, idvar, len = jjm)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inq_dimid (ncid, 'im', idvar)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!     stat = nf90_inquire_dimension (ncid, idvar, len = iim)
!   if (stat /= nf90_noerr) call netcdf_err(stat)
!
!   if(bio%nphy .ne. nchl) then
!     write(drv%dia,*)'Error: Number of phytoplankton types different than in the input file.'
!     call f_exit(10)
!   endif
!
!   if(bio%ncmp .ne. ncmps) then
!     write(drv%dia,*)'Error: Number of phytoplankton components different than in the input file.'
!   endif
!
!   if(grd%km .ne. nlevs) then
!     write(drv%dia,*)'Error: Number of levels with biological quotas is different than in the input file.'
!   endif
!
!   if(grd%im .ne. iim .or. grd%jm .ne. jjm) then
!     write(drv%dia,*)'Error: Grid dimensions are different than in the quotas input file.'
!   endif


   write(drv%dia,*)'Number of phytoplankton types is ', bio%nphy
   write(drv%dia,*)'Number of phytoplankton components is ', bio%ncmp

! Allocate quotas arrys
   ALLOCATE ( bio%pquot( grd%im, grd%jm, grd%km, bio%nphy))
   ALLOCATE ( bio%cquot( grd%im, grd%jm, grd%km, bio%nphy, bio%ncmp))

    bio%pquot(:,:,:,:) = dble(1)
    bio%cquot(:,:,:,:,:) = dble(1)

!    stat = nf90_inq_varid (ncid, 'pquot', idvar)
!  if (stat /= nf90_noerr) call netcdf_err(stat)
!    stat = nf90_get_var (ncid, idvar, bio%pquot)
!  if (stat /= nf90_noerr) call netcdf_err(stat)
!    stat = nf90_inq_varid (ncid, 'cquot', idvar)
!  if (stat /= nf90_noerr) call netcdf_err(stat)
!    stat = nf90_get_var (ncid, idvar, bio%cquot)
!  if (stat /= nf90_noerr) call netcdf_err(stat)
!
!   stat = nf90_close(ncid)

end subroutine rdquotas
