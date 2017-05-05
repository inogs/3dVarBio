subroutine readAnisotropy

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2015 Anna teruzzi, OGS Trieste                              !
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


  use set_knd
  use drv_str
  use netcdf
  use grd_str
  use cns_str
  use filenames
  use mpi_str

  implicit none

  INTEGER(i4)                    :: stat, ncid, idvar,imr,jmr

    stat = nf90_open(trim(ANIS_FILE), NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)

! Get dimensions 
      stat = nf90_inq_dimid (ncid, 'im', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = imr)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'jm', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = jmr)
    if (stat /= nf90_noerr) call netcdf_err(stat)

! Check on dimensions
    if ((imr .ne. GlobalRow).OR.(jmr.ne.GlobalCol)) then
       write(drv%dia,*)'Error: dimensions of rcorr different from grid ones'
       call f_exit(24)
    endif


!  Allocate rcorr arrays
     ALLOCATE ( rcf%rtx(GlobalRow,GlobalCol))
     ALLOCATE ( rcf%rty(GlobalRow,GlobalCol))


       stat = nf90_inq_varid (ncid, 'kx_n', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,rcf%rtx)
    if (stat /= nf90_noerr) call netcdf_err(stat)
       stat = nf90_inq_varid (ncid, 'ky_n', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,rcf%rty)
    if (stat /= nf90_noerr) call netcdf_err(stat)




end subroutine readAnisotropy