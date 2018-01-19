subroutine readNutCov

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


!---------------------------------------------------------------------------
! read covariances between Nitrate and Phsosphate                          !
! to update Phosphate with assimilation of Nitrate                         !
!---------------------------------------------------------------------------


  !use set_knd
  use drv_str
  use grd_str
  use bio_str
  !use cns_str
  use filenames
  !use rcfl

  use mpi_str
  use pnetcdf

  implicit none

  integer(i4)            :: stat, ncid, idvar
  integer(KIND=MPI_OFFSET_KIND)  :: GlobalStart(3), GlobalCount(3)
  real(r4), ALLOCATABLE  :: x3(:,:,:)

  !write(*,*)trim(RCORR_FILE)
  stat = nf90mpi_open(Var3DCommunicator, trim(NUTCOV_FILE), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open",stat)

  ALLOCATE ( x3(GlobalRow,GlobalCol,grd%km))
  GlobalStart(:) = 1
  GlobalCount(1) = GlobalRow
  GlobalCount(2) = GlobalCol
  GlobalCount(3) = grd%km

  stat = nf90mpi_inq_varid (ncid, 'covn3n_n1p', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid radius",stat)
  stat = nfmpi_get_vara_real_all (ncid, idvar, GlobalStart, GlobalCount, x3)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all radius",stat)

  bio%covn3n_n1p(:,:,:) = x3(:,:,:)


  stat = nf90mpi_close(ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_close", stat)

  DEALLOCATE(x3)

end subroutine readNutCov
