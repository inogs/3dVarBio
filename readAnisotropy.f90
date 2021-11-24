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
  use grd_str
  use cns_str
  use filenames

  use mpi_str
  use pnetcdf

  implicit none

  INTEGER(i4)                    :: stat, ncid, idvar
  integer(KIND=MPI_OFFSET_KIND)  :: GlobalStart(2), GlobalCount(2),imr,jmr
  real(r4), ALLOCATABLE  :: x2(:,:)

  stat = nf90mpi_open(Var3DCommunicator, trim(ANIS_FILE), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open "//trim(ANIS_FILE),stat)

  ALLOCATE ( x2(GlobalRow,GlobalCol))
  GlobalStart(:) = 1
  GlobalCount(1) = GlobalRow
  GlobalCount(2) = GlobalCol

  ! Get dimensions 
  stat = nf90mpi_inq_dimid (ncid, 'im', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid",stat)
  stat = nfmpi_inq_dimlen (ncid, idvar, imr)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen",stat)
  stat = nf90mpi_inq_dimid (ncid, 'jm', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid",stat)
  stat = nfmpi_inq_dimlen(ncid, idvar, jmr)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen",stat)

  ! Check on dimensions
  if ((imr .ne. GlobalRow).OR.(jmr.ne.GlobalCol)) then
    write(drv%dia,*)'Error: dimensions of rcorr differ from grid ones'
    call MPI_Abort(MPI_COMM_WORLD, -1, stat)
  endif


  !  Allocate rcorr arrays
  ALLOCATE ( rcf%rtx(GlobalRow,GlobalCol))
  ALLOCATE ( rcf%rty(GlobalRow,GlobalCol))


  stat = nf90mpi_inq_varid (ncid, 'kx_n', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid",stat)
  stat = nfmpi_get_vara_real_all (ncid, idvar, GlobalStart, GlobalCount, x2)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all",stat)
  rcf%rtx(:,:) = x2(:,:)

  stat = nf90mpi_inq_varid (ncid, 'ky_n', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid",stat)
  stat = nfmpi_get_vara_real_all (ncid, idvar, GlobalStart, GlobalCount, x2)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all",stat)
  rcf%rty(:,:) = x2(:,:)

  stat = nf90mpi_close(ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_close", stat)  
  DEALLOCATE(x2)

end subroutine readAnisotropy
