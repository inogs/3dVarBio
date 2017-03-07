subroutine parallel_rdrcorr

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
  use rcfl

  use mpi_str
  use pnetcdf

  implicit none

  integer(i4)            :: stat, ncid, idvar
  integer(KIND=MPI_OFFSET_KIND)  :: GlobalStart(3), GlobalCount(3)
  real(r4), ALLOCATABLE  :: x3(:,:,:)

  !write(*,*)trim(RCORR_FILE)
  stat = nf90mpi_open(Var3DCommunicator, trim(RCORR_FILE), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open",stat)

  ALLOCATE ( rcf%Lxyz(GlobalRow,GlobalCol,grd%km))
  ALLOCATE ( x3(GlobalRow,GlobalCol,grd%km))
  GlobalStart(:) = 1
  GlobalCount(1) = GlobalRow
  GlobalCount(2) = GlobalCol
  GlobalCount(3) = grd%km

  stat = nf90mpi_inq_varid (ncid, 'radius', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid radius",stat)
  stat = nfmpi_get_vara_real_all (ncid, idvar, GlobalStart, GlobalCount, x3)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all radius",stat)
  rcf%Lxyz(:,:,:) = x3(:,:,:)

  !laura from km to meter
  where (rcf%Lxyz<=0.0001) 
      rcf%Lxyz=rcf%Lxyz/1000
  end where
    rcf%Lxyz= rcf%Lxyz*1000  !from km to meter

  DEALLOCATE(x3)

end subroutine parallel_rdrcorr
