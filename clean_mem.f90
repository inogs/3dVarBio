subroutine clean_mem
  
  
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
  ! Save the result on the coarse grid                                   !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  use cns_str
  use rcfl
  use mpi_str
  
  implicit none

  integer :: ierr
  
  if (drv%argo_obs .eq. 1) then
     ! deallocate argo arrays
     DEALLOCATE ( arg%flc)
     DEALLOCATE ( arg%inc)
     DEALLOCATE ( arg%ib, arg%jb, arg%kb)
     DEALLOCATE ( arg%pq1, arg%pq2, arg%pq3, arg%pq4)
     DEALLOCATE ( arg%pq5, arg%pq6, arg%pq7, arg%pq8)
     DEALLOCATE ( arg%par)
     DEALLOCATE (grd%lon, grd%lat)
  endif

  ! chlorophyll structure
  if(drv%sat_obs .eq. 1) then
    DEALLOCATE ( sat%flc)
    DEALLOCATE ( sat%inc)
    DEALLOCATE ( sat%ib)
    DEALLOCATE ( sat%jb)
    DEALLOCATE ( sat%pq1)
    DEALLOCATE ( sat%pq2)
    DEALLOCATE ( sat%pq3)
    DEALLOCATE ( sat%pq4)
    DEALLOCATE ( sat%dzr)
  endif

  ! Constants structure
  DEALLOCATE ( rcf%al)
  DEALLOCATE ( rcf%sc)

  DEALLOCATE(SendCountX3D, SendDisplX3D)
  DEALLOCATE(RecCountX3D, RecDisplX3D)

  DEALLOCATE(ChlExtended)
  DEALLOCATE(SendBottom, RecTop)
  DEALLOCATE(SendTop, RecBottom)

  if(MyId .eq. 0) then
    write(*,*) ' ALL MEMORY CLEAN'
    write(*,*) ''
  endif

end subroutine clean_mem
