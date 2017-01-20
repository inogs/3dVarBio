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
  use mpi
  
  implicit none

  integer :: ierr
  
  ! Deallocate everithing related to the old grid
  DEALLOCATE ( drv%grid, drv%ratco, drv%ratio)
  DEALLOCATE ( drv%mask, drv%dda, drv%ddi)

  ! chlorophyll structure
  DEALLOCATE ( chl%flg)
  DEALLOCATE ( chl%flc)
  DEALLOCATE ( chl%inc)
  DEALLOCATE ( chl%err)
  DEALLOCATE ( chl%res)
  DEALLOCATE ( chl%ib)
  DEALLOCATE ( chl%pb)
  DEALLOCATE ( chl%jb)
  DEALLOCATE ( chl%qb)
  DEALLOCATE ( chl%pq1)
  DEALLOCATE ( chl%pq2)
  DEALLOCATE ( chl%pq3)
  DEALLOCATE ( chl%pq4)
  DEALLOCATE ( chl%dzr)

  ! Constants structure
  DEALLOCATE ( rcf%al)
  DEALLOCATE ( rcf%sc)

  DEALLOCATE(SendCountX2D, SendCountX4D)
  DEALLOCATE(SendDisplX2D, SendDisplX4D)
  DEALLOCATE(RecCountX2D, RecCountX4D)
  DEALLOCATE(RecDisplX2D, RecDisplX4D)

  ! DEALLOCATE(ChlExtended)
  ! DEALLOCATE(SendRight, RecLeft)
  ! DEALLOCATE(SendLeft, RecRight)
  ! DEALLOCATE(SendBottom, RecTop)
  ! DEALLOCATE(SendTop, RecBottom)

  ! DEALLOCATE(ChlExtended4D, ChlExtendedAD_4D)
  ! DEALLOCATE(SendLeft2D, RecRight2D)
  ! DEALLOCATE(SendRight2D, RecLeft2D)
  ! DEALLOCATE(SendTop2D, RecBottom2D)
  ! DEALLOCATE(SendBottom2D, RecTop2D)

  call MPI_Comm_free(CommSliceX, ierr)
  call MPI_Comm_free(CommSliceY, ierr)
  if(MyRank .eq. 0) then
    write(*,*) ' ALL MEMORY CLEAN'
    write(*,*) ''
  endif

end subroutine clean_mem
