subroutine sav_itr
  
  
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
  
  ! ---
  ! Save grid dimensions
  
  drv%im = grd%im
  drv%jm = grd%jm
  ! Save eigenvalues
  if (1.eq.0) then ! We do not know the reason of these lines
     ALLOCATE ( drv%ro(drv%im,drv%jm,ros%neof))    ; drv%ro   (:,:,:) = grd%ro   (:,:,:)
     ALLOCATE ( drv%msk(drv%im,drv%jm))            ; drv%msk  (:,:)   = grd%msr  (:,:,1)
  endif
  ! ---
  ! Grid structure
  DEALLOCATE ( grd%reg)
  DEALLOCATE ( grd%msk)
  DEALLOCATE ( grd%dep)
  DEALLOCATE ( grd%dx, grd%dy)
  DEALLOCATE ( grd%alx )
  DEALLOCATE ( grd%aly )
  DEALLOCATE ( grd%btx )
  DEALLOCATE ( grd%bty )
  DEALLOCATE ( grd%scx )
  DEALLOCATE ( grd%scy )
  DEALLOCATE ( grd%msr )
  DEALLOCATE ( grd%imx, grd%jmx)
  DEALLOCATE ( grd%istp, grd%jstp)
  DEALLOCATE ( grd%inx, grd%jnx)
  DEALLOCATE ( grd%aex)
  DEALLOCATE ( grd%aey)
  DEALLOCATE ( grd%bex)
  DEALLOCATE ( grd%bey)
  ! Biological vectors
  DEALLOCATE ( grd%chl)
  DEALLOCATE ( grd%chl_ad)
  ! Observational vector
  DEALLOCATE ( obs%inc, obs%amo, obs%res)
  DEALLOCATE ( obs%err, obs%gra)
  ! Covariances structure
  DEALLOCATE ( grd%ro)
  DEALLOCATE ( grd%ro_ad)
  DEALLOCATE ( ros%evc, ros%eva )

  if (drv%argo .eq. 1) then
     DEALLOCATE(grd%lon, grd%lat)
     ! deallocate argo arrays
     DEALLOCATE ( arg%flc)
     DEALLOCATE ( arg%inc)
     DEALLOCATE ( arg%err)
     DEALLOCATE ( arg%res)
     DEALLOCATE ( arg%ib, arg%jb, arg%kb)
     DEALLOCATE ( arg%pq1, arg%pq2, arg%pq3, arg%pq4)
     DEALLOCATE ( arg%pq5, arg%pq6, arg%pq7, arg%pq8)
  endif

  ! Control structure
  DEALLOCATE( ctl%x_c, ctl%g_c)

  DEALLOCATE (SurfaceWaterPoints)  
  
  DEALLOCATE ( a_rcx)
  DEALLOCATE ( b_rcx)
  DEALLOCATE ( c_rcx)
  DEALLOCATE ( a_rcy)
  DEALLOCATE ( b_rcy)
  DEALLOCATE ( c_rcy)
  DEALLOCATE ( alp_rcx)
  DEALLOCATE ( bta_rcx)
  DEALLOCATE ( alp_rcy)
  DEALLOCATE ( bta_rcy)
  DEALLOCATE (Dump_chl, Dump_vip, Dump_msk)
  
  if(MyRank .eq. 0) write(*,*) ' DEALLOCATION DONE'
  
end subroutine sav_itr
