subroutine sav_itr
  
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2018 Anna Teruzzi, OGS, Trieste                         !
  !                                                                          !
  !    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                          !
  !                                                                          !
  !    3DVarBio is  free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    3DVarBio is  distributed in the hope that it will be useful,           !
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
  ! Version 1: A. Teruzzi 2018                                           !
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
  use bio_str

  implicit none
  
  ! free MPI RMA Windows
  call FreeWindows

  ! ---
  ! Grid structure
  DEALLOCATE( grd%reg)
  DEALLOCATE( grd%msk)
  DEALLOCATE( grd%dep)
  DEALLOCATE( grd%dx, grd%dy)
  DEALLOCATE( grd%alx )
  DEALLOCATE( grd%aly )
  DEALLOCATE( grd%btx )
  DEALLOCATE( grd%bty )
  DEALLOCATE( grd%scx )
  DEALLOCATE( grd%scy )
  DEALLOCATE( grd%imx, grd%jmx)
  DEALLOCATE( grd%istp, grd%jstp)
  DEALLOCATE( grd%inx, grd%jnx)
  DEALLOCATE( grd%aex)
  DEALLOCATE( grd%aey)
  DEALLOCATE( grd%bex)
  DEALLOCATE( grd%bey)
 
  ! Biological vectors
  DEALLOCATE( grd%chl)
  DEALLOCATE( grd%chl_ad)
 
  ! Observational vector
  DEALLOCATE( obs%inc, obs%amo, obs%res)
  DEALLOCATE( obs%err, obs%gra)
 
  ! Covariances structure
  DEALLOCATE( grd%ro)
  DEALLOCATE( grd%ro_ad)
  DEALLOCATE( ros%evc, ros%eva )

  ! Control structure
  DEALLOCATE( ctl%x_c, ctl%g_c)

  ! Bio structure
  if(drv%bio_assim .eq. 1) then
    DEALLOCATE( bio%phy, bio%phy_ad)
    DEALLOCATE( bio%cquot, bio%pquot)
    DEALLOCATE( bio%InitialChl)
  endif

  DEALLOCATE(SurfaceWaterPoints)  
  
  DEALLOCATE( a_rcx)
  DEALLOCATE( b_rcx)
  DEALLOCATE( c_rcx)
  DEALLOCATE( a_rcy)
  DEALLOCATE( b_rcy)
  DEALLOCATE( c_rcy)
  DEALLOCATE( alp_rcx)
  DEALLOCATE( bta_rcx)
  DEALLOCATE( alp_rcy)
  DEALLOCATE( bta_rcy)
  
  if(MyId .eq. 0) write(*,*) ' DEALLOCATION DONE'
  
end subroutine sav_itr

subroutine FreeWindows

  use grd_str
  use mpi_str

  implicit none

  integer ierr

  call MPI_Win_free(MpiWinChl, ierr)
  call MPI_Win_free(MpiWinChlAd, ierr)

end subroutine FreeWindows