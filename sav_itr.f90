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
  use bio_str
  use da_params

  implicit none
  

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
 
  ! Chlorophyll vectors
  if(drv%multiv.eq.0) then

  if(drv%chl_assim .eq. 1) then
    DEALLOCATE( grd%chl)
    DEALLOCATE( grd%chl_ad)
  endif
  if(drv%nut .eq. 1) then
    if(bio%n3n .eq. 1) then
      DEALLOCATE( grd%n3n)
      DEALLOCATE( grd%n3n_ad)
    endif
    if(bio%o2o .eq. 1) then
      DEALLOCATE( grd%o2o)
      DEALLOCATE( grd%o2o_ad)
    endif    
  endif
  
  endif
  
  if(drv%multiv.eq.1) then
    DEALLOCATE( grd%chl)
    DEALLOCATE( grd%chl_ad)
    DEALLOCATE( grd%n3n)
    DEALLOCATE( grd%n3n_ad)
  endif
  
  ! Observational vector
  DEALLOCATE( obs%inc, obs%amo, obs%res)
  DEALLOCATE( obs%err, obs%gra)
 
  ! Covariances structure
  DEALLOCATE( grd%ro)
  DEALLOCATE( grd%ro_ad)
  if(drv%chl_assim .eq. 1) then
    DEALLOCATE( ros%evc_chl, ros%eva_chl )
  endif
  if(drv%nut .eq. 1) then
    if(bio%N3n .eq. 1) then
      DEALLOCATE( ros%evc_n3n, ros%eva_n3n )
    endif
    if(bio%O2o .eq. 1) then
      DEALLOCATE( ros%evc_o2o, ros%eva_o2o )
    endif
  endif

  ! Control structure
  DEALLOCATE( ctl%x_c, ctl%g_c)

  ! Bio structure
  if(drv%multiv.eq.0) then
    if(drv%chl_assim .eq. 1) then
      DEALLOCATE( bio%phy, bio%phy_ad)
      DEALLOCATE( bio%cquot, bio%pquot)
      DEALLOCATE( bio%InitialChl)
      if ((drv%nut .eq. 0) .and. (NNutVar .gt. 0)) then
        DEALLOCATE( bio%InitialNut)
        if (drv%chl_upnut.eq.1) then
          ! DEALLOCATE( bio%covn3n_n1p)
          DEALLOCATE( bio%covn1p_chl)
          DEALLOCATE( bio%covn3n_chl)
        endif  
      endif
    endif
    if(drv%nut .eq. 1) then
      DEALLOCATE( bio%InitialNut)
      if(bio%N3n.eq.1 .AND. bio%updateN1p.eq.1)  DEALLOCATE( bio%covn3n_n1p)
    endif
  endif

  if(drv%multiv.eq.1) then
    DEALLOCATE( bio%phy, bio%phy_ad)
    DEALLOCATE( bio%cquot, bio%pquot)
    DEALLOCATE( bio%InitialChl)
    DEALLOCATE( bio%InitialNut)
    if(bio%updateN1p.eq.1)  DEALLOCATE( bio%covn3n_n1p)
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
  DEALLOCATE( grd%global_msk)
  
  if(MyId .eq. 0) write(*,*) ' DEALLOCATION DONE'
  
end subroutine sav_itr


