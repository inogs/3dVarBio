subroutine cnv_inn
  
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
  ! Convert w to correction in physical space                            !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  use drv_str
  use bio_str
  
  implicit none
  
  ! --------
  ! Convert the control vector to v
  call cnv_ctv
  
  if(drv%chl_assim .eq. 1) then
    call ver_hor_chl
  endif
  if(drv%nut .eq. 1) then
    if(bio%N3n .eq. 1) then
      call ver_hor_nut(grd%n3n, grd%n3n_ad,'N')
    endif
    if(bio%O2o .eq. 1) then
      call ver_hor_nut(grd%o2o, grd%o2o_ad,'O')
    endif
  endif
  
  ! ---
  ! Apply biological repartition of the chlorophyll
  call bio_mod

end subroutine cnv_inn
