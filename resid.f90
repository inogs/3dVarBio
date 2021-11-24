subroutine resid
  
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
  ! Calculate analysis - observation                                     !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------

  
  use set_knd
  use obs_str
  use drv_str
  
  implicit none
  
  INTEGER(i4)   :: i, k
  
  k = 0
  
  ! ---
  ! ARGO observations
  if (drv%argo_obs .eq. 1) then
     do i=1,arg%no
        if(arg%flc(i).eq.1)then
           k = k + 1
           obs%inc(k) = arg%inc(i)
           obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
        endif
     enddo
  endif

  
  ! ---
  ! Observations of satellite chlorophyll
  if(drv%sat_obs .eq. 1) then
    do i=1,sat%no
     if(sat%flc(i).eq.1)then
        k = k + 1
        obs%inc(k) = sat%inc(i) 
        obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
     endif
    enddo
  endif
  
end subroutine resid
