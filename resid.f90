subroutine resid
  
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
  ! Calculate analysis - observation                                     !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
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
  ! Observations of chlorophyll
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
