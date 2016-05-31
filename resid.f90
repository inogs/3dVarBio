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

 implicit none

 INTEGER(i4)   :: i, k

  k = 0

! ---
! Satellite observations of SLA
 do i=1,sla%no
  if(sla%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = sla%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! ARGO observations 
 do i=1,arg%no
  if(arg%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = arg%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! XBT observations 
 do i=1,xbt%no
  if(xbt%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = xbt%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Glider observations 
 do i=1,gld%no
  if(gld%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = gld%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Observations of Argo float positions
 do i=1,tra%no
  if(tra%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = tra%inx(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo
 do i=1,tra%no
  if(tra%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = tra%iny(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Observations of positions of surface drifters
 do i=1,trd%no
  if(trd%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = trd%inx(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo
 do i=1,trd%no
  if(trd%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = trd%iny(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Observations of velocity by drifters
 do i=1,vdr%no
  if(vdr%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = vdr%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Observations of velocity by gliders
 do i=1,gvl%no
  if(gvl%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = gvl%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

! ---
! Observations of chlorophyll
 do i=1,chl%no
  if(chl%flc(i).eq.1)then
   k = k + 1
   obs%inc(k) = chl%inc(i) 
   obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
  endif
 enddo

end subroutine resid
