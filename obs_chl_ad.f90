subroutine obs_chl_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
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
! Apply observational operator for velocities from gliders  (adjoint)  !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk, l

 do kk = 1,chl%no

  if(chl%flc(kk).eq.1)then

   obs%k = obs%k + 1

   i=chl%ib(kk)
   j=chl%jb(kk)

  do l=1,grd%nchl
   do k=1,chl%kb(kk)
    grd%chl_ad(i  ,j  ,k,l) = grd%chl_ad(i  ,j  ,k,l) + chl%pq1(kk) * chl%dzr(k,kk) * obs%gra(obs%k)
    grd%chl_ad(i+1,j  ,k,l) = grd%chl_ad(i+1,j  ,k,l) + chl%pq2(kk) * chl%dzr(k,kk) * obs%gra(obs%k)
    grd%chl_ad(i  ,j+1,k,l) = grd%chl_ad(i  ,j+1,k,l) + chl%pq3(kk) * chl%dzr(k,kk) * obs%gra(obs%k)
    grd%chl_ad(i+1,j+1,k,l) = grd%chl_ad(i+1,j+1,k,l) + chl%pq4(kk) * chl%dzr(k,kk) * obs%gra(obs%k)
   enddo
  enddo

  endif

 enddo


end subroutine obs_chl_ad
