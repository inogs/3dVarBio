subroutine obs_vdr_ad

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
! Apply observational operator for velocities from drifters (adjoint)  !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk

 do kk = 1,vdr%no

  if(vdr%flc(kk).eq.1 .and. vdr%par(kk).eq.1 )then

    obs%k = obs%k + 1

    i=vdr%ib(kk)
    j=vdr%jb(kk)
    k=vdr%kb(kk)

    grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + vdr%pq1(kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + vdr%pq2(kk) * obs%gra(obs%k)
    grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + vdr%pq3(kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + vdr%pq4(kk) * obs%gra(obs%k)
    grd%uvl_ad(i  ,j  ,k+1) = grd%uvl_ad(i  ,j  ,k+1) + vdr%pq5(kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j  ,k+1) = grd%uvl_ad(i+1,j  ,k+1) + vdr%pq6(kk) * obs%gra(obs%k)
    grd%uvl_ad(i  ,j+1,k+1) = grd%uvl_ad(i  ,j+1,k+1) + vdr%pq7(kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j+1,k+1) = grd%uvl_ad(i+1,j+1,k+1) + vdr%pq8(kk) * obs%gra(obs%k)

  else if(vdr%flc(kk).eq.1 .and. vdr%par(kk).eq.2 )then

    obs%k = obs%k + 1

    i=vdr%ib(kk)
    j=vdr%jb(kk)
    k=vdr%kb(kk)

    grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + vdr%pq1(kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + vdr%pq2(kk) * obs%gra(obs%k)
    grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + vdr%pq3(kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + vdr%pq4(kk) * obs%gra(obs%k)
    grd%vvl_ad(i  ,j  ,k+1) = grd%vvl_ad(i  ,j  ,k+1) + vdr%pq5(kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j  ,k+1) = grd%vvl_ad(i+1,j  ,k+1) + vdr%pq6(kk) * obs%gra(obs%k)
    grd%vvl_ad(i  ,j+1,k+1) = grd%vvl_ad(i  ,j+1,k+1) + vdr%pq7(kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j+1,k+1) = grd%vvl_ad(i+1,j+1,k+1) + vdr%pq8(kk) * obs%gra(obs%k)

  endif

 enddo


end subroutine obs_vdr_ad
