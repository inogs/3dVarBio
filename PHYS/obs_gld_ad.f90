subroutine obs_gld_ad


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
! Apply observational operator for gliders (adjoint)                   !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk

 do kk = 1,gld%no

  if(gld%flc(kk).eq.1 .and. gld%par(kk).eq.1 )then

    obs%k = obs%k + 1

    i=gld%ib(kk)
    j=gld%jb(kk)
    k=gld%kb(kk)

    grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + gld%pq1(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + gld%pq2(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + gld%pq3(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + gld%pq4(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + gld%pq5(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + gld%pq6(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + gld%pq7(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + gld%pq8(kk) * obs%gra(obs%k)

  else if(gld%flc(kk).eq.1 .and. gld%par(kk).eq.2 )then

    obs%k = obs%k + 1

    i=gld%ib(kk)
    j=gld%jb(kk)
    k=gld%kb(kk)

    grd%sal_ad(i  ,j  ,k  ) = grd%sal_ad(i  ,j  ,k  ) + gld%pq1(kk) * obs%gra(obs%k)
    grd%sal_ad(i+1,j  ,k  ) = grd%sal_ad(i+1,j  ,k  ) + gld%pq2(kk) * obs%gra(obs%k)
    grd%sal_ad(i  ,j+1,k  ) = grd%sal_ad(i  ,j+1,k  ) + gld%pq3(kk) * obs%gra(obs%k)
    grd%sal_ad(i+1,j+1,k  ) = grd%sal_ad(i+1,j+1,k  ) + gld%pq4(kk) * obs%gra(obs%k)
    grd%sal_ad(i  ,j  ,k+1) = grd%sal_ad(i  ,j  ,k+1) + gld%pq5(kk) * obs%gra(obs%k)
    grd%sal_ad(i+1,j  ,k+1) = grd%sal_ad(i+1,j  ,k+1) + gld%pq6(kk) * obs%gra(obs%k)
    grd%sal_ad(i  ,j+1,k+1) = grd%sal_ad(i  ,j+1,k+1) + gld%pq7(kk) * obs%gra(obs%k)
    grd%sal_ad(i+1,j+1,k+1) = grd%sal_ad(i+1,j+1,k+1) + gld%pq8(kk) * obs%gra(obs%k)

  endif

 enddo


end subroutine obs_gld_ad
