subroutine obs_xbt_ad

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
! Apply observational operator for ARGO floats (adjoint)               !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk

 do kk = 1,xbt%no

  if(xbt%flc(kk).eq.1 .and. xbt%par(kk).eq.1 )then

    obs%k = obs%k + 1

    i=xbt%ib(kk)
    j=xbt%jb(kk)
    k=xbt%kb(kk)

    grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + xbt%pq1(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + xbt%pq2(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + xbt%pq3(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + xbt%pq4(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + xbt%pq5(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + xbt%pq6(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + xbt%pq7(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + xbt%pq8(kk) * obs%gra(obs%k)


  endif

 enddo


end subroutine obs_xbt_ad
