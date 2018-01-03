subroutine obs_sla_ad

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
! Apply observational operator for SLA (adjoint)                       !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k

 do k=1,sla%no

  if(sla%flc(k).eq.1)then

    obs%k = obs%k + 1

    i=sla%ib(k)
    j=sla%jb(k)

    grd%eta_ad(i  ,j  ) = grd%eta_ad(i  ,j  ) + sla%pq1(k) * obs%gra(obs%k)
    grd%eta_ad(i+1,j  ) = grd%eta_ad(i+1,j  ) + sla%pq2(k) * obs%gra(obs%k)
    grd%eta_ad(i  ,j+1) = grd%eta_ad(i  ,j+1) + sla%pq3(k) * obs%gra(obs%k)
    grd%eta_ad(i+1,j+1) = grd%eta_ad(i+1,j+1) + sla%pq4(k) * obs%gra(obs%k)

  endif

 enddo


end subroutine obs_sla_ad
