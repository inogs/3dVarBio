subroutine obs_gvl_ad

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
! Apply observational operator for velocities from gliders  (adjoint)  !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk

 do kk = 1,gvl%no

  if(gvl%flc(kk).eq.1 .and. gvl%par(kk).eq.1 )then

   obs%k = obs%k + 1

   i=gvl%ib(kk)
   j=gvl%jb(kk)

   do k=1,gvl%kb(kk)+1
    grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
   enddo

  else if(gvl%flc(kk).eq.1 .and. gvl%par(kk).eq.2 )then

    obs%k = obs%k + 1

    i=gvl%ib(kk)
    j=gvl%jb(kk)

   do k=1,gvl%kb(kk)+1
    grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
    grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
   enddo

  endif

 enddo


end subroutine obs_gvl_ad
