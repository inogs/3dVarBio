subroutine obs_vdr

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
! Apply observational operator for velocities from drifters            !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk

 do kk = 1,vdr%no

  if(vdr%flc(kk).eq.1 .and. vdr%par(kk).eq.1 )then

    i=vdr%ib(kk)
    j=vdr%jb(kk)
    k=vdr%kb(kk)

    vdr%inc(kk) = vdr%pq1(kk) * grd%uvl(i  ,j  ,k  ) +       &
                  vdr%pq2(kk) * grd%uvl(i+1,j  ,k  ) +       &
                  vdr%pq3(kk) * grd%uvl(i  ,j+1,k  ) +       &
                  vdr%pq4(kk) * grd%uvl(i+1,j+1,k  ) +       &
                  vdr%pq5(kk) * grd%uvl(i  ,j  ,k+1) +       &
                  vdr%pq6(kk) * grd%uvl(i+1,j  ,k+1) +       &
                  vdr%pq7(kk) * grd%uvl(i  ,j+1,k+1) +       &
                  vdr%pq8(kk) * grd%uvl(i+1,j+1,k+1)  

  else if(vdr%flc(kk).eq.1 .and. vdr%par(kk).eq.2 )then

    i=vdr%ib(kk)
    j=vdr%jb(kk)
    k=vdr%kb(kk)

    vdr%inc(kk) = vdr%pq1(kk) * grd%vvl(i  ,j  ,k  ) +       &
                  vdr%pq2(kk) * grd%vvl(i+1,j  ,k  ) +       &
                  vdr%pq3(kk) * grd%vvl(i  ,j+1,k  ) +       &
                  vdr%pq4(kk) * grd%vvl(i+1,j+1,k  ) +       &
                  vdr%pq5(kk) * grd%vvl(i  ,j  ,k+1) +       &
                  vdr%pq6(kk) * grd%vvl(i+1,j  ,k+1) +       &
                  vdr%pq7(kk) * grd%vvl(i  ,j+1,k+1) +       &
                  vdr%pq8(kk) * grd%vvl(i+1,j+1,k+1)  

  endif

 enddo

end subroutine obs_vdr
