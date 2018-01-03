subroutine obs_gvl

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

 do kk = 1,gvl%no

  if(gvl%flc(kk).eq.1 .and. gvl%par(kk).eq.1 )then

    i=gvl%ib(kk)
    j=gvl%jb(kk)

     gvl%inc(kk) = 0.0

    do k=1,gvl%kb(kk)+1
     gvl%inc(kk) = gvl%inc(kk) + (                        &
                   gvl%pq1(kk) * grd%uvl(i  ,j  ,k  ) +       &
                   gvl%pq2(kk) * grd%uvl(i+1,j  ,k  ) +       &
                   gvl%pq3(kk) * grd%uvl(i  ,j+1,k  ) +       &
                   gvl%pq4(kk) * grd%uvl(i+1,j+1,k  ) ) * gvl%dzr(k,kk)
    enddo

  else if(gvl%flc(kk).eq.1 .and. gvl%par(kk).eq.2 )then

    i=gvl%ib(kk)
    j=gvl%jb(kk)

     gvl%inc(kk) = 0.0

    do k=1,gvl%kb(kk)+1
     gvl%inc(kk) = gvl%inc(kk) + (                        &
                   gvl%pq1(kk) * grd%vvl(i  ,j  ,k  ) +       &
                   gvl%pq2(kk) * grd%vvl(i+1,j  ,k  ) +       &
                   gvl%pq3(kk) * grd%vvl(i  ,j+1,k  ) +       &
                   gvl%pq4(kk) * grd%vvl(i+1,j+1,k  ) ) * gvl%dzr(k,kk)
    enddo

  endif

 enddo

end subroutine obs_gvl
