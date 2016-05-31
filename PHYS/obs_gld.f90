subroutine obs_gld


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
! Apply observational operator for gliders                             !
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

    i=gld%ib(kk)
    j=gld%jb(kk)
    k=gld%kb(kk)

    gld%inc(kk) = gld%pq1(kk) * grd%tem(i  ,j  ,k  ) +       &
                  gld%pq2(kk) * grd%tem(i+1,j  ,k  ) +       &
                  gld%pq3(kk) * grd%tem(i  ,j+1,k  ) +       &
                  gld%pq4(kk) * grd%tem(i+1,j+1,k  ) +       &
                  gld%pq5(kk) * grd%tem(i  ,j  ,k+1) +       &
                  gld%pq6(kk) * grd%tem(i+1,j  ,k+1) +       &
                  gld%pq7(kk) * grd%tem(i  ,j+1,k+1) +       &
                  gld%pq8(kk) * grd%tem(i+1,j+1,k+1)  

  else if(gld%flc(kk).eq.1 .and. gld%par(kk).eq.2 )then

    i=gld%ib(kk)
    j=gld%jb(kk)
    k=gld%kb(kk)

    gld%inc(kk) = gld%pq1(kk) * grd%sal(i  ,j  ,k  ) +       &
                  gld%pq2(kk) * grd%sal(i+1,j  ,k  ) +       &
                  gld%pq3(kk) * grd%sal(i  ,j+1,k  ) +       &
                  gld%pq4(kk) * grd%sal(i+1,j+1,k  ) +       &
                  gld%pq5(kk) * grd%sal(i  ,j  ,k+1) +       &
                  gld%pq6(kk) * grd%sal(i+1,j  ,k+1) +       &
                  gld%pq7(kk) * grd%sal(i  ,j+1,k+1) +       &
                  gld%pq8(kk) * grd%sal(i+1,j+1,k+1)  

  endif

 enddo

end subroutine obs_gld
