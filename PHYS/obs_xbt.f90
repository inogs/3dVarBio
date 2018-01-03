subroutine obs_xbt

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
! Apply observational operator for XBT profiles                        !
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

    i=xbt%ib(kk)
    j=xbt%jb(kk)
    k=xbt%kb(kk)

    xbt%inc(kk) = xbt%pq1(kk) * grd%tem(i  ,j  ,k  ) +       &
                  xbt%pq2(kk) * grd%tem(i+1,j  ,k  ) +       &
                  xbt%pq3(kk) * grd%tem(i  ,j+1,k  ) +       &
                  xbt%pq4(kk) * grd%tem(i+1,j+1,k  ) +       &
                  xbt%pq5(kk) * grd%tem(i  ,j  ,k+1) +       &
                  xbt%pq6(kk) * grd%tem(i+1,j  ,k+1) +       &
                  xbt%pq7(kk) * grd%tem(i  ,j+1,k+1) +       &
                  xbt%pq8(kk) * grd%tem(i+1,j+1,k+1)

  endif

 enddo

end subroutine obs_xbt
