subroutine cnv_ctv_ad

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
! Convert from control to v - adjoint                                  !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


 use grd_str
 use ctl_str
 use eof_str

 implicit none

 INTEGER(i4)     :: i,j,k, kk
 INTEGER(i4)   :: jumpInd, indSupWP

   do k=1,ros%neof
     jumpInd =  (k -1 )* nSurfaceWaterPoints
        do indSupWP=1,nSurfaceWaterPoints
           i = SurfaceWaterPoints(1,indSupWP)
           j = SurfaceWaterPoints(2,indSupWP)
           kk = jumpInd + indSupWP
           ctl%g_c(kk) = grd%ro_ad(i,j,k)
        enddo
   enddo

end subroutine cnv_ctv_ad
