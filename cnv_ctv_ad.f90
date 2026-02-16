subroutine cnv_ctv_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Convert from control to v - adjoint                                  !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use grd_str
 use ctl_str
 use eof_str
 use drv_str
 use mpi_str

 implicit none

 INTEGER(i4)     :: i,j,k, kk
 INTEGER(i4)   :: jumpInd, indSupWP

   if (MyId .eq. 0) then
      print*, 'DIAG cnv_ctv_ad: before sum(grd%ro_ad)=', sum(grd%ro_ad), ' max=', maxval(grd%ro_ad)
      write(drv%dia,*) 'DIAG cnv_ctv_ad: before sum(grd%ro_ad)=', sum(grd%ro_ad), ' max=', maxval(grd%ro_ad)
   endif

   do k=1,ros%neof
     jumpInd =  (k -1 )* nSurfaceWaterPoints
        do indSupWP=1,nSurfaceWaterPoints
           i = SurfaceWaterPoints(1,indSupWP)
           j = SurfaceWaterPoints(2,indSupWP)
           kk = jumpInd + indSupWP
           ctl%g_c(kk) = grd%ro_ad(i,j,k)
        enddo

   enddo
   if (MyId .eq. 0) then
      print*, 'DIAG cnv_ctv_ad: after sum(ctl%g_c)=', sum(ctl%g_c), ' max=', maxval(ctl%g_c)
      write(drv%dia,*) 'DIAG cnv_ctv_ad: after sum(ctl%g_c)=', sum(ctl%g_c), ' max=', maxval(ctl%g_c)
   endif

end subroutine cnv_ctv_ad
