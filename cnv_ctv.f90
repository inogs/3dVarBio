subroutine cnv_ctv

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
! Convert from control to v                                            !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use ctl_str
 use eof_str

 implicit none

 INTEGER(i4)   :: i,j,k, kk
 INTEGER(i4)   :: jumpInd, indSupWP
! INTEGER(i4) mycounter
!       kk = 0
!   do k=1,ros%neof
!    do j=1,grd%jm
!     do i=1,grd%im
!       kk = kk+1
!       grd%ro(i,j,k) = ctl%x_c(kk)
!     enddo
!    enddo
!   enddo
!mycounter = 0


   do k=1,ros%neof
   jumpInd =  (k -1 )*nSurfaceWaterPoints
       do indSupWP = 1,nSurfaceWaterPoints
           i = SurfaceWaterPoints(1,indSupWP)
           j = SurfaceWaterPoints(2,indSupWP)
           kk = jumpInd + indSupWP
           grd%ro(i,j,k) = ctl%x_c(kk)
       enddo
   enddo




end subroutine cnv_ctv
