subroutine bio_conv_ad

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
  ! Biological model                                                     !
  !                                                                      !
  ! Version 1: A.Teruzzi 2012                                           !
  !-----------------------------------------------------------------------

  use grd_str
  use bio_str

  implicit none

  INTEGER(i4)   ::  i, j, k, l

  bio%phy_ad(:,:,:,:,:) = 0.0

  do l = 1,bio%nphy
    do k = 1,grd%km
      do j = 1,grd%jm
        do i = 1,grd%im
          bio%phy_ad(i,j,k,l,1) = bio%phy_ad(i,j,k,l,1) + grd%chl_ad(i,j,k)
        enddo
      enddo
    enddo
  enddo

  

end subroutine bio_conv_ad