subroutine obs_chl

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
! Apply observational operator for chlorophyll                         !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, k, l, kk

 do kk = 1,chl%no

  if(chl%flc(kk).eq.1 )then


    i=chl%ib(kk)
    j=chl%jb(kk)

     chl%inc(kk) = 0.0

   do l=1,grd%nchl
    do k=1,chl%kb(kk)
     chl%inc(kk) = chl%inc(kk) + (                        &
                   chl%pq1(kk) * grd%chl(i  ,j  ,k,l) +       &
                   chl%pq2(kk) * grd%chl(i+1,j  ,k,l) +       &
                   chl%pq3(kk) * grd%chl(i  ,j+1,k,l) +       &
                   chl%pq4(kk) * grd%chl(i+1,j+1,k,l) ) * chl%dzr(k,kk)
    enddo
   enddo

  endif

 enddo

end subroutine obs_chl
