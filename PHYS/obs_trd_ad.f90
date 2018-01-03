subroutine obs_trd_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2018 Anna Teruzzi, OGS, Trieste, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
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
! Call trajectory model of surface drifters (adjoint)                  !
!                                                                      !
! Version 1: V. Taillandier, S. Dobricic 2007                          !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, ii, jj, k

 do k=1,trd%no

  if(trd%flc(k).eq.1 )then

     obs%k = obs%k + 1

     trd%xtl_ad(k) = obs%gra(obs%k)

  endif

 enddo

 do k=1,trd%no

  if(trd%flc(k).eq.1 )then

     obs%k = obs%k + 1

     trd%ytl_ad(k) = obs%gra(obs%k)

  endif

 enddo

  call mod_trj_ad( trd%im,trd%jm,trd%umn,trd%vmn,trd%dx,trd%dy,trd%flg,   &
                   trd%nt,trd%no,trd%xmn,trd%ymn,trd%dtm,                 &
                   trd%uvl_ad,trd%vvl_ad,trd%xtl_ad,trd%ytl_ad )

    do jj =1,trd%jm
    do ii =1,trd%im

      i = trd%i1(ii,jj)
      j = trd%j1(ii,jj)

      grd%uvl_ad(i  ,j  ,trd%km) = grd%uvl_ad(i  ,j  ,trd%km) + trd%uvl_ad(ii,jj)*trd%pq1(i,j)  
      grd%uvl_ad(i+1,j  ,trd%km) = grd%uvl_ad(i+1,j  ,trd%km) + trd%uvl_ad(ii,jj)*trd%pq2(i,j)  
      grd%uvl_ad(i  ,j+1,trd%km) = grd%uvl_ad(i  ,j+1,trd%km) + trd%uvl_ad(ii,jj)*trd%pq3(i,j)  
      grd%uvl_ad(i+1,j+1,trd%km) = grd%uvl_ad(i+1,j+1,trd%km) + trd%uvl_ad(ii,jj)*trd%pq4(i,j)  

      grd%vvl_ad(i  ,j  ,trd%km) = grd%vvl_ad(i  ,j  ,trd%km) + trd%vvl_ad(ii,jj)*trd%pq1(i,j)  
      grd%vvl_ad(i+1,j  ,trd%km) = grd%vvl_ad(i+1,j  ,trd%km) + trd%vvl_ad(ii,jj)*trd%pq2(i,j)  
      grd%vvl_ad(i  ,j+1,trd%km) = grd%vvl_ad(i  ,j+1,trd%km) + trd%vvl_ad(ii,jj)*trd%pq3(i,j)  
      grd%vvl_ad(i+1,j+1,trd%km) = grd%vvl_ad(i+1,j+1,trd%km) + trd%vvl_ad(ii,jj)*trd%pq4(i,j)  

    enddo
    enddo


end subroutine obs_trd_ad
