subroutine obs_tra_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
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
! Call Argo trajectory model (adjoint)                                 !
!                                                                      !
! Version 1: V. Taillandier, S. Dobricic 2007                          !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, ii, jj, k

 do k=1,tra%no

  if(tra%flc(k).eq.1 )then

     obs%k = obs%k + 1

     tra%xtl_ad(k) = obs%gra(obs%k)

  endif

 enddo

 do k=1,tra%no

  if(tra%flc(k).eq.1 )then

     obs%k = obs%k + 1

     tra%ytl_ad(k) = obs%gra(obs%k)

  endif

 enddo

  call mod_trj_ad( tra%im,tra%jm,tra%umn,tra%vmn,tra%dx,tra%dy,tra%flg,  &
                   tra%nt,tra%no,tra%xmn,tra%ymn,tra%dtm,                &
                   tra%uvl_ad,tra%vvl_ad,tra%xtl_ad,tra%ytl_ad )

    do jj =1,tra%jm
    do ii =1,tra%im

      i = tra%i1(ii,jj)
      j = tra%j1(ii,jj)

      grd%uvl_ad(i  ,j  ,tra%km) = grd%uvl_ad(i  ,j  ,tra%km) + tra%uvl_ad(ii,jj)*tra%pq1(i,j)  
      grd%uvl_ad(i+1,j  ,tra%km) = grd%uvl_ad(i+1,j  ,tra%km) + tra%uvl_ad(ii,jj)*tra%pq2(i,j)  
      grd%uvl_ad(i  ,j+1,tra%km) = grd%uvl_ad(i  ,j+1,tra%km) + tra%uvl_ad(ii,jj)*tra%pq3(i,j)  
      grd%uvl_ad(i+1,j+1,tra%km) = grd%uvl_ad(i+1,j+1,tra%km) + tra%uvl_ad(ii,jj)*tra%pq4(i,j)  

      grd%vvl_ad(i  ,j  ,tra%km) = grd%vvl_ad(i  ,j  ,tra%km) + tra%vvl_ad(ii,jj)*tra%pq1(i,j)  
      grd%vvl_ad(i+1,j  ,tra%km) = grd%vvl_ad(i+1,j  ,tra%km) + tra%vvl_ad(ii,jj)*tra%pq2(i,j)  
      grd%vvl_ad(i  ,j+1,tra%km) = grd%vvl_ad(i  ,j+1,tra%km) + tra%vvl_ad(ii,jj)*tra%pq3(i,j)  
      grd%vvl_ad(i+1,j+1,tra%km) = grd%vvl_ad(i+1,j+1,tra%km) + tra%vvl_ad(ii,jj)*tra%pq4(i,j)  

    enddo
    enddo


end subroutine obs_tra_ad
