subroutine obs_tra

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2018 Anna Teruzzi, OGS, Trieste, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
!                                                                          !
!    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Call Argo trajectory model                                           !
!                                                                      !
! Version 1: V. Taillandier, S. Dobricic 2007                          !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, ii, jj, k

    do jj =1,tra%jm
    do ii =1,tra%im

      i = tra%i1(ii,jj)
      j = tra%j1(ii,jj)

      tra%uvl(ii,jj) = tra%pq1(i,j) * grd%uvl(i  ,j  ,tra%km) +     &
                       tra%pq2(i,j) * grd%uvl(i+1,j  ,tra%km) +     &
                       tra%pq3(i,j) * grd%uvl(i  ,j+1,tra%km) +     &
                       tra%pq4(i,j) * grd%uvl(i+1,j+1,tra%km) 

      tra%vvl(ii,jj) = tra%pq1(i,j) * grd%vvl(i  ,j  ,tra%km) +     &
                       tra%pq2(i,j) * grd%vvl(i+1,j  ,tra%km) +     &
                       tra%pq3(i,j) * grd%vvl(i  ,j+1,tra%km) +     &
                       tra%pq4(i,j) * grd%vvl(i+1,j+1,tra%km) 

    enddo
    enddo

  call mod_trj_tl( tra%im,tra%jm,tra%umn,tra%vmn,tra%dx,tra%dy,tra%flg,  &
                   tra%nt,tra%no,tra%xmn,tra%ymn,tra%dtm,                &
                   tra%uvl,tra%vvl,tra%xtl,tra%ytl )


 do k=1,tra%no

  if(tra%flc(k).eq.1 )then

     tra%inx(k) = tra%xtl(k)
     tra%iny(k) = tra%ytl(k)

     i=int(tra%xmn(tra%nt+1,k)+tra%xtl(k))
     j=int(tra%ymn(tra%nt+1,k)+tra%ytl(k))
     tra%loa(k) = tra%lon(i,j) +     &
                  (tra%xmn(tra%nt+1,k)+tra%xtl(k)-i)*(tra%lon(i+1,j)-tra%lon(i,j))
     tra%laa(k) = tra%lat(i,j) +     &
                  (tra%ymn(tra%nt+1,k)+tra%ytl(k)-j)*(tra%lat(i,j+1)-tra%lat(i,j))

  endif

 enddo

end subroutine obs_tra
