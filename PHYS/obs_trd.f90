subroutine obs_trd

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
! Call trajectory model of surface drifters                            !
!                                                                      !
! Version 1: V. Taillandier, S. Dobricic 2007                          !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, ii, jj, k

    do jj =1,trd%jm
    do ii =1,trd%im

      i = trd%i1(ii,jj)
      j = trd%j1(ii,jj)

      trd%uvl(ii,jj) = trd%pq1(i,j) * grd%uvl(i  ,j  ,trd%km) +     &
                       trd%pq2(i,j) * grd%uvl(i+1,j  ,trd%km) +     &
                       trd%pq3(i,j) * grd%uvl(i  ,j+1,trd%km) +     &
                       trd%pq4(i,j) * grd%uvl(i+1,j+1,trd%km) 

      trd%vvl(ii,jj) = trd%pq1(i,j) * grd%vvl(i  ,j  ,trd%km) +     &
                       trd%pq2(i,j) * grd%vvl(i+1,j  ,trd%km) +     &
                       trd%pq3(i,j) * grd%vvl(i  ,j+1,trd%km) +     &
                       trd%pq4(i,j) * grd%vvl(i+1,j+1,trd%km) 

    enddo
    enddo

  call mod_trj_tl( trd%im,trd%jm,trd%umn,trd%vmn,trd%dx,trd%dy,trd%flg,   &
                   trd%nt,trd%no,trd%xmn,trd%ymn,trd%dtm,                 &
                   trd%uvl,trd%vvl,trd%xtl,trd%ytl )


 do k=1,trd%no

  if(trd%flc(k).eq.1 )then

     trd%inx(k) = trd%xtl(k)
     trd%iny(k) = trd%ytl(k)

     i=int(trd%xmn(trd%nt+1,k)+trd%xtl(k))
     j=int(trd%ymn(trd%nt+1,k)+trd%ytl(k))
     trd%loa(k) = trd%lon(i,j) +     &
                  (trd%xmn(trd%nt+1,k)+trd%xtl(k)-i)*(trd%lon(i+1,j)-trd%lon(i,j))
     trd%laa(k) = trd%lat(i,j) +     &
                  (trd%ymn(trd%nt+1,k)+trd%ytl(k)-j)*(trd%lat(i,j+1)-trd%lat(i,j))

  endif

 enddo

end subroutine obs_trd
