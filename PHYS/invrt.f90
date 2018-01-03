subroutine invrt( im, jm, fld, msk, rgh, a1, a2, a3, a4, a0,           &
                  bnm, ovr, resem, ncnt, itr )

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
! Implicit solver - overrelaxation                                     !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


  use set_knd

  implicit none

 INTEGER(i4)    :: im, jm

 REAL(r8)       :: fld(im,jm), msk(im,jm)
 REAL(r8)       :: a0(im,jm), a1(im,jm), a2(im,jm), a3(im,jm), a4(im,jm)
 REAL(r8)       :: rgh(im,jm)
 REAL(r8)       :: res(im,jm)

 REAL(r8)       :: bnm, reser, ovr, resem
 INTEGER(i4)    :: i, j, icnt, ncnt, itr




            reser = 1.e20

       itr = 0

      do icnt=1,ncnt

       if(reser.gt.resem)then

         itr = itr + 1

         do j=2,jm-1,2
          do i=2,im-1,2
            res(i,j) = a1(i,j)*fld(i+1,j  )+a2(i,j)*fld(i-1,j  )+        &
                       a3(i,j)*fld(i  ,j+1)+a4(i,j)*fld(i  ,j-1)-        &
                       a0(i,j)*fld(i  ,j  ) - rgh(i,j)
               fld(i,j) = fld(i,j)+ovr*res(i,j)/a0(i,j)*msk(i,j)
          enddo
         enddo
         do j=3,jm-1,2
          do i=3,im-1,2
            res(i,j) = a1(i,j)*fld(i+1,j  )+a2(i,j)*fld(i-1,j  )+        &
                       a3(i,j)*fld(i  ,j+1)+a4(i,j)*fld(i  ,j-1)-        &
                       a0(i,j)*fld(i  ,j  ) - rgh(i,j)
               fld(i,j) = fld(i,j)+ovr*res(i,j)/a0(i,j)*msk(i,j)
          enddo
         enddo
         do j=3,jm-1,2
          do i=2,im-1,2
            res(i,j) = a1(i,j)*fld(i+1,j  )+a2(i,j)*fld(i-1,j  )+        &
                       a3(i,j)*fld(i  ,j+1)+a4(i,j)*fld(i  ,j-1)-        &
                       a0(i,j)*fld(i  ,j  ) - rgh(i,j)
               fld(i,j) = fld(i,j)+ovr*res(i,j)/a0(i,j)*msk(i,j)
          enddo
         enddo
         do j=2,jm-1,2
          do i=3,im-1,2
            res(i,j) = a1(i,j)*fld(i+1,j  )+a2(i,j)*fld(i-1,j  )+        &
                       a3(i,j)*fld(i  ,j+1)+a4(i,j)*fld(i  ,j-1)-        &
                       a0(i,j)*fld(i  ,j  ) - rgh(i,j)
               fld(i,j) = fld(i,j)+ovr*res(i,j)/a0(i,j)*msk(i,j)
          enddo
         enddo

            reser = 0.0
         do j=2,jm-1
          do i=2,im-1
            reser = reser + abs(res(i,j))
          enddo
         enddo
            reser = reser/bnm

       endif

      enddo



end subroutine invrt
