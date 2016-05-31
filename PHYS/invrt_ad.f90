subroutine invrt_ad( im, jm, fld, msk, rgh, a1, a2, a3, a4, a0,           &
                  bnm, ovr, resem, ncnt, itr )


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
! Implicit solver - overrelaxation (adjoint)                           !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
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


   res(:,:) = 0.0

      do icnt=itr,1,-1

         do j=2,jm-1,2
          do i=3,im-1,2
            res(i,j) = res(i,j) + fld(i,j)*ovr/a0(i,j)*msk(i,j)
            fld(i,j-1) = fld(i,j-1) + res(i,j)*a4(i,j)
            fld(i,j+1) = fld(i,j+1) + res(i,j)*a3(i,j)
            fld(i-1,j) = fld(i-1,j) + res(i,j)*a2(i,j)
            fld(i+1,j) = fld(i+1,j) + res(i,j)*a1(i,j)
            fld(i,j) = fld(i,j) - res(i,j)*a0(i,j)
            rgh(i,j) = rgh(i,j) - res(i,j)
            res(i,j) = 0.0
          enddo
         enddo
         do j=3,jm-1,2
          do i=2,im-1,2
            res(i,j) = res(i,j) + fld(i,j)*ovr/a0(i,j)*msk(i,j)
            fld(i,j-1) = fld(i,j-1) + res(i,j)*a4(i,j)
            fld(i,j+1) = fld(i,j+1) + res(i,j)*a3(i,j)
            fld(i-1,j) = fld(i-1,j) + res(i,j)*a2(i,j)
            fld(i+1,j) = fld(i+1,j) + res(i,j)*a1(i,j)
            fld(i,j) = fld(i,j) - res(i,j)*a0(i,j)
            rgh(i,j) = rgh(i,j) - res(i,j)
            res(i,j) = 0.0
          enddo
         enddo
         do j=3,jm-1,2
          do i=3,im-1,2
            res(i,j) = res(i,j) + fld(i,j)*ovr/a0(i,j)*msk(i,j)
            fld(i,j-1) = fld(i,j-1) + res(i,j)*a4(i,j)
            fld(i,j+1) = fld(i,j+1) + res(i,j)*a3(i,j)
            fld(i-1,j) = fld(i-1,j) + res(i,j)*a2(i,j)
            fld(i+1,j) = fld(i+1,j) + res(i,j)*a1(i,j)
            fld(i,j) = fld(i,j) - res(i,j)*a0(i,j)
            rgh(i,j) = rgh(i,j) - res(i,j)
            res(i,j) = 0.0
          enddo
         enddo
         do j=2,jm-1,2
          do i=2,im-1,2
            res(i,j) = res(i,j) + fld(i,j)*ovr/a0(i,j)*msk(i,j)
            fld(i,j-1) = fld(i,j-1) + res(i,j)*a4(i,j)
            fld(i,j+1) = fld(i,j+1) + res(i,j)*a3(i,j)
            fld(i-1,j) = fld(i-1,j) + res(i,j)*a2(i,j)
            fld(i+1,j) = fld(i+1,j) + res(i,j)*a1(i,j)
            fld(i,j) = fld(i,j) - res(i,j)*a0(i,j)
            rgh(i,j) = rgh(i,j) - res(i,j)
            res(i,j) = 0.0
          enddo
         enddo

            fld(:,:) = fld(:,:) * msk(:,:)

      enddo



end subroutine invrt_ad
