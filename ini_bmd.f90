subroutine ini_bmd

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
! Initialise the barotropic model                                      !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  INTEGER(i4)    :: i, j

       bmd%g     = 9.81
       bmd%nstp  = int( 24. * 3600. / bmd%dt )
       bmd%nstps = int( bmd%ndy * bmd%nstp )
       bmd%nstpa = int( bmd%ady * bmd%nstp )
       bmd%alp2  = 1.0 - bmd%alp1


       bmd%df1 = bmd%fc1 * grd%adxdy**2
       bmd%df2 = bmd%fc2 * grd%adxdy**2

     ALLOCATE ( bmd%itr(bmd%nstps) ) ; bmd%itr = huge(bmd%itr(1))
     ALLOCATE ( bmd%mst(grd%im,grd%jm)); bmd%mst = huge(bmd%mst(1,1))
     ALLOCATE ( bmd%msu(grd%im,grd%jm)); bmd%msu = huge(bmd%msu(1,1))
     ALLOCATE ( bmd%msv(grd%im,grd%jm)); bmd%msv = huge(bmd%msv(1,1))
     ALLOCATE ( bmd%hgt(grd%im,grd%jm)); bmd%hgt = huge(bmd%hgt(1,1))
     ALLOCATE ( bmd%hgu(grd%im,grd%jm)); bmd%hgu = huge(bmd%hgu(1,1))
     ALLOCATE ( bmd%hgv(grd%im,grd%jm)); bmd%hgv = huge(bmd%hgv(1,1))
     ALLOCATE ( bmd%dxu(grd%im,grd%jm)); bmd%dxu = huge(bmd%dxu(1,1))
     ALLOCATE ( bmd%dxv(grd%im,grd%jm)); bmd%dxv = huge(bmd%dxv(1,1))
     ALLOCATE ( bmd%dyu(grd%im,grd%jm)); bmd%dyu = huge(bmd%dyu(1,1))
     ALLOCATE ( bmd%dyv(grd%im,grd%jm)); bmd%dyv = huge(bmd%dyv(1,1))
     ALLOCATE ( bmd%a1(grd%im,grd%jm)) ; bmd%a1  = huge( bmd%a1(1,1))
     ALLOCATE ( bmd%a2(grd%im,grd%jm)) ; bmd%a2  = huge( bmd%a2(1,1))
     ALLOCATE ( bmd%a3(grd%im,grd%jm)) ; bmd%a3  = huge( bmd%a3(1,1))
     ALLOCATE ( bmd%a4(grd%im,grd%jm))  ; bmd%a4  = huge( bmd%a4(1,1))
     ALLOCATE ( bmd%a0(grd%im,grd%jm))  ; bmd%a0  = huge( bmd%a0(1,1))
     ALLOCATE ( bmd%a00(grd%im,grd%jm)); bmd%a00 = huge(bmd%a00(1,1))
     ALLOCATE ( bmd%bx(grd%im,grd%jm)) ; bmd%bx  = huge(bmd%bx(1,1))
     ALLOCATE ( bmd%by(grd%im,grd%jm)) ; bmd%by  = huge(bmd%by(1,1))
     ALLOCATE ( bmd%b_x(grd%im,grd%jm,grd%km)) ; bmd%b_x = huge(bmd%b_x(1,1,1))
     ALLOCATE ( bmd%b_y(grd%im,grd%jm,grd%km)) ; bmd%b_y = huge(bmd%b_y(1,1,1))
     ALLOCATE ( bmd%dns(grd%im,grd%jm,grd%km)) ; bmd%dns = huge(bmd%dns(1,1,1))
     ALLOCATE ( bmd%bxby(grd%im,grd%jm))  ; bmd%bxby = huge(bmd%bxby(1,1))
     ALLOCATE ( bmd%rgh(grd%im,grd%jm))   ; bmd%rgh  = huge(bmd%rgh(1,1))
     ALLOCATE ( bmd%etb(grd%im,grd%jm))    ; bmd%etb = huge(bmd%etb(1,1))
     ALLOCATE ( bmd%ub(grd%im,grd%jm))    ; bmd%ub  = huge(bmd%ub(1,1))
     ALLOCATE ( bmd%vb(grd%im,grd%jm))    ; bmd%vb  = huge(bmd%vb(1,1))
     ALLOCATE ( bmd%etn(grd%im,grd%jm))   ; bmd%etn = huge(bmd%etn(1,1))
     ALLOCATE ( bmd%un(grd%im,grd%jm))    ; bmd%un  = huge(bmd%un(1,1))
     ALLOCATE ( bmd%vn(grd%im,grd%jm))    ; bmd%vn  = huge(bmd%vn(1,1))
     ALLOCATE ( bmd%eta(grd%im,grd%jm))    ; bmd%eta = huge(bmd%eta(1,1))
     ALLOCATE ( bmd%ua(grd%im,grd%jm))    ; bmd%ua  = huge(bmd%ua(1,1))
     ALLOCATE ( bmd%va(grd%im,grd%jm))    ; bmd%va  = huge(bmd%va(1,1))
     ALLOCATE ( bmd%etm(grd%im,grd%jm))   ; bmd%etm = huge(bmd%etm(1,1))
     ALLOCATE ( bmd%um(grd%im,grd%jm))    ; bmd%um  = huge(bmd%um(1,1))
     ALLOCATE ( bmd%vm(grd%im,grd%jm))    ; bmd%vm  = huge(bmd%vm(1,1))
     ALLOCATE ( bmd%div(grd%im,grd%jm))   ; bmd%div = huge(bmd%div(1,1))
     ALLOCATE ( bmd%cu(grd%im,grd%jm))    ; bmd%cu  = huge(bmd%cu(1,1))
     ALLOCATE ( bmd%cv(grd%im,grd%jm))    ; bmd%cv  = huge(bmd%cv(1,1))
     ALLOCATE ( bmd%dux(grd%im,grd%jm))   ; bmd%dux = huge(bmd%dux(1,1))
     ALLOCATE ( bmd%duy(grd%im,grd%jm))   ; bmd%duy = huge(bmd%duy(1,1))
     ALLOCATE ( bmd%dvx(grd%im,grd%jm))   ; bmd%dvx = huge(bmd%dvx(1,1))
     ALLOCATE ( bmd%dvy(grd%im,grd%jm))   ; bmd%dvy = huge(bmd%dvy(1,1))
     ALLOCATE ( bmd%etx(grd%im,grd%jm))   ; bmd%etx = huge(bmd%etx(1,1))
     ALLOCATE ( bmd%ety(grd%im,grd%jm))   ; bmd%ety = huge(bmd%ety(1,1));


     bmd%itr(:) = 0

     do j=1,grd%jm
      do i=1,grd%im
       if(grd%hgt(i,j).gt.0.0) then
!         bmd%hgt(i,j) = max(15.,grd%hgt(i,j))
          bmd%hgt(i,j) = grd%hgt(i,j)
       else
          bmd%hgt(i,j) = 0.0
       endif
      enddo
     enddo
!          bmd%hgt(:,:) = grd%hgt(:,:)

     bmd%hgt(1,:) = 0.0
     bmd%hgt(:,1) = 0.0
     bmd%hgt(grd%im,:) = 0.0
     bmd%hgt(:,grd%jm) = 0.0

!-------------------------------------------------

     do j=1,grd%jm
      do i=1,grd%im
       if(bmd%hgt(i,j).gt.0.0) then
         bmd%mst(i,j) = 1.0
       else
         bmd%mst(i,j) = 0.0
       endif
      enddo
     enddo

     bmd%bnm = 0.0
     do j=2,grd%jm-1
      do i=2,grd%im-1
       if(bmd%mst(i,j).eq.1.) bmd%bnm = bmd%bnm + 1.
      enddo
     enddo

     do j=1,grd%jm
      do i=2,grd%im
!       bmd%hgu(i,j) = (bmd%hgt(i,j)+bmd%hgt(i-1,j))*0.5
       bmd%hgu(i,j) = min(bmd%hgt(i,j),bmd%hgt(i-1,j))
       bmd%dxu(i,j) = (grd%dx(i,j)+grd%dx(i-1,j))*0.5
       bmd%dyu(i,j) = (grd%dy(i,j)+grd%dy(i-1,j))*0.5
       bmd%msu(i,j) =  bmd%mst(i,j)*bmd%mst(i-1,j)
      enddo
     enddo
       bmd%dxu(1,:) = bmd%dxu(2,:)
       bmd%dyu(1,:) = bmd%dyu(2,:)
       bmd%msu(1,:) = 0.0
     do j=2,grd%jm
      do i=1,grd%im
!       bmd%hgv(i,j) = (bmd%hgt(i,j)+bmd%hgt(i,j-1))*0.5
       bmd%hgv(i,j) = min(bmd%hgt(i,j),bmd%hgt(i,j-1))
       bmd%dxv(i,j) = (grd%dx(i,j)+grd%dx(i,j-1))*0.5
       bmd%dyv(i,j) = (grd%dy(i,j)+grd%dy(i,j-1))*0.5
       bmd%msv(i,j) =  bmd%mst(i,j)*bmd%mst(i,j-1)
      enddo
     enddo
       bmd%dxv(:,1) = bmd%dxv(:,2)
       bmd%dyv(:,1) = bmd%dyv(:,2)
       bmd%msv(:,1) = 0.0

     do j= 2,grd%jm-1
      do i= 2,grd%im-1
        bmd%a1(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgu(i+1,j)/grd%dx(i,j)**2*bmd%msu(i+1,j)
        bmd%a2(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgu(i  ,j)/grd%dx(i,j)**2*bmd%msu(i  ,j)
        bmd%a3(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgv(i,j+1)/grd%dy(i,j)**2*bmd%msv(i,j+1)
        bmd%a4(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgv(i,j  )/grd%dy(i,j)**2*bmd%msv(i,j  )
        bmd%a0(i,j) = (bmd%a1(i,j)+bmd%a2(i,j)+bmd%a3(i,j)+bmd%a4(i,j) +1.0) 
       bmd%a00(i,j) = (bmd%a1(i,j)+bmd%a2(i,j)+bmd%a3(i,j)+bmd%a4(i,j)) *bmd%mst(i,j)
      enddo
     enddo


end subroutine ini_bmd
