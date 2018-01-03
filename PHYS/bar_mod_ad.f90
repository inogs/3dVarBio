subroutine bar_mod_ad

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
! Barotropic model  (adjoint)                                          !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  INTEGER(i4)    :: i, j, kstp


                 bmd%bxby(:,:) = 0.
                 bmd%cu(:,:) = 0.
                 bmd%cv(:,:) = 0.
                 bmd%div(:,:) = 0.
                 bmd%dux(:,:) = 0.
                 bmd%duy(:,:) = 0.
                 bmd%dvx(:,:) = 0.
                 bmd%dvy(:,:) = 0.
                 bmd%eta(:,:) = 0.
                 bmd%etb(:,:) = 0.
                 bmd%etx(:,:) = 0.
                 bmd%ety(:,:) = 0.
                 bmd%rgh(:,:) = 0.
                 bmd%ua(:,:) = 0.
                 bmd%ub(:,:) = 0.
                 bmd%un(:,:) = 0.
                 bmd%va(:,:) = 0.
                 bmd%vb(:,:) = 0.
                 bmd%vn(:,:) = 0.
                 grd%bx(:,:) = 0.
                 grd%by(:,:) = 0.

!---

     do kstp = bmd%nstps, 1, -1



        if (mod(kstp,bmd%nstpa) .eq. 0) then

          bmd%etm(:,:) = 0.0  
           bmd%vm(:,:) = 0.0
           bmd%um(:,:) = 0.0

          if(kstp.eq.bmd%nstpa) bmd%etm(:,:) = grd%eta_ad(:,:)

           bmd%vm(:,:) =  bmd%vm(:,:)/bmd%nstpa
           bmd%um(:,:) =  bmd%um(:,:)/bmd%nstpa
          bmd%etm(:,:) = bmd%etm(:,:)/bmd%nstpa

        endif

          bmd%etb(:,:) = bmd%etb(:,:) + bmd%etm(:,:)
           bmd%ub(:,:) =  bmd%ub(:,:) +  bmd%um(:,:)
           bmd%vb(:,:) =  bmd%vb(:,:) +  bmd%vm(:,:)

           bmd%ua(:,:) =  bmd%ua(:,:) +  bmd%un(:,:)
           bmd%va(:,:) =  bmd%va(:,:) +  bmd%vn(:,:)
           bmd%un(:,:) = 0.0
           bmd%vn(:,:) = 0.0

          bmd%eta(:,:) = bmd%eta(:,:) + bmd%etb(:,:)
           bmd%un(:,:) =  bmd%un(:,:) +  bmd%ub(:,:) 
           bmd%vn(:,:) =  bmd%vn(:,:) +  bmd%vb(:,:) 
          bmd%etb(:,:) = 0.0
           bmd%ub(:,:) = 0.0
           bmd%vb(:,:) = 0.0

           bmd%ua(:,:) =  bmd%ua(:,:) + bmd%un(:,:)*0.05
           bmd%ub(:,:) =  bmd%ub(:,:) + bmd%un(:,:)*0.05
           bmd%un(:,:) =  bmd%un(:,:) * (1.0 - 2.0*0.05)
           bmd%va(:,:) =  bmd%va(:,:) + bmd%vn(:,:)*0.05
           bmd%vb(:,:) =  bmd%vb(:,:) + bmd%vn(:,:)*0.05
           bmd%vn(:,:) =  bmd%vn(:,:) * (1.0 - 2.0*0.05)

       do j= 2,grd%jm-1
        do i= 2,grd%im-1
          bmd%dux(i+1,j) = bmd%dux(i+1,j) + bmd%df2*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%dux(i  ,j) = bmd%dux(i  ,j) - bmd%df2*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%duy(i,j+1) = bmd%duy(i,j+1) + bmd%df2*bmd%ua(i,j)/bmd%dyu(i,j)*bmd%msu(i,j)
          bmd%duy(i,j  ) = bmd%duy(i,j  ) - bmd%df2*bmd%ua(i,j)/bmd%dyu(i,j)*bmd%msu(i,j)
        enddo
       enddo
       do j= 2,grd%jm-1
        do i= 2,grd%im-1
          bmd%dvx(i+1,j) = bmd%dvx(i+1,j) + bmd%df2*bmd%va(i,j)/bmd%dxv(i,j)*bmd%msv(i,j)
          bmd%dvx(i  ,j) = bmd%dvx(i  ,j) - bmd%df2*bmd%va(i,j)/bmd%dxv(i,j)*bmd%msv(i,j)
          bmd%dvy(i,j+1) = bmd%dvy(i,j+1) + bmd%df2*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
          bmd%dvy(i,j  ) = bmd%dvy(i,j  ) - bmd%df2*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
        enddo
       enddo

        do j=1,grd%jm
         do i=2,grd%im
          bmd%ua(i  ,j) = bmd%ua(i  ,j) + bmd%dux(i,j)/bmd%dxu(i,j)
          bmd%ua(i-1,j) = bmd%ua(i-1,j) - bmd%dux(i,j)/bmd%dxu(i,j)
          bmd%va(i  ,j) = bmd%va(i  ,j) + bmd%dvx(i,j)/bmd%dxv(i,j)
          bmd%va(i-1,j) = bmd%va(i-1,j) - bmd%dvx(i,j)/bmd%dxv(i,j)
          bmd%dux(i,j) = 0.0
          bmd%dvx(i,j) = 0.0
         enddo
        enddo
        do j=2,grd%jm
         do i=1,grd%im
          bmd%ua(i,j  ) = bmd%ua(i,j  ) + bmd%duy(i,j)/bmd%dyu(i,j)
          bmd%ua(i,j-1) = bmd%ua(i,j-1) - bmd%duy(i,j)/bmd%dyu(i,j)
          bmd%va(i,j  ) = bmd%va(i,j  ) + bmd%dvy(i,j)/bmd%dyv(i,j)
          bmd%va(i,j-1) = bmd%va(i,j-1) - bmd%dvy(i,j)/bmd%dyv(i,j)
          bmd%duy(i,j) = 0.0
          bmd%dvy(i,j) = 0.0
         enddo
        enddo


       do j= 2,grd%jm-1
        do i= 2,grd%im
          grd%bx(i,j)    = grd%bx(i,j) - bmd%dt*bmd%ua(i,j)*bmd%msu(i,j)
          bmd%cu(i,j)    = bmd%cu(i,j) - bmd%dt*bmd%ua(i,j)*bmd%msu(i,j)
          bmd%dux(i+1,j) = bmd%dux(i+1,j) + bmd%df1*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%dux(i  ,j) = bmd%dux(i  ,j) - bmd%df1*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%duy(i,j+1) = bmd%duy(i,j+1) + bmd%df1*bmd%ua(i,j)/bmd%dyu(i,j)*bmd%msu(i,j)
          bmd%duy(i,j  ) = bmd%duy(i,j  ) - bmd%df1*bmd%ua(i,j)/bmd%dyu(i,j)*bmd%msu(i,j)
          bmd%eta(i  ,j) = bmd%eta(i  ,j) - bmd%dt*bmd%alp1*bmd%g*bmd%hgu(i,j)*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%eta(i-1,j) = bmd%eta(i-1,j) + bmd%dt*bmd%alp1*bmd%g*bmd%hgu(i,j)*bmd%ua(i,j)/bmd%dxu(i,j)*bmd%msu(i,j)
          bmd%etx(i,j)   = bmd%etx(i,j) - bmd%ua(i,j)
          bmd%ub(i,j)    = bmd%ub(i,j) + bmd%ua(i,j)
          bmd%ua(i,j)    = 0.0
        enddo
       enddo

       do j= 2,grd%jm
        do i= 2,grd%im-1
          grd%by(i,j)    = grd%by(i,j) - bmd%dt*bmd%va(i,j)*bmd%msv(i,j)
          bmd%cv(i,j)    = bmd%cv(i,j) - bmd%dt*bmd%va(i,j)*bmd%msv(i,j)
          bmd%dvx(i+1,j) = bmd%dvx(i+1,j) + bmd%df1*bmd%va(i,j)/bmd%dxv(i,j)*bmd%msv(i,j)
          bmd%dvx(i  ,j) = bmd%dvx(i  ,j) - bmd%df1*bmd%va(i,j)/bmd%dxv(i,j)*bmd%msv(i,j)
          bmd%dvy(i,j+1) = bmd%dvy(i,j+1) + bmd%df1*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
          bmd%dvy(i,j  ) = bmd%dvy(i,j  ) - bmd%df1*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
          bmd%eta(i,j  ) = bmd%eta(i,j  ) - bmd%dt*bmd%alp1*bmd%g*bmd%hgv(i,j)*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
          bmd%eta(i,j-1) = bmd%eta(i,j-1) + bmd%dt*bmd%alp1*bmd%g*bmd%hgv(i,j)*bmd%va(i,j)/bmd%dyv(i,j)*bmd%msv(i,j)
          bmd%ety(i,j)   = bmd%ety(i,j) - bmd%va(i,j)
          bmd%vb(i,j)    = bmd%vb(i,j) + bmd%va(i,j)
          bmd%va(i,j)    = 0.0
        enddo
       enddo

          bmd%ub(2:grd%im  ,:) = bmd%ub(2:grd%im  ,:) + bmd%dux(2:grd%im,:)/bmd%dxu(2:grd%im,:)
          bmd%ub(1:grd%im-1,:) = bmd%ub(1:grd%im-1,:) - bmd%dux(2:grd%im,:)/bmd%dxu(2:grd%im,:)
          bmd%dux(2:grd%im,:)  = 0.0
          bmd%vb(2:grd%im  ,:) = bmd%vb(2:grd%im  ,:) + bmd%dvx(2:grd%im,:)/bmd%dxv(2:grd%im,:)
          bmd%vb(1:grd%im-1,:) = bmd%vb(1:grd%im-1,:) - bmd%dvx(2:grd%im,:)/bmd%dxv(2:grd%im,:)
          bmd%dvx(2:grd%im,:)  = 0.0

          bmd%ub(:,2:grd%jm  ) = bmd%ub(:,2:grd%jm  ) + bmd%duy(:,2:grd%jm)/bmd%dyu(:,2:grd%jm)
          bmd%ub(:,1:grd%jm-1) = bmd%ub(:,1:grd%jm-1) - bmd%duy(:,2:grd%jm)/bmd%dyu(:,2:grd%jm)
          bmd%duy(2:grd%im,:)  = 0.0
          bmd%vb(:,2:grd%jm  ) = bmd%vb(:,2:grd%jm  ) + bmd%dvy(:,2:grd%jm)/bmd%dyv(:,2:grd%jm)
          bmd%vb(:,1:grd%jm-1) = bmd%vb(:,1:grd%jm-1) - bmd%dvy(:,2:grd%jm)/bmd%dyv(:,2:grd%jm)
          bmd%dvy(2:grd%im,:)  = 0.0

       do j=1,grd%jm-1
        do i=2,grd%im-1
         bmd%vn(i-1,j+1) = bmd%vn(i-1,j+1) - bmd%cu(i,j)*bmd%dxv(i-1,j+1)*grd%f(i,j+1)/bmd%dxu(i,j)*0.25
         bmd%vn(i-1,j  ) = bmd%vn(i-1,j  ) - bmd%cu(i,j)*bmd%dxv(i-1,j  )*grd%f(i,j  )/bmd%dxu(i,j)*0.25
         bmd%vn(i  ,j+1) = bmd%vn(i  ,j+1) - bmd%cu(i,j)*bmd%dxv(i  ,j+1)*grd%f(i,j+1)/bmd%dxu(i,j)*0.25
         bmd%vn(i  ,j  ) = bmd%vn(i  ,j  ) - bmd%cu(i,j)*bmd%dxv(i  ,j  )*grd%f(i,j  )/bmd%dxu(i,j)*0.25
         bmd%cu(i,j)     = 0.0
        enddo
       enddo
       do j=2,grd%jm-1
        do i=1,grd%im-1
         bmd%un(i+1,j-1) = bmd%un(i+1,j-1) + bmd%cv(i,j)*bmd%dyu(i+1,j-1)*grd%f(i+1,j)/bmd%dyv(i,j)*0.25
         bmd%un(i  ,j-1) = bmd%un(i  ,j-1) + bmd%cv(i,j)*bmd%dyu(i  ,j-1)*grd%f(i  ,j)/bmd%dyv(i,j)*0.25
         bmd%un(i+1,j  ) = bmd%un(i+1,j  ) + bmd%cv(i,j)*bmd%dyu(i+1,j  )*grd%f(i+1,j)/bmd%dyv(i,j)*0.25
         bmd%un(i  ,j  ) = bmd%un(i  ,j  ) + bmd%cv(i,j)*bmd%dyu(i  ,j  )*grd%f(i  ,j)/bmd%dyv(i,j)*0.25
         bmd%cv(i,j)     = 0.0
        enddo
       enddo

       call invrt_ad( grd%im, grd%jm, bmd%eta, bmd%mst, bmd%rgh, bmd%a1, bmd%a2, bmd%a3, bmd%a4, bmd%a0,           &
                      bmd%bnm, bmd%ovr, bmd%resem, bmd%ncnt, bmd%itr(kstp) )

       do j= 2,grd%jm-1
        do i= 2,grd%im-1
         bmd%bxby(i,j)  = bmd%bxby(i,j) + bmd%rgh(i,j)
         bmd%div(i,j)   = bmd%div(i,j) + bmd%dt*bmd%rgh(i,j)
         bmd%etb(i,j)   = bmd%etb(i,j) - bmd%rgh(i,j)
         bmd%etx(i+1,j) = bmd%etx(i+1,j) - bmd%alp1*bmd%dt*bmd%rgh(i,j)/grd%dx(i,j)
         bmd%etx(i  ,j) = bmd%etx(i  ,j) + bmd%alp1*bmd%dt*bmd%rgh(i,j)/grd%dx(i,j)
         bmd%ety(i,j+1) = bmd%ety(i,j+1) - bmd%alp1*bmd%dt*bmd%rgh(i,j)/grd%dy(i,j)
         bmd%ety(i,j  ) = bmd%ety(i,j  ) + bmd%alp1*bmd%dt*bmd%rgh(i,j)/grd%dy(i,j)
         bmd%rgh(i,j  ) = 0.0
        enddo
       enddo

       do j= 2,grd%jm
        do i= 2,grd%im-1
          bmd%etb(i,j  ) = bmd%etb(i,j  ) + bmd%alp2*bmd%dt*bmd%g*bmd%hgv(i,j)*bmd%ety(i,j)/bmd%dyv(i,j) * bmd%msv(i,j)
          bmd%etb(i,j-1) = bmd%etb(i,j-1) - bmd%alp2*bmd%dt*bmd%g*bmd%hgv(i,j)*bmd%ety(i,j)/bmd%dyv(i,j) * bmd%msv(i,j)
          bmd%ety(i,j)   = 0.0
        enddo
       enddo
       do j= 2,grd%jm-1
        do i= 2,grd%im
          bmd%etb(i  ,j) = bmd%etb(i  ,j) + bmd%alp2*bmd%dt*bmd%g*bmd%hgu(i,j)*bmd%etx(i,j)/bmd%dxu(i,j) * bmd%msu(i,j)
          bmd%etb(i-1,j) = bmd%etb(i-1,j) - bmd%alp2*bmd%dt*bmd%g*bmd%hgu(i,j)*bmd%etx(i,j)/bmd%dxu(i,j) * bmd%msu(i,j)
          bmd%etx(i,j)   = 0.0
        enddo
       enddo

       do j= 2,grd%jm-1
        do i= 2,grd%im-1
          bmd%ub(i+1,j) = bmd%ub(i+1,j) + bmd%div(i,j)/grd%dx(i,j) * bmd%mst(i,j)
          bmd%ub(i  ,j) = bmd%ub(i  ,j) - bmd%div(i,j)/grd%dx(i,j) * bmd%mst(i,j)
          bmd%vb(i,j+1) = bmd%vb(i,j+1) + bmd%div(i,j)/grd%dy(i,j) * bmd%mst(i,j)
          bmd%vb(i,j  ) = bmd%vb(i,j  ) - bmd%div(i,j)/grd%dy(i,j) * bmd%mst(i,j)
          bmd%div(i,j)  = 0.0
        enddo
       enddo

     enddo   ! kstp


      do j= 2,grd%jm-1
       do i= 2,grd%im-1
        bmd%bxby(i,j) = bmd%bxby(i,j) + bmd%rgh(i,j) * bmd%alp1**2
        bmd%rgh(i,j)  = 0.0
        grd%bx(i+1,j) = grd%bx(i+1,j) - (bmd%dt**2)*bmd%bxby(i,j)/grd%dx(i,j) * bmd%mst(i,j)
        grd%bx(i  ,j) = grd%bx(i  ,j) + (bmd%dt**2)*bmd%bxby(i,j)/grd%dx(i,j) * bmd%mst(i,j)
        grd%by(i,j+1) = grd%by(i,j+1) - (bmd%dt**2)*bmd%bxby(i,j)/grd%dy(i,j) * bmd%mst(i,j)
        grd%by(i,j  ) = grd%by(i,j  ) + (bmd%dt**2)*bmd%bxby(i,j)/grd%dy(i,j) * bmd%mst(i,j)
        bmd%bxby(i,j) = 0.0
       enddo
      enddo

    
end subroutine bar_mod_ad
