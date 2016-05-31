subroutine get_obs_vdr

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
! Load observations of drifter velocities                              !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, kk, i
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

   vdr%no = 0
   vdr%nc = 0


  open(511,file='vdr_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations

   read(511) vdr%no

   write(drv%dia,*) 'Number of velocity observations by gliders: ',vdr%no

   if(vdr%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( vdr%ino(vdr%no), vdr%flg(vdr%no), vdr%flc(vdr%no), vdr%par(vdr%no))
   ALLOCATE ( vdr%lon(vdr%no), vdr%lat(vdr%no), vdr%dpt(vdr%no), vdr%kdp(vdr%no))
   ALLOCATE ( vdr%tim(vdr%no), vdr%tms(vdr%no), vdr%tme(vdr%no))
   ALLOCATE ( vdr%val(vdr%no), vdr%bac(vdr%no), vdr%inc(vdr%no))
   ALLOCATE ( vdr%bia(vdr%no), vdr%err(vdr%no))
   ALLOCATE ( vdr%res(vdr%no), vdr%b_a(vdr%no))
   ALLOCATE ( vdr%ib(vdr%no), vdr%jb(vdr%no), vdr%kb(vdr%no))
   ALLOCATE ( vdr%pb(vdr%no), vdr%qb(vdr%no), vdr%rb(vdr%no))
   ALLOCATE ( vdr%pq1(vdr%no), vdr%pq2(vdr%no), vdr%pq3(vdr%no), vdr%pq4(vdr%no))
   ALLOCATE ( vdr%pq5(vdr%no), vdr%pq6(vdr%no), vdr%pq7(vdr%no), vdr%pq8(vdr%no))
   ALLOCATE ( vdr%nav(vdr%no))

   vdr%bia(:) = 0.0

       read(511)                                               &
       vdr%ino(1:vdr%no), vdr%flg(1:vdr%no), vdr%par(1:vdr%no) &
      ,vdr%lon(1:vdr%no), vdr%lat(1:vdr%no), vdr%dpt(1:vdr%no) &
      ,vdr%tim(1:vdr%no), vdr%tms(1:vdr%no), vdr%tme(1:vdr%no) &
      ,vdr%val(1:vdr%no), vdr%bac(1:vdr%no)                    &
      ,vdr%err(1:vdr%no), vdr%nav(1:vdr%no)                    &
      ,vdr%ib(1:vdr%no), vdr%jb(1:vdr%no)                      &
      ,vdr%pb(1:vdr%no), vdr%qb(1:vdr%no)

    close(511)

! ---
! Initialise quality flag
   if(obs%vdr.eq.0) vdr%flg(:) = -1

! ---
! Vertical interpolation parameters
    do k = 1,vdr%no
     if(vdr%flg(k).eq.1)then
       vdr%kb(k) = grd%km-1
     do kk = 1,grd%km-1
      if( vdr%dpt(k).ge.grd%dep(kk) .and. vdr%dpt(k).lt.grd%dep(kk+1) ) then
       vdr%kb(k) = kk
       vdr%rb(k) = (vdr%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
      endif
     enddo
     endif
    enddo

! calculate residuals
  do k=1,vdr%no
   if(vdr%nav(k).gt.0 ) then
     vdr%bac(k) = vdr%bac(k)/vdr%nav(k)
     vdr%res(k) = vdr%val(k)-vdr%bac(k)
   endif
  enddo

! residual check
  do k=1,vdr%no
   if(abs(vdr%res(k)).gt.1.0) vdr%flg(k) = 0
  enddo



! ---
! Count good observations
    vdr%nc = 0
  do k=1,vdr%no
   if(vdr%flg(k).eq.1)then
    vdr%nc = vdr%nc + 1
   else
    vdr%bia(k) = 0.
    vdr%res(k) = 0.
    vdr%inc(k) = 0.
    vdr%b_a(k) = 0.
    vdr%pq1(k) = 0.
    vdr%pq2(k) = 0.
    vdr%pq3(k) = 0.
    vdr%pq4(k) = 0.
    vdr%pq5(k) = 0.
    vdr%pq6(k) = 0.
    vdr%pq7(k) = 0.
    vdr%pq8(k) = 0.
   endif
  enddo

  vdr%flc(:) = vdr%flg(:)


1111 continue


end subroutine get_obs_vdr

subroutine int_par_vdr

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, j1, k1, i, idep
  REAL(r8)      ::  p1, q1
  REAL(r8)      ::  msk4, r1, div_x, div_y
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1


 if(vdr%no.gt.0) then


   vdr%flc(:) = vdr%flg(:)

! ---
! Horizontal interpolation parameters
    do k = 1,vdr%no
     if(vdr%par(k).eq.2)then
      q1 = (vdr%lat(k) - grd%lat(1,1)-(grd%lat(1,2)-grd%lat(1,1))*0.5) /    &
           (grd%lat(1,2)-grd%lat(1,1)) + 1.0
     else
      q1 = (vdr%lat(k) - grd%lat(1,1)) / (grd%lat(1,2)-grd%lat(1,1)) + 1.0
     endif
     j1 = int(q1)
     if(vdr%par(k).eq.1)then
      p1 = (vdr%lon(k) - grd%lon(1,1)-(grd%lon(2,1)-grd%lon(1,1))*0.5) /    &
           (grd%lon(2,1)-grd%lon(1,1)) + 1.0
     else
      p1 = (vdr%lon(k) - grd%lon(1,1)) / (grd%lon(2,1)-grd%lon(1,1)) + 1.0
     endif
     j1 = int(q1)
     i1 = int(p1)
     if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
       vdr%ib(k) = i1
       vdr%jb(k) = j1
       vdr%pb(k) = (p1-i1)
       vdr%qb(k) = (q1-j1)
     else
       vdr%flc(k) = 0
     endif
    enddo

! ---
! Undefine masked
    do k = 1,vdr%no
     if(vdr%flc(k).eq.1)then
      i1 = vdr%ib(k)
      j1 = vdr%jb(k)
      idep = vdr%kb(k)+1
       if(vdr%par(k).eq.1) then
         msk4 = grd%ums(i1,j1  ,idep) + grd%ums(i1+1,j1  ,idep) +      &
                grd%ums(i1,j1+1,idep) + grd%ums(i1+1,j1+1,idep)
       endif
       if(vdr%par(k).eq.2) then
         msk4 = grd%vms(i1,j1  ,idep) + grd%vms(i1+1,j1  ,idep) +      &
                grd%vms(i1,j1+1,idep) + grd%vms(i1+1,j1+1,idep)
       endif
      if(msk4.lt.2)vdr%flc(k) = 0
     endif
    enddo

! Horizontal interpolation parameters for each masked grid
       do k = 1,vdr%no
        if(vdr%flc(k) .eq. 1 .and. vdr%par(k).eq.1) then

         i1=vdr%ib(k)
         p1=vdr%pb(k)
         j1=vdr%jb(k)
         q1=vdr%qb(k)
         r1=vdr%rb(k)

         k1=vdr%kb(k)
         div_y =  (1.-q1) * max(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
          vdr%pq1(k) = grd%ums(i1,j1,k1)                                    &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                               &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq2(k) = grd%ums(i1+1,j1,k1)                                  &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
          vdr%pq3(k) = grd%ums(i1,j1+1,k1)                                  &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                 &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq4(k) = grd%ums(i1+1,j1+1,k1)                                &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                 &
                      /( div_x * div_y + 1.e-16 )

         k1=vdr%kb(k) + 1
         div_y =  (1.-q1) * max(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
          vdr%pq5(k) = grd%ums(i1,j1,k1)                                    &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                               &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq6(k) = grd%ums(i1+1,j1,k1)                                  &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
          vdr%pq7(k) = grd%ums(i1,j1+1,k1)                                  &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                 &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq8(k) = grd%ums(i1+1,j1+1,k1)                                &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                 &
                      /( div_x * div_y + 1.e-16 )

          vdr%pq1(k) = (1.-r1) * vdr%pq1(k)
          vdr%pq2(k) = (1.-r1) * vdr%pq2(k)
          vdr%pq3(k) = (1.-r1) * vdr%pq3(k)
          vdr%pq4(k) = (1.-r1) * vdr%pq4(k)
          vdr%pq5(k) =     r1  * vdr%pq5(k)
          vdr%pq6(k) =     r1  * vdr%pq6(k)
          vdr%pq7(k) =     r1  * vdr%pq7(k)
          vdr%pq8(k) =     r1  * vdr%pq8(k)

        else if(vdr%flc(k) .eq. 1 .and. vdr%par(k).eq.2) then

         i1=vdr%ib(k)
         p1=vdr%pb(k)
         j1=vdr%jb(k)
         q1=vdr%qb(k)
         r1=vdr%rb(k)


         k1=vdr%kb(k)
         div_y =  (1.-q1) * max(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
          vdr%pq1(k) = grd%vms(i1,j1,k1)                                       &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq2(k) = grd%vms(i1+1,j1,k1)                                     &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
          vdr%pq3(k) = grd%vms(i1,j1+1,k1)                                     &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq4(k) = grd%vms(i1+1,j1+1,k1)                                   &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

         k1=vdr%kb(k) + 1
         div_y =  (1.-q1) * max(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
          vdr%pq5(k) = grd%vms(i1,j1,k1)                                       &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq6(k) = grd%vms(i1+1,j1,k1)                                     &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
          vdr%pq7(k) = grd%vms(i1,j1+1,k1)                                     &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          vdr%pq8(k) = grd%vms(i1+1,j1+1,k1)                                   &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                 &
                      /( div_x * div_y + 1.e-16 )

          vdr%pq1(k) = (1.-r1) * vdr%pq1(k)
          vdr%pq2(k) = (1.-r1) * vdr%pq2(k)
          vdr%pq3(k) = (1.-r1) * vdr%pq3(k)
          vdr%pq4(k) = (1.-r1) * vdr%pq4(k)
          vdr%pq5(k) =     r1  * vdr%pq5(k)
          vdr%pq6(k) =     r1  * vdr%pq6(k)
          vdr%pq7(k) =     r1  * vdr%pq7(k)
          vdr%pq8(k) =     r1  * vdr%pq8(k)

        endif

       enddo


! ---
! Count good observations
    vdr%nc = 0
  do k=1,vdr%no
   if(vdr%flc(k).eq.1)then
    vdr%nc = vdr%nc + 1
   endif
  enddo
 
 endif




end subroutine int_par_vdr
