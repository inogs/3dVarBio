subroutine get_obs_gvl

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
! Load observations of glider velocities                               !
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
  REAL(r4)      ::  zbo, zbn
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

   gvl%no = 0
   gvl%nc = 0


  open(511,file='gvl_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations

   read(511) gvl%no

   write(drv%dia,*) 'Number of velocity observations by gliders: ',gvl%no

   if(gvl%no.eq.0)then
      close(511)
      return
   endif


   ALLOCATE ( gvl%ino(gvl%no), gvl%flg(gvl%no), gvl%flc(gvl%no), gvl%par(gvl%no))
   ALLOCATE ( gvl%lon(gvl%no), gvl%lat(gvl%no), gvl%dpt(gvl%no), gvl%kdp(gvl%no))
   ALLOCATE ( gvl%tim(gvl%no), gvl%tms(gvl%no), gvl%tme(gvl%no))
   ALLOCATE ( gvl%val(gvl%no), gvl%bac(gvl%no), gvl%inc(gvl%no))
   ALLOCATE ( gvl%bia(gvl%no), gvl%err(gvl%no))
   ALLOCATE ( gvl%res(gvl%no), gvl%b_a(gvl%no))
   ALLOCATE ( gvl%ib(gvl%no), gvl%jb(gvl%no), gvl%kb(gvl%no))
   ALLOCATE ( gvl%pb(gvl%no), gvl%qb(gvl%no), gvl%rb(gvl%no))
   ALLOCATE ( gvl%pq1(gvl%no), gvl%pq2(gvl%no), gvl%pq3(gvl%no), gvl%pq4(gvl%no))
   ALLOCATE ( gvl%pq5(gvl%no), gvl%pq6(gvl%no), gvl%pq7(gvl%no), gvl%pq8(gvl%no))
   ALLOCATE ( gvl%nav(gvl%no))
   ALLOCATE ( gvl%dzr(grd%km,gvl%no))

   gvl%bia(:) = 0.0
   gvl%b_a(:) = 0.0

       read(511)                                               &
       gvl%ino(1:gvl%no), gvl%flg(1:gvl%no), gvl%par(1:gvl%no) &
      ,gvl%lon(1:gvl%no), gvl%lat(1:gvl%no), gvl%dpt(1:gvl%no) &
      ,gvl%tim(1:gvl%no), gvl%tms(1:gvl%no), gvl%tme(1:gvl%no) &
      ,gvl%val(1:gvl%no), gvl%bac(1:gvl%no), gvl%res(1:gvl%no) &
      ,gvl%err(1:gvl%no), gvl%nav(1:gvl%no)                    &
      ,gvl%ib(1:gvl%no), gvl%jb(1:gvl%no)                      &
      ,gvl%pb(1:gvl%no), gvl%qb(1:gvl%no)

    close(511)

! ---
! Initialise quality flag
   if(obs%gvl.eq.0) gvl%flg(:) = -1

! ---
! Vertical interpolation parameters
    do k = 1,gvl%no

     if(gvl%flg(k).eq.1)then
       gvl%kb(k) = grd%km-1
      zbo = 0.0
     do kk = 1,grd%km-1
      zbn = grd%dep(kk)+grd%dz(kk)*0.5
      if( gvl%dpt(k) .gt. zbn) then
       gvl%dzr(kk,k) = grd%dz(kk)/gvl%dpt(k)
      else if( gvl%dpt(k) .ge. zbo) then
       gvl%kb(k) = kk
       gvl%dzr(kk,k) = (gvl%dpt(k)-zbo)/gvl%dpt(k)
      else 
       gvl%dzr(kk,k) = 0.0
      endif
      zbo=zbn
     enddo
     endif

    enddo

! residual check
  do k=1,gvl%no
   if(abs(gvl%res(k)).gt.1.0) gvl%flg(k) = 0
  enddo



! ---
! Count good observations
    gvl%nc = 0
  do k=1,gvl%no
   if(gvl%flg(k).eq.1)then
    gvl%nc = gvl%nc + 1
   else
    gvl%bia(k) = 0.
    gvl%res(k) = 0.
    gvl%inc(k) = 0.
    gvl%b_a(k) = 0.
    gvl%pq1(k) = 0.
    gvl%pq2(k) = 0.
    gvl%pq3(k) = 0.
    gvl%pq4(k) = 0.
    gvl%pq5(k) = 0.
    gvl%pq6(k) = 0.
    gvl%pq7(k) = 0.
    gvl%pq8(k) = 0.
   endif
  enddo

  gvl%flc(:) = gvl%flg(:)


1111 continue


end subroutine get_obs_gvl

subroutine int_par_gvl

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
  REAL(r8)      ::  msk4u, msk4v, r1, div_x, div_y
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1


 if(gvl%no.gt.0) then


   gvl%flc(:) = gvl%flg(:)

! ---
! Horizontal interpolation parameters
    do k = 1,gvl%no
     if(gvl%par(k).eq.2)then
      q1 = (gvl%lat(k) - grd%lat(1,1)-(grd%lat(1,2)-grd%lat(1,1))*0.5) /    &
           (grd%lat(1,2)-grd%lat(1,1)) + 1.0
     else
      q1 = (gvl%lat(k) - grd%lat(1,1)) / (grd%lat(1,2)-grd%lat(1,1)) + 1.0
     endif
     j1 = int(q1)
     if(gvl%par(k).eq.1)then
      p1 = (gvl%lon(k) - grd%lon(1,1)-(grd%lon(2,1)-grd%lon(1,1))*0.5) /    &
           (grd%lon(2,1)-grd%lon(1,1)) + 1.0
     else
      p1 = (gvl%lon(k) - grd%lon(1,1)) / (grd%lon(2,1)-grd%lon(1,1)) + 1.0
     endif
     j1 = int(q1)
     i1 = int(p1)
     if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
       gvl%ib(k) = i1
       gvl%jb(k) = j1
       gvl%pb(k) = (p1-i1)
       gvl%qb(k) = (q1-j1)
     else
       gvl%flc(k) = 0
     endif
    enddo

! ---
! Undefine masked
    do k = 1,gvl%no
     if(gvl%flc(k).eq.1)then
      i1 = gvl%ib(k)
      j1 = gvl%jb(k)
      idep = gvl%kb(k)+1
        msk4u = grd%ums(i1,j1  ,idep) + grd%ums(i1+1,j1  ,idep) +      &
                grd%ums(i1,j1+1,idep) + grd%ums(i1+1,j1+1,idep)
        msk4v = grd%vms(i1,j1  ,idep) + grd%vms(i1+1,j1  ,idep) +      &
                grd%vms(i1,j1+1,idep) + grd%vms(i1+1,j1+1,idep)
      if(msk4u.lt.4 .or. msk4v.lt.4)gvl%flc(k) = 0
     endif
    enddo

! Horizontal interpolation parameters for each masked grid
       do k = 1,gvl%no
        if(gvl%flc(k) .eq. 1 .and. gvl%par(k).eq.1) then

         i1=gvl%ib(k)
         p1=gvl%pb(k)
         j1=gvl%jb(k)
         q1=gvl%qb(k)

         k1=gvl%kb(k)+1
         div_y =  (1.-q1) * max(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
          gvl%pq1(k) = grd%ums(i1,j1,k1)                                    &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                               &
                      /( div_x * div_y + 1.e-16 )
          gvl%pq2(k) = grd%ums(i1+1,j1,k1)                                  &
                      * max(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
          gvl%pq3(k) = grd%ums(i1,j1+1,k1)                                  &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                 &
                      /( div_x * div_y + 1.e-16 )
          gvl%pq4(k) = grd%ums(i1+1,j1+1,k1)                                &
                      * max(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                 &
                      /( div_x * div_y + 1.e-16 )

        else if(gvl%flc(k) .eq. 1 .and. gvl%par(k).eq.2) then

         i1=gvl%ib(k)
         p1=gvl%pb(k)
         j1=gvl%jb(k)
         q1=gvl%qb(k)


         k1=gvl%kb(k)+1
         div_y =  (1.-q1) * max(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
          gvl%pq1(k) = grd%vms(i1,j1,k1)                                       &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          gvl%pq2(k) = grd%vms(i1+1,j1,k1)                                     &
                      * max(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
          gvl%pq3(k) = grd%vms(i1,j1+1,k1)                                     &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          gvl%pq4(k) = grd%vms(i1+1,j1+1,k1)                                   &
                      * max(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

        endif

       enddo


! ---
! Count good observations
    gvl%nc = 0
  do k=1,gvl%no
   if(gvl%flc(k).eq.1)then
    gvl%nc = gvl%nc + 1
   endif
  enddo
 
 endif




end subroutine int_par_gvl
