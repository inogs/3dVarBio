subroutine get_obs_xbt

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
! Load SLA observations                                                !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i8)   ::  i1, kk, i
  REAL(r8)      ::  cnti
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

   xbt%no = 0
   xbt%nc = 0

   
  open(511,file='xbt_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations 

  read(511) xbt%no

   write(drv%dia,*) 'No xbt: ',xbt%no

   if(xbt%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( xbt%ino(xbt%no), xbt%flg(xbt%no), xbt%flc(xbt%no), xbt%par(xbt%no))
   ALLOCATE ( xbt%lon(xbt%no), xbt%lat(xbt%no), xbt%dpt(xbt%no), xbt%tim(xbt%no))
   ALLOCATE ( xbt%val(xbt%no), xbt%bac(xbt%no), xbt%inc(xbt%no))
   ALLOCATE ( xbt%bia(xbt%no), xbt%err(xbt%no))
   ALLOCATE ( xbt%res(xbt%no), xbt%b_a(xbt%no))
   ALLOCATE ( xbt%ib(xbt%no), xbt%jb(xbt%no), xbt%kb(xbt%no))
   ALLOCATE ( xbt%pb(xbt%no), xbt%qb(xbt%no), xbt%rb(xbt%no))
   ALLOCATE ( xbt%pq1(xbt%no), xbt%pq2(xbt%no), xbt%pq3(xbt%no), xbt%pq4(xbt%no))
   ALLOCATE ( xbt%pq5(xbt%no), xbt%pq6(xbt%no), xbt%pq7(xbt%no), xbt%pq8(xbt%no))


   xbt%bia(:) = 0.0

       read (511)                                               &
        xbt%ino(1:xbt%no), xbt%flg(1:xbt%no), xbt%par(1:xbt%no) &
       ,xbt%lon(1:xbt%no), xbt%lat(1:xbt%no)                    &
       ,xbt%dpt(1:xbt%no), xbt%tim(1:xbt%no)                    &
       ,xbt%val(1:xbt%no), xbt%bac(1:xbt%no)                    &
       ,xbt%err(1:xbt%no), xbt%res(1:xbt%no)                    &
       ,xbt%ib(1:xbt%no), xbt%jb(1:xbt%no), xbt%kb(1:xbt%no)    &
       ,xbt%pb(1:xbt%no), xbt%qb(1:xbt%no), xbt%rb(1:xbt%no)
      close (511)

! ---
! Initialise quality flag
   if(obs%xbt.eq.0) xbt%flg(:) = -1

! ---
! Vertical interpolation parameters
    do k = 1,xbt%no
     if(xbt%flg(k).eq.1)then
       xbt%kb(k) = grd%km-1
     do kk = 1,grd%km-1
      if( xbt%dpt(k).ge.grd%dep(kk) .and. xbt%dpt(k).lt.grd%dep(kk+1) ) then
       xbt%kb(k) = kk
       xbt%rb(k) = (xbt%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
      endif
     enddo
     endif
    enddo

! ---
! Thin observations
    do k = 1,xbt%no-1
     if(xbt%flg(k).eq.1)then
         kk = k + 1
         cnti = 1.
       do while(kk.le.xbt%no .and. xbt%kb(k).eq.xbt%kb(min(kk,xbt%no)) .and. xbt%flg(min(kk,xbt%no)).eq.1)
         xbt%val(k) = xbt%val(k) + xbt%val(kk)
         xbt%bac(k) = xbt%bac(k) + xbt%bac(kk)
         xbt%res(k) = xbt%res(k) + xbt%res(kk)
         xbt%dpt(k) = xbt%dpt(k) + xbt%dpt(kk)
         xbt%flg(kk) = 0
         cnti = cnti + 1.
         kk = kk + 1
       enddo
         xbt%val(k) = xbt%val(k)/cnti
         xbt%bac(k) = xbt%bac(k)/cnti
         xbt%res(k) = xbt%res(k)/cnti
         xbt%dpt(k) = xbt%dpt(k)/cnti
         xbt%rb(k) = (xbt%dpt(k) - grd%dep(xbt%kb(k))) / (grd%dep(xbt%kb(k)+1) - grd%dep(xbt%kb(k)))
     endif
    enddo


! residual check
  do k=1,xbt%no
   if(xbt%par(k).eq.1 .and. abs(xbt%res(k)).gt.5.0) xbt%flg(k) = 0
  enddo


! ---
! Count good observations
    xbt%nc = 0
  do k=1,xbt%no
   if(xbt%flg(k).eq.1)then
    xbt%nc = xbt%nc + 1
   else
    xbt%bia(k) = 0.
    xbt%res(k) = 0.
    xbt%inc(k) = 0.
    xbt%b_a(k) = 0.
    xbt%pq1(k) = 0.
    xbt%pq2(k) = 0.
    xbt%pq3(k) = 0.
    xbt%pq4(k) = 0.
    xbt%pq5(k) = 0.
    xbt%pq6(k) = 0.
    xbt%pq7(k) = 0.
    xbt%pq8(k) = 0.
   endif
  enddo


   xbt%flc(:) = xbt%flg(:)


1111 continue


end subroutine get_obs_xbt



subroutine int_par_xbt

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k, k1
  INTEGER(i4)   ::  i1, j1, i, idep
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

 if(xbt%no.gt.0) then

    xbt%flc(:) = xbt%flg(:)

! ---
! Horizontal interpolation parameters
    do k = 1,xbt%no
     q1 = (xbt%lat(k) - grd%lat(1,1)) / grd%dlt + 1.0
     j1 = int(q1)
     p1 = (xbt%lon(k) - grd%lon(1,1)) / grd%dln + 1.0
     i1 = int(p1)
     if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
       xbt%ib(k) = i1
       xbt%jb(k) = j1
       xbt%pb(k) = (p1-i1)
       xbt%qb(k) = (q1-j1)
     else
       xbt%flc(k) = 0
     endif
    enddo

! ---
! Undefine masked for multigrid
    do k = 1,xbt%no
     if(xbt%flc(k).eq.1)then
      i1 = xbt%ib(k)
      j1 = xbt%jb(k)
      idep = xbt%kb(k)+1
       msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
      if(msk4.lt.1.0)xbt%flc(k) = 0
     endif
    enddo

! ---
! Horizontal interpolation parameters for each masked grid
       do k = 1,xbt%no
        if(xbt%flc(k) .eq. 1) then

         i1=xbt%ib(k)
         p1=xbt%pb(k)
         j1=xbt%jb(k)
         q1=xbt%qb(k)


         k1=xbt%kb(k)
         div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
          xbt%pq1(k) = grd%msk(i1,j1,k1)                                       &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          xbt%pq2(k) = grd%msk(i1+1,j1,k1)                                     &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
          xbt%pq3(k) = grd%msk(i1,j1+1,k1)                                     &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          xbt%pq4(k) = grd%msk(i1+1,j1+1,k1)                                   &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

         k1=xbt%kb(k) + 1
         div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
          xbt%pq5(k) = grd%msk(i1,j1,k1)                                       &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          xbt%pq6(k) = grd%msk(i1+1,j1,k1)                                     &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
          xbt%pq7(k) = grd%msk(i1,j1+1,k1)                                     &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          xbt%pq8(k) = grd%msk(i1+1,j1+1,k1)                                   &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

         r1=xbt%rb(k)
          xbt%pq1(k) = (1.-r1) * xbt%pq1(k)
          xbt%pq2(k) = (1.-r1) * xbt%pq2(k)
          xbt%pq3(k) = (1.-r1) * xbt%pq3(k)
          xbt%pq4(k) = (1.-r1) * xbt%pq4(k)
          xbt%pq5(k) =     r1  * xbt%pq5(k)
          xbt%pq6(k) =     r1  * xbt%pq6(k)
          xbt%pq7(k) =     r1  * xbt%pq7(k)
          xbt%pq8(k) =     r1  * xbt%pq8(k)

        endif
       enddo


! ---
! Count good observations
    xbt%nc = 0
  do k=1,xbt%no
   if(xbt%flc(k).eq.1)then
    xbt%nc = xbt%nc + 1
   endif
  enddo


 endif


end subroutine int_par_xbt
