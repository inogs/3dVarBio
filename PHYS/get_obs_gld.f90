subroutine get_obs_gld

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
! Load GLIDER observations                                             !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, kk, i, cnti
  REAL(r8)      ::  houa
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

    gld%no = 0
    gld%nc = 0

   
  open(511,file='gld_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations 

   read(511) gld%no

   write(drv%dia,*) 'Number of glider observations: ',gld%no

   if(gld%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( gld%ino(gld%no), gld%flg(gld%no), gld%flc(gld%no), gld%par(gld%no))
   ALLOCATE ( gld%lon(gld%no), gld%lat(gld%no), gld%dpt(gld%no), gld%tim(gld%no))
   ALLOCATE ( gld%val(gld%no), gld%bac(gld%no), gld%inc(gld%no))
   ALLOCATE ( gld%bia(gld%no), gld%err(gld%no))
   ALLOCATE ( gld%res(gld%no), gld%b_a(gld%no))
   ALLOCATE ( gld%ib(gld%no), gld%jb(gld%no), gld%kb(gld%no))
   ALLOCATE ( gld%pb(gld%no), gld%qb(gld%no), gld%rb(gld%no))
   ALLOCATE ( gld%pq1(gld%no), gld%pq2(gld%no), gld%pq3(gld%no), gld%pq4(gld%no))
   ALLOCATE ( gld%pq5(gld%no), gld%pq6(gld%no), gld%pq7(gld%no), gld%pq8(gld%no))

       read (511)                                               &
        gld%ino(1:gld%no), gld%flg(1:gld%no), gld%par(1:gld%no) &
       ,gld%lon(1:gld%no), gld%lat(1:gld%no)                    &
       ,gld%dpt(1:gld%no), gld%tim(1:gld%no)                    &
       ,gld%val(1:gld%no), gld%bac(1:gld%no)                    &
       ,gld%err(1:gld%no), gld%res(1:gld%no)                    &
       ,gld%ib(1:gld%no), gld%jb(1:gld%no), gld%kb(1:gld%no)    &
       ,gld%pb(1:gld%no), gld%qb(1:gld%no), gld%rb(1:gld%no)
    close(511)

! ---
! Initialise quality flag
   if(obs%gld.eq.0) gld%flg(:) = -1

! ---
! Vertical interpolation parameters
    do k = 1,gld%no
     if(gld%flg(k).eq.1)then
       gld%kb(k) = grd%km-1
     do kk = 1,grd%km-1
      if( gld%dpt(k).ge.grd%dep(kk) .and. gld%dpt(k).lt.grd%dep(kk+1) ) then
       gld%kb(k) = kk
       gld%rb(k) = (gld%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
      endif
     enddo
     endif
    enddo

! residual check
  do k=1,gld%no
   if(gld%par(k).eq.1 .and. abs(gld%res(k)).gt.5.0) gld%flg(k) = 0
   if(gld%par(k).eq.2 .and. abs(gld%res(k)).gt.2.0) gld%flg(k) = 0
  enddo

! ---
! Thin observations
        houa = 12.0
    do k = 1,gld%no-1
     if(gld%flg(k).eq.1)then
         kk = k + 1
         cnti = 1.
       do kk=k+1,gld%no
        if( gld%par(k).eq.gld%par(kk) .and.                         &
            abs(gld%tim(k)-gld%tim(kk)).le.houa/24. .and.           &
            gld%kb(k).eq.gld%kb(kk) .and. gld%flg(kk).eq.1) then
         gld%lon(k) = gld%lon(k) + gld%lon(kk)
         gld%lat(k) = gld%lat(k) + gld%lat(kk)
         gld%tim(k) = gld%tim(k) + gld%tim(kk)
         gld%val(k) = gld%val(k) + gld%val(kk)
         gld%bac(k) = gld%bac(k) + gld%bac(kk)
         gld%res(k) = gld%res(k) + gld%res(kk)
         gld%dpt(k) = gld%dpt(k) + gld%dpt(kk)
         gld%flg(kk) = 0
         cnti = cnti + 1.
        endif
       enddo
         gld%lon(k) = gld%lon(k)/cnti
         gld%lat(k) = gld%lat(k)/cnti
         gld%tim(k) = gld%tim(k)/cnti
         gld%val(k) = gld%val(k)/cnti
         gld%bac(k) = gld%bac(k)/cnti
         gld%res(k) = gld%res(k)/cnti
         gld%dpt(k) = gld%dpt(k)/cnti
         gld%rb(k) = (gld%dpt(k) - grd%dep(gld%kb(k))) / (grd%dep(gld%kb(k)+1) - grd%dep(gld%kb(k)))
     endif
    enddo


! ---
! Count good observations
    gld%nc = 0
  do k=1,gld%no
   if(gld%flg(k).eq.1)then
    gld%nc = gld%nc + 1
   else
    gld%bia(k) = 0.
    gld%res(k) = 0.
    gld%inc(k) = 0.
    gld%b_a(k) = 0.
    gld%pq1(k) = 0.
    gld%pq2(k) = 0.
    gld%pq3(k) = 0.
    gld%pq4(k) = 0.
    gld%pq5(k) = 0.
    gld%pq6(k) = 0.
    gld%pq7(k) = 0.
    gld%pq8(k) = 0.
   endif
  enddo

  gld%flc(:) = gld%flg(:)


1111 continue


end subroutine get_obs_gld



subroutine int_par_gld

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

  INTEGER(i4)   ::  i, k
  INTEGER(i4)   ::  i1, j1, k1, idep
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

 if(gld%no.gt.0) then

   gld%flc(:) = gld%flg(:)

! ---
! Horizontal interpolation parameters
    do k = 1,gld%no
     q1 = (gld%lat(k) - grd%lat(1,1)) / grd%dlt + 1.0
     j1 = int(q1)
     p1 = (gld%lon(k) - grd%lon(1,1)) / grd%dln + 1.0
     i1 = int(p1)
     if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
       gld%ib(k) = i1
       gld%jb(k) = j1
       gld%pb(k) = (p1-i1)
       gld%qb(k) = (q1-j1)
     else
       gld%flc(k) = 0
     endif
    enddo

! ---
! Undefine masked for multigrid
    do k = 1,gld%no
     if(gld%flc(k).eq.1)then
      i1 = gld%ib(k)
      j1 = gld%jb(k)
      idep = gld%kb(k)+1
      msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
      if(msk4.lt.1.) gld%flc(k) = 0
     endif
    enddo

! ---
! Horizontal interpolation parameters for each masked grid
       do k = 1,gld%no
        if(gld%flc(k) .eq. 1) then

         i1=gld%ib(k)
         p1=gld%pb(k)
         j1=gld%jb(k)
         q1=gld%qb(k)


         k1=gld%kb(k)
         div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
          gld%pq1(k) = grd%msk(i1,j1,k1)                                       &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          gld%pq2(k) = grd%msk(i1+1,j1,k1)                                     &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
          gld%pq3(k) = grd%msk(i1,j1+1,k1)                                     &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          gld%pq4(k) = grd%msk(i1+1,j1+1,k1)                                   &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

         k1=gld%kb(k) + 1
         div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                 +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
         div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
          gld%pq5(k) = grd%msk(i1,j1,k1)                                       &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                       * (1.-p1) * (1.-q1)                                     &
                      /( div_x * div_y + 1.e-16 )
          gld%pq6(k) = grd%msk(i1+1,j1,k1)                                     &
                      * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                      *     p1  * (1.-q1)                                      &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
          gld%pq7(k) = grd%msk(i1,j1+1,k1)                                     &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      * (1.-p1) *     q1                                       &
                      /( div_x * div_y + 1.e-16 )
          gld%pq8(k) = grd%msk(i1+1,j1+1,k1)                                   &
                      * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                      *     p1  *     q1                                       &
                      /( div_x * div_y + 1.e-16 )

         r1=gld%rb(k)
          gld%pq1(k) = (1.-r1) * gld%pq1(k)
          gld%pq2(k) = (1.-r1) * gld%pq2(k)
          gld%pq3(k) = (1.-r1) * gld%pq3(k)
          gld%pq4(k) = (1.-r1) * gld%pq4(k)
          gld%pq5(k) =     r1  * gld%pq5(k)
          gld%pq6(k) =     r1  * gld%pq6(k)
          gld%pq7(k) =     r1  * gld%pq7(k)
          gld%pq8(k) =     r1  * gld%pq8(k)

        endif
       enddo


! ---
! Count good observations
    gld%nc = 0
  do k=1,gld%no
   if(gld%flc(k).eq.1)then
    gld%nc = gld%nc + 1
   endif
  enddo

 endif

end subroutine int_par_gld
