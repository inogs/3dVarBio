subroutine get_obs_sla

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Load SLA observations                                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, i, iter
  REAL(r8)      ::  sumt, sumi, timp, dxx, dyy, dsm
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

  sla%no = 0
  sla%nc = 0

   
  open(511,file='sla_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations 

   read(511) sla%no
   write(drv%dia,*) 'Number of SLA observations: ',  sla%no, obs%sla


   if(sla%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( sla%ino(sla%no), sla%flg(sla%no), sla%flc(sla%no))
   ALLOCATE ( sla%lon(sla%no), sla%lat(sla%no), sla%tim(sla%no))
   ALLOCATE ( sla%val(sla%no), sla%bac(sla%no), sla%inc(sla%no))
   ALLOCATE ( sla%bia(sla%no), sla%err(sla%no))
   ALLOCATE ( sla%res(sla%no), sla%b_a(sla%no))
   ALLOCATE ( sla%ib(sla%no), sla%jb(sla%no))
   ALLOCATE ( sla%pb(sla%no), sla%qb(sla%no))
   ALLOCATE ( sla%pq1(sla%no), sla%pq2(sla%no), sla%pq3(sla%no), sla%pq4(sla%no))
   ALLOCATE ( sla%dpt(sla%no))

! ---
! Initialise quality flag
   sla%flc(:) = 1

! ---
! Level corresponding to the minimum depth
   sla%kdp=grd%km
  do k=grd%km, 1, -1
   if(grd%dep(k).ge.sla%dep) sla%kdp = k
  enddo

       read (511)                                              &
        sla%ino(1:sla%no), sla%flg(1:sla%no)                   &
       ,sla%lon(1:sla%no), sla%lat(1:sla%no), sla%tim(1:sla%no)&
       ,sla%val(1:sla%no), sla%bac(1:sla%no)                   &
       ,sla%err(1:sla%no), sla%res(1:sla%no)                   &
       ,sla%ib(1:sla%no), sla%jb(1:sla%no)                     &
       ,sla%pb(1:sla%no), sla%qb(1:sla%no)
    close(511)


   if(obs%sla.eq.0) sla%flg(:) = -1

! ---
! Remove bias along each track and obseravtaions with large residuals

 do iter = 1,3

!bias
   sla%bia(:) = 0.0
   timp = sla%tim(1)
   dsm = 100.
   i1 = 1
 do k=2,sla%no

      dxx = 6371.*3.14/180. * (sla%lon(k)-sla%lon(k-1)) * cos(sla%lat(k)*3.14/180.)
      dyy = 6371.*3.14/180. * (sla%lat(k)-sla%lat(k-1))

  if((sla%tim(k).ne.timp .or. sqrt(dxx**2+dyy**2).gt.dsm) .and. k.gt.i1)then
     sumt = 0.0
     sumi = 0.0
    do i=i1,k-1
     if(sla%flg(i).eq.1)then
      sumt = sumt + sla%res(i)
      sumi = sumi + 1.0
     endif
    enddo
     if(sumi.gt.0.) sumt = sumt/sumi
    do i=i1,k-1
     sla%res(i) = sla%res(i) - sumt
     sla%bia(i) = sumt
    enddo
     timp = sla%tim(k)
     i1 = k
  else if(k.eq.sla%no .and. k.ge.i1)then
     sumt = 0.0
     sumi = 0.0
    do i=i1,k
     if(sla%flg(i).eq.1)then
      sumt = sumt + sla%res(i)
      sumi = sumi + 1.0
     endif
    enddo
     if(sumi.gt.0.) sumt = sumt/sumi
    do i=i1,k
     sla%res(i) = sla%res(i) - sumt
     sla%bia(i) = sumt
    enddo
  endif
 enddo

 enddo ! iter

! residual check
  do k=1,sla%no
   if(abs(sla%res(k)).gt.0.3) sla%flg(k) = 0
  enddo


! ---
! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flg(k).eq.1)then
    sla%nc = sla%nc + 1
   else
    sla%flc(k) = 0
!    sla%bia(k) = 0.
    sla%res(k) = 0.
    sla%inc(k) = 0.
    sla%b_a(k) = 0.
    sla%pq1(k) = 0.
    sla%pq2(k) = 0.
    sla%pq3(k) = 0.
    sla%pq4(k) = 0.
   endif
  enddo

  sla%flc(:) = sla%flg(:)



1111 continue


end subroutine get_obs_sla

subroutine int_par_sla

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  integer(i4)   ::  k
  integer(i4)   ::  i1, kk, i, j1
  real(r8)      ::  p1, q1
  real(r8)      ::  msk4, div_x, div_y
  logical       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

  write(drv%dia,*) 'Number of SLA observations:  >>>>>>>>>>>>>',sla%no

 if(sla%no.gt.0) then
 
       sla%flc(:) = sla%flg(:)

       sla%dpt(:) = 0.0

! ---
! Interpolation parameters
    do kk = 1,sla%no
     q1 = (sla%lat(kk) - grd%lat(1,1)) / grd%dlt + 1.0
     j1 = int(q1)
     p1 = (sla%lon(kk) - grd%lon(1,1)) / grd%dln + 1.0
     i1 = int(p1)
     if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
       sla%dpt(kk) = max(grd%hgt(i1,j1),max(grd%hgt(i1+1,j1),max(grd%hgt(i1,j1+1),grd%hgt(i1+1,j1+1))))
       msk4 = grd%msk(i1,j1,sla%kdp) + grd%msk(i1+1,j1,sla%kdp) + grd%msk(i1,j1+1,sla%kdp) + grd%msk(i1+1,j1+1,sla%kdp)
      if(msk4.ge.4.0)then
       sla%ib(kk) = i1
       sla%jb(kk) = j1
       sla%pb(kk) = (p1-i1)
       sla%qb(kk) = (q1-j1)
      else
       sla%flc(kk) = 0
      endif
     else
       sla%flc(kk) = 0
     endif
    enddo


! ---
! Horizontal interpolation parameters for each masked grid
       do k = 1,sla%no
        if(sla%flc(k) .eq. 1) then

         i1=sla%ib(k)
         p1=sla%pb(k)
         j1=sla%jb(k)
         q1=sla%qb(k)

         div_y =  (1.-q1) * max(grd%msk(i1,j1  ,1),grd%msk(i1+1,j1  ,1))     &
                 +    q1  * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))
         div_x =  (1.-p1) * grd%msk(i1  ,j1,1) + p1 * grd%msk(i1+1,j1,1)
          sla%pq1(k) = grd%msk(i1,j1,1)                                      &
                      * max(grd%msk(i1,j1,1),grd%msk(i1+1,j1,1))             &
                       * (1.-p1) * (1.-q1)                                   &
                      /( div_x * div_y + 1.e-16 )
          sla%pq2(k) = grd%msk(i1+1,j1,1)                                    &
                      * max(grd%msk(i1,j1,1),grd%msk(i1+1,j1,1))             &
                      *     p1  * (1.-q1)                                    &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p1) * grd%msk(i1  ,j1+1,1) + p1 * grd%msk(i1+1,j1+1,1)
          sla%pq3(k) = grd%msk(i1,j1+1,1)                                    &
                      * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))         &
                      * (1.-p1) *     q1                                     &
                      /( div_x * div_y + 1.e-16 )
          sla%pq4(k) = grd%msk(i1+1,j1+1,1)                                  &
                      * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))         &
                      *     p1  *     q1                                     &
                      /( div_x * div_y + 1.e-16 )

        endif
       enddo
   
! ---
! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flc(k).eq.1)then
    sla%nc = sla%nc + 1
   endif
  enddo

  print*,'Good sla observations: ',sla%nc


 endif




end subroutine int_par_sla
