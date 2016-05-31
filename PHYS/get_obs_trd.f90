subroutine get_obs_trd

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
!                                                                          !
!    This file is part of OceanVar.                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Load trajectory observations by surface drifters                     !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, j1, i, j
  LOGICAL       ::  ins

  ins(i,i1) = i.ge.1 .and. i.lt.i1

   trd%no = 0
   trd%nc = 0

   
    open(511,file='trd_obs.dat',form='unformatted',status='old',err=1111)

    read(511) trd%no, trd%nc, trd%nt, trd%im, trd%jm, trd%km, trd%dpt

! ---
! Allocate memory for observations 

   write(drv%dia,*) 'Number of trajectory observations by surface drifters: ',  trd%no


   if(trd%no.eq.0)then
      close(511)
      return
   endif

   allocate ( trd%ino(trd%no), trd%flg(trd%no), trd%flc(trd%no))
   allocate ( trd%loi(trd%no), trd%lai(trd%no), trd%tim(trd%no), trd%dtm(trd%no))
   allocate ( trd%lof(trd%no), trd%laf(trd%no))
   allocate ( trd%err(trd%no))
   allocate ( trd%lob(trd%nt+1,trd%no), trd%lab(trd%nt+1,trd%no) )
   allocate ( trd%loa(trd%no), trd%laa(trd%no) )
   allocate ( trd%xob(trd%no), trd%xmn(trd%nt+1,trd%no), trd%erx(trd%no) )
   allocate ( trd%yob(trd%no), trd%ymn(trd%nt+1,trd%no), trd%ery(trd%no) )
   allocate ( trd%umn(trd%im,trd%jm), trd%vmn(trd%im,trd%jm) )
   allocate ( trd%dx(trd%im,trd%jm), trd%dy(trd%im,trd%jm) )
   allocate ( trd%lon(trd%im,trd%jm), trd%lat(trd%im,trd%jm) )

   allocate ( trd%rex(trd%no), trd%inx(trd%no))
   allocate ( trd%rey(trd%no), trd%iny(trd%no))
   allocate ( trd%xtl(trd%no), trd%ytl(trd%no) )
   allocate ( trd%xtl_ad(trd%no), trd%ytl_ad(trd%no) )

   allocate(  trd%i1(trd%im,trd%jm),  trd%j1(trd%im,trd%jm) )
   allocate( trd%pq1(trd%im,trd%jm), trd%pq2(trd%im,trd%jm) )
   allocate( trd%pq3(trd%im,trd%jm), trd%pq4(trd%im,trd%jm) )
   allocate( trd%uvl(trd%im,trd%jm), trd%vvl(trd%im,trd%jm) )
   allocate( trd%uvl_ad(trd%im,trd%jm), trd%vvl_ad(trd%im,trd%jm) )


   read(511)  trd%ino, trd%flg, trd%tim, trd%dtm,      &
              trd%loi, trd%lai, trd%lof, trd%laf,      &
              trd%err, trd%lob, trd%lab,               &
              trd%xob, trd%xmn, trd%erx,               &
              trd%yob, trd%ymn, trd%ery,               &
              trd%umn, trd%vmn, trd%dx, trd%dy,        &
              trd%lon, trd%lat

   close(511)

! ---
! Initialise quality flag
   if(obs%trd.eq.0)then
    do j=1,trd%no
     if(trd%flg(j).ne.0)trd%flg(j)=-1
    enddo
   endif
   trd%flc(:) = trd%flg(:)

    do j=1,trd%no
     trd%rex(j) = trd%xob(j) - trd%xmn(trd%nt+1,j)
     trd%rey(j) = trd%yob(j) - trd%ymn(trd%nt+1,j)
     trd%loa(j) = trd%lob(trd%nt+1,j)
     trd%laa(j) = trd%lab(trd%nt+1,j)
    enddo


! residual check
  do j=1,trd%no
   if(abs(trd%lof(j)-trd%lob(trd%nt+1,j)).gt.0.5 .or.           &
      abs(trd%laf(j)-trd%lab(trd%nt+1,j)).gt.0.5) trd%flg(j) = 0
  enddo

! ---
! Count good observations
    trd%nc = 0
  do k=1,trd%no
   if(trd%flg(k).eq.1)then
    trd%nc = trd%nc + 1
   endif
  enddo

  trd%flc(:) = trd%flg(:)


1111 continue
   

end subroutine get_obs_trd

subroutine int_par_trd

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str
 use obs_str
 use drv_str

 implicit none

  INTEGER(i4)   ::  ii, jj, i, j
  REAL(r8)      ::  ri, rj, p, q
  LOGICAL       ::  ins

 if(trd%no.gt.0) then

! ---
! Interpolate between grids

     do jj=1,trd%jm
     do ii=1,trd%im
       ri=max(1.,min(real(grd%im-1),real(ii-1)/real(drv%ratco(drv%ktr)) + 1.))
        i=int(ri)
        p=ri-i
       rj=max(1.,min(real(grd%jm-1),real(jj-1)/real(drv%ratco(drv%ktr)) + 1.))
        j=int(rj)
        q=rj-j

           trd%i1(ii,jj) = i
           trd%j1(ii,jj) = j

          trd%pq1(ii,jj) = (1.-p) * (1.-q) 
          trd%pq2(ii,jj) =     p  * (1.-q)
          trd%pq3(ii,jj) = (1.-p) *     q   
          trd%pq4(ii,jj) =     p  *     q

     enddo
     enddo


 endif



end subroutine int_par_trd
