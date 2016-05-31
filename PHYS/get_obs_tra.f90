subroutine get_obs_tra

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
! Load Argo trajectory observations                                    !
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

   tra%no = 0
   tra%nc = 0

   
    open(511,file='tra_obs.dat',form='unformatted',status='old',err=1111)

    read(511) tra%no, tra%nc, tra%nt, tra%im, tra%jm, tra%km, tra%dpt

! ---
! Allocate memory for observations 

   write(drv%dia,*) 'Number of Argo trajectory observations: ',  tra%no,tra%nc


   if(tra%no.eq.0)then
      close(511)
      return
   endif

   allocate ( tra%ino(tra%no), tra%flg(tra%no), tra%flc(tra%no))
   allocate ( tra%loi(tra%no), tra%lai(tra%no), tra%tim(tra%no), tra%dtm(tra%no))
   allocate ( tra%lof(tra%no), tra%laf(tra%no))
   allocate ( tra%err(tra%no))
   allocate ( tra%lob(tra%nt+1,tra%no), tra%lab(tra%nt+1,tra%no) )
   allocate ( tra%loa(tra%no), tra%laa(tra%no) )
   allocate ( tra%xob(tra%no), tra%xmn(tra%nt+1,tra%no), tra%erx(tra%no) )
   allocate ( tra%yob(tra%no), tra%ymn(tra%nt+1,tra%no), tra%ery(tra%no) )
   allocate ( tra%umn(tra%im,tra%jm), tra%vmn(tra%im,tra%jm) )
   allocate ( tra%dx(tra%im,tra%jm), tra%dy(tra%im,tra%jm) )
   allocate ( tra%lon(tra%im,tra%jm), tra%lat(tra%im,tra%jm) )

   allocate ( tra%rex(tra%no), tra%inx(tra%no))
   allocate ( tra%rey(tra%no), tra%iny(tra%no))
   allocate ( tra%xtl(tra%no), tra%ytl(tra%no) )
   allocate ( tra%xtl_ad(tra%no), tra%ytl_ad(tra%no) )

   allocate(  tra%i1(tra%im,tra%jm),  tra%j1(tra%im,tra%jm) )
   allocate( tra%pq1(tra%im,tra%jm), tra%pq2(tra%im,tra%jm) )
   allocate( tra%pq3(tra%im,tra%jm), tra%pq4(tra%im,tra%jm) )
   allocate( tra%uvl(tra%im,tra%jm), tra%vvl(tra%im,tra%jm) )
   allocate( tra%uvl_ad(tra%im,tra%jm), tra%vvl_ad(tra%im,tra%jm) )


   read(511)  tra%ino, tra%flg, tra%tim, tra%dtm,      &
              tra%loi, tra%lai, tra%lof, tra%laf,      &
              tra%err, tra%lob, tra%lab,               &
              tra%xob, tra%xmn, tra%erx,               &
              tra%yob, tra%ymn, tra%ery,               &
              tra%umn, tra%vmn, tra%dx, tra%dy,        &
              tra%lon, tra%lat

   close(511)

! ---
! Initialise quality flag
   if(obs%tra.eq.0)then
    do j=1,tra%no
     if(tra%flg(j).ne.0)tra%flg(j)=-1
    enddo
   endif
   tra%flc(:) = tra%flg(:)

    do j=1,tra%no
     tra%rex(j) = tra%xob(j) - tra%xmn(tra%nt+1,j)
     tra%rey(j) = tra%yob(j) - tra%ymn(tra%nt+1,j)
     tra%loa(j) = tra%lob(tra%nt+1,j)
     tra%laa(j) = tra%lab(tra%nt+1,j)
    enddo


! residual check
  do j=1,tra%no
   if(abs(tra%rex(j)).gt.10. .or. abs(tra%rey(j)).gt.10.) tra%flg(j) = 0
  enddo

! ---
! Count good observations
    tra%nc = 0
  do k=1,tra%no
   if(tra%flg(k).eq.1)then
    tra%nc = tra%nc + 1
   endif
  enddo

  tra%flc(:) = tra%flg(:)

   write(drv%dia,*) 'Number of good Argo trajectory observations: ',  tra%nc

1111 continue


end subroutine get_obs_tra

subroutine int_par_tra

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

 if(tra%no.gt.0) then

! ---
! Interpolate between grids

     do jj=1,tra%jm
     do ii=1,tra%im
       ri=max(1.,min(real(grd%im-1),real(ii-1)/real(drv%ratco(drv%ktr)) + 1.))
        i=int(ri)
        p=ri-i
       rj=max(1.,min(real(grd%jm-1),real(jj-1)/real(drv%ratco(drv%ktr)) + 1.))
        j=int(rj)
        q=rj-j

           tra%i1(ii,jj) = i
           tra%j1(ii,jj) = j

          tra%pq1(ii,jj) = (1.-p) * (1.-q) 
          tra%pq2(ii,jj) =     p  * (1.-q)
          tra%pq3(ii,jj) = (1.-p) *     q   
          tra%pq4(ii,jj) =     p  *     q

     enddo
     enddo


 endif



end subroutine int_par_tra
