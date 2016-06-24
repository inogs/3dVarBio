subroutine ver_hor
  
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
  ! Apply horizontal filter                                              !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  ! Version 2: S.Dobricic 2007                                           !
  !     Symmetric calculation in presence of coastal boundaries          !
  !     eta_ad, tem_ad, and sal_ad are here temporary arrays             !
  ! Version 3: A. Teruzzi 2013                                           !
  !     Attenuation of correction near the cost where d<200m             !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use eof_str
  use cns_str
  use drv_str
  use obs_str
  implicit none
  
  INTEGER(i4)    :: i,j,k, ione, l
  INTEGER        :: jp,nestr
  REAL(r8)          :: chlapp(8),chlsum
  
  ione = 1
  
  ! ---
  ! Vertical EOFs           
  call veof
  !return
  !goto 103 !No Vh
  
  ! ---
  ! Load temporary arrays
  if(drv%mask(drv%ktr) .gt. 1) then
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl_ad(:,:,k,l) = grd%chl(:,:,k,l) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
  endif
  
  !
  ! Attenuation of the correction near the cost and where d<200m 
  do l=1,grd%nchl
     do j=2,grd%jm-1
        do i=2,grd%im-1
           if ((grd%msk(i,j,chl%kdp).eq.0).and.  &
                (grd%msk(i,j,1).eq.1)) then
              do k=1,grd%km
                 if(grd%msk(i,j,k).eq.1) then
                    chlapp(1)=grd%chl(i+1,j,  k,l)
                    chlapp(2)=grd%chl(i-1,j,  k,l)
                    chlapp(3)=grd%chl(i,  j+1,k,l)
                    chlapp(4)=grd%chl(i,  j-1,k,l)
                    chlapp(5)=grd%chl(i+1,j+1,k,l)
                    chlapp(6)=grd%chl(i+1,j-1,k,l)
                    chlapp(7)=grd%chl(i-1,j+1,k,l)
                    chlapp(8)=grd%chl(i-1,j-1,k,l)
                    nestr=0
                    chlsum=0.
                    do jp=1,8
                       if ((chlapp(jp).ne.0).and.(chlapp(jp)/chlapp(jp).eq.1)) then
                          nestr=nestr+1;
                          chlsum=chlsum+chlapp(jp)
                       endif
                    enddo ! do on jp
                    if (nestr.ne.0) then
                       grd%chl(i,j,k,l)=.1*chlsum/nestr
                    endif
                 endif !if on k
              enddo ! do on k
           endif ! if on grd%chl(i,j,1,l)
        enddo ! do on i
     enddo ! do on j
  enddo ! do on l
  ! ---
  ! x direction
  call rcfl_x( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
  
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scx(:,:) 
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  enddo
  
  ! ---
  ! y direction
  call rcfl_y( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
  
  ! ---
  ! Scale by the scaling factor
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scy(:,:) 
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  enddo
  
  ! ---
  ! Transpose calculation in the presense of coastal boundaries
  if(drv%mask(drv%ktr) .gt. 1) then
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scy(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
     ! ---
     ! y direction
     call rcfl_y_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)

     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scx(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo

     ! ---
     ! x direction
     call rcfl_x_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
     
     ! ---
     ! Average
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl(:,:,k,l)   = (grd%chl(:,:,k,l) + grd%chl_ad(:,:,k,l) ) * 0.5 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
  endif
  
  
  ! ---
  ! Scale for boundaries
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl(:,:,k,l)   = grd%chl(:,:,k,l) * grd%fct(:,:,k)  
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  enddo
  
  
  !103 continue
  ! Correction is zero out of mask (for correction near the coast)
  do k=1,grd%km
     do j=1,grd%jm
        do i=1,grd%im
           if (grd%msk(i,j,k).eq.0) then
              grd%chl(i,j,k,:) = 0.
           endif
        enddo  !i
     enddo  !j
  enddo  !k
  

end subroutine ver_hor
