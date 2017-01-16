subroutine parallel_ver_hor

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
  use mpi_str

  implicit none

  INTEGER(i4)    :: i,j,k, ione, l
  INTEGER        :: jp, SurfaceIndex
  REAL(r8)          :: chlapp(8),chlsum
  INTEGER(i4)    :: iProc, ierr
  ! type(DoubleGrid), allocatable, dimension(:,:,:,:) :: SendBuf4D
  ! type(DoubleGrid), allocatable, dimension(:)       :: RecBuf1D(:)
  ! REAL(r8), allocatable, dimension(:,:,:,:) :: DefBufChl, DefBufChlAd
  
  ione = 1
  
  ! ---
  ! Vertical EOFs
  if(MyRank .eq. 0) then
    call veof
    !return
    !goto 103 !No Vh
  
    ! ---
    ! Load temporary arrays
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
  
    !********** APPLY RECURSIVE FILTERS ********** !
    ! ---
    ! y direction
  
    ! ---
    ! Scale by the scaling factor
  !   do l=1,grd%nchl
  !     !$OMP PARALLEL  &
  !     !$OMP PRIVATE(k)
  !     !$OMP DO
  !     do k=1,grd%km
  !         grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scy(:,:)
  !     enddo
  !     !$OMP END DO
  !     !$OMP END PARALLEL
  !   enddo

  endif
  
  ! Apply recursive filter in y direction
  call rcfl_y_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
  
  ! ---
  ! x direction

  ! if(MyRank .eq. 0) then 
  !    ! ---
  !    ! Scale by the scaling factor
  !    do l=1,grd%nchl
  !       !$OMP PARALLEL  &
  !       !$OMP PRIVATE(k)
  !       !$OMP DO
  !       do k=1,grd%km
  !          grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scx(:,:)
  !       enddo
  !       !$OMP END DO
  !       !$OMP END PARALLEL
  !    enddo
  ! endif
  
  call rcfl_x_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
     
  
  
  ! ---
  ! x direction
  call rcfl_x( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)

  ! if(MyRank .eq. 0) then     
  !    do l=1,grd%nchl
  !       !$OMP PARALLEL  &
  !       !$OMP PRIVATE(k)
  !       !$OMP DO
  !       do k=1,grd%km
  !          grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scx(:,:)
  !       enddo
  !       !$OMP END DO
  !       !$OMP END PARALLEL
  !    enddo
     
  ! endif
  
  ! ---
  ! y direction
  ! Apply recursive filter in y direction
  call rcfl_y( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
  
  ! ---
  ! Scale by the scaling factor
  if(MyRank .eq. 0) then
    
    ! do l=1,grd%nchl
    !   !$OMP PARALLEL  &
    !   !$OMP PRIVATE(k)
    !   !$OMP DO
    !   do k=1,grd%km
    !     grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scy(:,:)
    !   enddo
    !   !$OMP END DO
    !   !$OMP END PARALLEL
    ! enddo
  
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
  
    ! ---
    ! Scale for boundaries
    do l=1,grd%nchl
      !$OMP PARALLEL  &
      !$OMP PRIVATE(k)
      !$OMP DO
      do k=1,grd%km
         grd%chl(:,:,k,l)   = grd%chl(:,:,k,l) * grd%msk(:,:,k)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
    
  
    !103 continue
    ! Correction is zero out of mask (for correction near the coast)
    ! THE FOLLOWING CYCLE IS REPLACED WITH THE PREVIOUS ONE 
    ! (BY MEANS OF THE PRODUCT WITH THE GRID MASK)
    ! do k=1,grd%km
    !   do j=1,grd%jm
    !       do i=1,grd%im
    !         if (grd%msk(i,j,k).eq.0) then
    !             grd%chl(i,j,k,:) = 0.
    !         endif
    !       enddo  !i
    !   enddo  !j
    ! enddo  !k
  
  endif ! MyRank .eq. 0

end subroutine parallel_ver_hor

subroutine parallel_ver_hor_ad
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Transformation from physical to control space                        !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  ! Version 2: S.Dobricic 2007                                           !
  !     Symmetric calculation in presence of coastal boundaries          !
  !     eta, tem, and sal are here temporary arrays                      !
  ! Version 3: A.Teruzzi 2013                                            !
  !     Smoothing of the solution at d<200m                              !
  !-----------------------------------------------------------------------


  use set_knd
  use grd_str
  use eof_str
  use cns_str
  use drv_str
  use obs_str
  use mpi_str

  implicit none

  INTEGER(i4)    :: i,j,k, ione, l
  INTEGER        :: jp, SurfaceIndex
  REAL(r8)       :: chlapp(8),chlsum
  INTEGER(i4)    :: iProc, ierr
  ! type(DoubleGrid), allocatable, dimension(:,:,:,:) :: SendBuf4D
  ! type(DoubleGrid), allocatable, dimension(:)       :: RecBuf1D
  ! REAL(r8), allocatable, dimension(:,:,:,:) :: DefBufChl, DefBufChlAd

  if(MyRank .eq. 0) then
    ! ---
    ! Correction is zero out of mask (for correction near the coast)
    ! THE FOLLOWING CYCLE IS REPLACED WITH THE PREVIOUS ONE 
    ! (BY MEANS OF THE PRODUCT WITH THE GRID MASK)
    ! do k=1,grd%km
    !   do j=1,grd%jm
    !       do i=1,grd%im
    !         if (grd%msk(i,j,k).eq.0) then
    !             grd%chl_ad(i,j,k,:) = 0.
    !         endif
    !       enddo  !i
    !   enddo  !j
    ! enddo  !k
  
  
    !goto 103 ! No Vh
    ione = 1
  
    ! ---
    ! Scale for boundaries
    do l=1,grd%nchl
      !$OMP PARALLEL  &
      !$OMP PRIVATE(k)
      !$OMP DO
      do k=1,grd%km
         grd%chl_ad(:,:,k,l)   = grd%chl_ad(:,:,k,l) * grd%msk(:,:,k)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
  
  
    ! ---
    ! Load temporary arrays
    do l=1,grd%nchl
      !$OMP PARALLEL  &
      !$OMP PRIVATE(k)
      !$OMP DO
      do k=1,grd%km
          grd%chl(:,:,k,l)    = grd%chl_ad(:,:,k,l)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
  
    ! ---
    ! y direction
    ! ---
    ! Scale by the scaling factor
    ! do l=1,grd%nchl
    !   !$OMP PARALLEL  &
    !   !$OMP PRIVATE(k)
    !   !$OMP DO
    !   do k=1,grd%km
    !       grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scy(:,:)
    !   enddo
    !   !$OMP END DO
    !   !$OMP END PARALLEL
    ! enddo
  
  endif

  ! Apply recursive filter in y direction
  call rcfl_y_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
  
  ! ---
  ! x direction
  ! if(MyRank .eq. 0) then
  !    ! ---
  !    ! Scale by the scaling factor
  !    do l=1,grd%nchl
  !       !$OMP PARALLEL  &
  !       !$OMP PRIVATE(k)
  !       !$OMP DO
  !       do k=1,grd%km
  !          grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scx(:,:)
  !       enddo
  !       !$OMP END DO
  !       !$OMP END PARALLEL
  !    enddo
     
  ! endif
    
  call rcfl_x_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)

  
  ! ---
  ! x direction
  call rcfl_x( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
     
  ! ---
  ! Scale by the scaling factor
  ! if(MyRank .eq. 0) then
  !    do l=1,grd%nchl
  !       !$OMP PARALLEL  &
  !       !$OMP PRIVATE(k)
  !       !$OMP DO
  !       do k=1,grd%km
  !          grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scx(:,:)
  !       enddo
  !       !$OMP END DO
  !       !$OMP END PARALLEL
  !    enddo
  ! end if
  
  
  ! ! ---
  ! ! y direction
  ! Apply recursive filter in y direction
  call rcfl_y( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
  
  ! ---
  ! Scale by the scaling factor
  if(MyRank .eq. 0) then
  !   do l=1,grd%nchl
  !     !$OMP PARALLEL  &
  !     !$OMP PRIVATE(k)
  !     !$OMP DO
  !     do k=1,grd%km
  !         grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scy(:,:)
  !     enddo
  !     !$OMP END DO
  !     !$OMP END PARALLEL
  !   enddo
  
  
    ! ---
    ! Average
    do l=1,grd%nchl
      !$OMP PARALLEL  &
      !$OMP PRIVATE(k)
      !$OMP DO
      do k=1,grd%km
          grd%chl_ad(:,:,k,l)  = (grd%chl_ad(:,:,k,l) + grd%chl(:,:,k,l) ) * 0.5
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo
  
  
    !103 continue
    ! ---
    ! Vertical EOFs
    call veof_ad

  endif
  
end subroutine parallel_ver_hor_ad
