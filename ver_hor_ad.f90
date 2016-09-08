subroutine ver_hor_ad
  
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
#ifdef _USE_MPI
  use mpi_str
#endif
  implicit none
  
  INTEGER(i4)    :: i,j,k, ione, l
  INTEGER        :: jp,nestr
  REAL(r8)       :: chlapp(8),chlsum
#ifdef _USE_MPI
  INTEGER(i4)    :: iProc, ierr
  REAL(r8), allocatable :: SendBuf4D(:,:,:,:), RecBuf4D(:,:,:,:), DefBuf4D(:,:,:,:)
#endif
  
  ! ---
  ! Correction is zero out of mask (for correction near the coast)
  do k=1,grd%km
     do j=1,grd%jm
        do i=1,grd%im
           if (grd%msk(i,j,k).eq.0) then
              grd%chl_ad(i,j,k,:) = 0.
           endif
        enddo  !i
     enddo  !j
  enddo  !k
  
  
  !goto 103 ! No Vh
  ione = 1
    
  ! ---
  ! Scale for boundaries
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl_ad(:,:,k,l)   = grd%chl_ad(:,:,k,l) * grd%fct(:,:,k) 
     enddo
     !$OMP END DO
     !$OMP END PARALLEL  
  enddo
  
  
  if(drv%mask(drv%ktr) .gt. 1) then
     
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
     ! x direction
#ifdef _USE_MPI
     call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
#else
     call rcfl_x( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
#endif

     ! ---
     ! Scale by the scaling factor
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
#ifdef _USE_MPI
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
     ALLOCATE( RecBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
     ALLOCATE( DefBuf4D(localRow, GlobalCol, grd%km, grd%nchl))
     
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,j,i) = grd%chl(i,j,k,l)
              end do
           end do
        end do
     end do
     
     call MPI_Alltoall(SendBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, &
          RecBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, MPI_COMM_WORLD, ierr)
     
     do i=1,localRow
        do iProc=0, Size-1
           do j=1,grd%jm
              do k=1,grd%km
                 DefBuf4D(i,j + iProc*localCol,k,1) = RecBuf4D(1,k,j,i + iProc*localRow)
              end do
           end do
        end do
     end do
     
     ! Apply recursive filter in y direction
     call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)
     
     ! Reordering data to send back
     DEALLOCATE(SendBuf4D, RecBuf4D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
     ALLOCATE( RecBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
     
     do j=1,GlobalCol
        do i=1,localRow
           do k=1,grd%km
              SendBuf4D(1,k,i,j) = DefBuf4D(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoall(SendBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, &
          RecBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, MPI_COMM_WORLD, ierr)
     
     do i=1,localRow
        do iProc=0, Size-1
           do j=1,grd%jm
              do k=1,grd%km
                 grd%chl(i + iProc*localRow,j,k,1) = RecBuf4D(1,k,i,j + iProc*localCol)
              end do
           end do
        end do
     end do
     
#else
     call rcfl_y( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
#endif
     
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
     
  endif
  
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
#ifdef _USE_MPI
  DEALLOCATE(SendBuf4D, RecBuf4D)
  ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
  ALLOCATE( RecBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
  ! ALLOCATE( DefBuf4D(localRow, GlobalCol, grd%km, grd%nchl))
  
  do l=1,grd%nchl
     do k=1,grd%km
        do j=1,grd%jm
           do i=1,grd%im
              SendBuf4D(l,k,j,i) = grd%chl_ad(i,j,k,l)
           end do
        end do
     end do
  end do
  
  call MPI_Alltoall(SendBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, &
       RecBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, MPI_COMM_WORLD, ierr)
  
  do i=1,localRow
     do iProc=0, Size-1
        do j=1,grd%jm
           do k=1,grd%km
              DefBuf4D(i,j + iProc*localCol,k,1) = RecBuf4D(1,k,j,i + iProc*localRow)
           end do
        end do
     end do
  end do

  ! Apply recursive filter in y direction
  call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)

  ! Reordering data to send back
  DEALLOCATE(SendBuf4D, RecBuf4D)
  ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
  ALLOCATE( RecBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
  
  do j=1,GlobalCol
     do i=1,localRow
        do k=1,grd%km
           SendBuf4D(1,k,i,j) = DefBuf4D(i,j,k,1)
        end do
     end do
  end do
  
  call MPI_Alltoall(SendBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, &
       RecBuf4D, grd%nchl*grd%km*grd%jm*grd%im/size, MPI_REAL8, MPI_COMM_WORLD, ierr)
  
  do i=1,localRow
     do iProc=0, Size-1
        do j=1,grd%jm
           do k=1,grd%km
              grd%chl_ad(i + iProc*localRow,j,k,1) = RecBuf4D(1,k,i,j + iProc*localCol)
           end do
        end do
     end do
  end do
  
#else
  call rcfl_y_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
#endif
  !    write(*,*) 'AFTER rcfl_y_ad', grd%chl_ad(71,12,1,1)
  
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
#ifdef _USE_MPI
  call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
#else
  call rcfl_x_ad( grd%im, grd%jm, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
#endif
  
  
  ! ---
  ! Average
  if(drv%mask(drv%ktr) .gt. 1) then
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
  endif
  
  !anna sreduction of correction d<200m
  do l=1,grd%nchl
     do j=2,grd%jm-1  ! OMP
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
                    do jp=1,8
                       if ((chlapp(jp).ne.0).and.(chlapp(jp)/chlapp(jp).eq.1)) then
                          nestr=nestr+1;
                       endif
                    enddo ! do on jp
                    if (nestr.ne.0) then
                       grd%chl_ad(i+1,j,  k,l)=grd%chl_ad(i+1,j,  k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i-1,j,  k,l)=grd%chl_ad(i-1,j,  k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i,  j+1,k,l)=grd%chl_ad(i  ,j+1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i,  j-1,k,l)=grd%chl_ad(i  ,j-1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i+1,j+1,k,l)=grd%chl_ad(i+1,j+1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i+1,j-1,k,l)=grd%chl_ad(i+1,j-1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i-1,j+1,k,l)=grd%chl_ad(i-1,j+1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i-1,j-1,k,l)=grd%chl_ad(i-1,j-1,k,l)+  &
                            .1*grd%chl_ad(i,j,k,l)/nestr
                       grd%chl_ad(i,j,k,l)=0.
                    endif
                 endif !if on k
              enddo ! do on k
           endif ! if on grd%chl(i,j,1,l)
        enddo ! do on i
     enddo ! do on j
  enddo ! do on l
  
  
  
  
  !103 continue
  ! ---
  ! Vertical EOFs           
  call veof_ad
  ! print*, "DONE WITH VER_HOR_AD"
#ifdef _USE_MPI
  DEALLOCATE(SendBuf4D, RecBuf4D, DefBuf4D)
#endif

  
end subroutine ver_hor_ad
