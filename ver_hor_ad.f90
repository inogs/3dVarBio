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

subroutine ver_hor_ad
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Transformation from physical to control space                        !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !     Symmetric calculation in presence of coastal boundaries          !
  !     eta, tem, and sal are here temporary arrays                      !
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
  INTEGER        :: jp, SurfaceIndex, TmpOffset, LinearIndex
  INTEGER(i4)    :: iProc, ierr
  type(DoubleGrid), allocatable, dimension(:,:,:) :: SendBuf3D
  type(DoubleGrid), allocatable, dimension(:)       :: RecBuf1D
  REAL(r8), allocatable, dimension(:,:,:) :: DefBufChl, DefBufChlAd  
  
  ione = 1
  
  if(drv%bio_assim .eq. 1) &
    call bio_mod_ad
  
  ! goto 103 ! No Vh

  ! ---
  ! Scale for boundaries
  do k=1,grd%km
    grd%chl_ad(:,:,k)   = grd%chl_ad(:,:,k) * grd%msk(:,:,k)
  enddo
  
  
  ! ---
  ! Load temporary arrays
  do k=1,grd%km
    grd%chl(:,:,k)    = grd%chl_ad(:,:,k)
  enddo
  
  ! ---
  ! y direction
  ! ---
  ! Scale by the scaling factor
  do k=1,grd%km
    grd%chl_ad(:,:,k) = grd%chl_ad(:,:,k) * grd%scy(:,:,k)
  enddo
  
  ! Apply recursive filter in y direction
  call rcfl_y_ad( localRow, GlobalCol, grd%km, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
  
  ! ---
  ! x direction
  if(NumProcI .gt. 1) then
     ALLOCATE(SendBuf3D(grd%km, grd%im, grd%jm))
     ALLOCATE( RecBuf1D(grd%km*localCol*GlobalRow))
     ALLOCATE(DefBufChl(GlobalRow, localCol, grd%km))
     ALLOCATE(DefBufChlAd(GlobalRow, localCol, grd%km))
     
     do k=1,grd%km
        do j=1,grd%jm
           do i=1,grd%im
              SendBuf3D(k,i,j)%chl = grd%chl(i,j,k)
           end do
        end do
     end do
     do k=1,grd%km
        do j=1,grd%jm
           do i=1,grd%im
              SendBuf3D(k,i,j)%chl_ad = grd%chl_ad(i,j,k)
           end do
        end do
     end do

     call MPI_Alltoallv(SendBuf3D, SendCountX3D, SendDisplX3D, MyPair, &
          RecBuf1D, RecCountX3D, RecDisplX3D, MyPair, Var3DCommunicator, ierr)
     
     SurfaceIndex = localCol*grd%km
     do j=1,localCol
        do iProc=0, NumProcI-1
           TmpOffset = RecDisplX3D(iProc+1)/SurfaceIndex
           do i=1,RecCountX3D(iProc+1)/SurfaceIndex
              LinearIndex = (i-1)*grd%km + (j-1)*RecCountX3D(iProc+1)/localCol + RecDisplX3D(iProc+1)
              do k=1,grd%km
                 DefBufChl(i + TmpOffset,j,k) = RecBuf1D(k + LinearIndex)%chl
              end do
           end do
        end do
     end do
     do j=1,localCol
        do iProc=0, NumProcI-1
           TmpOffset = RecDisplX3D(iProc+1)/SurfaceIndex
           do i=1,RecCountX3D(iProc+1)/SurfaceIndex
              LinearIndex = (i-1)*grd%km + (j-1)*RecCountX3D(iProc+1)/localCol + RecDisplX3D(iProc+1)
              do k=1,grd%km
                 DefBufChlAd(i + TmpOffset,j,k) = RecBuf1D(k + LinearIndex)%chl_ad
              end do
           end do
        end do
     end do
     
     ! ---
     ! Scale by the scaling factor
     do k=1,grd%km
        DefBufChlAd(:,:,k) = DefBufChlAd(:,:,k) * grd%scx(:,:,k)
     enddo
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km, grd%imax, grd%aex, grd%bex, DefBufChlAd, grd%inx, grd%imx)
     
  else ! NumProcI .eq. 1
     ! ---
     ! Scale by the scaling factor
     do k=1,grd%km
        grd%chl_ad(:,:,k) = grd%chl_ad(:,:,k) * grd%scx(:,:,k)
     enddo
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
  end if
  
  
  ! ---
  ! x direction
  if(NumProcI .gt. 1) then
     
     call rcfl_x( GlobalRow, localCol, grd%km, grd%imax, grd%aex, grd%bex, DefBufChl, grd%inx, grd%imx)
     
     ! ---
     ! Scale by the scaling factor
     do k=1,grd%km
        DefBufChl(:,:,k) = DefBufChl(:,:,k) * grd%scx(:,:,k)
     enddo
     
     ! Reordering data to send back
     DEALLOCATE(SendBuf3D, RecBuf1D)
     ALLOCATE(SendBuf3D(grd%km, localCol, GlobalRow))
     ALLOCATE( RecBuf1D(grd%km*grd%jm*grd%im))
     
     do k=1,grd%km
        do j=1,localCol
           do i=1,GlobalRow
              SendBuf3D(k,j,i)%chl = DefBufChl(i,j,k)
           end do
        end do
     end do
     do k=1,grd%km
        do j=1,localCol
           do i=1,GlobalRow
              SendBuf3D(k,j,i)%chl_ad = DefBufChlAd(i,j,k)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf3D, RecCountX3D, RecDisplX3D, MyPair, &
          RecBuf1D, SendCountX3D, SendDisplX3D, MyPair, Var3DCommunicator, ierr)
     
     SurfaceIndex = grd%im*grd%km
     do i=1,grd%im
        do iProc=0, NumProcI-1
           TmpOffset = SendDisplX3D(iProc+1)/SurfaceIndex
           do j=1,SendCountX3D(iProc+1)/SurfaceIndex
              LinearIndex = (j-1)*grd%km +(i-1)*SendCountX3D(iProc+1)/grd%im + SendDisplX3D(iProc+1)
              do k=1,grd%km
                 grd%chl(i, j + TmpOffset,k) = RecBuf1D(k + LinearIndex)%chl
              end do
           end do
        end do
     end do
     do i=1,grd%im
        do iProc=0, NumProcI-1
           TmpOffset = SendDisplX3D(iProc+1)/SurfaceIndex
           do j=1,SendCountX3D(iProc+1)/SurfaceIndex
              LinearIndex = (j-1)*grd%km +(i-1)*SendCountX3D(iProc+1)/grd%im + SendDisplX3D(iProc+1)
              do k=1,grd%km
                 grd%chl_ad(i, j + TmpOffset,k) = RecBuf1D(k + LinearIndex)%chl_ad
              end do
           end do
        end do
     end do

     DEALLOCATE(SendBuf3D, RecBuf1D, DefBufChl, DefBufChlAd)
     
  else ! NumProcI .eq. 1
     call rcfl_x( GlobalRow, localCol, grd%km, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
     
     ! ---
     ! Scale by the scaling factor
     do k=1,grd%km
        grd%chl(:,:,k) = grd%chl(:,:,k) * grd%scx(:,:,k)
     enddo
  end if
  
  
  ! ! ---
  ! ! y direction
  ! Apply recursive filter in y direction
  call rcfl_y( localRow, GlobalCol, grd%km, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
  
  ! ---
  ! Scale by the scaling factor
  do k=1,grd%km
     grd%chl(:,:,k) = grd%chl(:,:,k) * grd%scy(:,:,k)
  enddo
  
  
  ! ---
  ! Average
  do k=1,grd%km
     grd%chl_ad(:,:,k)  = (grd%chl_ad(:,:,k) + grd%chl(:,:,k) ) * 0.5
  enddo
  
  
  ! 103 continue
  ! ---
  ! Vertical EOFs
  call veof_ad
  
end subroutine ver_hor_ad
