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
  INTEGER        :: jp,nestr
  REAL(r8)          :: chlapp(8),chlsum
  INTEGER(i4)    :: iProc, ierr
  REAL(r8), allocatable, dimension(:,:,:,:) :: SendBuf4D, DefBuf4D
  REAL(r8), allocatable, dimension(:)       :: RecBuf1D(:)
  
  INTEGER   :: ReqRecvRight, ReqSendRight, ReqSendLeft, ReqRecvLeft
  INTEGER   :: ReqRecvTop, ReqSendTop, ReqSendBottom, ReqRecvBottom
  INTEGER   :: StatRight(MPI_STATUS_SIZE), StatLeft(MPI_STATUS_SIZE)
  INTEGER   :: StatTop(MPI_STATUS_SIZE), StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag  

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
  
  !********** APPLY RECURSIVE FILTERS ********** !

  ! ---
  ! x direction

  if(NumProcI .gt. 1) then
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%im, grd%jm))
     ALLOCATE( RecBuf1D(GlobalRow*localCol*grd%km*grd%nchl))
     ALLOCATE( DefBuf4D(GlobalRow, localCol, grd%km, grd%nchl))
     
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,i,j) = grd%chl(i,j,k,l)
              end do
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, SendCountX4D, SendDisplX4D, MPI_REAL8, &
          RecBuf1D, RecCountX4D, RecDisplX4D, MPI_REAL8, CommSliceX, ierr)
     
     do j=1,localCol
        do iProc=0, NumProcI-1
           do i=1,RecCountX4D(iProc+1)/(localCol*grd%km)
              do k=1,grd%km
                 DefBuf4D(i + RecDisplX4D(iProc+1)/(localCol*grd%km),j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km +(j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))
              end do
           end do
        end do
     end do
     
     call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBuf4D, grd%inx, grd%imx)
     
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scx(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
     ! Reordering data to send back
     DEALLOCATE(SendBuf4D, RecBuf1D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, localCol, GlobalRow))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
     
     do k=1,grd%km
        do j=1,localCol
           do i=1,GlobalRow
              SendBuf4D(1,k,j,i) = DefBuf4D(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, RecCountX4D, RecDisplX4D, MPI_REAL8, &
          RecBuf1D, SendCountX4D, SendDisplX4D, MPI_REAL8, CommSliceX, ierr)
     
     do i=1,grd%im
        do iProc=0, NumProcI-1
           do j=1,SendCountX4D(iProc+1)/(grd%im*grd%km)
              do k=1,grd%km
                 grd%chl(i, j + SendDisplX4D(iProc+1)/(grd%im*grd%km),k,1) = &
                      RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))
              end do
           end do
        end do
     end do
     DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)

  else ! NumProcI .eq. 1

     call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
     
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

  end if

  ! ---
  ! y direction
  if(NumProcJ .gt. 1) then
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*localRow*GlobalCol))
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
     
     call MPI_Alltoallv(SendBuf4D, SendCountY4D, SendDisplY4D, MPI_REAL8, &
          RecBuf1D, RecCountY4D, RecDisplY4D, MPI_REAL8, CommSliceY, ierr)
  
     do i=1,localRow
        do iProc=0, NumProcJ-1
           do j=1,RecCountY4D(iProc+1)/(localRow*grd%km)
              do k=1,grd%km
                 DefBuf4D(i,j+RecDisplY4D(iProc+1)/(localRow*grd%km),k,1) = &
                      RecBuf1D(k + (j-1)*grd%km + (i-1)*RecCountY4D(iProc+1)/localRow + RecDisplY4D(iProc+1))
              end do
           end do
        end do
     end do
     
     ! Apply recursive filter in y direction
     call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scy(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
     ! Reordering data to send back
     DEALLOCATE(SendBuf4D, RecBuf1D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
     
     do j=1,GlobalCol
        do i=1,localRow
           do k=1,grd%km
              SendBuf4D(1,k,i,j) = DefBuf4D(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, RecCountY4D, RecDisplY4D, MPI_REAL8, &
          RecBuf1D, SendCountY4D, SendDisplY4D, MPI_REAL8, CommSliceY, ierr)
     
     do j=1,grd%jm
        do iProc=0, NumProcJ-1
           do i=1, SendCountY4D(iProc+1)/(grd%jm*grd%km)
              do k=1,grd%km
                 grd%chl(i + SendDisplY4D(iProc+1)/(grd%jm*grd%km),j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km + (j-1)*SendCountY4D(iProc+1)/grd%jm + SendDisplY4D(iProc+1))
              end do
           end do
        end do
     end do
     DEALLOCATE(SendBuf4D, RecBuf1D)
  else
     ! Apply recursive filter in y direction
     call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
     
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
  end if
  
  ! ---
  ! Transpose calculation in the presense of coastal boundaries
  if(drv%mask(drv%ktr) .gt. 1) then
          
     ! ---
     ! y direction
     if(NumProcJ .gt. 1) then
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*localRow*GlobalCol))
     
        do l=1,grd%nchl
           do k=1,grd%km
              do j=1,grd%jm
                 do i=1,grd%im
                    SendBuf4D(l,k,j,i) = grd%chl_ad(i,j,k,l)
                 end do
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, SendCountY4D, SendDisplY4D, MPI_REAL8, &
             RecBuf1D, RecCountY4D, RecDisplY4D, MPI_REAL8, CommSliceY, ierr)
        
        do i=1,localRow
           do iProc=0, NumProcJ-1
              do j=1,RecCountY4D(iProc+1)/(localRow*grd%km)
                 do k=1,grd%km
                    DefBuf4D(i,j+RecDisplY4D(iProc+1)/(localRow*grd%km),k,1) = &
                         RecBuf1D(k + (j-1)*grd%km + (i-1)*RecCountY4D(iProc+1)/localRow + RecDisplY4D(iProc+1))
                 end do
              end do
           end do
        end do
        
        ! ---
        ! Scale by the scaling factor
        do l=1,grd%nchl
           !$OMP PARALLEL  &
           !$OMP PRIVATE(k)
           !$OMP DO
           do k=1,grd%km
              DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scy(:,:) 
           enddo
           !$OMP END DO
           !$OMP END PARALLEL
        enddo

        ! Apply recursive filter in y direction
        call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)
        
        DEALLOCATE(SendBuf4D, RecBuf1D)
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
        
        do l=1,grd%nchl
           do j=1,GlobalCol
              do i=1,localRow
                 do k=1,grd%km
                    SendBuf4D(l,k,i,j) = DefBuf4D(i,j,k,1)
                 end do
              end do
           end do
        end do
  
        call MPI_Alltoallv(SendBuf4D, RecCountY4D, RecDisplY4D, MPI_REAL8, &
             RecBuf1D, SendCountY4D, SendDisplY4D, MPI_REAL8, CommSliceY, ierr)
        
        do j=1,grd%jm
           do iProc=0, NumProcJ-1
              do i=1, SendCountY4D(iProc+1)/(grd%jm*grd%km)
                 do k=1,grd%km
                    grd%chl_ad(i + SendDisplY4D(iProc+1)/(grd%jm*grd%km),j,k,1) = &
                         RecBuf1D(k + (i-1)*grd%km + (j-1)*SendCountY4D(iProc+1)/grd%jm + SendDisplY4D(iProc+1))
                 end do
              end do
           end do
        end do
        DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)
     else
        
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

        ! Apply recursive filter in y direction
        call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
     end if
     
     ! ---
     ! x direction
     if(NumProcI .gt. 1) then
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%im, grd%jm))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*GlobalRow*localCol))
        ALLOCATE( DefBuf4D(GlobalRow, localCol, grd%km, grd%nchl))
        
        do l=1,grd%nchl
           do k=1,grd%km
              do j=1,grd%jm
                 do i=1,grd%im
                    SendBuf4D(l,k,i,j) = grd%chl_ad(i,j,k,l)
                 end do
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, SendCountX4D, SendDisplX4D, MPI_REAL8, &
             RecBuf1D, RecCountX4D, RecDisplX4D, MPI_REAL8, CommSliceX, ierr)
        
        do j=1,localCol
           do iProc=0, NumProcI-1
              do i=1,RecCountX4D(iProc+1)/(localCol*grd%km)
                 do k=1,grd%km
                    DefBuf4D(i + RecDisplX4D(iProc+1)/(localCol*grd%km),j,k,1) = &
                         RecBuf1D(k + (i-1)*grd%km + (j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))
                 end do
              end do
           end do
        end do
        
        ! ---
        ! Scale by the scaling factor
        do l=1,grd%nchl
           !$OMP PARALLEL  &
           !$OMP PRIVATE(k)
           !$OMP DO
           do k=1,grd%km
              DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scx(:,:) 
           enddo
           !$OMP END DO
           !$OMP END PARALLEL
        enddo
        
        call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBuf4D, grd%inx, grd%imx)
        
        ! Reordering data to send back
        DEALLOCATE(SendBuf4D, RecBuf1D)
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, localCol, GlobalRow))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
        
        do k=1,grd%km
           do j=1,localCol
              do i=1,GlobalRow
                 SendBuf4D(1,k,j,i) = DefBuf4D(i,j,k,1)
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, RecCountX4D, RecDisplX4D, MPI_REAL8, &
             RecBuf1D, SendCountX4D, SendDisplX4D, MPI_REAL8, CommSliceX, ierr)
        
        do i=1,grd%im
           do iProc=0, NumProcI-1
              do j=1,SendCountX4D(iProc+1)/(grd%im*grd%km)
                 do k=1,grd%km
                    grd%chl_ad(i, j + SendDisplX4D(iProc+1)/(grd%im*grd%km),k,1) = &
                         RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))
                 end do
              end do
           end do
        end do
        DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)

     else
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
        
        call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)

     end if
     
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
  INTEGER        :: jp,nestr
  REAL(r8)       :: chlapp(8),chlsum
  INTEGER(i4)    :: iProc, ierr
  REAL(r8), allocatable, dimension(:,:,:,:) :: SendBuf4D, DefBuf4D
  REAL(r8), allocatable, dimension(:)       :: RecBuf1D
  INTEGER   :: ReqRecvRight, ReqSendRight, ReqSendLeft, ReqRecvLeft
  INTEGER   :: ReqRecvTop, ReqSendTop, ReqSendBottom, ReqRecvBottom
  INTEGER   :: StatRight(MPI_STATUS_SIZE), StatLeft(MPI_STATUS_SIZE)
  INTEGER   :: StatTop(MPI_STATUS_SIZE), StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag
  
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
     if(NumProcI .gt. 1) then
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%im, grd%jm))
        ALLOCATE( RecBuf1D(GlobalRow*localCol*grd%km*grd%nchl))
        ALLOCATE( DefBuf4D(GlobalRow, localCol, grd%km, grd%nchl))
        
        do l=1,grd%nchl
           do k=1,grd%km
              do j=1,grd%jm
                 do i=1,grd%im
                    SendBuf4D(l,k,i,j) = grd%chl(i,j,k,l)
                 end do
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, SendCountX4D, SendDisplX4D, MPI_REAL8, &
             RecBuf1D, RecCountX4D, RecDisplX4D, MPI_REAL8, CommSliceX, ierr)
        
        do j=1,localCol
           do iProc=0, NumProcI-1
              do i=1,RecCountX4D(iProc+1)/(localCol*grd%km)
                 do k=1,grd%km
                    DefBuf4D(i + RecDisplX4D(iProc+1)/(localCol*grd%km),j,k,1) = &
                         RecBuf1D(k + (i-1)*grd%km +(j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))
                 end do
              end do
           end do
        end do
     
        call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBuf4D, grd%inx, grd%imx)
        
        ! ---
        ! Scale by the scaling factor
        do l=1,grd%nchl
           !$OMP PARALLEL  &
           !$OMP PRIVATE(k)
           !$OMP DO
           do k=1,grd%km
              DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scx(:,:) 
           enddo
           !$OMP END DO
           !$OMP END PARALLEL  
        enddo
        
        ! Reordering data to send back
        DEALLOCATE(SendBuf4D, RecBuf1D)
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, localCol, GlobalRow))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
        
        do k=1,grd%km
           do j=1,localCol
              do i=1,GlobalRow
                 SendBuf4D(1,k,j,i) = DefBuf4D(i,j,k,1)
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, RecCountX4D, RecDisplX4D, MPI_REAL8, &
             RecBuf1D, SendCountX4D, SendDisplX4D, MPI_REAL8, CommSliceX, ierr)
        
        do i=1,grd%im
           do iProc=0, NumProcI-1
              do j=1,SendCountX4D(iProc+1)/(grd%im*grd%km)
                 do k=1,grd%km
                    grd%chl(i, j + SendDisplX4D(iProc+1)/(grd%im*grd%km),k,1) = &
                         RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))
                 end do
              end do
           end do
        end do
        DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)
        
     else ! NumProcI .eq. 1
        call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
        
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
     end if


     ! ! ---
     ! ! y direction
     if(NumProcJ .gt. 1) then
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*localRow*GlobalCol))     
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
        
        call MPI_Alltoallv(SendBuf4D, SendCountY4D, SendDisplY4D, MPI_REAL8, &
             RecBuf1D, RecCountY4D, RecDisplY4D, MPI_REAL8, CommSliceY, ierr)
        
        do i=1,localRow
           do iProc=0, NumProcJ-1
              do j=1,RecCountY4D(iProc+1)/(localRow*grd%km)
                 do k=1,grd%km
                    DefBuf4D(i,j+RecDisplY4D(iProc+1)/(localRow*grd%km),k,1) = &
                         RecBuf1D(k + (j-1)*grd%km + (i-1)*RecCountY4D(iProc+1)/localRow + RecDisplY4D(iProc+1))
                 end do
              end do
           end do
        end do
        
        ! Apply recursive filter in y direction
        call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)
        
        ! ---
        ! Scale by the scaling factor
        do l=1,grd%nchl
           !$OMP PARALLEL  &
           !$OMP PRIVATE(k)
           !$OMP DO
           do k=1,grd%km
              DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scy(:,:) 
           enddo
           !$OMP END DO
           !$OMP END PARALLEL  
        enddo
        
        ! Reordering data to send back
        DEALLOCATE(SendBuf4D, RecBuf1D)
        ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
        ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
        
        do j=1,GlobalCol
           do i=1,localRow
              do k=1,grd%km
                 SendBuf4D(1,k,i,j) = DefBuf4D(i,j,k,1)
              end do
           end do
        end do
        
        call MPI_Alltoallv(SendBuf4D, RecCountY4D, RecDisplY4D, MPI_REAL8, &
             RecBuf1D, SendCountY4D, SendDisplY4D, MPI_REAL8, CommSliceY, ierr)
        
        do j=1,grd%jm
           do iProc=0, NumProcJ-1
              do i=1, SendCountY4D(iProc+1)/(grd%jm*grd%km)
                 do k=1,grd%km
                    grd%chl(i + SendDisplY4D(iProc+1)/(grd%jm*grd%km),j,k,1) = &
                         RecBuf1D(k + (i-1)*grd%km + (j-1)*SendCountY4D(iProc+1)/grd%jm + SendDisplY4D(iProc+1))
                 end do
              end do
           end do
        end do
        
     else ! NumProcJ .eq. 1
        ! Apply recursive filter in y direction
        call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
        
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
     end if

  endif
  
  ! ---
  ! y direction
  if(NumProcJ .gt. 1) then
     DEALLOCATE(SendBuf4D, RecBuf1D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%jm, grd%im))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*localRow*GlobalCol))
     
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,j,i) = grd%chl_ad(i,j,k,l)
              end do
           end do
        end do
     end do

     call MPI_Alltoallv(SendBuf4D, SendCountY4D, SendDisplY4D, MPI_REAL8, &
          RecBuf1D, RecCountY4D, RecDisplY4D, MPI_REAL8, CommSliceY, ierr)
     
     do i=1,localRow
        do iProc=0, NumProcJ-1
           do j=1,RecCountY4D(iProc+1)/(localRow*grd%km)
              do k=1,grd%km
                 DefBuf4D(i,j+RecDisplY4D(iProc+1)/(localRow*grd%km),k,1) = &
                      RecBuf1D(k + (j-1)*grd%km + (i-1)*RecCountY4D(iProc+1)/localRow + RecDisplY4D(iProc+1))
              end do
           end do
        end do
     end do
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scy(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL  
     enddo
     
     ! Apply recursive filter in y direction
     call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, DefBuf4D, grd%jnx, grd%jmx)
     
     ! Reordering data to send back
     DEALLOCATE(SendBuf4D, RecBuf1D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, localRow, GlobalCol))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
     
     do j=1,GlobalCol
        do i=1,localRow
           do k=1,grd%km
              SendBuf4D(1,k,i,j) = DefBuf4D(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, RecCountY4D, RecDisplY4D, MPI_REAL8, &
          RecBuf1D, SendCountY4D, SendDisplY4D, MPI_REAL8, CommSliceY, ierr)
     
     do j=1,grd%jm
        do iProc=0, NumProcJ-1
           do i=1, SendCountY4D(iProc+1)/(grd%jm*grd%km)
              do k=1,grd%km
                 grd%chl_ad(i + SendDisplY4D(iProc+1)/(grd%jm*grd%km),j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km + (j-1)*SendCountY4D(iProc+1)/grd%jm + SendDisplY4D(iProc+1))
              end do
           end do
        end do
     end do
     DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)
     
  else ! NumProcJ .eq. 1
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
     
     ! Apply recursive filter in y direction
     call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
  end if
  
  ! ---
  ! x direction
  if(NumProcI .gt. 1) then
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%im, grd%jm))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*localCol*GlobalRow))
     ALLOCATE( DefBuf4D(GlobalRow, localCol, grd%km, grd%nchl))
     
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,i,j) = grd%chl_ad(i,j,k,l)
              end do
           end do
        end do
     end do
  
     call MPI_Alltoallv(SendBuf4D, SendCountX4D, SendDisplX4D, MPI_REAL8, &
          RecBuf1D, RecCountX4D, RecDisplX4D, MPI_REAL8, CommSliceX, ierr)
     
     do j=1,localCol
        do iProc=0, NumProcI-1
           do i=1,RecCountX4D(iProc+1)/(localCol*grd%km)
              do k=1,grd%km
                 DefBuf4D(i + RecDisplX4D(iProc+1)/(localCol*grd%km),j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km + (j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))
              end do
           end do
        end do
     end do
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           DefBuf4D(:,:,k,l) = DefBuf4D(:,:,k,l) * grd%scx(:,:) 
        enddo
        !$OMP END DO
        !$OMP END PARALLEL  
     enddo
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBuf4D, grd%inx, grd%imx)
  
     ! Reordering data to send back
     DEALLOCATE(SendBuf4D, RecBuf1D)
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, localCol, GlobalRow))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*grd%jm*grd%im))
     
     do k=1,grd%km
        do j=1,localCol
           do i=1,GlobalRow
              SendBuf4D(1,k,j,i) = DefBuf4D(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, RecCountX4D, RecDisplX4D, MPI_REAL8, &
          RecBuf1D, SendCountX4D, SendDisplX4D, MPI_REAL8, CommSliceX, ierr)
     
     do i=1,grd%im
        do iProc=0, NumProcI-1
           do j=1,SendCountX4D(iProc+1)/(grd%im*grd%km)
              do k=1,grd%km
                 grd%chl_ad(i, j + SendDisplX4D(iProc+1)/(grd%im*grd%km),k,1) = &
                      RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))
              end do
           end do
        end do
     end do
     
     DEALLOCATE(SendBuf4D, RecBuf1D, DefBuf4D)
  
  else ! NumProcI .eq. 1
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
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
  end if
  
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
  
    

  !103 continue
  ! ---
  ! Vertical EOFs           
  call veof_ad
  
end subroutine parallel_ver_hor_ad
