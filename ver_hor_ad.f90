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
  type(DoubleGrid), allocatable, dimension(:,:,:,:) :: SendBuf4D
  type(DoubleGrid), allocatable, dimension(:)       :: RecBuf1D
  REAL(r8), allocatable, dimension(:,:,:,:) :: DefBufChl, DefBufChlAd  
  
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
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scy(:,:,k)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  enddo
  
  ! Apply recursive filter in y direction
  call rcfl_y_ad( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl_ad, grd%jnx, grd%jmx)
  
  ! ---
  ! x direction
  if(NumProcI .gt. 1) then
     ALLOCATE(SendBuf4D(grd%nchl, grd%km, grd%im, grd%jm))
     ALLOCATE( RecBuf1D(grd%nchl*grd%km*localCol*GlobalRow))
     ALLOCATE(DefBufChl(GlobalRow, localCol, grd%km, grd%nchl))
     ALLOCATE(DefBufChlAd(GlobalRow, localCol, grd%km, grd%nchl))
     
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,i,j)%chl = grd%chl(i,j,k,l)
              end do
           end do
        end do
     end do
     do l=1,grd%nchl
        do k=1,grd%km
           do j=1,grd%jm
              do i=1,grd%im
                 SendBuf4D(l,k,i,j)%chl_ad = grd%chl_ad(i,j,k,l)
              end do
           end do
        end do
     end do

     call MPI_Alltoallv(SendBuf4D, SendCountX4D, SendDisplX4D, MyPair, &
          RecBuf1D, RecCountX4D, RecDisplX4D, MyPair, Var3DCommunicator, ierr)
     
     SurfaceIndex = localCol*grd%km
     do j=1,localCol
        do iProc=0, NumProcI-1
           do i=1,RecCountX4D(iProc+1)/SurfaceIndex
              do k=1,grd%km
                 DefBufChl(i + RecDisplX4D(iProc+1)/SurfaceIndex,j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km + (j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))%chl
              end do
           end do
        end do
     end do
     do j=1,localCol
        do iProc=0, NumProcI-1
           do i=1,RecCountX4D(iProc+1)/SurfaceIndex
              do k=1,grd%km
                 DefBufChlAd(i + RecDisplX4D(iProc+1)/SurfaceIndex,j,k,1) = &
                      RecBuf1D(k + (i-1)*grd%km + (j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))%chl_ad
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
           DefBufChlAd(:,:,k,l) = DefBufChlAd(:,:,k,l) * grd%scx(:,:,k)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBufChlAd, grd%inx, grd%imx)
     
  else ! NumProcI .eq. 1
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl_ad(:,:,k,l) = grd%chl_ad(:,:,k,l) * grd%scx(:,:,k)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
     
     call rcfl_x_ad( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl_ad, grd%inx, grd%imx)
  end if
  
  
  ! ---
  ! x direction
  if(NumProcI .gt. 1) then
     
     call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, DefBufChl, grd%inx, grd%imx)
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           DefBufChl(:,:,k,l) = DefBufChl(:,:,k,l) * grd%scx(:,:,k)
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
              SendBuf4D(1,k,j,i)%chl = DefBufChl(i,j,k,1)
           end do
        end do
     end do
     do k=1,grd%km
        do j=1,localCol
           do i=1,GlobalRow
              SendBuf4D(1,k,j,i)%chl_ad = DefBufChlAd(i,j,k,1)
           end do
        end do
     end do
     
     call MPI_Alltoallv(SendBuf4D, RecCountX4D, RecDisplX4D, MyPair, &
          RecBuf1D, SendCountX4D, SendDisplX4D, MyPair, Var3DCommunicator, ierr)
     
     SurfaceIndex = grd%im*grd%km
     do i=1,grd%im
        do iProc=0, NumProcI-1
           do j=1,SendCountX4D(iProc+1)/SurfaceIndex
              do k=1,grd%km
                 grd%chl(i, j + SendDisplX4D(iProc+1)/SurfaceIndex,k,1) = &
                      RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))%chl
              end do
           end do
        end do
     end do
     do i=1,grd%im
        do iProc=0, NumProcI-1
           do j=1,SendCountX4D(iProc+1)/SurfaceIndex
              do k=1,grd%km
                 grd%chl_ad(i, j + SendDisplX4D(iProc+1)/SurfaceIndex,k,1) = &
                      RecBuf1D(k + (j-1)*grd%km +(i-1)*SendCountX4D(iProc+1)/grd%im + SendDisplX4D(iProc+1))%chl_ad
              end do
           end do
        end do
     end do

     DEALLOCATE(SendBuf4D, RecBuf1D, DefBufChl, DefBufChlAd)
     
  else ! NumProcI .eq. 1
     call rcfl_x( GlobalRow, localCol, grd%km*grd%nchl, grd%imax, grd%aex, grd%bex, grd%chl, grd%inx, grd%imx)
     
     ! ---
     ! Scale by the scaling factor
     do l=1,grd%nchl
        !$OMP PARALLEL  &
        !$OMP PRIVATE(k)
        !$OMP DO
        do k=1,grd%km
           grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scx(:,:,k)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
     enddo
  end if
  
  
  ! ! ---
  ! ! y direction
  ! Apply recursive filter in y direction
  call rcfl_y( localRow, GlobalCol, grd%km*grd%nchl, grd%jmax, grd%aey, grd%bey, grd%chl, grd%jnx, grd%jmx)
  
  ! ---
  ! Scale by the scaling factor
  do l=1,grd%nchl
     !$OMP PARALLEL  &
     !$OMP PRIVATE(k)
     !$OMP DO
     do k=1,grd%km
        grd%chl(:,:,k,l) = grd%chl(:,:,k,l) * grd%scy(:,:,k)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  enddo
  
  
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
  
end subroutine parallel_ver_hor_ad
