subroutine obs_sat_ad

  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply observational operator for velocities from gliders  (adjoint)  !
  !                                                                      !
  ! Version 1: S.Dobricic 2007                                           !
  !-----------------------------------------------------------------------


  use set_knd
  use grd_str
  use obs_str
  use mpi_str
  use bio_str

  implicit none

  INTEGER(i4)   ::  i, j, k, kk, l
  INTEGER   :: ReqBottom, ReqTop, ierr
  INTEGER   :: StatTop(MPI_STATUS_SIZE), StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag

  ! Filling array to send
  do j=1,grd%jm
     SendTop(j)  = grd%chl_ad(1,j,1)
  end do

  MyTag = 42
  RecBottom(:) = 0

  call MPI_Isend(SendTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       Var3DCommunicator, ReqTop, ierr)
  call MPI_Irecv(RecBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       Var3DCommunicator, ReqBottom, ierr)

  do j=1,grd%jm
     do i=1,grd%im
        ChlExtended(i,j) = grd%chl_ad(i,j,1)
     end do
  end do

  call MPI_Wait(ReqBottom, StatBottom, ierr)
  do j=1,grd%jm
     ChlExtended(grd%im+1,j) = RecBottom(j)
  end do

  do kk = 1,sat%no

     if(sat%flc(kk).eq.1)then

        obs%k = obs%k + 1

        i=sat%ib(kk)
        j=sat%jb(kk)

        ChlExtended(i  ,j  ) = ChlExtended(i  ,j  ) + sat%pq1(kk) * sat%dzr(1,kk) * obs%gra(obs%k)
        ChlExtended(i+1,j  ) = ChlExtended(i+1,j  ) + sat%pq2(kk) * sat%dzr(1,kk) * obs%gra(obs%k)
        ChlExtended(i  ,j+1) = ChlExtended(i  ,j+1) + sat%pq3(kk) * sat%dzr(1,kk) * obs%gra(obs%k)
        ChlExtended(i+1,j+1) = ChlExtended(i+1,j+1) + sat%pq4(kk) * sat%dzr(1,kk) * obs%gra(obs%k)
     endif
  enddo

  do j=1,grd%jm
     SendBottom(j) = ChlExtended(grd%im+1,j)
  end do

  RecTop(:)  = SendTop(:)

  call MPI_Isend(SendBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       Var3DCommunicator, ReqBottom, ierr)
  call MPI_Irecv(RecTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       Var3DCommunicator, ReqTop, ierr)

  do j=1,grd%jm
     do i=1,grd%im
        grd%chl_ad(i,j,1) = ChlExtended(i,j)
     end do
  end do

  call MPI_Wait(ReqTop, StatTop, ierr)
  do j=1,grd%jm
     grd%chl_ad(1,j,1) = grd%chl_ad(1,j,1) + RecTop(j) - SendTop(j)
  end do

  bio%phy_ad(:,:,:,:,:) = 0.0

  do l = 1,grd%nchl
    do k = 1,grd%km
      do j = 1,grd%jm
        do i = 1,grd%im
          bio%phy_ad(i,j,k,l,1) = bio%phy_ad(i,j,k,l,1) + grd%chl_ad(i,j,k)
        enddo
      enddo
    enddo
  enddo


end subroutine obs_sat_ad
