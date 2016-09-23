subroutine parallel_obs_chl
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
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
  ! Apply observational operator for chlorophyll                         !
  !                                                                      !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use obs_str
  use mpi_str
  use mpi
  
  implicit none
  
  INTEGER(i4)   ::  i, j, l, kk
  ! REAL(r8), POINTER    ::  ChlExtended(:,:,:)
  ! REAL(r8), POINTER    ::  SendLeft(:), RecRight(:), SendTop(:), RecBottom(:)
  INTEGER   :: ReqLeft, ReqRight, ReqTop, ReqBottom, ierr
  INTEGER   :: StatRight(MPI_STATUS_SIZE), StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag
  
  ! ALLOCATE(ChlExtended(grd%im+1, grd%jm+1, grd%nchl))
  ! ALLOCATE(SendLeft(grd%im), RecRight(grd%im))
  ! ALLOCATE(SendTop(grd%im), RecBottom(grd%im))

  ! Filling array to send
  do i=1,grd%im
     SendLeft(i) = grd%chl(i,1,1,1)
  end do
  do j=1,grd%jm
     SendTop(j)  = grd%chl(1,j,1,1)
  end do
  
  MyTag = 42
  RecRight(:)  = 0
  RecBottom(:) = 0
  
  call MPI_Isend(SendLeft, grd%im, MPI_REAL8, ProcLeft, MyTag, &
       MPI_COMM_WORLD, ReqLeft, ierr)
  call MPI_Irecv(RecRight, grd%im, MPI_REAL8, ProcRight, MyTag, &
       MPI_COMM_WORLD, ReqRight, ierr)

  call MPI_Isend(SendTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       MPI_COMM_WORLD, ReqTop, ierr)
  call MPI_Irecv(RecBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       MPI_COMM_WORLD, ReqBottom, ierr)
  
  do j=1,grd%jm
     do i=1,grd%im
        ChlExtended(i,j,1) = grd%chl(i,j,1,1)
     end do
  end do
  
  call MPI_Wait(ReqRight, StatRight, ierr)
  call MPI_Wait(ReqBottom, StatBottom, ierr)
  do i=1,grd%im
     ChlExtended(i,grd%jm+1,1) = RecRight(i)
  end do
  do j=1,grd%jm
     ChlExtended(grd%im+1,j,1) = RecBottom(j)
  end do
  
  do kk = 1,chl%no
     
     if(chl%flc(kk).eq.1 )then
        
        
        i=chl%ib(kk)
        j=chl%jb(kk)
        
        chl%inc(kk) = 0.0
        
        do l=1,grd%nchl
           chl%inc(kk) = chl%inc(kk) + (                        &
                chl%pq1(kk) * ChlExtended(i  ,j  ,1) +       &
                chl%pq2(kk) * ChlExtended(i+1,j  ,1) +       &
                chl%pq3(kk) * ChlExtended(i  ,j+1,1) +       &
                chl%pq4(kk) * ChlExtended(i+1,j+1,1) ) * chl%dzr(1,kk)
        enddo
        
     endif
     
  enddo

  ! DEALLOCATE(ChlExtended)
  ! DEALLOCATE(SendLeft, RecRight)
  ! DEALLOCATE(SendTop, RecBottom)
  
end subroutine parallel_obs_chl

subroutine parallel_obs_chl_ad  
  
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
  use mpi
  
  implicit none
  
  INTEGER(i4)   ::  i, j, kk, l
  INTEGER   :: ReqRight, ReqBottom, ReqLeft, ReqTop, ierr
  INTEGER   :: StatLeft(MPI_STATUS_SIZE), StatRight(MPI_STATUS_SIZE)
  INTEGER   :: StatTop(MPI_STATUS_SIZE), StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag
  
  ! Filling array to send
  do i=1,grd%im
     SendLeft(i) = grd%chl_ad(i,1,1,1)
  end do
  do j=1,grd%jm
     SendTop(j)  = grd%chl_ad(1,j,1,1)
  end do

  MyTag = 42
  RecRight(:)  = 0
  RecBottom(:) = 0
  
  call MPI_Isend(SendLeft, grd%im, MPI_REAL8, ProcLeft, MyTag, &
       MPI_COMM_WORLD, ReqLeft, ierr)
  call MPI_Irecv(RecRight, grd%im, MPI_REAL8, ProcRight, MyTag, &
       MPI_COMM_WORLD, ReqRight, ierr)

  call MPI_Isend(SendTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       MPI_COMM_WORLD, ReqTop, ierr)
  call MPI_Irecv(RecBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       MPI_COMM_WORLD, ReqBottom, ierr)
  
  do j=1,grd%jm
     do i=1,grd%im
        ChlExtended(i,j,1) = grd%chl_ad(i,j,1,1)
     end do
  end do
  
  call MPI_Wait(ReqRight, StatRight, ierr)
  call MPI_Wait(ReqBottom, StatBottom, ierr)
  do i=1,grd%im
     ChlExtended(i,grd%jm+1,1) = RecRight(i)
  end do
  do j=1,grd%jm
     ChlExtended(grd%im+1,j,1) = RecBottom(j)
  end do
  
  do kk = 1,chl%no
     
     if(chl%flc(kk).eq.1)then
        
        obs%k = obs%k + 1
        
        i=chl%ib(kk)
        j=chl%jb(kk)
        
        do l=1,grd%nchl
           ChlExtended(i  ,j  ,1) = ChlExtended(i  ,j  ,1) + chl%pq1(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
           ChlExtended(i+1,j  ,1) = ChlExtended(i+1,j  ,1) + chl%pq2(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
           ChlExtended(i  ,j+1,1) = ChlExtended(i  ,j+1,1) + chl%pq3(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
           ChlExtended(i+1,j+1,1) = ChlExtended(i+1,j+1,1) + chl%pq4(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
        enddo
     endif   
  enddo

  do i=1,grd%im
     SendRight(i)  = ChlExtended(i,grd%jm+1,1)
  end do
  do j=1,grd%jm
     SendBottom(j) = ChlExtended(grd%im+1,j,1)
  end do
  
  RecLeft(:) = SendLeft(:)
  RecTop(:)  = SendTop(:)
  
  call MPI_Isend(SendRight, grd%im, MPI_REAL8, ProcRight, MyTag, &
       MPI_COMM_WORLD, ReqRight, ierr)
  call MPI_Irecv(RecLeft, grd%im, MPI_REAL8, ProcLeft, MyTag, &
       MPI_COMM_WORLD, ReqLeft, ierr)

  call MPI_Isend(SendBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       MPI_COMM_WORLD, ReqBottom, ierr)
  call MPI_Irecv(RecTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       MPI_COMM_WORLD, ReqTop, ierr)
  
  do j=1,grd%jm
     do i=1,grd%im
        grd%chl_ad(i,j,1,1) = ChlExtended(i,j,1)
     end do
  end do

  call MPI_Wait(ReqLeft, StatLeft, ierr)  
  do i=1,grd%im
     grd%chl_ad(i,1,1,1) = grd%chl_ad(i,1,1,1) + RecLeft(i) - SendLeft(i)
  end do

  call MPI_Wait(ReqTop, StatTop, ierr)  
  do j=1,grd%jm
     grd%chl_ad(1,j,1,1) = grd%chl_ad(1,j,1,1) + RecTop(j) - SendTop(j)
  end do

end subroutine parallel_obs_chl_ad
