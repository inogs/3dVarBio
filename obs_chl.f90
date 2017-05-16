subroutine obs_chl
  
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
  
  implicit none
  
  INTEGER(i4)   ::  i, j, k, l, kk
  INTEGER   :: ReqTop, ReqBottom, ierr
  INTEGER   :: StatBottom(MPI_STATUS_SIZE)
  INTEGER   :: MyTag
  
  grd%chl(:,:,:) = 0.0

  do l = 1,grd%nchl
    do k = 1,grd%km
      do j = 1,grd%jm
        do i = 1,grd%im
          grd%chl(i,j,k) = grd%chl(i,j,k) + bio%phy(i,j,k,l,1)
        enddo
      enddo
    enddo
  enddo


  ! Filling array to send
  do j=1,grd%jm
    SendTop(j)  = grd%chl(1,j,1)
  end do
  
  MyTag = 42
  RecBottom(:) = 0
  
  call MPI_Isend(SendTop, grd%jm, MPI_REAL8, ProcTop, MyTag, &
       Var3DCommunicator, ReqTop, ierr)
  call MPI_Irecv(RecBottom, grd%jm, MPI_REAL8, ProcBottom, MyTag, &
       Var3DCommunicator, ReqBottom, ierr)
  
  do j=1,grd%jm
     do i=1,grd%im
        ChlExtended(i,j) = grd%chl(i,j,1)
     end do
  end do
  
  call MPI_Wait(ReqBottom, StatBottom, ierr)
  do j=1,grd%jm
     ChlExtended(grd%im+1,j) = RecBottom(j)
  end do
  
  do kk = 1,chl%no
     
     if(chl%flc(kk).eq.1 )then
        
        
        i=chl%ib(kk)
        j=chl%jb(kk)
        
        chl%inc(kk) = 0.0
        
        chl%inc(kk) = chl%inc(kk) + (                        &
            chl%pq1(kk) * ChlExtended(i  ,j  ) +       &
            chl%pq2(kk) * ChlExtended(i+1,j  ) +       &
            chl%pq3(kk) * ChlExtended(i  ,j+1) +       &
            chl%pq4(kk) * ChlExtended(i+1,j+1) ) * chl%dzr(1,kk)
        
     endif
     
  enddo

end subroutine obs_chl

