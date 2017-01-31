subroutine onesided_obs_chl
  
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
  
  INTEGER(i4)   ::  i, j, l, kk, ierr
  REAL(r8)      :: FirstData, LastData
  INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset

  do kk = 1,chl%no
     
     if(chl%flc(kk).eq.1 )then
        
        
        i=chl%ib(kk)
        j=chl%jb(kk)
        
        chl%inc(kk) = 0.0
        
        if(i .lt. grd%im) then
          do l=1,grd%nchl
            chl%inc(kk) = chl%inc(kk) + (                        &
                  chl%pq1(kk) * grd%chl(i  ,j  ,1,l) +       &
                  chl%pq2(kk) * grd%chl(i+1,j  ,1,l) +       &
                  chl%pq3(kk) * grd%chl(i  ,j+1,1,l) +       &
                  chl%pq4(kk) * grd%chl(i+1,j+1,1,l) ) * chl%dzr(1,kk)
          enddo
        else
          call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChl, ierr )

          TargetOffset = j-1
          call MPI_Get (FirstData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)
          TargetOffset = j
          call MPI_Get (LastData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)

          call MPI_Win_unlock(ProcBottom, MpiWinChl, ierr)


          do l=1,grd%nchl
            chl%inc(kk) = chl%inc(kk) + (                        &
                  chl%pq1(kk) * grd%chl(i  ,j  ,1,l) +       &
                  chl%pq2(kk) * FirstData +       &
                  chl%pq3(kk) * grd%chl(i  ,j+1,1,l) +       &
                  chl%pq4(kk) * LastData ) * chl%dzr(1,kk)
          enddo

        endif

     endif
     
  enddo

end subroutine onesided_obs_chl

subroutine onesided_obs_chl_ad
  
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
  ! Apply observational operator for velocities from gliders  (adjoint)  !
  !                                                                      !
  ! Version 1: S.Dobricic 2007                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use obs_str
  use mpi_str

  implicit none
  
  INTEGER(i4)   ::  i, j, kk, l, ierr
  REAL(r8)      :: ToSum
  INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset
  
  do kk = 1,chl%no
     
     if(chl%flc(kk).eq.1)then
        
        obs%k = obs%k + 1
        
        i=chl%ib(kk)
        j=chl%jb(kk)
        
        if(i .lt. grd%im) then
          do l=1,grd%nchl
            grd%chl_ad(i  ,j  ,1,l) = grd%chl_ad(i  ,j  ,1,l) + chl%pq1(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            grd%chl_ad(i+1,j  ,1,l) = grd%chl_ad(i+1,j  ,1,l) + chl%pq2(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            grd%chl_ad(i  ,j+1,1,l) = grd%chl_ad(i  ,j+1,1,l) + chl%pq3(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            grd%chl_ad(i+1,j+1,1,l) = grd%chl_ad(i+1,j+1,1,l) + chl%pq4(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
          enddo
        else

          do l=1,grd%nchl
            grd%chl_ad(i  ,j  ,1,l) = grd%chl_ad(i  ,j  ,1,l) + chl%pq1(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            ! grd%chl_ad(i+1,j  ,1,l) = grd%chl_ad(i+1,j  ,1,l) + chl%pq2(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            grd%chl_ad(i  ,j+1,1,l) = grd%chl_ad(i  ,j+1,1,l) + chl%pq3(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
            ! grd%chl_ad(i+1,j+1,1,l) = grd%chl_ad(i+1,j+1,1,l) + chl%pq4(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
          enddo

          call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChlAd, ierr )

          TargetOffset = j-1
          ToSum = chl%pq2(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
          call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)
          TargetOffset = j
          ToSum = chl%pq4(kk) * chl%dzr(1,kk) * obs%gra(obs%k)
          call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)

          call MPI_Win_unlock(ProcBottom, MpiWinChlAd, ierr)

        endif


     endif   
  enddo
  
  call MPI_Barrier(MyCommWorld, ierr)
  
end subroutine onesided_obs_chl_ad
