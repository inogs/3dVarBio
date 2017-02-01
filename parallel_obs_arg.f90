subroutine parallel_obs_arg


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
! Apply observational operator for ARGO floats                         !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str
 use mpi_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk
 REAL(r8)      :: FirstData, SecData, ThirdData, LastData
 INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset
 INTEGER       :: ierr

 do kk = 1,arg%no

  if(arg%flc(kk).eq.1)then

    i=arg%ib(kk)
    j=arg%jb(kk)
    k=arg%kb(kk)

    if(i .lt. grd%im) then
      arg%inc(kk) = arg%pq1(kk) * grd%chl(i  ,j  ,k  ,1) +       &
                    arg%pq2(kk) * grd%chl(i+1,j  ,k  ,1) +       &
                    arg%pq3(kk) * grd%chl(i  ,j+1,k  ,1) +       &
                    arg%pq4(kk) * grd%chl(i+1,j+1,k  ,1) +       &
                    arg%pq5(kk) * grd%chl(i  ,j  ,k+1,1) +       &
                    arg%pq6(kk) * grd%chl(i+1,j  ,k+1,1) +       &
                    arg%pq7(kk) * grd%chl(i  ,j+1,k+1,1) +       &
                    arg%pq8(kk) * grd%chl(i+1,j+1,k+1,1)  

    else
      call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChl, ierr )

      TargetOffset = j-1 + (k-1)*grd%jm*NextLocalRow
      call MPI_Get (FirstData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)
      TargetOffset = j + (k-1)*grd%jm*NextLocalRow
      call MPI_Get ( SecData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)
      TargetOffset = j-1 + k*grd%jm*NextLocalRow
      call MPI_Get (ThirdData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)
      TargetOffset = j + k*grd%jm*NextLocalRow
      call MPI_Get (LastData, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MpiWinChl, ierr)
      

      call MPI_Win_unlock(ProcBottom, MpiWinChl, ierr)

      arg%inc(kk) = arg%pq1(kk) * grd%chl(i  ,j  ,k  ,1) +       &
                    arg%pq2(kk) * FirstData +       &
                    arg%pq3(kk) * grd%chl(i  ,j+1,k  ,1) +       &
                    arg%pq4(kk) * SecData +       &
                    arg%pq5(kk) * grd%chl(i  ,j  ,k+1,1) +       &
                    arg%pq6(kk) * ThirdData +       &
                    arg%pq7(kk) * grd%chl(i  ,j+1,k+1,1) +       &
                    arg%pq8(kk) * LastData  

    endif

  endif

 enddo

end subroutine parallel_obs_arg

subroutine parallel_obs_arg_ad


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
! Apply observational operator for ARGO floats (adjoint)  
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str
 use mpi_str

 implicit none

 INTEGER(i4)   ::  i, j, k, kk
 REAL(r8)      :: ToSum
 INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset
 INTEGER       :: ierr

 do kk = 1,arg%no

  if(arg%flc(kk).eq.1)then

    obs%k = obs%k + 1

    i=arg%ib(kk)
    j=arg%jb(kk)
    k=arg%kb(kk)

    if(i .lt. grd%im) then
  
      grd%chl_ad(i  ,j  ,k  ,1) = grd%chl_ad(i  ,j  ,k  ,1) + arg%pq1(kk) * obs%gra(obs%k)
      grd%chl_ad(i+1,j  ,k  ,1) = grd%chl_ad(i+1,j  ,k  ,1) + arg%pq2(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j+1,k  ,1) = grd%chl_ad(i  ,j+1,k  ,1) + arg%pq3(kk) * obs%gra(obs%k)
      grd%chl_ad(i+1,j+1,k  ,1) = grd%chl_ad(i+1,j+1,k  ,1) + arg%pq4(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j  ,k+1,1) = grd%chl_ad(i  ,j  ,k+1,1) + arg%pq5(kk) * obs%gra(obs%k)
      grd%chl_ad(i+1,j  ,k+1,1) = grd%chl_ad(i+1,j  ,k+1,1) + arg%pq6(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j+1,k+1,1) = grd%chl_ad(i  ,j+1,k+1,1) + arg%pq7(kk) * obs%gra(obs%k)
      grd%chl_ad(i+1,j+1,k+1,1) = grd%chl_ad(i+1,j+1,k+1,1) + arg%pq8(kk) * obs%gra(obs%k)

    else
      
      grd%chl_ad(i  ,j  ,k  ,1) = grd%chl_ad(i  ,j  ,k  ,1) + arg%pq1(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j+1,k  ,1) = grd%chl_ad(i  ,j+1,k  ,1) + arg%pq3(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j  ,k+1,1) = grd%chl_ad(i  ,j  ,k+1,1) + arg%pq5(kk) * obs%gra(obs%k)
      grd%chl_ad(i  ,j+1,k+1,1) = grd%chl_ad(i  ,j+1,k+1,1) + arg%pq7(kk) * obs%gra(obs%k)

      call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChlAd, ierr )

      TargetOffset = j-1 + (k-1)*grd%jm*NextLocalRow
      ToSum = arg%pq2(kk) * obs%gra(obs%k)
      call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)
      
      TargetOffset = j + (k-1)*grd%jm*NextLocalRow
      ToSum = arg%pq4(kk) * obs%gra(obs%k)
      call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)
      
      TargetOffset = j-1 + k*grd%jm*NextLocalRow
      ToSum = arg%pq6(kk) * obs%gra(obs%k)
      call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)
      
      TargetOffset = j-1 + k*grd%jm*NextLocalRow
      ToSum = arg%pq8(kk) * obs%gra(obs%k)
      call MPI_Accumulate (ToSum, 1, MPI_REAL8, ProcBottom, TargetOffset, 1, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)

      call MPI_Win_unlock(ProcBottom, MpiWinChlAd, ierr)

    endif

  endif

 enddo

 call MPI_Barrier(MyCommWorld, ierr)
 
end subroutine parallel_obs_arg_ad
