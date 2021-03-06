subroutine obs_arg_ad
  
  
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
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply observational operator for ARGO floats (adjoint)               !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use obs_str
  use mpi_str
  use filenames
  use drv_str

  implicit none
  
  INTEGER(i4)   ::  i, j, k, kk
  REAL(r8)      :: ToSum
  INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset
  INTEGER       :: ierr, NData
  real(r8), pointer, dimension(:,:,:) :: MatrixToSum
  
  do kk = 1,arg%no
     
     if(arg%flc(kk).eq.1)then
        
        obs%k = obs%k + 1
        
        i=arg%ib(kk)
        j=arg%jb(kk)
        k=arg%kb(kk)
        
        if(i .lt. grd%im) then
           
           grd%chl_ad(i  ,j  ,k  ) = grd%chl_ad(i  ,j  ,k  ) + arg%pq1(kk) * obs%gra(obs%k)
           grd%chl_ad(i+1,j  ,k  ) = grd%chl_ad(i+1,j  ,k  ) + arg%pq2(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j+1,k  ) = grd%chl_ad(i  ,j+1,k  ) + arg%pq3(kk) * obs%gra(obs%k)
           grd%chl_ad(i+1,j+1,k  ) = grd%chl_ad(i+1,j+1,k  ) + arg%pq4(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j  ,k+1) = grd%chl_ad(i  ,j  ,k+1) + arg%pq5(kk) * obs%gra(obs%k)
           grd%chl_ad(i+1,j  ,k+1) = grd%chl_ad(i+1,j  ,k+1) + arg%pq6(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j+1,k+1) = grd%chl_ad(i  ,j+1,k+1) + arg%pq7(kk) * obs%gra(obs%k)
           grd%chl_ad(i+1,j+1,k+1) = grd%chl_ad(i+1,j+1,k+1) + arg%pq8(kk) * obs%gra(obs%k)
           
        else

           ALLOCATE(MatrixToSum(NextLocalRow,grd%jm,2))
           MatrixToSum(:,:,:) = dble(0)
           
           grd%chl_ad(i  ,j  ,k  ) = grd%chl_ad(i  ,j  ,k  ) + arg%pq1(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j+1,k  ) = grd%chl_ad(i  ,j+1,k  ) + arg%pq3(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j  ,k+1) = grd%chl_ad(i  ,j  ,k+1) + arg%pq5(kk) * obs%gra(obs%k)
           grd%chl_ad(i  ,j+1,k+1) = grd%chl_ad(i  ,j+1,k+1) + arg%pq7(kk) * obs%gra(obs%k)
           
           MatrixToSum(1,j  ,1) = arg%pq2(kk) * obs%gra(obs%k)
           MatrixToSum(1,j+1,1) = arg%pq4(kk) * obs%gra(obs%k)
           MatrixToSum(1,j  ,2) = arg%pq6(kk) * obs%gra(obs%k)
           MatrixToSum(1,j+1,2) = arg%pq8(kk) * obs%gra(obs%k)
           
           call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChlAd, ierr )
           NData = NextLocalRow*grd%jm*2
           TargetOffset = (k-1)*grd%jm*NextLocalRow
           call MPI_Accumulate (MatrixToSum, NData, MPI_REAL8, ProcBottom, TargetOffset, NData, MPI_REAL8, MPI_SUM, MpiWinChlAd, ierr)
           
           call MPI_Win_unlock(ProcBottom, MpiWinChlAd, ierr)
           DEALLOCATE(MatrixToSum)

        endif
        
     endif
     
  enddo
  
end subroutine obs_arg_ad
