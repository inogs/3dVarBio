subroutine obs_arg
  
  
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
  REAL(r8)      :: FirstData, SecData, ThirdData, LastData, Test
  INTEGER(kind=MPI_ADDRESS_KIND) :: TargetOffset
  INTEGER       :: ierr, NData
  real(r8), pointer, dimension(:,:,:)   :: GetData
  
  do kk = 1,arg%no
     
    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.0)then
        
      i=arg%ib(kk)
      j=arg%jb(kk)
      k=arg%kb(kk)
        
      if(i .lt. grd%im) then
        arg%inc(kk) = arg%pq1(kk) * grd%chl(i  ,j  ,k) +       &
          arg%pq2(kk) * grd%chl(i+1,j  ,k  ) +       &
          arg%pq3(kk) * grd%chl(i  ,j+1,k  ) +       &
          arg%pq4(kk) * grd%chl(i+1,j+1,k  ) +       &
          arg%pq5(kk) * grd%chl(i  ,j  ,k+1) +       &
          arg%pq6(kk) * grd%chl(i+1,j  ,k+1) +       &
          arg%pq7(kk) * grd%chl(i  ,j+1,k+1) +       &
          arg%pq8(kk) * grd%chl(i+1,j+1,k+1)  
          
      else
        ALLOCATE(GetData(NextLocalRow,grd%jm,2))
          
        NData = NextLocalRow*grd%jm*2
        call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinChl, ierr )
        TargetOffset = (k-1)*grd%jm*NextLocalRow
        call MPI_Get (GetData, NData, MPI_REAL8, ProcBottom, TargetOffset, NData, MPI_REAL8, MpiWinChl, ierr)
        call MPI_Win_unlock(ProcBottom, MpiWinChl, ierr)
        
        arg%inc(kk) = arg%pq1(kk) * grd%chl(i  ,j  ,k) +       &
          arg%pq2(kk) * GetData(1  ,j  ,1  ) +       &
          arg%pq3(kk) * grd%chl(i  ,j+1,k  ) +       &
          arg%pq4(kk) * GetData(1  ,j+1,1  ) +       &
          arg%pq5(kk) * grd%chl(i  ,j  ,k+1) +       &
          arg%pq6(kk) * GetData(1  ,j  ,2  ) +       &
          arg%pq7(kk) * grd%chl(i  ,j+1,k+1) +       &
          arg%pq8(kk) * GetData(1  ,j+1,2  )  
          
          
          DEALLOCATE(GetData)
      endif

    else if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.1) then
        
      i=arg%ib(kk)
      j=arg%jb(kk)
      k=arg%kb(kk)
        
      if(i .lt. grd%im) then
        arg%inc(kk) = arg%pq1(kk) * grd%n3n(i  ,j  ,k) +       &
          arg%pq2(kk) * grd%n3n(i+1,j  ,k  ) +       &
          arg%pq3(kk) * grd%n3n(i  ,j+1,k  ) +       &
          arg%pq4(kk) * grd%n3n(i+1,j+1,k  ) +       &
          arg%pq5(kk) * grd%n3n(i  ,j  ,k+1) +       &
          arg%pq6(kk) * grd%n3n(i+1,j  ,k+1) +       &
          arg%pq7(kk) * grd%n3n(i  ,j+1,k+1) +       &
          arg%pq8(kk) * grd%n3n(i+1,j+1,k+1)  
          
      else
        ALLOCATE(GetData(NextLocalRow,grd%jm,2))
          
        NData = NextLocalRow*grd%jm*2
        call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinN3n, ierr )
        TargetOffset = (k-1)*grd%jm*NextLocalRow
        call MPI_Get (GetData, NData, MPI_REAL8, ProcBottom, TargetOffset, NData, MPI_REAL8, MpiWinN3n, ierr)
        call MPI_Win_unlock(ProcBottom, MpiWinN3n, ierr)
        
        arg%inc(kk) = arg%pq1(kk) * grd%n3n(i  ,j  ,k) +       &
          arg%pq2(kk) * GetData(1  ,j  ,1  ) +       &
          arg%pq3(kk) * grd%n3n(i  ,j+1,k  ) +       &
          arg%pq4(kk) * GetData(1  ,j+1,1  ) +       &
          arg%pq5(kk) * grd%n3n(i  ,j  ,k+1) +       &
          arg%pq6(kk) * GetData(1  ,j  ,2  ) +       &
          arg%pq7(kk) * grd%n3n(i  ,j+1,k+1) +       &
          arg%pq8(kk) * GetData(1  ,j+1,2  )  
          
          
          DEALLOCATE(GetData)
      endif
        
    else if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.2) then
        
      i=arg%ib(kk)
      j=arg%jb(kk)
      k=arg%kb(kk)
        
      if(i .lt. grd%im) then
        arg%inc(kk) = arg%pq1(kk) * grd%o2o(i  ,j  ,k) +       &
          arg%pq2(kk) * grd%o2o(i+1,j  ,k  ) +       &
          arg%pq3(kk) * grd%o2o(i  ,j+1,k  ) +       &
          arg%pq4(kk) * grd%o2o(i+1,j+1,k  ) +       &
          arg%pq5(kk) * grd%o2o(i  ,j  ,k+1) +       &
          arg%pq6(kk) * grd%o2o(i+1,j  ,k+1) +       &
          arg%pq7(kk) * grd%o2o(i  ,j+1,k+1) +       &
          arg%pq8(kk) * grd%o2o(i+1,j+1,k+1)  
          
      else
        ALLOCATE(GetData(NextLocalRow,grd%jm,2))
          
        NData = NextLocalRow*grd%jm*2
        call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, ProcBottom, 0, MpiWinO2o, ierr )
        TargetOffset = (k-1)*grd%jm*NextLocalRow
        call MPI_Get (GetData, NData, MPI_REAL8, ProcBottom, TargetOffset, NData, MPI_REAL8, MpiWinO2o, ierr)
        call MPI_Win_unlock(ProcBottom, MpiWinO2o, ierr)
        
        arg%inc(kk) = arg%pq1(kk) * grd%o2o(i  ,j  ,k) +       &
          arg%pq2(kk) * GetData(1  ,j  ,1  ) +       &
          arg%pq3(kk) * grd%o2o(i  ,j+1,k  ) +       &
          arg%pq4(kk) * GetData(1  ,j+1,1  ) +       &
          arg%pq5(kk) * grd%o2o(i  ,j  ,k+1) +       &
          arg%pq6(kk) * GetData(1  ,j  ,2  ) +       &
          arg%pq7(kk) * grd%o2o(i  ,j+1,k+1) +       &
          arg%pq8(kk) * GetData(1  ,j+1,2  )  
          
          
          DEALLOCATE(GetData)
      endif
        
    endif
     
  enddo

end subroutine obs_arg
