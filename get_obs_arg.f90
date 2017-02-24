subroutine get_obs_arg
  
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
  ! Load ARGO observations                                                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use mpi_str
  
  implicit none
  
  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, kk, i
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: TmpFlc, TmpPar, TmpLon, TmpLat
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: TmpDpt, TmpTim, TmpRes, TmpErr, TmpIno 
  INTEGER(i4)   :: GlobalArgNum, Counter, ierr
  character(len=1024) :: filename
  
  arg%no  = 0
  arg%nc  = 0
  Counter = 0
  
  
  ! ---
  ! Allocate memory for observations
  if(MyRank .eq. 0) then
    ! open(511,file='arg_datnew.dat',form='formatted')
    open(511,file='arg_mis.dat')
    read(511,'(I4)') GlobalArgNum
    write(drv%dia,*)'Number of ARGO observations: ', GlobalArgNum
  endif

  call MPI_Bcast(GlobalArgNum, 1, MPI_INT, 0, MyCommWorld, ierr)

  if(GlobalArgNum .eq. 0)then
    if(MyRank .eq. o) &
      close(511)
    return
  endif

  ALLOCATE( TmpFlc(GlobalArgNum), TmpPar(GlobalArgNum))
  ALLOCATE( TmpLon(GlobalArgNum), TmpLat(GlobalArgNum))
  ALLOCATE( TmpDpt(GlobalArgNum), TmpTim(GlobalArgNum))
  ALLOCATE( TmpRes(GlobalArgNum), TmpErr(GlobalArgNum))
  ALLOCATE( TmpIno(GlobalArgNum))

  if(MyRank .eq. 0) then
    ! each process reads all the argo observations
    do k=1,GlobalArgNum
      read (511,'(I5,I5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,I8)') &
      ! read (511,*) &
            TmpFlc(k), TmpPar(k), &
            TmpLon(k), TmpLat(k), &
            TmpDpt(k), TmpTim(k), &
            TmpRes(k), TmpErr(k), TmpIno(k)
    end do
    close (511)
  endif

  call MPI_Bcast(TmpFlc, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpPar, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpLon, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpLat, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpDpt, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpTim, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpRes, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpErr, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)
  call MPI_Bcast(TmpIno, GlobalArgNum, MPI_REAL8, 0, MyCommWorld, ierr)

  ! Counting the number of observations that falls in the domain
  do k=1,GlobalArgNum
    if( TmpLon(k) .ge. grd%lon(1,1) .and. TmpLon(k) .lt. grd%NextLongitude .and. &
        TmpLat(k) .ge. grd%lat(1,1) .and. TmpLat(k) .lt. grd%lat(grd%im,grd%jm) ) then
      Counter = Counter + 1
    endif
  enddo

  if(drv%Verbose .eq. 1) &
       print*, "MyRank", MyRank, "has",Counter,"ARGO observations"

  arg%no  = Counter

  ALLOCATE ( arg%ino(arg%no), arg%flg(arg%no), arg%flc(arg%no), arg%par(arg%no))
  ALLOCATE ( arg%lon(arg%no), arg%lat(arg%no), arg%dpt(arg%no), arg%tim(arg%no))
  ALLOCATE ( arg%inc(arg%no))
  ALLOCATE ( arg%err(arg%no))
  ALLOCATE ( arg%res(arg%no))
  ALLOCATE ( arg%ib(arg%no), arg%jb(arg%no), arg%kb(arg%no))
  ALLOCATE ( arg%pb(arg%no), arg%qb(arg%no), arg%rb(arg%no))
  ALLOCATE ( arg%pq1(arg%no), arg%pq2(arg%no), arg%pq3(arg%no), arg%pq4(arg%no))
  ALLOCATE ( arg%pq5(arg%no), arg%pq6(arg%no), arg%pq7(arg%no), arg%pq8(arg%no))  

  Counter = 0
  do k=1,GlobalArgNum
    if( TmpLon(k) .ge. grd%lon(1,1) .and. TmpLon(k) .lt. grd%NextLongitude .and. &
        TmpLat(k) .ge. grd%lat(1,1) .and. TmpLat(k) .lt. grd%lat(grd%im,grd%jm) ) then
      Counter = Counter + 1
      arg%flc(Counter) = TmpFlc(k)
      arg%par(Counter) = TmpPar(k)
      arg%lon(Counter) = TmpLon(k)
      arg%lat(Counter) = TmpLat(k)
      arg%dpt(Counter) = TmpDpt(k)
      arg%res(Counter) = TmpRes(k)
      arg%err(Counter) = TmpErr(k)
      arg%ino(Counter) = TmpIno(k)
    endif

  enddo  
  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  ! arg%res(:) = 1
  ! arg%err(:) = 1.d-2

  
  ! ---
  ! Initialise quality flag
  arg%flg(:) = 1
  arg%rb(:) = 0
  
  ! ---
! Vertical interpolation parameters
  do k = 1,arg%no
     if(arg%flg(k).eq.1)then
        arg%kb(k) = grd%km-1
        do kk = 1,grd%km-1
           if( arg%dpt(k).ge.grd%dep(kk) .and. arg%dpt(k).lt.grd%dep(kk+1) ) then
              arg%kb(k) = kk
              arg%rb(k) = (arg%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
           endif
        enddo
     endif
  enddo
  
  
  ! ---
  ! Count good observations
  arg%nc = 0
  do k=1,arg%no
     if(arg%flg(k).eq.1)then
        arg%nc = arg%nc + 1
     else
        arg%res(k) = 0.
        arg%inc(k) = 0.
        arg%pq1(k) = 0.
        arg%pq2(k) = 0.
        arg%pq3(k) = 0.
        arg%pq4(k) = 0.
        arg%pq5(k) = 0.
        arg%pq6(k) = 0.
        arg%pq7(k) = 0.
        arg%pq8(k) = 0.
     endif
  enddo
  arg%flc(:) = arg%flg(:)

  DEALLOCATE( TmpFlc, TmpPar)
  DEALLOCATE( TmpLon, TmpLat)
  DEALLOCATE( TmpDpt, TmpTim)
  DEALLOCATE( TmpRes, TmpErr)
  DEALLOCATE( TmpIno)  
  
end subroutine get_obs_arg



subroutine int_par_arg
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Get interpolation parameters for a grid                              !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use mpi_str

  implicit none
  
  INTEGER(i4)   ::  i, j, k, ierr
  INTEGER(i4)   ::  i1, j1, k1, idep
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y
  LOGICAL       ::  ins
  
  ins(i,i1) = i.ge.1 .and. i.le.i1
  
  if(arg%no.gt.0) then
     
     arg%flc(:) = arg%flg(:)

     ! ---
     ! Horizontal interpolation parameters
     do k = 1,arg%no
        do j=1,grd%jm-1
           do i=1,grd%im-1
              if( grd%lat(i,j).le.arg%lat(k) .and. grd%lat(i,j+1).gt.arg%lat(k) .and.   &
                   grd%lon(i,j).le.arg%lon(k) .and. grd%lon(i+1,j).gt.arg%lon(k) ) then
                 j1 = j
                 i1 = i
                 q1 = j1 + (arg%lat(k) - grd%lat(i,j)) / (grd%lat(i,j+1) - grd%lat(i,j))
                 p1 = i1 + (arg%lon(k) - grd%lon(i,j)) / (grd%lon(i+1,j) - grd%lon(i,j))
              else if( i .eq. grd%im-1 .and. grd%lat(i,j).le.arg%lat(k) .and. grd%lat(i,j+1).gt.arg%lat(k) .and.   &
                   grd%lon(grd%im,j).le.arg%lon(k) .and. grd%NextLongitude.gt.arg%lon(k) ) then
                 j1 = j
                 i1 = grd%im
                 q1 = j1 + (arg%lat(k) - grd%lat(i,j)) / (grd%lat(i,j+1) - grd%lat(i,j))
                 p1 = i1 + (arg%lon(k) - grd%lon(grd%im,j)) / (grd%NextLongitude - grd%lon(grd%im,j))
              endif
           enddo
        enddo
        
        !     q1 = (arg%lat(k) - grd%lat(1,1)) / grd%dlt + 1.0
        !     j1 = int(q1)
        !     p1 = (arg%lon(k) - grd%lon(1,1)) / grd%dln + 1.0
        !     i1 = int(p1)
        if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
           arg%ib(k) = i1
           arg%jb(k) = j1
           arg%pb(k) = (p1-i1)
           arg%qb(k) = (q1-j1)
        else
           arg%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Undefine masked for multigrid
     do k = 1,arg%no
        if(arg%flc(k).eq.1)then
           i1 = arg%ib(k)
           j1 = arg%jb(k)
           idep = arg%kb(k)+1
           msk4 = grd%global_msk(GlobalRowOffset+i1,j1,idep) + grd%global_msk(GlobalRowOffset+i1+1,j1,idep) + &
            grd%global_msk(GlobalRowOffset+i1,j1+1,idep) + grd%global_msk(GlobalRowOffset+i1+1,j1+1,idep)
           if(msk4.lt.1.) arg%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Horizontal interpolation parameters for each masked grid
     do k = 1,arg%no
        if(arg%flc(k) .eq. 1) then
           
           i1=arg%ib(k)
           p1=arg%pb(k)
           j1=arg%jb(k)
           q1=arg%qb(k)
           
           
           k1=arg%kb(k)
           div_y =  (1.-q1) * max(grd%global_msk(GlobalRowOffset+i1,j1  ,k1),grd%global_msk(GlobalRowOffset+i1+1,j1  ,k1))     &
                +    q1  * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1,k1)
           arg%pq1(k) = grd%global_msk(GlobalRowOffset+i1,j1,k1)                                       &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           arg%pq2(k) = grd%global_msk(GlobalRowOffset+i1+1,j1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1+1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)
           arg%pq3(k) = grd%global_msk(GlobalRowOffset+i1,j1+1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           arg%pq4(k) = grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)                                   &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           k1=arg%kb(k) + 1
           div_y =  (1.-q1) * max(grd%global_msk(GlobalRowOffset+i1,j1  ,k1),grd%global_msk(GlobalRowOffset+i1+1,j1  ,k1))     &
                +    q1  * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1,k1)
           arg%pq5(k) = grd%global_msk(GlobalRowOffset+i1,j1,k1)                                       &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           arg%pq6(k) = grd%global_msk(GlobalRowOffset+i1+1,j1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1+1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)
           arg%pq7(k) = grd%global_msk(GlobalRowOffset+i1,j1+1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           arg%pq8(k) = grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)                                   &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           r1=arg%rb(k)
           arg%pq1(k) = (1.-r1) * arg%pq1(k)
           arg%pq2(k) = (1.-r1) * arg%pq2(k)
           arg%pq3(k) = (1.-r1) * arg%pq3(k)
           arg%pq4(k) = (1.-r1) * arg%pq4(k)
           arg%pq5(k) =     r1  * arg%pq5(k)
           arg%pq6(k) =     r1  * arg%pq6(k)
           arg%pq7(k) =     r1  * arg%pq7(k)
           arg%pq8(k) =     r1  * arg%pq8(k)

           if(arg%pq1(k) .lt. 1.E-16) arg%pq1(k) = dble(0)
           if(arg%pq2(k) .lt. 1.E-16) arg%pq2(k) = dble(0)
           if(arg%pq3(k) .lt. 1.E-16) arg%pq3(k) = dble(0)
           if(arg%pq4(k) .lt. 1.E-16) arg%pq4(k) = dble(0)
           if(arg%pq5(k) .lt. 1.E-16) arg%pq5(k) = dble(0)
           if(arg%pq6(k) .lt. 1.E-16) arg%pq6(k) = dble(0)
           if(arg%pq7(k) .lt. 1.E-16) arg%pq7(k) = dble(0)
           if(arg%pq8(k) .lt. 1.E-16) arg%pq8(k) = dble(0)


        endif
     enddo
     
     
     ! ---
     ! Count good observations
     arg%nc = 0
     do k=1,arg%no
        if(arg%flc(k).eq.1)then
           arg%nc = arg%nc + 1
        endif
     enddo
     
  endif
  
  arg%nc_global = 0
  call MPI_Allreduce(arg%nc, arg%nc_global, 1, MPI_INT, MPI_SUM, MyCommWorld, ierr)

  if(MyRank .eq. 0) then
     write(drv%dia,*)'Real number of ARGO observations: ',arg%nc_global
     print*,'Good argo observations: ',arg%nc_global
  end if

  DEALLOCATE ( arg%ino, arg%flg, arg%par)  
  DEALLOCATE ( arg%lon, arg%lat, arg%dpt, arg%tim)
  DEALLOCATE ( arg%pb, arg%qb, arg%rb)

end subroutine int_par_arg
