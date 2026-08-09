subroutine get_densincr_arg
  
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
  ! Load Density increments                                                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use dnc_str
!  use obs_str
  use mpi_str
  use filenames
!  use bio_str
  
  implicit none
  
  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, kk, i
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: TmpFlc, TmpLon, TmpLat
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: TmpDpt, TmpErr!, TmpStd
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: TmpInc!, TmpCorr 
  INTEGER(i4)   :: GlobalDncNum, Counter, ierr
  character(len=1024) :: filename
  
  dnc%no  = 0
  dnc%nc  = 0
  
  
  ! ---
  ! Allocate memory for observations
  if(MyId .eq. 0) then
    open(511,file=trim(INCR_FILE))
    read(511,'(I4)') GlobalDncNum
    write(drv%dia,*)'Number of density increments: ', GlobalDncNum
  endif

  call MPI_Bcast(GlobalDncNum, 1, MPI_INT, 0, Var3DCommunicator, ierr)

  if(GlobalDncNum .eq. 0)then
    if(MyId .eq. 0) &
      close(511)
    return
  endif

  ALLOCATE( TmpFlc(GlobalDncNum)) !, TmpPar(GlobalDncNum))
  ALLOCATE( TmpLon(GlobalDncNum), TmpLat(GlobalDncNum))
  ALLOCATE( TmpDpt(GlobalDncNum))!, TmpTim(GlobalDncNum))
  ALLOCATE( TmpInc(GlobalDncNum))!, TmpCorr(GlobalDncNum))
  ALLOCATE( TmpErr(GlobalDncNum))!, TmpStd(GlobalDncNum))

  if(MyId .eq. 0) then
    ! process 0 reads all the density increments
    do k=1,GlobalDncNum
      read (511,*) &
            TmpFlc(k), & !TmpPar(k), &
            TmpLon(k), TmpLat(k), &
            TmpDpt(k), &!TmpTim(k), &
            TmpInc(k), &!TmpCorr(k), &
            TmpErr(k) !, TmpStd(k)
    end do
    close (511)
  endif

  
  call MPI_Bcast(TmpFlc, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
!   call MPI_Bcast(TmpPar, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  call MPI_Bcast(TmpLon, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  call MPI_Bcast(TmpLat, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  call MPI_Bcast(TmpDpt, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  !   call MPI_Bcast(TmpTim, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  call MPI_Bcast(TmpInc, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
!   call MPI_Bcast(TmpCorr, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  call MPI_Bcast(TmpErr, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
!   call MPI_Bcast(TmpStd, GlobalDncNum, MPI_REAL8, 0, Var3DCommunicator, ierr)
  
  ! Counting the number of observations that falls in the domain
  Counter = 0
  do k=1,GlobalDncNum
   if( TmpLon(k) .ge. grd%lon(1,1) .and. TmpLon(k) .lt. grd%NextLongitude .and. &
   TmpLat(k) .ge. grd%lat(1,1) .and. TmpLat(k) .lt. grd%lat(grd%im,grd%jm) ) then
      ! if(drv%dnc.eq.1) then
         Counter = Counter + 1
      ! endif
   endif
enddo


  if(drv%Verbose .eq. 1) &
       print*, "MyId", MyId, "has",Counter,"Density increments"

  dnc%no  = Counter

  ALLOCATE ( dnc%flg(dnc%no), dnc%flc(dnc%no)) !, dnc%par(dnc%no))
  ALLOCATE ( dnc%lon(dnc%no), dnc%lat(dnc%no), dnc%dpt(dnc%no)) !, dnc%tim(dnc%no))
  ALLOCATE ( dnc%inc(dnc%no))
!   ALLOCATE ( dnc%corr(dnc%no))
  ALLOCATE ( dnc%err(dnc%no))
  ALLOCATE ( dnc%res(dnc%no))
!   ALLOCATE ( dnc%std(dnc%no))
  ALLOCATE ( dnc%ib(dnc%no), dnc%jb(dnc%no), dnc%kb(dnc%no))
  ALLOCATE ( dnc%pb(dnc%no), dnc%qb(dnc%no), dnc%rb(dnc%no))
  ALLOCATE ( dnc%pq1(dnc%no), dnc%pq2(dnc%no), dnc%pq3(dnc%no), dnc%pq4(dnc%no))
  ALLOCATE ( dnc%pq5(dnc%no), dnc%pq6(dnc%no), dnc%pq7(dnc%no), dnc%pq8(dnc%no))

  Counter = 0
  do k=1,GlobalDncNum
    if( TmpLon(k) .ge. grd%lon(1,1) .and. TmpLon(k) .lt. grd%NextLongitude .and. &
        TmpLat(k) .ge. grd%lat(1,1) .and. TmpLat(k) .lt. grd%lat(grd%im,grd%jm) ) then
      !   if(drv%dnc.eq.1) then
            Counter = Counter + 1
            dnc%flc(Counter) = TmpFlc(k)
            ! dnc%par(Counter) = TmpPar(k)
            dnc%lon(Counter) = TmpLon(k)
            dnc%lat(Counter) = TmpLat(k)
            dnc%dpt(Counter) = TmpDpt(k)
            dnc%res(Counter) = TmpInc(k) !called res for analogy with get_obs_arg
            ! dnc%corr(Counter) = TmpCorr(k)
            dnc%err(Counter) = TmpErr(k)
            ! dnc%std(Counter) = TmpStd(k)
            ! dnc%ino(Counter) = TmpIno(k)
      !   endif
    endif

  enddo
  
  
  ! ---
  ! Initialise quality flag
  dnc%flg(:) = 1
  dnc%rb(:) = 0
  
  ! ---
! Vertical interpolation parameters
  do k = 1,dnc%no
     if(dnc%flg(k).eq.1)then
        dnc%kb(k) = grd%km-1
        do kk = 1,grd%km-1
           if( dnc%dpt(k).ge.grd%dep(kk) .and. dnc%dpt(k).lt.grd%dep(kk+1) ) then
              dnc%kb(k) = kk
              dnc%rb(k) = (dnc%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
           else if ( dnc%dpt(k).ge.0 .and. dnc%dpt(k).lt.grd%dep(1)) then
              dnc%kb(k) = 1
           endif
        enddo
     endif
  enddo
  
  
  ! ---
  ! Count good observations
  dnc%nc = 0
  do k=1,dnc%no
     if(dnc%flg(k).eq.1)then
        dnc%nc = dnc%nc + 1
     else
        dnc%res(k) = 0.
        dnc%inc(k) = 0.
        dnc%pq1(k) = 0.
        dnc%pq2(k) = 0.
        dnc%pq3(k) = 0.
        dnc%pq4(k) = 0.
        dnc%pq5(k) = 0.
        dnc%pq6(k) = 0.
        dnc%pq7(k) = 0.
        dnc%pq8(k) = 0.
     endif
  enddo
  dnc%flc(:) = dnc%flg(:)

  DEALLOCATE( TmpFlc)!, TmpPar)
  DEALLOCATE( TmpLon, TmpLat)
  DEALLOCATE( TmpDpt)!, TmpTim)
  DEALLOCATE( TmpInc)!, TmpTim)
!   DEALLOCATE( TmpCorr)!, TmpTim)
  DEALLOCATE( TmpErr)
!   DEALLOCATE( TMPStd)
!   DEALLOCATE( TmpIno)
  
end subroutine get_densincr_arg



subroutine int_par_dnc
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Get interpolation parameters for a grid                              !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
!   use eof_str
!   use obs_str
  use dnc_str
  use mpi_str

  implicit none
  
  INTEGER(i4)   ::  i, j, k, ierr, kind, kk
  INTEGER(i4)   ::  i1, j1, k1, idep
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y
  LOGICAL       ::  ins
  
  ins(i,i1) = i.ge.1 .and. i.le.i1
  
  if(dnc%no.gt.0) then
     
     dnc%flc(:) = dnc%flg(:)

     ! ---
     ! Horizontal interpolation parameters
     do k = 1,dnc%no
        do j=1,grd%jm-1
           do i=1,grd%im-1
              if( grd%lat(i,j).le.dnc%lat(k) .and. grd%lat(i,j+1).gt.dnc%lat(k) .and.   &
                   grd%lon(i,j).le.dnc%lon(k) .and. grd%lon(i+1,j).gt.dnc%lon(k) ) then
                 j1 = j
                 i1 = i
                 q1 = j1 + (dnc%lat(k) - grd%lat(i,j)) / (grd%lat(i,j+1) - grd%lat(i,j))
                 p1 = i1 + (dnc%lon(k) - grd%lon(i,j)) / (grd%lon(i+1,j) - grd%lon(i,j))
              else if( i .eq. grd%im-1 .and. grd%lat(i,j).le.dnc%lat(k) .and. grd%lat(i,j+1).gt.dnc%lat(k) .and.   &
                   grd%lon(grd%im,j).le.dnc%lon(k) .and. grd%NextLongitude.gt.dnc%lon(k) ) then
                 j1 = j
                 i1 = grd%im
                 q1 = j1 + (dnc%lat(k) - grd%lat(i,j)) / (grd%lat(i,j+1) - grd%lat(i,j))
                 p1 = i1 + (dnc%lon(k) - grd%lon(grd%im,j)) / (grd%NextLongitude - grd%lon(grd%im,j))
              endif
           enddo
        enddo
        
        !     q1 = (dnc%lat(k) - grd%lat(1,1)) / grd%dlt + 1.0
        !     j1 = int(q1)
        !     p1 = (dnc%lon(k) - grd%lon(1,1)) / grd%dln + 1.0
        !     i1 = int(p1)
        if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
           dnc%ib(k) = i1
           dnc%jb(k) = j1
           dnc%pb(k) = (p1-i1)
           dnc%qb(k) = (q1-j1)
        else
           dnc%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Undefine masked for multigrid
     do k = 1,dnc%no
        if(dnc%flc(k).eq.1)then
           i1 = dnc%ib(k)
           j1 = dnc%jb(k)
           idep = dnc%kb(k)+1
           msk4 = grd%global_msk(GlobalRowOffset+i1,j1,idep) + grd%global_msk(GlobalRowOffset+i1+1,j1,idep) + &
            grd%global_msk(GlobalRowOffset+i1,j1+1,idep) + grd%global_msk(GlobalRowOffset+i1+1,j1+1,idep)
           if(msk4.lt.1.) dnc%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Horizontal interpolation parameters for each masked grid
     do k = 1,dnc%no
        if(dnc%flg(k) .eq. 1) then
        ! if(dnc%flg(k) .eq. 1) then ! to verify that it works also in this case
           
           i1=dnc%ib(k)
           p1=dnc%pb(k)
           j1=dnc%jb(k)
           q1=dnc%qb(k)
           
           
           k1=dnc%kb(k)
           div_y =  (1.-q1) * max(grd%global_msk(GlobalRowOffset+i1,j1  ,k1),grd%global_msk(GlobalRowOffset+i1+1,j1  ,k1))     &
                +    q1  * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1,k1)
           dnc%pq1(k) = grd%global_msk(GlobalRowOffset+i1,j1,k1)                                       &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           dnc%pq2(k) = grd%global_msk(GlobalRowOffset+i1+1,j1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1+1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)
           dnc%pq3(k) = grd%global_msk(GlobalRowOffset+i1,j1+1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           dnc%pq4(k) = grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)                                   &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           k1=dnc%kb(k) + 1
           div_y =  (1.-q1) * max(grd%global_msk(GlobalRowOffset+i1,j1  ,k1),grd%global_msk(GlobalRowOffset+i1+1,j1  ,k1))     &
                +    q1  * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1,k1)
           dnc%pq5(k) = grd%global_msk(GlobalRowOffset+i1,j1,k1)                                       &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           dnc%pq6(k) = grd%global_msk(GlobalRowOffset+i1+1,j1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%global_msk(GlobalRowOffset+i1  ,j1+1,k1) + p1 * grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)
           dnc%pq7(k) = grd%global_msk(GlobalRowOffset+i1,j1+1,k1)                                     &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           dnc%pq8(k) = grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1)                                   &
                * max(grd%global_msk(GlobalRowOffset+i1,j1+1,k1),grd%global_msk(GlobalRowOffset+i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           r1=dnc%rb(k)
           dnc%pq1(k) = (1.-r1) * dnc%pq1(k)
           dnc%pq2(k) = (1.-r1) * dnc%pq2(k)
           dnc%pq3(k) = (1.-r1) * dnc%pq3(k)
           dnc%pq4(k) = (1.-r1) * dnc%pq4(k)
           dnc%pq5(k) =     r1  * dnc%pq5(k)
           dnc%pq6(k) =     r1  * dnc%pq6(k)
           dnc%pq7(k) =     r1  * dnc%pq7(k)
           dnc%pq8(k) =     r1  * dnc%pq8(k)

           if(dnc%pq1(k) .lt. 1.E-16) dnc%pq1(k) = dble(0)
           if(dnc%pq2(k) .lt. 1.E-16) dnc%pq2(k) = dble(0)
           if(dnc%pq3(k) .lt. 1.E-16) dnc%pq3(k) = dble(0)
           if(dnc%pq4(k) .lt. 1.E-16) dnc%pq4(k) = dble(0)
           if(dnc%pq5(k) .lt. 1.E-16) dnc%pq5(k) = dble(0)
           if(dnc%pq6(k) .lt. 1.E-16) dnc%pq6(k) = dble(0)
           if(dnc%pq7(k) .lt. 1.E-16) dnc%pq7(k) = dble(0)
           if(dnc%pq8(k) .lt. 1.E-16) dnc%pq8(k) = dble(0)


        endif
     enddo
     

     
     ! ---
     ! Count good observations
     dnc%nc = 0
     do k=1,dnc%no
        if(dnc%flc(k).eq.1)then
           dnc%nc = dnc%nc + 1
        endif
     enddo
     
  endif
  
  dnc%nc_global = 0
  call MPI_Allreduce(dnc%nc, dnc%nc_global, 1, MPI_INT, MPI_SUM, Var3DCommunicator, ierr)
  
  if(MyId .eq. 0) then
     write(drv%dia,*)'Real number of density increments: ',dnc%nc_global
     print*,'Good density increments: ',dnc%nc_global
  end if

  DEALLOCATE ( dnc%flg)  
  DEALLOCATE ( dnc%lon, dnc%lat, dnc%dpt) !, dnc%tim)
  DEALLOCATE ( dnc%pb, dnc%qb, dnc%rb)

end subroutine int_par_dnc
