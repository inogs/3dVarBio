subroutine parallel_def_cov

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
  ! Define filter constants, EOFs, etc.                                  !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use eof_str
  use cns_str
  use rcfl
  use mpi_str
  use mpi
  
  implicit none
  
  INTEGER(i4)                 :: k, nspl, i, j, kk, l
  REAL(r8)                    :: E, dst, Lmean, mean_rad
  REAL(r8)    , ALLOCATABLE   :: sfct(:), al(:), bt(:)
  INTEGER(i4) , ALLOCATABLE   :: jnxx(:)
  INTEGER nthreads, threadid
  INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  INTEGER(i8) :: ierr, iProc
  REAL(r8), allocatable :: SendBuf2D(:,:), RecBuf2D(:,:), DefBuf2D(:,:)
  REAL(r8), allocatable :: SendBuf1D(:), RecBuf1D(:), TmpBuf1D(:), DefBuf1D(:)
  REAL(r8), allocatable :: SendBuf3D(:,:,:), RecBuf3D(:,:,:), DefBuf3D(:,:,:)
  REAL(r8), allocatable :: ColBuf3D(:,:,:)
  
  call parallel_rdrcorr
  
  ! decomment the next line to use one single correlation radius
  ! rcf%Lxyz = rcf%L

  nthreads = 1
  threadid = 0
  !$OMP PARALLEL
  !$ nthreads = OMP_GET_NUM_THREADS()
  !$ threadid = OMP_GET_THREAD_NUM()
  if(threadid.eq.0) then
     if(MyRank .eq. 0) &
          write(*,*) "OMP version with threads = ", nthreads
  endif
  !$OMP END PARALLEL
  ! ---
  ! Recursive filter constants
  !---------
  ! Create table
  
  !nspl = max(grd%jm,grd%im)
  nspl = max(GlobalRow,GlobalCol)
  ALLOCATE ( sfct(nspl)) ; sfct = huge(sfct(1))
  ALLOCATE ( jnxx(nspl)) ; jnxx = huge(jnxx(1))
  ALLOCATE ( al(nspl))   ; al   = huge(al(1))
  ALLOCATE ( bt(nspl))   ; bt   = huge(bt(1))
  
  
  rcf%ntb = min(20,min(grd%jm,grd%im))
  
  ! KB grid problem (chl assimilation)
  rcf%ntb = 1000
  !
  
  ALLOCATE ( rcf%al(rcf%ntb)) ; rcf%al = huge(rcf%al(1))
  ALLOCATE ( rcf%sc(grd%km,rcf%ntb)) ; rcf%sc = huge(rcf%sc(1,1))
  
  rcf%dsmn =  1.e20
  rcf%dsmx = -1.e20
  do j=1,grd%jm
     do i=1,grd%im
        rcf%dsmn = min(rcf%dsmn,min(grd%dx(i,j),grd%dy(i,j)))
        rcf%dsmx = max(rcf%dsmx,max(grd%dx(i,j),grd%dy(i,j)))
     enddo
  enddo

  ! Computes the global maximum and minimum
  call MPI_Allreduce(MPI_IN_PLACE, rcf%dsmx, 1, MPI_REAL8, MPI_MAX, MyCommWorld, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, rcf%dsmn, 1, MPI_REAL8, MPI_MIN, MyCommWorld, ierr)
  
  rcf%dsmx = rcf%dsmx + max(1.d0,(rcf%dsmx-rcf%dsmn)/(rcf%ntb-2.))
  
  rcf%dsl = (rcf%dsmx-rcf%dsmn) / (rcf%ntb-1.)
  
  do k=1,grd%km
     Lmean = mean_rad(k,rcf%Lxyz(:,:,k)) !for each level, mean of Lxyz
     do l=1,rcf%ntb
        dst = rcf%dsmn + (l-1.) * rcf%dsl
        !E = (2. * rcf%ntr) * dst**2 / (4. * rcf%L**2)  !da capire bene
        E   = (2. * rcf%ntr) * dst**2 / (4. * Lmean**2)
        rcf%al(l) = 1. + E - sqrt(E*(E+2.))
        rcf%alp   = rcf%al(l)
        sfct(:) = 0.
        al(:) = rcf%al(l)
        bt(:) = rcf%al(l)
        do j=1,nspl
                jnxx(j) = j
        enddo
        sfct(nspl/2+1) = 1.
        call rcfl_y_init    ( 1, nspl, 1, nspl, al, bt, sfct, jnxx, nspl)
        call rcfl_y_ad_init ( 1, nspl, 1, nspl, al, bt, sfct, jnxx, nspl)
        rcf%sc(k,l) = sfct(nspl/2+1)
    enddo
  enddo

  DEALLOCATE ( sfct, jnxx, al, bt ) 
  
  ! preparing the call to MPI_ALLTOALL
  ALLOCATE(RecBuf1D(GlobalRow*localCol))
  ALLOCATE(DefBuf2D(GlobalRow, localCol))

  call MPI_Alltoallv(grd%dx, SendCountX2D, SendDisplX2D, MPI_REAL8, &
       RecBuf1D, RecCountX2D, RecDisplX2D, MPI_REAL8, MyCommWorld, ierr)

  do j=1, localCol
     do iProc=0, NumProcI-1
        do i=1, RecCountX2D(iProc+1)/localCol
           DefBuf2D(i + RecDisplX2D(iProc+1)/localCol, j) = &
                RecBuf1D(i + (j-1)*RecCountX2D(iProc+1)/localCol + RecDisplX2D(iProc+1))
        end do
     end do
  end do

  do k=1,grd%km
     do j=1,localCol
        do i=1,GlobalRow
           dst = ( DefBuf2D(i,j) - rcf%dsmn )/rcf%dsl
           l = int(dst) + 1
           dst = dst - real(l-1)
           grd%scx(i,j,k) = sqrt( 1./ (rcf%sc(k,l)*(1.-dst) + rcf%sc(k,l+1)*dst) ) 
        enddo
     enddo
  enddo

  do k=1,grd%km
     do j=1,localCol
        do i=2,GlobalRow
           dst = (DefBuf2D(i-1,j) + DefBuf2D(i,j)) * 0.5 
           E   = (2. * rcf%ntr) * dst**2 / (4. * rcf%Lxyz(i,j,k)**2)
           grd%alx(i,j,k) = 1. + E - sqrt(E*(E+2.))
        enddo
        do i=1,GlobalRow-1
           dst = (DefBuf2D(i,j) + DefBuf2D(i+1,j)) * 0.5 
           E   = (2. * rcf%ntr) * dst**2 / (4. * rcf%Lxyz(i,j,k)**2)
           grd%btx(i,j,k) = 1. + E - sqrt(E*(E+2.))
        enddo
     enddo
  enddo

  do k=1,grd%km
     do j=1,localCol
        do i=1,GlobalRow
           grd%istp(i,j,k) = int( rcf%Lxyz(i,j,k) * rcf%efc / DefBuf2D(i,j) )+1
        enddo
     enddo
  enddo



  ! MPI_ALLTOALL SECTION
  DEALLOCATE(RecBuf1D, DefBuf2D)
  ALLOCATE(SendBuf2D(grd%jm, grd%im))
  ALLOCATE( RecBuf1D(localRow*GlobalCol))
  ALLOCATE( DefBuf2D(localRow, GlobalCol))

  do j=1,grd%jm
     do i=1,grd%im
        SendBuf2D(j,i) = grd%dy(i,j)
     end do
  end do

  call MPI_Alltoallv(SendBuf2D, SendCountY2D, SendDisplY2D, MPI_REAL8, &
       RecBuf1D, RecCountY2D, RecDisplY2D, MPI_REAL8, CommSliceY, ierr)

  do i=1, localRow
     do iProc=0, NumProcJ-1
        do j=1, RecCountY2D(iProc+1)/localRow
           DefBuf2D(i, j + RecDisplY2D(iProc+1)/localRow) = &
                RecBuf1D(j + (i-1)*RecCountY2D(iProc+1)/localRow + RecDisplY2D(iProc+1))
        end do
     end do
  end do

  do k=1,grd%km
     do j=1,GlobalCol
        do i=1,localRow
           dst = ( DefBuf2D(i,j) - rcf%dsmn )/rcf%dsl
           l = int(dst) + 1
           dst = dst - real(l-1)
           grd%scy(i,j,k) = sqrt( 1./ (rcf%sc(k,l)*(1.-dst) + rcf%sc(k,l+1)*dst) ) 
        enddo
     enddo
  enddo

  do k=1,grd%km
     do j=2,GlobalCol
        do i=1,localRow
           dst = (DefBuf2D(i,j-1) + DefBuf2D(i,j)) * 0.5 
           E   = (2. * rcf%ntr) * dst**2 / (4. * rcf%Lxyz(i,j,k)**2)
           grd%aly(i,j,k) = 1. + E - sqrt(E*(E+2.))
        enddo
     enddo
     do j=1,GlobalCol-1
        do i=1,localRow
           dst = (DefBuf2D(i,j) + DefBuf2D(i,j+1)) * 0.5 
           E   = (2. * rcf%ntr) * dst**2 / (4. * rcf%Lxyz(i,j,k)**2)
           grd%bty(i,j,k) = 1. + E - sqrt(E*(E+2.))
        enddo
     enddo
  enddo
  
  do k=1,grd%km
     grd%alx(     1,   :,   k) = grd%alx(       2,   :,    k)
     grd%btx(GlobalRow,:,   k) = grd%btx(GlobalRow-1,:,    k)
     grd%aly(   :     ,1,   k) = grd%aly(  :  ,   2   ,    k)
     grd%bty(   :,GlobalCol,k) = grd%bty(  :  ,GlobalCol-1,k)
  enddo
  
  !---
  ! Define extended grids
  do k=1,grd%km
     do j=1,GlobalCol
        do i=1,localRow
           grd%jstp(i,j,k) = int( rcf%Lxyz(i,j,k) * rcf%efc / DefBuf2D(i,j) )+1
        enddo
     enddo
  enddo
  grd%imax   = 0
  grd%jmax   = 0

  DEALLOCATE(RecBuf1D)
  ALLOCATE(SendBuf3D(grd%km,grd%jm,grd%im))
  ALLOCATE( RecBuf1D(localRow*GlobalCol*grd%km))
  ALLOCATE( DefBuf3D(localRow, GlobalCol, grd%km))

  !************* HORIZONTAL SLICING *************!
  do k=1,grd%km
     do j=1,grd%jm
        do i=1,grd%im
           SendBuf3D(k,j,i) = grd%msr(i,j,k)
        end do
     end do
  end do

  call MPI_Alltoallv(SendBuf3D, SendCountY4D, SendDisplY4D, MPI_REAL8, &
       RecBuf1D, RecCountY4D, RecDisplY4D, MPI_REAL8, CommSliceY, ierr)
  
  ! Reordering data
  if(size .eq. 1) then
     do k=1,grd%km
        do j=1,grd%jm
           do i=1,grd%im
              DefBuf3D(i,j,k) = grd%msr(i,j,k)
           end do
        end do
     end do
  else
     do i=1, localRow
        do iProc=0, NumProcJ-1
           do j=1, RecCountY4D(iProc+1)/(localRow*grd%km)
              do k=1, grd%km
                 DefBuf3D(i,j+RecDisplY4D(iProc+1)/(localRow*grd%km),k) = &
                        RecBuf1D(k + (j-1)*grd%km + (i-1)*RecCountY4D(iProc+1)/localRow + RecDisplY4D(iProc+1))
              end do
           end do
        end do
     end do

  end if

  !************* VERTICAL SLICING *************!
  DEALLOCATE(SendBuf3D, RecBuf1D)
  ALLOCATE(SendBuf3D(grd%km, grd%im, grd%jm))
  ALLOCATE( RecBuf1D(GlobalRow*localCol*grd%km))
  ALLOCATE( ColBuf3D(GlobalRow, localCol, grd%km))

  do k=1,grd%km
     do j=1,grd%jm
        do i=1,grd%im
           SendBuf3D(k,i,j) = grd%msr(i,j,k)
        end do
     end do
  end do

  call MPI_Alltoallv(SendBuf3D, SendCountX4D, SendDisplX4D, MPI_REAL8, &
       RecBuf1D, RecCountX4D, RecDisplX4D, MPI_REAL8, MyCommWorld, ierr)

  do j=1, localCol
     do iProc=0, NumProcI-1
        do i=1, RecCountX4D(iProc+1)/(localCol*grd%km)
           do k=1, grd%km
              ColBuf3D(i + RecDisplX4D(iProc+1)/(localCol*grd%km), j, k) = &
                   RecBuf1D(k + (i-1)*grd%km + (j-1)*RecCountX4D(iProc+1)/localCol + RecDisplX4D(iProc+1))
           end do
        end do
     end do
  end do

  do k = 1, grd%km
     
     grd%imx(k) = 0
     do j = 1, localCol
        kk = grd%istp(1,j,k)
        if( ColBuf3D(1,j,k).eq.1. ) kk = kk + 1
        grd%inx(1,j,k) = kk
        do i = 2, GlobalRow
           if( ColBuf3D(i,j,k).eq.0. .and. ColBuf3D(i-1,j,k).eq.1. ) then
              kk = kk + grd%istp(i,j,k)
           else if( ColBuf3D(i,j,k).eq.1. .and. ColBuf3D(i-1,j,k).eq.0. ) then
              kk = kk + grd%istp(i,j,k) + 1
           else if( ColBuf3D(i,j,k).eq.1. ) then
              kk = kk + 1
           endif
           grd%inx(i,j,k) = kk
        enddo
        grd%imx(k) = max( grd%imx(k), kk+grd%istp(GlobalRow,j,k))
     enddo
     grd%imax   = max( grd%imax, grd%imx(k))
     
     grd%jmx(k) = 0
     do i = 1, localRow
        kk = grd%jstp(i,1,k)
        if( DefBuf3D(i,1,k).eq.1. ) kk = kk + 1
        grd%jnx(i,1,k) = kk
        do j = 2, GlobalCol
           if( DefBuf3D(i,j,k).eq.0. .and. DefBuf3D(i,j-1,k).eq.1. ) then
              kk = kk + grd%jstp(i,j,k)
           else if( DefBuf3D(i,j,k).eq.1. .and. DefBuf3D(i,j-1,k).eq.0. ) then
              kk = kk + grd%jstp(i,j,k) + 1
           else if( DefBuf3D(i,j,k).eq.1. ) then
              kk = kk + 1
           endif
           grd%jnx(i,j,k) = kk
        enddo
        grd%jmx(k) = max( grd%jmx(k), kk+grd%jstp(i,GlobalCol,k))
     enddo
     grd%jmax   = max( grd%jmax, grd%jmx(k))
     
  enddo

  call MPI_Allreduce(MPI_IN_PLACE, grd%imax, 1, MPI_INT, MPI_MAX, MyCommWorld, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, grd%jmax, 1, MPI_INT, MPI_MAX, MyCommWorld, ierr)
  
  ALLOCATE( grd%aex(localCol,grd%imax,grd%km)) ; grd%aex(:,:,:) = 0.0
  ALLOCATE( grd%bex(localCol,grd%imax,grd%km)) ; grd%bex(:,:,:) = 0.0
  ALLOCATE( grd%aey(localRow,grd%jmax,grd%km)) ; grd%aey(:,:,:) = 0.0
  ALLOCATE( grd%bey(localRow,grd%jmax,grd%km)) ; grd%bey(:,:,:) = 0.0
  
  
  do k = 1, grd%km
        
     do j = 1, localCol
        kk = grd%istp(1,j,k)
        if( ColBuf3D(1,j,k).eq.1. ) then
           kk = kk + 1
           grd%aex(j,1:kk,k) = grd%alx(1,j,k)
           grd%bex(j,1:kk,k) = grd%btx(1,j,k)
        endif
        do i = 2, GlobalRow
           if( ColBuf3D(i,j,k).eq.0. .and. ColBuf3D(i-1,j,k).eq.1. ) then
              grd%aex(j,kk+1:kk+grd%istp(i,j,k),k) = grd%alx(i,j,k)
              grd%bex(j,kk+1:kk+grd%istp(i,j,k),k) = grd%btx(i,j,k)
              kk = kk + grd%istp(i,j,k)
           else if( ColBuf3D(i,j,k).eq.1. .and. ColBuf3D(i-1,j,k).eq.0. ) then
              grd%aex(j,kk+1:kk+grd%istp(i,j,k)+1,k) = grd%alx(i,j,k)
              grd%bex(j,kk+1:kk+grd%istp(i,j,k)+1,k) = grd%btx(i,j,k)
              kk = kk + grd%istp(i,j,k) + 1
           else if( ColBuf3D(i,j,k).eq.1. ) then
              grd%aex(j,kk+1,k) = grd%alx(i,j,k)
              grd%bex(j,kk+1,k) = grd%btx(i,j,k)
              kk = kk + 1
           endif
        enddo
     enddo
     
     do i = 1, localRow
        kk = grd%jstp(i,1,k)
        if( DefBuf3D(i,1,k).eq.1. ) then
           kk = kk + 1
           grd%aey(i,1:kk,k) = grd%aly(i,1,k)
           grd%bey(i,1:kk,k) = grd%bty(i,1,k)
        endif
        do j = 2, GlobalCol
           if( DefBuf3D(i,j,k).eq.0. .and. DefBuf3D(i,j-1,k).eq.1. ) then
              grd%aey(i,kk+1:kk+grd%jstp(i,j,k),k) = grd%aly(i,j,k)
              grd%bey(i,kk+1:kk+grd%jstp(i,j,k),k) = grd%bty(i,j,k)
              kk = kk + grd%jstp(i,j,k)
           else if( DefBuf3D(i,j,k).eq.1. .and. DefBuf3D(i,j-1,k).eq.0. ) then
              grd%aey(i,kk+1:kk+grd%jstp(i,j,k)+1,k) = grd%aly(i,j,k)
              grd%bey(i,kk+1:kk+grd%jstp(i,j,k)+1,k) = grd%bty(i,j,k)
              kk = kk + grd%jstp(i,j,k) + 1
           else if( DefBuf3D(i,j,k).eq.1. ) then
              grd%aey(i,kk+1,k) = grd%aly(i,j,k)
              grd%bey(i,kk+1,k) = grd%bty(i,j,k)
              kk = kk + 1
           endif
        enddo
     enddo
     
  enddo

  ! ---
  ! Vertical EOFs
           
  ros%kmt = grd%km * grd%nchl 

  call parallel_rdeofs
  
  ALLOCATE ( grd%ro(    grd%im, grd%jm, ros%neof))   ; grd%ro    = 0.0
  ALLOCATE ( grd%ro_ad( grd%im, grd%jm, ros%neof))   ; grd%ro_ad = 0.0
  ALLOCATE ( Dump_vip ( grd%im, grd%jm, ros%neof))   ; Dump_vip  = 0.0
  
  if(MyRank .eq. 0) &
       write(*,*) 'rcfl allocation :', grd%jmax, grd%imax, nthreads
  ALLOCATE ( a_rcx(localCol,grd%imax,nthreads)) ; a_rcx = huge(a_rcx(1,1,1))
  ALLOCATE ( b_rcx(localCol,grd%imax,nthreads)) ; b_rcx = huge(b_rcx(1,1,1))
  ALLOCATE ( c_rcx(localCol,grd%imax,nthreads)) ; c_rcx = huge(c_rcx(1,1,1))
  
  ALLOCATE ( a_rcy(localRow,grd%jmax,nthreads)) ; a_rcy = huge(a_rcy(1,1,1))
  ALLOCATE ( b_rcy(localRow,grd%jmax,nthreads)) ; b_rcy = huge(b_rcy(1,1,1))
  ALLOCATE ( c_rcy(localRow,grd%jmax,nthreads)) ; c_rcy = huge(c_rcy(1,1,1))
  
  
  ALLOCATE ( alp_rcx(localCol,grd%imax,nthreads)) ; alp_rcx = huge(alp_rcx(1,1,1))
  ALLOCATE ( bta_rcx(localCol,grd%imax,nthreads)) ; bta_rcx = huge(bta_rcx(1,1,1))
  
  ALLOCATE ( alp_rcy(localRow,grd%jmax,nthreads)) ; alp_rcy = huge(alp_rcy(1,1,1))
  ALLOCATE ( bta_rcy(localRow,grd%jmax,nthreads)) ; bta_rcy = huge(bta_rcy(1,1,1))

  DEALLOCATE(SendBuf2D, RecBuf1D, DefBuf2D)
  DEALLOCATE(SendBuf3D)
  DEALLOCATE(ColBuf3D,  DefBuf3D)
  
end subroutine parallel_def_cov
