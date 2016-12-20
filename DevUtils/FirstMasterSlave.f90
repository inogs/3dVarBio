program CommunicationTime
  
  use mpi

  implicit none

  integer :: ierr, MyRank, size, i, j, k
  integer :: im, jm, km, nlev
  real(8), pointer, dimension(:,:,:) :: ToSend !, ToRec

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr) 

  ! im = 722
  ! jm = 255
  ! km = 26
  ! nlev = 1

  im = 10
  jm = 10
  km = 10
  nlev = 1
  
  if(MyRank .eq. 0) then
     ALLOCATE(ToSend(im,jm,km))
     
     do k=1,km
        do j=1,jm
           do i=1,im
              ToSend(i,j,k) = 0. !dble(i + (j-1)*jm + (k-1)*km)
           end do
        end do
     end do

     call Master(ToSend, im, jm, km, nlev, size)
     
  else
     call Slave(im, jm, km, nlev, MyRank)
  end if

  print*, "process ", MyRank, " ending"
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if(MyRank .eq. 0) then
     do k=1,km
        print*, "Printing level", k
        ! do j=1,jm
           ! do i=1,im
           print*, ToSend(:,:,k)
           ! end do
        ! end do
     end do
  end if

  call MPI_Finalize(ierr)

end program CommunicationTime

subroutine ReadySlave(WhoIs, Res, ComputedLevel, im, jm, nlev)

  use mpi

  implicit none

  integer :: WhoIs, ComputedLevel, im, jm, nlev, ierr
  integer :: MyStatus(MPI_STATUS_SIZE)
  real(8), dimension(im, jm, nlev) :: Res

  ! print*, "Rank 0 within ReadySlave"
  call MPI_Recv(Res, im*jm*nlev, MPI_REAL8, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MyStatus, ierr)
  WhoIs = MyStatus(MPI_SOURCE)
  ComputedLevel = MyStatus(MPI_TAG)
  ! print*, "Rank 0 received from ", WhoIs, "with tag ", ComputedLevel

end subroutine ReadySlave

subroutine Master(ToSend, im, jm, km, nlev, size)

  use mpi
  implicit none
  
  integer :: im, jm, km, nlev, size
  integer :: i, j, k, ReadyProc, ierr, nlev_tmp
  integer :: ComputedLevel, tmpk, RealCounter
  real(8), dimension(im, jm, km) :: ToSend
  real(8), allocatable, dimension(:,:,:) :: RecArr, TmpBuf

  ALLOCATE(RecArr(im, jm, nlev))
  ALLOCATE(TmpBuf(im, jm, nlev))

  k = 1
  nlev_tmp = nlev
  RealCounter = 1

  ! do while(k .le. km)
  do while(RealCounter .le. km)
     call ReadySlave(ReadyProc, RecArr, ComputedLevel, im, jm, nlev_tmp)

     if(ComputedLevel .gt. 0) then
        RealCounter = RealCounter + 1
        do tmpk = ComputedLevel, ComputedLevel+nlev_tmp-1
           do j=1,jm
              do i=1,im
                 ToSend(i,j,tmpk) = RecArr(i,j,1)
              end do
           end do
        end do
     endif
     
     if(k .le. km) then
        do j=1,jm
           do i=1,im
              TmpBuf(i,j,1) = ToSend(i,j,k) !:k+nlev)
           end do
        end do
        
        ! print*, "Rank 0 sending level ", k, " to rank", ReadyProc
        print*, "Sending level ", k, "to process", ReadyProc
        call MPI_Send(TmpBuf, im*jm*nlev, MPI_REAL8, ReadyProc, k, MPI_COMM_WORLD, ierr)
        k = k + nlev
     endif

     ! if(k + nlev .ge. km) then
     !    nlev_tmp = km-k
     ! else        
  end do

  do i=1,size-1
     print*, "killing process ", i
     call MPI_Send(ToSend(:,:,1:nlev), im*jm*nlev, MPI_REAL8, i, km+1, MPI_COMM_WORLD, ierr)
  end do

  DEALLOCATE(RecArr, TmpBuf)

end subroutine Master

subroutine Slave(im, jm, km, nlev, MyRank)
  
  use mpi

  implicit none

  integer :: i, j
  integer :: im, jm, km, nlev, MyRank
  integer :: MyLevel, ierr
  integer :: MyStatus(MPI_STATUS_SIZE)
  real(8), allocatable, dimension(:,:,:) :: MyArr

  ALLOCATE(MyArr(im, jm, nlev))
  MyLevel = 0

  call MPI_Send(MyArr, im*jm*nlev, MPI_REAL8, 0, MyLevel, MPI_COMM_WORLD, ierr)
  
  do while(.true.)
     call MPI_Recv(MyArr, im*jm*nlev, MPI_REAL8, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MyStatus, ierr)
     
     MyLevel = MyStatus(MPI_TAG)

     if(MyLevel .le. km) then
        !
        ! OPERATIONS TO PERFORM
        !
        ! print*, "MyRank slave", MyRank
        do j=1, jm
           do i=1,im
              MyArr(i,j,nlev) = real(MyRank, 8)
           end do
        end do
        
        call MPI_Send(MyArr, im*jm*nlev, MPI_REAL8, 0, MyLevel, MPI_COMM_WORLD, ierr)

     else
        exit
     end if
  end do

  DEALLOCATE(MyArr)

end subroutine Slave
