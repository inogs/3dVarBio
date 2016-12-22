program CommunicationTime

  use mpi

  implicit none

  integer :: ierr, MyRank, size, i, j, k
  integer :: im, jm, km, nlev, NRep, Counter
  real(8), pointer, dimension(:,:,:) :: ToSend

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

  NRep = 1000

  im = 722
  jm = 255
  km = 26
  nlev = 1

  ! im = 10
  ! jm = 10
  ! km = 10
  ! nlev = 1
  if(MyRank .eq. 0) then
     ALLOCATE(ToSend(im,jm,km))
     do k=1,km
        do j=1,jm
           do i=1,im
              ToSend(i,j,k) = dble(i + (j-1)*jm + (k-1)*km)
           end do
        end do
     end do


  end if

  ! Master-Slave execution
  do Counter=1, NRep

     if(MyRank .eq. 0) then
        ! print*, "Repetition", Counter
        call Master(ToSend, im, jm, km, nlev, size)
     else
        call Slave(im, jm, km, nlev, MyRank)
     end if

     ! print*, "process ", MyRank, " ending"
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

  end do

  ! Printing (if needed)
  ! if(MyRank .eq. 0) then
  !    do k=1,km
  !       print*, "Printing level", k
  !       ! do j=1,jm
  !          ! do i=1,im
  !          print*, ToSend(:,:,k)
  !          ! end do
  !       ! end do
  !    end do
  ! end if

  call MPI_Finalize(ierr)

end program CommunicationTime

subroutine ReadySlave(WhoIs, Res, ComputedLevel, im, jm, nlev)

  use mpi

  implicit none

  integer :: WhoIs, ComputedLevel, im, jm, nlev, ierr
  integer :: MyStatus(MPI_STATUS_SIZE)
  real(8), dimension(im, jm, nlev) :: Res

  call MPI_Recv(Res, im*jm*nlev, MPI_REAL8, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MyStatus, ierr)
  WhoIs = MyStatus(MPI_SOURCE)
  ComputedLevel = MyStatus(MPI_TAG)

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

        ! print*, "Sending level ", k, "to process", ReadyProc
        call MPI_Send(TmpBuf, im*jm*nlev, MPI_REAL8, ReadyProc, k, MPI_COMM_WORLD, ierr)
        k = k + nlev
      else
        call MPI_Send(ToSend(:,:,1:nlev), im*jm*nlev, MPI_REAL8, ReadyProc, km+1, MPI_COMM_WORLD, ierr)
     endif

     ! if(k + nlev .ge. km) then
     !    nlev_tmp = km-k
     ! else
  end do

  ! do i=1,size-1
  !    ! print*, "killing process ", i
  !    call MPI_Send(ToSend(:,:,1:nlev), im*jm*nlev, MPI_REAL8, i, km+1, MPI_COMM_WORLD, ierr)
  ! end do

  DEALLOCATE(RecArr, TmpBuf)

end subroutine Master

subroutine Slave(im, jm, km, nlev, MyRank)

  use mpi

  implicit none

  integer :: i, j, k
  integer :: im, jm, km, nlev, MyRank
  integer :: MyLevel, ierr
  integer :: MyStatus(MPI_STATUS_SIZE)
  real(8), allocatable, dimension(:,:,:) :: MyArr
  real(8), allocatable, dimension(:,:)   :: Vector

  ALLOCATE(MyArr(im, jm, nlev))
  ALLOCATE(Vector(im, jm))
  MyLevel = 0

  do j=1,jm
     do i=1,im
        Vector(i, j) = real(jm, 8)*0.5-real(j, 8) + 10.*(real(im, 8)*0.5-real(i, 8))
     end do
  end do

  call MPI_Send(MyArr, im*jm*nlev, MPI_REAL8, 0, MyLevel, MPI_COMM_WORLD, ierr)

  do while(.true.)
     call MPI_Recv(MyArr, im*jm*nlev, MPI_REAL8, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MyStatus, ierr)

     MyLevel = MyStatus(MPI_TAG)

     if(MyLevel .le. km) then

        !
        ! OPERATIONS TO PERFORM
        !
        ! print*, "MyRank slave", MyRank
        ! do j=1, jm
        !    do i=1,im
        !       MyArr(i,j,nlev) = real(MyRank, 8)
        !    end do
        ! end do
        do i=1,4
           do k=1, nlev
              MyArr(:,:,k) = matmul(MyArr(:,:,k), Vector)
           end do
        end do

        call MPI_Send(MyArr, im*jm*nlev, MPI_REAL8, 0, MyLevel, MPI_COMM_WORLD, ierr)

     else

        exit

     end if

  end do

  DEALLOCATE(MyArr, Vector)

end subroutine Slave
