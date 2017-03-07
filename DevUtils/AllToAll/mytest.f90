program myalltoall

  use mpi
  
  implicit none

  integer :: MyId, NPE, ierr, i, j, NData, iProc, blockNPE
  integer :: GlobalRow, localRow, localCol
  real, allocatable :: Buffer(:,:), TmpBuf(:,:), RecBuf(:,:), DefBuf(:,:), LastBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_NPE(MPI_COMM_WORLD, NPE, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyId, ierr)

  write(*,*) "Hello world from process ", MyId, " of ", NPE

  localRow = 6
  localCol = 3
  GlobalRow = 9

  ! localRow = 4
  ! localCol = 5
  ! GlobalRow = 10

  if(mod(localRow, NPE) .ne. 0) then
     if(MyId .eq. 0) then
        print*, ''
        print*, 'Warning! localRow % 2 == ', mod(localRow, 2)
        print*, ''
     end if
  end if

  blockNPE = localRow/NPE

  ALLOCATE(Buffer(localRow, localCol))
  ALLOCATE(TmpBuf(localCol, localRow))
  ALLOCATE(RecBuf(localCol, localRow))
  ALLOCATE(DefBuf(blockNPE, localCol*NPE))
  ALLOCATE(LastBuf(localRow, localCol))

  ! Store initial data
  do i=1,localRow
     do j=1,localCol
        Buffer(i, j) = j + (i - 1) * GlobalRow + MyId*localCol
     end do
  end do

  TmpBuf = TRANSPOSE(Buffer)
  NData = localCol * localRow / NPE

  ! Perform AllToAll communication
  call MPI_Alltoall(TmpBuf, NData, MPI_FLOAT, RecBuf, NData, MPI_FLOAT, MPI_COMM_WORLD, ierr)

  write(*,*) "MyId = ", MyId, " is starting with:"
  ! print*, ""
  ! print*, Buffer
  ! print*, ""
  ! print*, RecBuf
  ! print*, ""

  ! Reorder data
  do j = 1,blockNPE
     do iProc = 0, NPE-1
        
        do i=1,localCol
           DefBuf(j,i + iProc*localCol) = RecBuf(i, j + iProc*blockNPE)
        end do

     end do
  end do

  print*, ""
  print*, DefBuf
  print*, ""

  NData = blockNPE * localCol
  RecBuf = RESHAPE(RecBuf, (/blockNPE, localCol*NPE/))
  call MPI_Alltoall(DefBuf, NData, MPI_FLOAT, RecBuf, NData, MPI_FLOAT, MPI_COMM_WORLD, ierr)

  ! print*, ""
  ! print*, RecBuf
  ! print*, SHAPE(RecBuf)
  ! print*, ""

  do j=1,localCol
     do iProc=0, NPE-1
        do i=1,blockNPE
           LastBuf(i + iProc*blockNPE, j) = RecBuf(i, iProc*localCol + j)
        end do
     end do
  end do

  print*, ""
  print*, LastBuf
  print*, ""

  
  
  DEALLOCATE(Buffer, TmpBuf, RecBuf, DefBuf, LastBuf)

  call MPI_Finalize(ierr)

end program myalltoall
