program myalltoall

  use mpi
  
  implicit none

  integer :: MyRank, Size, ierr, i, j, NData, iProc, blockSize
  integer :: GlobalRow, localRow, localCol
  real, allocatable :: Buffer(:,:), TmpBuf(:,:), RecBuf(:,:), DefBuf(:,:), LastBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)

  write(*,*) "Hello world from process ", MyRank, " of ", Size

  localRow = 6
  localCol = 3
  GlobalRow = 9

  ! localRow = 4
  ! localCol = 5
  ! GlobalRow = 10

  if(mod(localRow, size) .ne. 0) then
     if(MyRank .eq. 0) then
        print*, ''
        print*, 'Warning! localRow % 2 == ', mod(localRow, 2)
        print*, ''
     end if
  end if

  blockSize = localRow/Size

  ALLOCATE(Buffer(localRow, localCol))
  ALLOCATE(TmpBuf(localCol, localRow))
  ALLOCATE(RecBuf(localCol, localRow))
  ALLOCATE(DefBuf(blockSize, localCol*Size))
  ALLOCATE(LastBuf(localRow, localCol))

  ! Store initial data
  do i=1,localRow
     do j=1,localCol
        Buffer(i, j) = j + (i - 1) * GlobalRow + MyRank*localCol
     end do
  end do

  TmpBuf = TRANSPOSE(Buffer)
  NData = localCol * localRow / Size

  ! Perform AllToAll communication
  call MPI_Alltoall(TmpBuf, NData, MPI_FLOAT, RecBuf, NData, MPI_FLOAT, MPI_COMM_WORLD, ierr)

  write(*,*) "MyRank = ", MyRank, " is starting with:"
  ! print*, ""
  ! print*, Buffer
  ! print*, ""
  ! print*, RecBuf
  ! print*, ""

  ! Reorder data
  do j = 1,blockSize
     do iProc = 0, Size-1
        
        do i=1,localCol
           DefBuf(j,i + iProc*localCol) = RecBuf(i, j + iProc*blockSize)
        end do

     end do
  end do

  print*, ""
  print*, DefBuf
  print*, ""

  NData = blockSize * localCol
  RecBuf = RESHAPE(RecBuf, (/blockSize, localCol*Size/))
  call MPI_Alltoall(DefBuf, NData, MPI_FLOAT, RecBuf, NData, MPI_FLOAT, MPI_COMM_WORLD, ierr)

  ! print*, ""
  ! print*, RecBuf
  ! print*, SHAPE(RecBuf)
  ! print*, ""

  do j=1,localCol
     do iProc=0, Size-1
        do i=1,blockSize
           LastBuf(i + iProc*blockSize, j) = RecBuf(i, iProc*localCol + j)
        end do
     end do
  end do

  print*, ""
  print*, LastBuf
  print*, ""

  
  
  DEALLOCATE(Buffer, TmpBuf, RecBuf, DefBuf, LastBuf)

  call MPI_Finalize(ierr)

end program myalltoall
