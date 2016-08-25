program myalltoall

  use mpi
  
  implicit none

  integer :: MyRank, Size, ierr, i, j, NData, iProc, blockSize
  integer :: GlobalRow, localRow, localCol !, offset
  real, allocatable :: Buffer(:,:), TmpBuf(:,:), RecBuf(:,:), DefBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)

  write(*,*) "Hello world from process ", MyRank, " of ", Size

  localRow = 4
  localCol = 5 !4
  GlobalRow = 10 !8

  if(mod(localRow, 2) .ne. 0) then
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

  ! offset = localCol !4

  ! Store initial data
  do i=1,localRow
     do j=1,localCol
        Buffer(i, j) = j + (i - 1) * GlobalRow + MyRank*localCol ! offset
     end do
  end do

  ! Set up array for AllToAll 
  ! do i=1,localCol
  !    do j=1,localRow
  !       TmpBuf(i,j) = Buffer(j,i)
  !    end do
  ! end do
  TmpBuf = TRANSPOSE(Buffer)
  
  NData = localCol * localRow / 2
  call MPI_Alltoall(TmpBuf, NData, MPI_FLOAT, RecBuf, NData, MPI_FLOAT, MPI_COMM_WORLD, ierr)
  
  write(*,*) "MyRank = ", MyRank, " is starting with:"
  print*, ""
  print*, Buffer
  print*, ""
  print*, RecBuf
  print*, ""

  do j = 1,blockSize
     do iProc = 0, Size-1
        
        do i=1,localCol
           DefBuf(j,i + iProc*localCol) = RecBuf(i, j + iProc*blockSize)
        end do

     end do
  end do

  print*, DefBuf
  print*, ""
  print*, "Have a nice day from proc ", MyRank
  print*, ""

  ! write(*,*) "MyRank = ", MyRank
  ! print*, ""
  ! print*, Buffer
  ! print*, ""
  ! print*, TmpBuf
  ! print*, ""
  
  DEALLOCATE(Buffer, TmpBuf, RecBuf, DefBuf)

  call MPI_Finalize(ierr)

end program myalltoall
