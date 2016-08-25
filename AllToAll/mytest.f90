program myalltoall

  use mpi
  
  implicit none

  integer :: MyRank, Size, ierr, i, j
  integer :: GlobalRow, localRow, localCol !, offset
  real, allocatable :: Buffer(:,:), TmpBuf(:,:), RecBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)

  write(*,*) "Hello world from process ", MyRank, " of ", Size

  localRow = 4
  localCol = 5 !4
  ALLOCATE(Buffer(localRow, localCol))
  ! ALLOCATE(TmpBuf(localRow, localCol))
  ALLOCATE(TmpBuf(localCol, localRow))
  ALLOCATE(RecBuf(localRow, localCol))

  GlobalRow = 10 !8
  ! offset = localCol !4

  ! Store initial data
  do i=1,localRow
     do j=1,localCol
        Buffer(i, j) = j + (i - 1) * GlobalRow + MyRank*localCol ! offset
     end do
  end do
  
  do i=1,localCol
     do j=1,localRow
        TmpBuf(i,j) = Buffer(j,i)
     end do
  end do

  write(*,*) "MyRank = ", MyRank
  print*, ""
  print*, Buffer
  print*, ""
  print*, TmpBuf
  print*, ""
  
  call MPI_Finalize(ierr)

  DEALLOCATE(Buffer, TmpBuf)
  
end program myalltoall
