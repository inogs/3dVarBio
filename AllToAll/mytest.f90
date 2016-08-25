program myalltoall

  use mpi
  
  implicit none

  integer :: MyRank, Size, ierr, i, j
  integer :: GlobalRow, offset, localRow, localCol
  real, allocatable :: Buffer(:,:), TmpBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)

  write(*,*) "Hello world from process ", MyRank, " of ", Size

  localRow = 4
  localCol = 4
  ALLOCATE(Buffer(localRow, localCol))
  ALLOCATE(TmpBuf(localRow, localCol))

  GlobalRow = 8
  offset = 4
  
  do i=1,localCol
     do j=1,localRow
        Buffer(i, j) = j + (i - 1) * GlobalRow + MyRank*offset
     end do
  end do
  
  ! do i=1,localCol
  !    do j=1,localRow
  !       Buffer(j, i) = j !+ (i - 1) * GlobalRow + MyRank*offset
  !       TmpBuf(j, i) = 0.
  !    end do
  ! end do

  TmpBuf(1,1) = 1.
  
  write(*,*) "MyRank = ", MyRank
  print*, ""
  print*, Buffer
  print*, ""
  
  call MPI_Finalize(ierr)

  DEALLOCATE(Buffer, TmpBuf)
  
end program myalltoall
