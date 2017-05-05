program myalltoall

  use mpi
  
  implicit none

  integer :: MyId, NPE, ierr, i, j
  integer :: GlobalRow, offset, localRow, localCol
  real, allocatable :: Buffer(:,:), TmpBuf(:) !,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_NPE(MPI_COMM_WORLD, NPE, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyId, ierr)

  write(*,*) "Hello world from process ", MyId, " of ", NPE

  localRow = 4
  localCol = 4
  ALLOCATE(Buffer(localRow, localCol))
  ALLOCATE(TmpBuf(localRow))

  GlobalRow = 8
  offset = 4
  
  do i=1,localCol
     do j=1,localRow
        Buffer(j, i) = j
     end do
  end do
  
  TmpBuf(:) = 0.
  TmpBuf(1) = 1.
  
  write(*,*) "MyId = ", MyId
  print*, ""
  print*, Buffer
  print*, ""
  print*, TmpBuf
  print*, ""
  print*, "Buffer*TmpBuf = ", MatMul(Buffer,TmpBuf)
  
  call MPI_Finalize(ierr)

  DEALLOCATE(Buffer, TmpBuf)
  
end program myalltoall
