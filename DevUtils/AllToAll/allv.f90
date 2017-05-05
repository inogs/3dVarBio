program myalltoall

  use mpi
  
  implicit none

  integer :: MyId, NPE, ierr, i, j, iProc, blockNPE
  integer :: GlobalCol, localRow, localCol, RestRow, RestCol, OffsetRow, OffsetCol
  integer, allocatable :: SendCount(:), SendDispl(:), RecCount(:), RecDispl(:)
  real, allocatable :: Buffer(:,:), TmpBuf(:,:), RecBuf(:), DefBuf(:,:), LastBuf(:,:)
  
  call MPI_Init(ierr)
  call MPI_Comm_NPE(MPI_COMM_WORLD, NPE, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyId, ierr)

  write(*,*) "Hello world from process ", MyId, " of ", NPE

  ALLOCATE(SendCount(NPE))
  ALLOCATE(SendDispl(NPE))
  ALLOCATE(RecCount(NPE))
  ALLOCATE(RecDispl(NPE))
  
  ! localRow = 6
  ! localCol = 3
  ! GlobalCol = 9

  ! localRow = 3 !4
  ! localCol = 5
  ! GlobalCol = 10

  localRow = 6
  GlobalCol = 10

  localCol = GlobalCol / NPE
  RestCol  = mod(GlobalCol, NPE)

  if(mod(localRow, NPE) .ne. 0) then
     if(MyId .eq. 0) then
        print*, ''
        print*, 'Warning! localRow % ', NPE,' == ', mod(localRow, NPE)
        print*, ''
     end if
  end if

  RestRow = mod(localRow, NPE)
  blockNPE = localRow/NPE

  if(RestRow .ne. 0) then
     if(MyId .lt. RestRow) then
        blockNPE = blockNPE + 1
     end if
  end if
        
  if(MyId .lt. RestCol) then
     localCol = localCol + 1
  end if

  ! construct SendCount array
  SendDispl(1) = 0
  RecDispl(1)  = 0
  
  do i = 1, NPE
     if(i-1 .lt. RestRow) then
        OffsetRow = 1
     else
        OffsetRow = 0
     end if

     if(i-1 .lt. RestCol) then
        OffsetCol = 1
     else
        OffsetCol = 0
     end if
     
     SendCount(i) = (localRow/NPE + OffsetRow) * localCol
     RecCount(i)  = blockNPE * (GlobalCol/NPE + OffsetCol)
     
     if(i .lt. NPE) then
        SendDispl(i+1) = SendDispl(i) + SendCount(i)
        RecDispl(i+1)  = RecDispl(i)  + RecCount(i)
     end if
  end do

  ! if(MyId .eq. 3) then
  !    print*, SendCount
  !    print*, SendDispl
  !    print*, RecCount
  !    print*, RecDispl
  ! end if

  ALLOCATE(Buffer(localRow, localCol))
  ALLOCATE(TmpBuf(localCol, localRow))
  ALLOCATE(RecBuf(blockNPE*GlobalCol))
  ALLOCATE(DefBuf(blockNPE, GlobalCol))
  ALLOCATE(LastBuf(localRow, localCol))

  OffsetCol = 0
  if(MyId .ge. RestCol) then
     OffsetCol = RestCol
  end if
  
  ! Store initial data
  do i=1,localRow
     do j=1,localCol
        Buffer(i, j) = j + (i - 1) * GlobalCol + MyId*localCol + OffsetCol
     end do
  end do

  TmpBuf = TRANSPOSE(Buffer)
  call MPI_Alltoallv(TmpBuf, SendCount, SendDispl, MPI_FLOAT, RecBuf, RecCount, RecDispl, MPI_FLOAT, MPI_COMM_WORLD, ierr)
  
  ! write(*,*) "MyId = ", MyId, " is starting with:"

  ! print*, ""
  ! print*, Buffer
  ! print*, ""
  ! print*, RecBuf, SHAPE(RecBuf)
  ! print*, ""

  ! ! Reorder data
  do j = 1,blockNPE
     do iProc = 0, NPE-1        
        do i=1, recCount(iProc+1)/blockNPE
           DefBuf(j,i + recDispl(iProc+1)/blockNPE) = RecBuf( RecDispl(iProc+1) + i + (j-1) * recCount(iProc+1)/blockNPE)
        end do
     end do
  end do

  print*, ""
  print*, DefBuf
  print*, ""
  

  DEALLOCATE(RecBuf)
  ALLOCATE(RecBuf(localRow*localCol))
  call MPI_Alltoallv(DefBuf, RecCount, RecDispl, MPI_FLOAT, RecBuf, SendCount, SendDispl, MPI_FLOAT, MPI_COMM_WORLD, ierr)

  ! print*, ""
  ! print*, RecBuf
  ! print*, SHAPE(RecBuf)
  ! print*, ""
  
  do j=1,localCol
     do iProc=0, NPE-1
        do i=1, SendCount(iProc+1)/localCol
           LastBuf(i + SendDispl(iProc+1)/localCol, j) = RecBuf(i + SendDispl(iProc+1) + (j-1) * SendCount(iProc+1)/localCol)
        end do
     end do
  end do
  
  print*, ""
  print*, LastBuf, SHAPE(LastBuf)
  print*, ""

  
  
  DEALLOCATE(Buffer, TmpBuf, RecBuf, DefBuf, LastBuf)
  DEALLOCATE(SendCount, RecCount, SendDispl, RecDispl)
  call MPI_Finalize(ierr)

end program myalltoall
