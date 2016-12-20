program vectest

#include "petsc/finclude/petscvecdef.h"
  use petscvec
  use mpi

  implicit none

  Vec MyState
  PetscErrorCode :: ierr
  PetscInt :: GlobalStart, MyEnd 
  integer :: state, size, MyRank, n, M, j
  integer, allocatable :: loc(:)
  PetscScalar, allocatable :: MyValues(:)
  PetscScalar, pointer, dimension(:)  :: TmpPtr
  real(8) :: zero

  call MPI_Init(state)
  call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, state)
  call MPI_Comm_size(MPI_COMM_WORLD, size, state)

  print*, "MPI_Init done by rank ", MyRank

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  print*, "PetscInitialize done by rank ", MyRank

  n = 5
  M = n*size

  ALLOCATE(loc(n), MyValues(n))

  call VecCreateMPI(MPI_COMM_WORLD, n, M, MyState, ierr)
  call VecGetOwnershipRange(MyState, GlobalStart, MyEnd, ierr)

  zero = 0.
  zero = zero + MyRank*n

  do j = 1, n
     zero = zero + 1
     loc(j) = GlobalStart + j - 1
     MyValues(j) = zero
  end do

  call VecSetValues(MyState, n, loc, MyValues, INSERT_VALUES, ierr)
  call VecAssemblyBegin(MyState, ierr)
  call VecAssemblyEnd(MyState, ierr)

  call VecView(MyState, PETSC_VIEWER_STDOUT_SELF, ierr)

  call VecGetOwnershipRange(MyState, GlobalStart, MyEnd, ierr)
  call VecGetArrayF90(MyState, TmpPtr, ierr)

  do j=1, n
     TmpPtr(j) = TmpPtr(j) - 1
  end do

  call VecRestoreArrayF90(MyState, TmpPtr, ierr)

  call VecView(MyState, PETSC_VIEWER_STDOUT_SELF, ierr)

  call PetscFinalize(ierr)
  call MPI_Finalize(state)
  
end program vectest
