subroutine tao_minimizer

#include "petsc/finclude/petscvecdef.h"

  use drv_str
  use ctl_str
  use mpi_str
  use petscvec

  implicit none

#include "petsc/finclude/petsctao.h"

  PetscErrorCode     ::   ierr
  Tao                ::   tao
  Vec                ::   MyState    ! array that stores the (temporary) state
  PetscInt           ::   n, M, GlobalStart, MyEnd, iter
  PetscReal          ::   fval, gnorm, cnorm, xdiff
  PetscScalar        ::   MyTolerance
  TaoConvergedReason ::   reason
  integer(i4)        ::   j
  real(8)            ::   MaxGrad

  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: MyValues
  PetscScalar, pointer                    :: xtmp(:)

  external MyFuncAndGradient

  if(MyId .eq. 0) then
     print*,'PETSc-TAO lmvm minimizer configuration'
     print*, ''
  endif

  if(MyId .eq. 0) then
     write(drv%dia,*) ''
     write(drv%dia,*) "Within tao_minimizer subroutine!"
  endif

  ! Allocate working arrays
  ! TAO needs to know both n (the local size of the array)
  ! and M (global size of the array)
  n = ctl%n
  M = ctl%n_global

  ALLOCATE(loc(n), MyValues(n))

  ! Create MyState array and fill it
  ! VecGetOwnergshipRange returns the indices related to the local
  ! section of the PETSc Vector (required in order to fill the initial state vector,
  ! further informations can be found in the PETSc user manual, Section 2.1 and 2.2)
  call VecCreateMPI(Var3DCommunicator, n, M, MyState, ierr)
  call VecGetOwnershipRange(MyState, GlobalStart, MyEnd, ierr)

  if(drv%Verbose .eq. 1) &
       print*, "MyState initialization by MyId ", MyId, "with indices: ", GlobalStart, MyEnd

  if( ctl%n .ne. MyEnd - GlobalStart ) then
     print*, ""
     print*, "WARNING!!"
     print*, "ctl%n .ne. GlobalStart - MyEnd"
     print*, "ctl%n = ", ctl%n
     print*, "GlobalStart = ", GlobalStart
     print*, "MyEnd = ", MyEnd
     print*, ""
  endif

  ! Setting initial state vector to zero
  do j = 1, ctl%n
     loc(j) = GlobalStart + j - 1
     MyValues(j) = 0.
  end do

  ! Setting only local values (since each process can access at all entries of MyState)
  call VecSetValues(MyState, ctl%n, loc, MyValues, INSERT_VALUES, ierr)
  call VecAssemblyBegin(MyState, ierr)
  call VecAssemblyEnd(MyState, ierr)

  ! Iteration counter initialization
  drv%MyCounter = 0

  ! Create Tao object and set type LMVM (ones that use BFGS minimization algorithm)
  call TaoCreate(Var3DCommunicator, tao, ierr)
  CHKERRQ(ierr)
  call TaoSetType(tao,"lmvm",ierr)
  CHKERRQ(ierr)

  ! Set initial solution array, MyBounds and MyFuncAndGradient routines
  call TaoSetInitialVector(tao, MyState, ierr)
  CHKERRQ(ierr)
  call TaoSetObjectiveAndGradientRoutine(tao, MyFuncAndGradient, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)

  ! Calling costf in order to compute
  ! the initial gradient. This will be used to
  ! set MyTolerance
  call costf
  MaxGrad = 0
  do j=1,ctl%n
     MaxGrad = max(MaxGrad, abs(ctl%g_c(j)))
  end do

  ! since MaxGrad is the local maximum value, the MPI_Allreduce call
  ! is required to find the global maximum value
  call MPI_Allreduce(MPI_IN_PLACE, MaxGrad, 1, MPI_REAL8, MPI_MAX, Var3DCommunicator, ierr)
  MyTolerance = ctl%pgper * MaxGrad

  if(MyId .eq. 0) then
     print*, "Setting MyTolerance", MyTolerance
     print*, ""
     write(drv%dia,*) "Setting MyTolerance", MyTolerance
  endif

  ! setting tolerances
  call TaoSetTolerances(tao, MyTolerance, 1.d-4, ctl%pgper, ierr)
  CHKERRQ(ierr)

  ! calling the solver to minimize the problem
  call TaoSolve(tao, ierr)
  CHKERRQ(ierr)

  ! printing solver info
  if(MyId .eq. 0) then
     print*, ''
     print*, 'Tao Solver Info:'
     print*, ''
  endif

  call TaoView(tao, PETSC_VIEWER_STDOUT_WORLD, ierr)

  call TaoGetSolutionStatus(tao, iter, fval, gnorm, cnorm, xdiff, reason, ierr)
  if(reason .lt. 0) then
    if(MyId .eq. 0) then
      print*, "TAO failed to find a solution"
      print*, "Aborting.."
    endif
    call MPI_Barrier(Var3DCommunicator, ierr)
    call MPI_Abort(Var3DCommunicator, -1, ierr)
  endif

  ! Get the solution and copy into ctl%x_c array
  call TaoGetSolutionVector(tao, MyState, ierr)
  CHKERRQ(ierr)
  call VecGetArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  do j = 1, ctl%n
     ctl%x_c(j) = xtmp(j)
  end do

  call VecRestoreArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  ! Deallocating variables
  DEALLOCATE(loc, MyValues)

  ! free memory
  call TaoDestroy(tao, ierr)
  CHKERRQ(ierr)

  call VecDestroy(MyState, ierr)
  CHKERRQ(ierr)

  if(MyId .eq. 0) then
     write(drv%dia,*) 'Minimization done with ', drv%MyCounter
     write(drv%dia,*) 'iterations'
     write(drv%dia,*) ''

     print*, ""
     print*, "Minimization done with ", drv%MyCounter, "iterations"
     print*, ""
  endif

end subroutine tao_minimizer

!-------------------------------------------------!
! subroutine that, given a state                  !
! MyState (provided by Tao solver),               !
! computes the value of cost function             !
! in MyState and its gradient                     !
!-------------------------------------------------!

subroutine MyFuncAndGradient(tao, MyState, CostFunc, Grad, dummy, ierr)

#include "petsc/finclude/petscvecdef.h"

  use set_knd
  use drv_str
  use ctl_str
  use petscvec
  use mpi_str

  implicit none

#include "petsc/finclude/petsctao.h"

  Tao             ::   tao
  Vec             ::   MyState, Grad
  PetscScalar     ::   CostFunc
  PetscErrorCode  ::   ierr
  integer(i4)     ::   dummy, j

  ! Working arrays
  PetscScalar, pointer, dimension(:)  :: my_grad
  PetscScalar, pointer, dimension(:)  :: xtmp

  ! read temporary state provided by Tao Solver
  ! and set it in ctl%x_c array in order to compute
  ! the actual value of Cost Function and the gradient
  ! with costf subroutine
  ! VecGetArrayReadF90 function puts MyState (provided by TAO)
  ! into xtmp pointer -> Now we can access to MyState values :)
  call VecGetArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  ! access to MyState values
  do j=1,ctl%n
     ctl%x_c(j) = xtmp(j)
  end do

  call VecRestoreArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  ! compute function and gradient
  call costf

  ! assign the Cost Function value computed by costf to CostFunc
  CostFunc = ctl%f_c

  call VecGetArrayF90(Grad, my_grad, ierr)

  do j = 1, ctl%n
     my_grad(j) = ctl%g_c(j)
  end do

  call VecRestoreArrayF90(Grad, my_grad, ierr)

  ! Update counter
  drv%MyCounter = drv%MyCounter + 1

  ! Exit without errors
  ierr = 0

end subroutine MyFuncAndGradient
