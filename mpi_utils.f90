subroutine var3d_mpi_init()

#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include "petsc/finclude/petscvec.h"
#else
#include "petsc/finclude/petscvecdef.h"
#endif

  use mpi_str
  use drv_str
  use petscvec

  implicit none

  integer :: ierr, zero
  PetscErrorCode  ::   stat


  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, MyId,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, NPE,ierr)

  zero = 0
  call MPI_Comm_split(MPI_COMM_WORLD, zero, MyId, Var3DCommunicator, ierr)

  ! initialize PETSc environment
  PETSC_COMM_WORLD = Var3DCommunicator
  call PetscInitialize(PETSC_NULL_CHARACTER,stat)
  CHKERRQ(stat)

  if(drv%Verbose .eq. 1) &
       print*, 'PetscInitialize(...) done by MyId ', MyId  

end subroutine var3d_mpi_init

subroutine my_3dvar_node()
  ! ---------------------------------------------------------------------
  !
  !                        routine mynode
  !                      ******************
  !
  !   Purpose :
  !   ---------
  !      Massively parallel processors
  !      Find processor unit
  !
  !    Input :
  !    -----
  !       argument                :
  !
  !    Modifications:
  !    --------------
  !        original  : 93-09 (M. Imbard)
  !        additions : 96-05 (j. Escobar)
  !        additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
  !                           SHMEM and MPI versions
  ! -----------------------------------------------------------------------


  use mpi_str

  ! -----------------------------------------------------------------------

  implicit none

  !
  !  MPI VERSION
  !
  !         -------------
  !         Enroll in MPI
  !         -------------
  !

  INTEGER ierr

  CALL mpi_comm_rank(Var3DCommunicator, MyId,ierr)
  CALL mpi_comm_size(Var3DCommunicator, NPE,ierr)  

  NumProcI = NPE
  NumProcJ = 1

  MyPosI = mod(MyId, NumProcI)
  MyPosJ = MyId / NumProcI

  if(NumProcI .gt. 1) then
     ProcTop  = MyId - 1
     if(ProcTop .lt. 0) ProcTop = MPI_PROC_NULL
     ProcBottom = MyId + 1
     if(ProcBottom .ge. NumProcI) ProcBottom = MPI_PROC_NULL
  else
     print*, ""
     print*, "You are using a single MPI Process!"
     ProcTop    = MPI_PROC_NULL
     ProcBottom = MPI_PROC_NULL
  end if

  call MPI_TYPE_CONTIGUOUS(2, MPI_REAL8, MyPair, ierr)
  call MPI_TYPE_COMMIT(MyPair, ierr)

  ALLOCATE(SendCountX2D(NumProcI), SendCountX3D(NumProcI))
  ALLOCATE(SendCountX3D_chl(NumProcI))
  ALLOCATE(SendDisplX2D(NumProcI), SendDisplX3D(NumProcI))
  ALLOCATE(SendDisplX3D_chl(NumProcI))
  ALLOCATE(RecCountX2D(NumProcI), RecCountX3D(NumProcI))
  ALLOCATE(RecCountX3D_chl(NumProcI))
  ALLOCATE(RecDisplX2D(NumProcI), RecDisplX3D(NumProcI))
  ALLOCATE(RecDisplX3D_chl(NumProcI))

  ! print for debug 
  ! write(*,*) "MyId", MyId, "PosI", MyPosI, "PosJ", MyPosJ, "Left", ProcLeft, "Right", ProcRight, "Top", ProcTop, "Bottom", ProcBottom

  if(NumProcI * NumProcJ .ne. NPE) then
     if(MyId .eq. 0) then
        WRITE(*,*) ""
        WRITE(*,*) " Error: gridX * gridY != nproc "
        WRITE(*,*) " Exit "
        WRITE(*,*) ""
     end if
     call MPI_Abort(Var3DCommunicator, -1, ierr)
  end if

  if(MyId .eq. 0) then
     WRITE(*,*) ' '
     WRITE(*,*) 'Dom_size'
     WRITE(*,*) ' '
     WRITE(*,*) ' number of processors following i : NumProcI   = ', NumProcI
     WRITE(*,*) ' number of processors following j : NumProcJ   = ', NumProcJ
     WRITE(*,*) ' '
  endif

end subroutine my_3dvar_node

subroutine mpi_sync

  use mpi_str

  implicit none

  INTEGER :: ierror

  CALL mpi_barrier(Var3DCommunicator, ierror)

end subroutine mpi_sync


subroutine mpi_stop

#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include "petsc/finclude/petscvec.h"
#else
#include "petsc/finclude/petscvecdef.h"
#endif

  use mpi_str
  use petscvec

  implicit none

  INTEGER info

  integer :: ierr
  PetscErrorCode  ::   stat



  call PetscFinalize(stat)

  call MPI_Comm_free(Var3DCommunicator, ierr)
  call mpi_finalize(info)

end subroutine mpi_stop
