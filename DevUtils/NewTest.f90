PROGRAM mpitest

IMPLICIT NONE
include "mpif.h"

INTEGER :: ierr, npe, rnk, win
INTEGER (KIND=MPI_ADDRESS_KIND) lowerbound, NPEofreal, MyDispl
REAL :: val = 1.0, oval = 2.0

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rnk, ierr )
CALL MPI_COMM_NPE( MPI_COMM_WORLD, npe, ierr )

CALL MPI_TYPE_GET_EXTENT(MPI_REAL, lowerbound, NPEofreal, ierr)

CALL MPI_WIN_CREATE(val, NPEofreal, NPEofreal, MPI_INFO_NULL, MPI_COMM_WORLD, win, ierr)
call MPI_Win_fence(0, win, ierr)
IF( rnk .EQ. 1 ) THEN
    MyDispl = 0
    lowerbound = 0
   print*,oval
   CALL MPI_WIN_LOCK( MPI_LOCK_SHARED, 0, 0, win, ierr )
   CALL MPI_GET( oval, 1, MPI_REAL, lowerbound, 0, 1, MPI_REAL, win, ierr )
   CALL MPI_WIN_UNLOCK( 0, win, ierr )
   print*,oval
END IF

CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
CALL MPI_WIN_FREE(win, ierr)
CALL MPI_FINALIZE(ierr)

END PROGRAM mpitest