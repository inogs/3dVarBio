MODULE myalloc_mpi
  
  IMPLICIT NONE
  
  public
  
  
#include <mpif.h>
  
  !-----------------------------------------------!
  !     MPI vaiables
  !
  !     size : number of process
  !     rank : process number  [ 0 - size-1 ]
  !-----------------------------------------------!
  
  INTEGER  :: size,rank
  INTEGER  :: mpi_pack_size =18
  
  
END MODULE myalloc_mpi
