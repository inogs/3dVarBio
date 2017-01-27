program OneSideComm

    ! use mpi

    implicit none
    include "mpif.h"

    INTEGER  :: Size, MyRank, ierr
    INTEGER  :: im, jm, km, i, j, k
    INTEGER  :: mpi_win_obj, nbytes, lenreal
    real*8   :: Test
    real*8, pointer, dimension(:,:,:) :: a

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)

    im = 3
    jm = 4
    km = 1
    lenreal = 8

    nbytes = lenreal*im*jm*km

    print*, "Hello World by", MyRank, "allocating matrix", im,"x",jm,"x",km
    ALLOCATE(a(im, jm, km))

    do k=1,km
        do j=1,jm
            do  i=1,im
                a(i,j,k) = dble(MyRank*100 + j-1 + (i-1)*jm + (k-1)*im*jm)
            enddo
        enddo
    enddo

    call MPI_Win_create(a, nbytes, lenreal, MPI_INFO_NULL, MPI_COMM_WORLD, mpi_win_obj, ierr)

    print*, "MyRank", MyRank, "print matrix:", a


    call MPI_Win_fence(0, mpi_win_obj, ierr)
    if(MyRank .eq. 0) then
        call MPI_Get(Test, 1, MPI_REAL8, MyRank+1, 0, 1, MPI_REAL8, mpi_win_obj, ierr)
        print*, ""
        print*, ""
        print*, "Process", MyRank, "get", Test
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call MPI_Win_free(mpi_win_obj, ierr)
    DEALLOCATE(a)
    call MPI_Finalize(ierr)

end program OneSideComm