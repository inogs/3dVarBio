program OneSideComm

    ! use mpi

    implicit none
    include "mpif.h"

    INTEGER  :: Size, MyRank, ierr
    INTEGER  :: im, jm, km, i, j, k
    INTEGER  :: mpi_win_obj, MyCount
    INTEGER(kind=MPI_ADDRESS_KIND) :: nbytes, lenreal, target_displ
    real*8   :: Test
    
    real*8, pointer, dimension(:,:,:) :: a
    ! real a
    ! pointer (P, a(3,4,1))

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, Size, ierr)

    im = 3
    jm = 4
    km = 1
    lenreal = 8

    MyCount = im*jm*km
    nbytes = lenreal*MyCount

    print*, "Hello World by", MyRank, "allocating matrix", im,"x",jm,"x",km

    ALLOCATE(a(im, jm, km))
    ! call MPI_Alloc_mem(nbytes, MPI_INFO_NULL, P, ierr)

    do k=1,km
        do j=1,jm
            do  i=1,im
                a(i,j,k) = dble(MyRank*100 + j-1 + (i-1)*jm + (k-1)*im*jm)
            enddo
        enddo
    enddo

    call MPI_Win_create(a, nbytes, lenreal, MPI_INFO_NULL, MPI_COMM_WORLD, mpi_win_obj, ierr)
    call MPI_Win_fence(0, mpi_win_obj, ierr)

    print*, "MyRank", MyRank, "print matrix:", a
    
    if(MyRank .eq. 0) then
        target_displ = 0
        CALL MPI_WIN_LOCK( MPI_LOCK_SHARED, 0, 0, mpi_win_obj, ierr )

        call MPI_Get(Test, 1, MPI_REAL8, MyRank+1, target_displ, 1, MPI_REAL8, mpi_win_obj, ierr)
        ! call MPI_Accumulate(a, MyCount, MPI_REAL8, 1, target_displ, MyCount, MPI_REAL8, MPI_SUM, mpi_win_obj, ierr)
        
        call MPI_WIN_UNLOCK(0, mpi_win_obj, ierr)
        print*, ""
        print*, ""
        print*, ""
        print*, "Process", MyRank, "get", Test
        print*, ""
        print*, ""
        print*, ""
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    print*, "MyRank", MyRank, "print matrix:", a
    
    call MPI_Win_free(mpi_win_obj, ierr)

    ! call MPI_Free_mem(a, ierr)
    DEALLOCATE(a)

    call MPI_Finalize(ierr)

end program OneSideComm