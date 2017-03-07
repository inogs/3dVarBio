program OneSideComm

    ! use mpi

    implicit none
    include "mpif.h"

    INTEGER  :: NPE, MyId, ierr
    INTEGER  :: im, jm, km, i, j, k
    INTEGER  :: mpi_win_obj, MyCount
    INTEGER(kind=MPI_ADDRESS_KIND) :: nbytes, lenreal, target_displ
    real*8   :: Test
    
    real*8, pointer, dimension(:,:,:) :: a

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MyId, ierr)
    call MPI_Comm_NPE(MPI_COMM_WORLD, NPE, ierr)

    im = 3
    jm = 4
    km = 1
    lenreal = 8

    MyCount = im*jm*km
    nbytes = lenreal*MyCount

    print*, "Hello World by", MyId, "allocating matrix", im,"x",jm,"x",km

    ALLOCATE(a(im, jm, km))

    do k=1,km
        do j=1,jm
            do  i=1,im
                a(i,j,k) = dble(MyId*100 + j-1 + (i-1)*jm + (k-1)*im*jm)
            enddo
        enddo
    enddo

    call MPI_Win_create(a, nbytes, lenreal, MPI_INFO_NULL, MPI_COMM_WORLD, mpi_win_obj, ierr)

    print*, "MyId", MyId, "print matrix:", a
    
    if(MyId .eq. 0) then
        target_displ = 6
        CALL MPI_WIN_LOCK( MPI_LOCK_EXCLUSIVE, 1, 0, mpi_win_obj, ierr )

        call MPI_Get(Test, 1, MPI_REAL8, MyId+1, target_displ, 1, MPI_REAL8, mpi_win_obj, ierr)

        ! now perform "+=" operation on the array
        target_displ = 0
        call MPI_Accumulate(a, MyCount, MPI_REAL8, 1, target_displ, MyCount, MPI_REAL8, MPI_SUM, mpi_win_obj, ierr)
        
        call MPI_WIN_UNLOCK(1, mpi_win_obj, ierr)
        print*, ""
        print*, ""
        print*, ""
        print*, "Process", MyId, "get", Test
        print*, ""
        print*, ""
        print*, ""
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    print*, "MyId", MyId, "print matrix:", a
    call MPI_Win_free(mpi_win_obj, ierr)

    DEALLOCATE(a)

    call MPI_Finalize(ierr)

end program OneSideComm
