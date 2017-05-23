subroutine wrt_bio_stat

  use set_knd
  use grd_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use filenames

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l, m
  INTEGER(i4)        :: idP, iVar
  INTEGER(I4)        :: xid,yid,depid,timeId

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  CHARACTER(LEN=37)  :: BioRestart
  CHARACTER(LEN=6)   :: MyVarName

  real(r4), allocatable, dimension(:,:,:) :: DumpBio
  
  ALLOCATE(DumpBio(grd%im,grd%jm,72)); DumpBio(:,:,:) = 0.

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing bio structure'     
     write(*,*) 'writing bio structure'     
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = 72 ! grd%km
  MyTime = 1

  do m=1,bio%ncmp
    do l=1,bio%nphy
      iVar = l + bio%nphy*(m-1)

      if(iVar .gt. 17) CYCLE

      BioRestart = 'RESTARTS/RST.'//DA_DATE//'.'//bio%DA_VarList(iVar)//'.nc'

      if(drv%Verbose .eq. 1 .and. MyId .eq. 0) &
        print*, "Writing BioRestart ", BioRestart
      
      ierr = nf90mpi_create(Var3DCommunicator, BioRestart, NF90_CLOBBER, MPI_INFO_NULL, ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_create '//BioRestart, ierr)

      ierr = nf90mpi_def_dim(ncid,'z'    ,global_km, depid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim depth ', ierr)
      ierr = nf90mpi_def_dim(ncid,'y' ,global_jm ,yid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)
      ierr = nf90mpi_def_dim(ncid,'x',global_im ,xid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
      ierr = nf90mpi_def_dim(ncid,'time',MyTime ,timeId)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim time ', ierr)
      
      MyVarName='TRN'//bio%DA_VarList(iVar)

      ierr = nf90mpi_def_var(ncid, MyVarName, nf90_float, (/xid,yid,depid,timeId/), idP )
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)

      ierr = nf90mpi_enddef(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var'//bio%DA_VarList(iVar), ierr)

      do k=1,grd%km
        do j=1,grd%jm
          do i=1,grd%im
            DumpBio(i,j,k) = REAL(bio%phy(i,j,k,l,m), 4)
          enddo
        enddo
      enddo

      ierr = nf90mpi_put_var_all(ncid,idP,DumpBio,MyStart,MyCount)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)

      ierr = nf90mpi_close(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close '//BioRestart, ierr)
      
    enddo ! l
  enddo ! m

  DEALLOCATE(DumpBio)

end subroutine wrt_bio_stat