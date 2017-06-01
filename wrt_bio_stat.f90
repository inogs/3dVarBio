subroutine wrt_bio_stat

  use set_knd
  use grd_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use da_params

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l, m
  INTEGER(i4)        :: idP, iVar
  INTEGER(I4)        :: xid,yid,depid,timeId, idTim

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=37)  :: BioRestart
  CHARACTER(LEN=6)   :: MyVarName

  real(r8)           :: TmpVal
  real(r8), allocatable, dimension(:,:,:) :: DumpBio
  real(r8) :: TimeArr(1)
  
  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km)); DumpBio(:,:,:) = 1.e20

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing bio structure'     
     write(*,*) 'writing bio structure'     
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1

  MyCountSingle(1) = 1
  MyStartSingle(1) = 1
  TimeArr(1) = DA_JulianDate

  do m=1,bio%ncmp
    do l=1,bio%nphy
      iVar = l + bio%nphy*(m-1)

      if(iVar .gt. NBioVar) CYCLE

      BioRestart = 'RESTARTS/RST.'//DA_DATE//'.'//DA_VarList(iVar)//'.nc'

      if(drv%Verbose .eq. 1 .and. MyId .eq. 0) &
        print*, "Writing BioRestart ", BioRestart
      
      ierr = nf90mpi_create(Var3DCommunicator, BioRestart, NF90_CLOBBER, MPI_INFO_NULL, ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_create '//BioRestart, ierr)

      ierr = nf90mpi_def_dim(ncid,'x',global_im ,xid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
      ierr = nf90mpi_def_dim(ncid,'y' ,global_jm ,yid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)
      ierr = nf90mpi_def_dim(ncid,'z'    ,global_km, depid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim depth ', ierr)
      ierr = nf90mpi_def_dim(ncid,'time',MyTime ,timeId)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim time ', ierr)

      MyVarName='TRN'//DA_VarList(iVar)

      ierr = nf90mpi_def_var(ncid, MyVarName, nf90_double, (/xid,yid,depid,timeId/), idP )
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
      
      ierr = nf90mpi_def_var(ncid,'time'   , nf90_double,  (/timeid/)  , idTim)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
      ierr = nf90mpi_put_att(ncid,idP   , 'missing_value',1.e+20)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', ierr)

      ierr = nf90mpi_enddef(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef'//DA_VarList(iVar), ierr)

      do k=1,grd%km
        do j=1,grd%jm
          do i=1,grd%im

            if(bio%InitialChl(i,j,k) .lt. 1.e20) then

              if(grd%msk(i,j,k).eq.1) then

                ! check obtained values and eventually
                ! correct them in order to avoid negative concentrations
                ! if the correction is negative, the correction must be reduced
                TmpVal = bio%InitialChl(i,j,k) + grd%chl(i,j,k)
                if(TmpVal .lt. 0 .and. m .eq. 1) then
                  TmpVal = 0.01*bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
                  DumpBio(i,j,k) = TmpVal
                  bio%phy(i,j,k,l,1) = TmpVal - bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
                else
                  TmpVal = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*(bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,1))
                  DumpBio(i,j,k) = TmpVal
                endif
              else
                DumpBio(i,j,k) = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
              endif

            endif
          enddo
        enddo
      enddo

      ierr = nf90mpi_put_var_all(ncid,idP,DumpBio,MyStart,MyCount)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

      ierr = nf90mpi_put_var_all(ncid,idTim,TimeArr,MyStartSingle,MyCountSingle)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

      ierr = nf90mpi_close(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close '//BioRestart, ierr)
      
    enddo ! l
  enddo ! m

  DEALLOCATE(DumpBio)

end subroutine wrt_bio_stat