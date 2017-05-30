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
  INTEGER(i4)        :: idLon, idLat, idLev, idTim
  INTEGER(i4)        :: xaid, yaid, zaid

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=37)  :: BioRestart
  CHARACTER(LEN=6)   :: MyVarName

  real(r8)           :: TmpVal
  real(r8), allocatable, dimension(:,:,:) :: DumpBio
  real(r8), allocatable, dimension(:)     :: VoidArr
  
  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km)); DumpBio(:,:,:) = 1.e20
  ALLOCATE(VoidArr(72)); VoidArr(:) = 42.;

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing bio structure'     
     write(*,*) 'writing bio structure'     
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1

  MyCountSingle(1) = grd%km
  MyStartSingle(1) = 1

  ! check obtained values and eventually
  ! correct them in order to avoid negative concentrations
  do k=1,grd%km
    do j=1,grd%jm
      do i=1,grd%im

        ! value obtained with the assimilation algorithm
        TmpVal = bio%InitialChl(i,j,k) + grd%chl(i,j,k)

        if(TmpVal .lt. 0) then
          ! negative values are not allowed
          ! therefore the correction must be reduced
          do l=1,bio%nphy
            bio%phy(i,j,k,l,1) = 0.01*bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
          enddo
        endif

      enddo
    enddo
  enddo

  do m=1,bio%ncmp
    do l=1,bio%nphy
      iVar = l + bio%nphy*(m-1)

      if(iVar .gt. 17) CYCLE

      BioRestart = 'RESTARTS/RST.'//DA_DATE//'.'//bio%DA_VarList(iVar)//'.nc'

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

      ierr = nf90mpi_def_dim(ncid,'x_a',MyTime ,xaid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
      ierr = nf90mpi_def_dim(ncid,'y_a' ,MyTime ,yaid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)
      MyTime = 1
      ierr = nf90mpi_def_dim(ncid,'z_a'    ,MyTime, zaid)
      
      MyVarName='TRN'//bio%DA_VarList(iVar)

      ierr = nf90mpi_def_var(ncid, MyVarName, nf90_double, (/xid,yid,depid,timeId/), idP )
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
      
      ierr = nf90mpi_def_var(ncid,'nav_lon', nf90_double,  (/xid,yid/), idLon)
      ierr = nf90mpi_def_var(ncid,'nav_lat', nf90_double,  (/xid,yid/), idLat)
      ierr = nf90mpi_def_var(ncid,'nav_lev', nf90_double,  (/depid/)  , idLev)
      ierr = nf90mpi_def_var(ncid,'time'   , nf90_double,  (/timeid/)  , idTim)
      
      ierr = nf90mpi_put_att(ncid,idP   , 'missing_value',1.e+20)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', ierr)

      ierr = nf90mpi_enddef(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef'//bio%DA_VarList(iVar), ierr)

      do k=1,grd%km
        do j=1,grd%jm
          do i=1,grd%im

            ! if(bio%cquot(i,j,k,l,2).gt.MAX_N_CHL .or. bio%cquot(i,j,k,l,4) .and. grd%chl(i,j,k) .gt. 0.) then
            ! endif
            
            if(grd%msk(i,j,k).eq.1) then


              TmpVal = bio%InitialChl(i,j,k) + grd%chl(i,j,k)
              if(TmpVal .lt. 0 .and. m .eq. 1) then
                TmpVal = 0.01*bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
                DumpBio(i,j,k) = TmpVal
                bio%phy(i,j,k,l,1) = TmpVal - bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
              else
                TmpVal = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*(bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,1))
                DumpBio(i,j,k) = TmpVal
              endif

            endif
          enddo
        enddo
      enddo

      ! MyCount(3) = 72

      ierr = nf90mpi_put_var_all(ncid,idP,DumpBio,MyStart,MyCount)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)

      ierr = nf90mpi_put_var_all(ncid,idLon,grd%dx(:,:),MyStart,MyCount)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)
      ierr = nf90mpi_put_var_all(ncid,idLat,grd%dy(:,:),MyStart,MyCount)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)
      ierr = nf90mpi_put_var_all(ncid,idLev,VoidArr,MyStartSingle,MyCountSingle)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)

      MyCountSingle(1) = 1
      ierr = nf90mpi_put_var_all(ncid,idTim,VoidArr,MyStartSingle,MyCountSingle)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//bio%DA_VarList(iVar), ierr)

      ierr = nf90mpi_close(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close '//BioRestart, ierr)
      
    enddo ! l
  enddo ! m

  DEALLOCATE(DumpBio)
  DEALLOCATE(VoidArr)

end subroutine wrt_bio_stat