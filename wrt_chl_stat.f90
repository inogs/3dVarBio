subroutine wrt_chl_stat

  use set_knd
  use grd_str
  use eof_str
  use drv_str
  use mpi_str
  use bio_str
  use pnetcdf
  use da_params

  implicit none

  INTEGER(i4)        :: ncid, ierr, i, j, k, l, m, mm, my_km
  INTEGER(i4)        :: idP, iVar, idL
  INTEGER(I4)        :: xid,yid,depid,timeId, idTim
  INTEGER            :: system, SysErr

  INTEGER(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km, MyTime
  INTEGER(KIND=MPI_OFFSET_KIND) :: MyCountSingle(1), MyStartSingle(1)
  CHARACTER(LEN=45)    :: BioRestart
  CHARACTER(LEN=47)    :: BioRestartLong
  CHARACTER(LEN=6)     :: MyVarName
  CHARACTER(LEN=16)    :: LimVarName
  CHARACTER(LEN=49)    :: LimCorrfile
  LOGICAL, ALLOCATABLE :: MyConditions(:,:,:,:)
  INTEGER, ALLOCATABLE :: LimitCorr(:,:,:)


  real(r8)           :: TmpVal, MyCorr, MyRatio, SMALL
  real(r4), allocatable, dimension(:,:,:) :: ValuesToTest

  ! bug fix Intel 2018
  real(r4), allocatable, dimension(:,:,:,:) :: DumpBio
  integer(KIND=MPI_OFFSET_KIND) :: MyStart_4d(4), MyCount_4d(4)

  real(r8) :: TimeArr(1)
  real(r4) :: MAX_N_CHL, MAX_P_CHL, MAX_P_C, MAX_N_C
  real(r4) :: OPT_N_C, OPT_P_C, OPT_S_C, LIM_THETA

  MAX_N_CHL =  150.        ! Derived from max chl:c=0.02 (BFMconsortium)
  MAX_P_CHL =  10.
  MAX_P_C   =  7.86e-4*2   ! values from BFMconsortium parametrs document (P.Lazzari)
  OPT_P_C   =  7.86e-4
  MAX_N_C   =  1.26e-2*2   ! values from BFMconsortium parametrs document (P.Lazzari)
  OPT_N_C   =  1.26e-2
  OPT_S_C   =  0.01        ! values from BFMconsortium parametrs document (P.Lazzari)
  LIM_THETA =  0.001
  SMALL     =  1.e-5

  MyStart_4d(1:3) = MyStart(:)
  MyStart_4d(4) = 1
  MyCount_4d(1:3) = MyCount(:)
  MyCount_4d(4) = 1

  my_km = grd%km
  if(drv%multiv.eq.1) &
    my_km = ros%kmchl

  ALLOCATE(DumpBio(grd%im,grd%jm,grd%km,1)); DumpBio(:,:,:,:) = 1.e20
  ALLOCATE(ValuesToTest(grd%im,grd%jm,my_km)); ValuesToTest(:,:,:) = dble(0.)
  ALLOCATE(MyConditions(grd%im,grd%jm,my_km,bio%nphy))
  ALLOCATE(LimitCorr(grd%im,grd%jm,grd%km)); LimitCorr(:,:,:) = -99

  if(MyId .eq. 0) then
     write(drv%dia,*) 'writing chl structure'     
     write(*,*) 'writing chl structure'     
  endif

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  MyTime = 1

  MyCountSingle(1) = 1
  MyStartSingle(1) = 1
  TimeArr(1) = DA_JulianDate

  LimitCorr(:,:,1:my_km) = 0

  do k=1,my_km
    do j=1,grd%jm
      do i=1,grd%im

        if(bio%InitialChl(i,j,k) .lt. 1.e20) then
          ! check obtained values and eventually
          ! correct them in order to avoid negative concentrations
          ! if the correction is negative, the correction must be reduced
          ValuesToTest(i,j,k) = bio%InitialChl(i,j,k) + grd%chl(i,j,k)
          if(bio%ApplyConditions) then
            if(ValuesToTest(i,j,k) .gt. 10*bio%InitialChl(i,j,k)) then
              if(ValuesToTest(i,j,1) .gt. 0) then 
                LimitCorr(i,j,k) = 1
                ! print*, ValuesToTest(i,j,1), bio%InitialChl(i,j,1), grd%chl(i,j,1),i,j
              endif
 
              do m=1,bio%ncmp
                do l=1,bio%nphy
                  bio%phy(i,j,k,l,m) = 9.*bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
                enddo
              enddo

            endif
            
            ! limitations in case of high nutrient contents
            do l=1,bio%nphy
              MyConditions(i,j,k,l) = bio%cquot(i,j,k,l,3) .gt. MAX_N_CHL
              MyConditions(i,j,k,l) = MyConditions(i,j,k,l) .or. (bio%cquot(i,j,k,l,4) .gt. MAX_P_CHL)
              MyConditions(i,j,k,l) = MyConditions(i,j,k,l) .or. (bio%cquot(i,j,k,l,3)/bio%cquot(i,j,k,l,2) .gt. (4*MAX_N_C))
              MyConditions(i,j,k,l) = MyConditions(i,j,k,l) .or. (bio%cquot(i,j,k,l,4)/bio%cquot(i,j,k,l,2) .gt. (4*MAX_P_C))
            enddo
          endif
        endif
      enddo
    enddo
  enddo



  do m=1,bio%ncmp
    do l=1,bio%nphy
      iVar = l + bio%nphy*(m-1)

      if(iVar .gt. NPhytoVar) CYCLE

      BioRestart = 'DA__FREQ_1/RST_after.'//ShortDate//'.'//DA_VarList(iVar)//'.nc'
      BioRestartLong = 'DA__FREQ_1/RST_after.'//DA_DATE//'.'//DA_VarList(iVar)//'.nc'

      if(drv%Verbose .eq. 1 .and. MyId .eq. 0) &
        print*, "Writing Phyto Restart ", BioRestart
      
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

      ierr = nf90mpi_def_var(ncid, MyVarName, nf90_float, (/xid,yid,depid,timeId/), idP )
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
      
      ierr = nf90mpi_def_var(ncid,'time'   , nf90_double,  (/timeid/)  , idTim)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)
      ierr = nf90mpi_put_att(ncid,idP   , 'missing_value',1.e+20)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', ierr)

      ierr = nf90mpi_enddef(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef'//DA_VarList(iVar), ierr)

      do k=1,my_km
        do j=1,grd%jm
          do i=1,grd%im

            if(bio%InitialChl(i,j,k) .lt. 1.e20) then

              if(grd%msk(i,j,k).eq.1) then

                if(ValuesToTest(i,j,k) .lt. 0 .and. m .eq. 1) then
                  ! Excluding negative concentrations
                  ! This correction must be the first
                  ! condition applied (before apply corrections
                  ! on the other components)
                  TmpVal = 0.01*bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
                  if(TmpVal.gt.SMALL) then
                    TmpVal = SMALL
                  endif
                  DumpBio(i,j,k,1) = TmpVal
                  LimitCorr(i,j,k) = 2

                  ! the positiveness is applied to
                  ! the other components
                  bio%phy(i,j,k,l,1) = TmpVal - bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k)
                  do mm=2,bio%ncmp
                    bio%phy(i,j,k,l,mm) = bio%cquot(i,j,k,l,mm)*bio%phy(i,j,k,l,1)
                  enddo

                else

                  if(bio%ApplyConditions) then

                    if(bio%phy(i,j,k,l,m) .gt. 0 .and. MyConditions(i,j,k,l)) then
                      bio%phy(i,j,k,l,m) = 0.
                    endif

                    ! limitation on Carbon corrections
                    ! when chl/Carbon ratio is small
                    if(m .eq. 2) then
                      MyRatio = 1./bio%cquot(i,j,k,l,m)
                      if(MyRatio .lt. LIM_THETA .and. bio%phy(i,j,k,l,m) .gt. 0) then
                        MyCorr = bio%pquot(i,j,k,l)*bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,1)
                        MyCorr = MyCorr/LIM_THETA - bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
                        bio%phy(i,j,k,l,m) = max(0., MyCorr)
                        LimitCorr(i,j,k) = 3
                      endif
                    endif

                    ! limitation on Nitrogen corrections
                    ! to the optimal N/C ratio
                    if(m .eq. 3) then
                      ! compute N/C fraction
                      MyRatio = bio%cquot(i,j,k,l,m)/bio%cquot(i,j,k,l,2)
                      if(MyRatio .gt. OPT_N_C .and. bio%phy(i,j,k,l,m) .gt. 0) then
                        MyCorr = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,2)*bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,2)
                        MyCorr = MyCorr*OPT_N_C - bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
                        bio%phy(i,j,k,l,m) = max(0., MyCorr)
                        LimitCorr(i,j,k) = 4
                      endif

                    endif

                    ! limitation on Phosphorus corrections
                    ! to the optimal P/C ratio
                    if(m .eq. 4) then
                      ! compute P/C fraction
                      MyRatio = bio%cquot(i,j,k,l,m)/bio%cquot(i,j,k,l,2)
                      if(MyRatio .gt. OPT_P_C .and. bio%phy(i,j,k,l,m) .gt. 0) then
                        MyCorr = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,2)*bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,2)
                        MyCorr = MyCorr*OPT_P_C - bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
                        bio%phy(i,j,k,l,m) = max(0., MyCorr)
                        LimitCorr(i,j,k) = 5
                      endif

                    endif

                    ! limitation on Silicon corrections
                    ! to the optimal Si/C ratio
                    if(m .eq. 5) then
                      ! compute Si/C fraction
                      MyRatio = bio%cquot(i,j,k,l,m)/bio%cquot(i,j,k,l,2)
                      if(MyRatio .gt. OPT_S_C .and. bio%phy(i,j,k,l,m) .gt. 0) then
                        MyCorr = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,2)*bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,2)
                        MyCorr = MyCorr*OPT_S_C - bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
                        bio%phy(i,j,k,l,m) = max(0., MyCorr)
                        LimitCorr(i,j,k) = 6
                      endif

                    endif

                  endif ! ApplyConditions

                  DumpBio(i,j,k,1) = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k) + bio%phy(i,j,k,l,m)
                endif
              else
                DumpBio(i,j,k,1) = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
              endif

            endif
          enddo
        enddo
      enddo

      do k=my_km+1,grd%km
        do j=1,grd%jm
          do i=1,grd%im
            DumpBio(i,j,k,1) = bio%pquot(i,j,k,l)*bio%cquot(i,j,k,l,m)*bio%InitialChl(i,j,k)
          enddo
        enddo
      enddo

      ierr = nf90mpi_put_var_all(ncid,idP,DumpBio,MyStart_4d,MyCount_4d)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

      ierr = nf90mpi_put_var_all(ncid,idTim,TimeArr,MyStartSingle,MyCountSingle)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all '//DA_VarList(iVar), ierr)

      ierr = nf90mpi_close(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close '//BioRestart, ierr)

      call MPI_Barrier(Var3DCommunicator, ierr)
      ! only process 0 creates link to restart files
      if(MyId .eq. 0) then
        SysErr = system("ln -sf $PWD/"//BioRestart//" "//BioRestartLong)
        if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
      endif
    enddo ! l
  enddo ! m

! File for post check DA
! plus check variables
  LimCorrfile = 'DA__FREQ_1/limcorr.'//ShortDate//'.nc'
  ierr = nf90mpi_create(Var3DCommunicator, LimCorrfile, NF90_CLOBBER, MPI_INFO_NULL, ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('LimCorrfile ', ierr)

  ierr = nf90mpi_def_dim(ncid,'x',global_im ,xid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
  ierr = nf90mpi_def_dim(ncid,'y' ,global_jm ,yid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)
  ierr = nf90mpi_def_dim(ncid,'z' ,global_km, depid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim depth ', ierr)

  LimVarName='lim_to_corr_FLAG'

  ierr = nf90mpi_def_var(ncid, LimVarName, nf90_int, (/xid,yid,depid/), idL )
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)

  ! ierr = nf90mpi_def_var(ncid, 'valtest', nf90_float, (/xid,yid/), iVar1 )
  ! if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)

  ! ierr = nf90mpi_def_var(ncid, 'initchl', nf90_float, (/xid,yid/), iVar2 )
  ! if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)

  ierr = nf90mpi_enddef(ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef', ierr)

  ierr = nf90mpi_put_var_all(ncid,idL,LimitCorr,MyStart,MyCount)
  ! ierr = nf90mpi_put_var_all(ncid,idL,LimitCorr,MyStartLim,MyCountLim)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all ', ierr)

  ! ierr =
  ! nf90mpi_put_var_all(ncid,iVar1,ValuesToTest(:,:,1),MyStartLim,MyCountLim)
  ! if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all ', ierr)

  ! ierr =
  ! nf90mpi_put_var_all(ncid,iVar2,bio%InitialChl(:,:,1),MyStartLim,MyCountLim)
  ! if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all ', ierr)

  ierr = nf90mpi_close(ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('LimCorrfile ', ierr)




  DEALLOCATE(DumpBio, ValuesToTest, MyConditions)


end subroutine wrt_chl_stat
