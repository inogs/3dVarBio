subroutine readBioStat

!---------------------------------------------------------------------------
!anna                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!--------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!                                                                      !
! READ quotas for biological variables                                 !
!                                                                      !
! Version 1: A.Teruzzi 2012                                            !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use filenames
  use bio_str
  use grd_str
  use drv_str
  use mpi_str
  use pnetcdf

  implicit none

  INTEGER(i4)        :: ncid, VarId, ierr, iVar
  INTEGER(i4)        :: i,j,k,l,m
   
  CHARACTER(LEN=51)  :: RstFileName
  CHARACTER(LEN=3)   :: MyVarName
  REAL(4), ALLOCATABLE :: x3(:,:,:)

  ALLOCATE(x3(grd%im, grd%jm, grd%km))
  ALLOCATE(bio%InitialChl(grd%im, grd%jm, grd%km)) ; bio%InitialChl(:,:,:) = 0.0
  ! Allocate quotas arrys
  ALLOCATE(bio%pquot( grd%im, grd%jm, grd%km, bio%nphy)) ; bio%pquot(:,:,:,:) = 0.0
  ALLOCATE(bio%cquot( grd%im, grd%jm, grd%km, bio%nphy, bio%ncmp)) ; bio%cquot(:,:,:,:,:) = 0.0

  x3(:,:,:)     = 0.0

  bio%DA_VarList( 1)='P1l'
  bio%DA_VarList( 2)='P2l'
  bio%DA_VarList( 3)='P3l'
  bio%DA_VarList( 4)='P4l'

  bio%DA_VarList( 5)='P1n'
  bio%DA_VarList( 6)='P2n'
  bio%DA_VarList( 7)='P3n'
  bio%DA_VarList( 8)='P4n'

  bio%DA_VarList( 9)='P1c'
  bio%DA_VarList(10)='P2c'
  bio%DA_VarList(11)='P3c'
  bio%DA_VarList(12)='P4c'

  bio%DA_VarList(13)='P1p'
  bio%DA_VarList(14)='P2p'
  bio%DA_VarList(15)='P3p'
  bio%DA_VarList(16)='P4p'

  bio%DA_VarList(17)='P1s'

  do m=1,bio%ncmp
    do l=1,bio%nphy
      iVar = l + bio%nphy*(m-1)

      if(iVar .gt. 17) cycle

      MyVarName = bio%DA_VarList(iVar)
      RstFileName = 'DA__FREQ_1/RST.'//DA_DATE//'.'//MyVarName//'.nc'

      if(drv%Verbose .eq. 1) then
        if(MyId .eq. 0) &
          write(*,*) "Reading ", RstFileName, " date: ", DA_DATE
      endif

      ierr = nf90mpi_open(Var3DCommunicator, trim(RstFileName), NF90_NOWRITE, MPI_INFO_NULL, ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open RST', ierr)

      ierr = nf90mpi_inq_varid (ncid, bio%DA_VarList(iVar), VarId)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
      ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x3)
      if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all RST', ierr)

      if(m .eq. 1) then
        ! preparing pquot array
        do k=1,grd%km
          do j=1,grd%jm
            do i=1,grd%im
              ! bio%pquot(i,j,k,l) = grd%msk(i,j,k) * x3(i,j,k)
              bio%pquot(i,j,k,l) = x3(i,j,k)
              bio%InitialChl(i,j,k) = bio%InitialChl(i,j,k) + x3(i,j,k)

              ! debug print
              ! if((l.eq.3) .and. (i.eq.181) .and. (j.eq.158) .and. (k.eq.22)) &
              !   print*, "pquot ", bio%pquot(181,158,22,3)
            enddo
          enddo
        enddo
      endif

      do k=1,grd%km
        do j=1,grd%jm
          do i=1,grd%im
            if(bio%pquot(i,j,k,l).ne.0) &
              bio%cquot(i,j,k,l,m) = grd%msk(i,j,k) * x3(i,j,k) / bio%pquot(i,j,k,l)

            ! debug print
            ! if(bio%cquot(i,j,k,l,m) .gt. 1.e6) write(*,*) "i ", i, " j ", j, " k ", k, " l ", l, " m ", m, " bio%cquot ", bio%cquot(i,j,k,l,m)
            ! if((l.eq.3) .and. (i.eq.181) .and. (j.eq.158) .and. (k.eq.22) .and. (m.eq.4)) &
            !   print*, "cquot ", bio%cquot(181,158,22,3,4), " x3 ", x3(181,158,22)
          enddo
        enddo
      enddo
      

      ierr = nf90mpi_close(ncid)
      if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_close RST', ierr)
    
      if(drv%Verbose .eq. 1) then
        if(MyId .eq. 0) &
          write(*,*) "Restart ", RstFileName, " read"
      endif

    enddo
  enddo

  do l=1,bio%nphy
    do k=1,grd%km
      do j=1,grd%jm
        do i=1,grd%im
          if(bio%InitialChl(i,j,k).ne.0) &
            bio%pquot(i,j,k,l) = bio%pquot(i,j,k,l) / bio%InitialChl(i,j,k)
        enddo
      enddo
    enddo
  enddo
  

  if(MyId .eq. 0) then
    write(drv%dia,*)'Number of phytoplankton types is ', bio%nphy
    write(drv%dia,*)'Number of phytoplankton components is ', bio%ncmp
  endif

  DEALLOCATE(x3)

end subroutine readBioStat
