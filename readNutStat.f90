subroutine readNutStat

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!--------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!                                                                      !
! READ quotas for nutrient variables                                   !
!                                                                      !
! Version 1: A.Teruzzi 2012                                            !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  ! use filenames
  use da_params
  use bio_str
  use grd_str
  use drv_str
  use mpi_str
  use pnetcdf

  implicit none

  INTEGER(i4)        :: ncid, VarId, ierr, iVar
  INTEGER(i4)        :: i,j,k,l
   
  CHARACTER(LEN=47)  :: RstFileName
  CHARACTER(LEN=3)   :: MyVarName
  REAL(4), ALLOCATABLE :: x3(:,:,:)

  ALLOCATE(x3(grd%im, grd%jm, grd%km))
  ALLOCATE(bio%InitialNut(grd%im, grd%jm, grd%km, NNutVar)) ; bio%InitialNut(:,:,:,:) = 0.0
  
  x3(:,:,:)     = 0.0



  do l=1,NNutVar
    iVar = NPhytoVar + l

    if(iVar .gt. NBioVar) then
      if(MyId .eq. 0) &
        write(*,*) "Warning: Reading a variable not in the DA_VarList!"
    endif

    MyVarName = DA_VarList(iVar)
    RstFileName = 'DA__FREQ_1/RSTbefore.'//ShortDate//'.'//MyVarName//'.nc'

    if(drv%Verbose .eq. 1) then
      if(MyId .eq. 0) &
        write(*,*) "Reading ", RstFileName, " date: ", DA_DATE
    endif

    ierr = nf90mpi_open(Var3DCommunicator, trim(RstFileName), NF90_NOWRITE, MPI_INFO_NULL, ncid)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_open RST', ierr)

    ierr = nf90mpi_inq_varid (ncid, DA_VarList(iVar), VarId)
    if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', ierr)
    ierr = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, x3)
    if (ierr .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all RST', ierr)

    do k=1,grd%km
      do j=1,grd%jm
        do i=1,grd%im

          if(x3(i,j,k) .lt. 1.e20) then

            bio%InitialNut(i,j,k,l) = bio%InitialNut(i,j,k,l) + x3(i,j,k)
            
          else
            bio%InitialNut(i,j,k,l) = x3(i,j,k)
            if(grd%msk(i,j,k) .eq. 1) then
              write(*,*) "Warning!! Bad mask point in bio structure!"
              write(*,*) "i=",i," j=",j," k=",k
              write(*,*) "grd%msk(i,j,k)=",grd%msk(i,j,k)
              write(*,*) "bio%InitialNut(i,j,k)=",bio%InitialNut(i,j,k,l)
              write(*,*) "Aborting.."
              call MPI_Abort(Var3DCommunicator, -1, ierr)
            endif
          endif

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

  

  if(MyId .eq. 0) then
    write(drv%dia,*)'Number of Nutrients is ', NNutVar
  endif

  DEALLOCATE(x3)

end subroutine readNutStat
