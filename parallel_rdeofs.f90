subroutine parallel_rdeofs
  
  !---------------------------------------------------------------------------
  !                                                                          !
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
  ! READ parameters of the MFS_16_72 grid                                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  ! This routine will have effect only if compiled with netcdf library.  !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use eof_str
  use grd_str
  use filenames

  use mpi
  use mpi_str
  use pnetcdf
  
  implicit none
  
  INTEGER(i4)                    :: stat, ncid, idvar
  integer(8)                     :: neofs, nlevs, nregs
  integer(KIND=MPI_OFFSET_KIND)  :: GlobalStart(3), GlobalCount(3)
  real(4), allocatable           :: x3(:,:,:), x2(:,:)
  
  ! stat = nf90_open(trim(EOF_FILE), NF90_NOWRITE, ncid)
  stat = nf90mpi_open(MyCommWorld, trim(EOF_FILE), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open", stat)
  
  ! Get dimensions 
  stat = nf90mpi_inq_dimid (ncid, 'nreg', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid nreg", stat)
  stat = nfmpi_inq_dimlen (ncid, idvar, nregs)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen nregs", stat)
  stat = nf90mpi_inq_dimid (ncid, 'nlev', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid nlev", stat)
  stat = nfmpi_inq_dimlen (ncid, idvar, nlevs)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen nlevs", stat)
  stat = nf90mpi_inq_dimid (ncid, 'neof', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid neof", stat)
  stat = nfmpi_inq_dimlen (ncid, idvar, len = neofs)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen neofs", stat)

  if(MyRank .eq. 0) then
     write(drv%dia,*)'Eof dimensions are: ',ros%nreg, ros%kmt, neofs
     write(drv%dia,*)'Uses ',ros%neof,' eofs.'
  endif
  
  if(ros%nreg .ne. nregs) then

     if(MyRank .eq. 0) &
          write(drv%dia,*)'Error: ros%nreg differs from nregs'

     call MPI_Abort(MyCommWorld, -1, stat)
     
  endif
  
  if(ros%neof .gt. neofs) then

     if(MyRank .eq. 0) &
          write(drv%dia,*)'Error: Requires more Eofs than available in the input file.'
     call MPI_Abort(MyCommWorld, -1, stat)
     
  else if(ros%neof .lt. neofs) then
     
     if(MyRank .eq. 0) then
        write(drv%dia,*)'Warning: ros%neof < neofs!'
        write(drv%dia,*)'ros%neof =', ros%neof
        write(drv%dia,*)'neofs =', neofs
        write(drv%dia,*)'continue using ros%neof'
     endif
  endif
  
  if(ros%kmt .ne. nlevs) then
     if(MyRank .eq. 0) &
          write(drv%dia,*)'Error: Vertical dimension different than in the input file.'

     call MPI_Abort(MyCommWorld, -1, stat)
  endif
  
  !  Allocate eof arrays and get data
  ALLOCATE ( ros%evc( ros%nreg, ros%kmt, ros%neof) )  ; ros%evc = huge(ros%evc(1,1,1))
  ALLOCATE ( ros%eva( ros%nreg, ros%neof) )           ; ros%eva = huge(ros%eva(1,1))
  ALLOCATE ( x3( ros%nreg, ros%kmt, ros%neof) )
  ALLOCATE ( x2( ros%nreg, ros%neof) )
  GlobalStart(:) = 1
  GlobalCount(1) = ros%nreg
  GlobalCount(2) = ros%kmt
  GlobalCount(3) = ros%neof
  
  stat = nf90mpi_inq_varid(ncid, 'evc', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid evc", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart, GlobalCount, x3)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all eva", stat)

  ros%evc(:,:,:) = x3(:,:,:)
  
  GlobalCount(1) = ros%nreg
  GlobalCount(2) = ros%neof

  stat = nf90mpi_inq_varid(ncid, 'eva', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid eva", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart(1:2), GlobalCount(1:2), x2)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all", stat)
  ros%eva(:,:) = x2(:,:)
  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  ! ros%evc(:,:,:) = 1.
  ! ros%eva(:,:) = 1.
  
  stat = nf90mpi_close(ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_close", stat)

  DEALLOCATE(x3, x2)
  
end subroutine parallel_rdeofs


