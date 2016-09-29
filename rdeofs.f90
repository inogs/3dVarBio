subroutine rdeofs
  
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
  use netcdf
  use grd_str
  use filenames
  
  implicit none
  
  INTEGER(i4)                    :: stat, ncid, idvar, neofs, nlevs, nregs
  
  
  stat = nf90_open(trim(EOF_FILE), NF90_NOWRITE, ncid)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  
  ! Get dimensions 
  stat = nf90_inq_dimid (ncid, 'nreg', idvar)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inquire_dimension (ncid, idvar, len = nregs)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inq_dimid (ncid, 'nlev', idvar)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inquire_dimension (ncid, idvar, len = nlevs)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inq_dimid (ncid, 'neof', idvar)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inquire_dimension (ncid, idvar, len = neofs)
  if (stat /= nf90_noerr) call netcdf_err(stat)

  write(drv%dia,*)'Eof dimensions are: ',ros%nreg, ros%kmt, neofs
  write(drv%dia,*)'Uses ',ros%neof,' eofs.'

  if(ros%nreg .ne. nregs) then
     write(drv%dia,*)'Error: ros%nreg differs from nregs'
     !stop
     call f_exit(22)
  endif
  
  if(ros%neof .gt. neofs) then
     write(drv%dia,*)'Error: Requires more Eofs than available in the input file.'
     !stop
     call f_exit(22)
  else if(ros%neof .lt. neofs) then
     write(drv%dia,*)'Warning: ros%neof < neofs!'
     write(drv%dia,*)'ros%neof =', ros%neof
     write(drv%dia,*)'neofs =', neofs
     write(drv%dia,*)'continue using ros%neof'
  endif
  
  if(ros%kmt .ne. nlevs) then
     write(drv%dia,*)'Error: Vertical dimension different than in the input file.'
     !stop
     call f_exit(23)
  endif
  
  !  Allocate eof arrays and get data
  ALLOCATE ( ros%evc( ros%nreg, ros%kmt, ros%neof) )  ; ros%evc = huge(ros%evc(1,1,1))
  ALLOCATE ( ros%eva( ros%nreg, ros%neof) )           ; ros%eva = huge(ros%eva(1,1))
  
  stat = nf90_inq_varid (ncid, 'eva', idvar)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_get_var (ncid,idvar,ros%eva, &
       start = (/1,1/), count = (/ros%nreg, ros%neof/))
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_inq_varid (ncid, 'evc', idvar)
  if (stat /= nf90_noerr) call netcdf_err(stat)
  stat = nf90_get_var (ncid,idvar,ros%evc, &
       start = (/1,1,1/), count = (/ros%nreg, ros%kmt, ros%neof/))
  if (stat /= nf90_noerr) call netcdf_err(stat)
  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  !ros%evc(:,:,:) = 1.
  !ros%eva(:,:) = 1.
  
  stat = nf90_close(ncid)  
  
end subroutine rdeofs


