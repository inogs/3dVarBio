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

  INTEGER(i4)                    :: stat, ncid, idvar, neofs, nlevs


    stat = nf90_open(trim(EOF_FILE), NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)

! Get dimensions 
      stat = nf90_inq_dimid (ncid, 'nreg', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = ros%nreg)
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

    if(ros%neof .gt. neofs) then
      write(drv%dia,*)'Error: Requires more Eofs than available in the input file.'
      !stop
      call f_exit(22)
    endif
    if(ros%kmt .ne. nlevs) then
      write(drv%dia,*)'Error: Vertical dimension different than in the input file.'
      !stop
      call f_exit(23)
    endif
      
!  Allocate eof arrays
     ALLOCATE ( ros%evc( ros%nreg, ros%kmt, neofs))   ; ros%evc = huge(ros%evc(1,1,1))
     ALLOCATE (ros%eva( ros%nreg, neofs) )            ; ros%eva = huge(ros%eva(1,1))
     ALLOCATE ( ros%cor( ros%nreg, ros%neof, neofs) ) ; ros%cor = huge(ros%cor(1,1,1))

      stat = nf90_inq_varid (ncid, 'eva', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,ros%eva)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'evc', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,ros%evc)
    if (stat /= nf90_noerr) call netcdf_err(stat)

    stat = nf90_close(ncid)

!    ros%eva(:,:) = 0.1
!    ros%evc(:,:,:,:) = 0.1
    

!    do nec=1,ros%neof
!     do k=1,ros%neof
!      do nrg=1,ros%nreg
!       if(k.eq.nec)then
!!        ros%cor(nrg,k,nec) = ros%eva(nrg,nec)
!       else
!!        ros%cor(nrg,k,nec) = 0.0
!       endif
!      enddo
!     enddo
!    enddo


end subroutine rdeofs


