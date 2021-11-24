subroutine rdeofs_multi
  
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

  use mpi_str
  use pnetcdf
  
  implicit none
  
  INTEGER(i4)                    :: stat, ncid, idvar
  INTEGER(8)                     :: ik,ireg
  integer(8)                     :: neofs, nlevs, nregs, nregsstd
  integer(KIND=MPI_OFFSET_KIND)  :: GlobalStart(3), GlobalCount(3)
  real(4), allocatable           :: std_chl(:),std_n3n(:)
  real(4), allocatable           :: x3(:,:,:), x2(:,:), x1(:)
  

  ! Read std file for normalization
  stat = nf90mpi_open(Var3DCommunicator, trim(STD_FILE_MULTI), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open "//trim(STD_FILE_MULTI), stat)
  
  ! Get dimensions 
  stat = nf90mpi_inq_dimid (ncid, 'nreg', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_dimid nreg", stat)
  stat = nfmpi_inq_dimlen (ncid, idvar, nregsstd)
  if (stat /= nf90_noerr) call handle_err("nfmpi_inq_dimlen nregsstd", stat)

  if(MyId .eq. 0) then
     write(drv%dia,*)'Number of reg in std file : ',nregsstd
     write(drv%dia,*)'Uses ',ros%neof_multi,' eofs.'
  endif
  
  if(ros%nreg .ne. nregsstd) then

     if(MyId .eq. 0) &
          write(drv%dia,*)'Error: nregsstd differs from nregs'

     call MPI_Abort(Var3DCommunicator, -1, stat)
     
  endif
  
  !  Allocate eof arrays and get data
  ALLOCATE ( std_chl( ros%nreg ) )  ; std_chl = huge(std_chl(1))
  ALLOCATE ( std_n3n( ros%nreg ) )  ; std_n3n = huge(std_n3n(1))
  ALLOCATE ( x1( ros%nreg ) )
  GlobalStart(:) = 1
  GlobalCount(1) = ros%nreg
  GlobalCount(2) = ros%kmchl+ros%kmnit
  GlobalCount(3) = ros%neof_multi
  
  stat = nf90mpi_inq_varid(ncid, 'stdP_l', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid evc", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart(1), GlobalCount(1), x1)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all eva", stat)

  std_chl(:) = x1(:)
  
  stat = nf90mpi_inq_varid(ncid, 'stdN3n', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid eva", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart(1), GlobalCount(1), x1)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all", stat)

  std_n3n(:) = x1(:)
  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  ! ros%evc_multi(:,:,:) = 1.
  ! ros%eva_multi(:,:) = 1.
  
  stat = nf90mpi_close(ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_close", stat)




  stat = nf90mpi_open(Var3DCommunicator, trim(EOF_FILE_MULTI), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_open "//trim(EOF_FILE_MULTI), stat)
  
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

  if(MyId .eq. 0) then
     write(drv%dia,*)'Eof dimensions for multivariate are: ',ros%nreg, ros%kmchl+ros%kmnit, neofs
     write(drv%dia,*)'Uses ',ros%neof_multi,' eofs.'
  endif
  
  if(ros%nreg .ne. nregs) then

     if(MyId .eq. 0) &
          write(drv%dia,*)'Error: ros%nreg differs from nregs'

     call MPI_Abort(Var3DCommunicator, -1, stat)
     
  endif
  
  if(ros%neof_multi .gt. neofs) then

     if(MyId .eq. 0) &
          write(drv%dia,*)'Error: Requires more Eofs than available in the input file.'
     call MPI_Abort(Var3DCommunicator, -1, stat)
     
  else if(ros%neof_multi .lt. neofs) then
     
     if(MyId .eq. 0) then
        write(drv%dia,*)'Warning: ros%neof_multi < neofs!'
        write(drv%dia,*)'ros%neof_multi =', ros%neof_multi
        write(drv%dia,*)'neofs =', neofs
        write(drv%dia,*)'continue using ros%neof_multi'
        write(*,*)'Warning: ros%neof_multi < neofs!'
        write(*,*)'ros%neof_multi =', ros%neof_multi
        write(*,*)'neofs =', neofs
        write(*,*)'continue using ros%neof_multi'
     endif
  endif
  
  if((ros%kmchl+ros%kmnit) .ne. nlevs) then
     if(MyId .eq. 0) &
          write(drv%dia,*)'Error: Vertical dimension different than in the input file.'

     call MPI_Abort(Var3DCommunicator, -1, stat)
  endif
  
  !  Allocate eof arrays and get data
  ALLOCATE ( ros%evc_multi( ros%nreg, ros%kmchl+ros%kmnit, ros%neof_multi) )  ; ros%evc_multi = huge(ros%evc_multi(1,1,1))
  ALLOCATE ( ros%eva_multi( ros%nreg, ros%neof_multi) )           ; ros%eva_multi = huge(ros%eva_multi(1,1))
  ALLOCATE ( x3( ros%nreg, ros%kmchl+ros%kmnit, ros%neof_multi) )
  ALLOCATE ( x2( ros%nreg, ros%neof_multi) )
  GlobalStart(:) = 1
  GlobalCount(1) = ros%nreg
  GlobalCount(2) = ros%kmchl+ros%kmnit
  GlobalCount(3) = ros%neof_multi
  
  stat = nf90mpi_inq_varid(ncid, 'evc', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid evc", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart, GlobalCount, x3)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all eva", stat)

  ros%evc_multi(:,:,:) = x3(:,:,:)

  do ireg=1,ros%nreg
    do ik=1,ros%kmchl
      ros%evc_multi(ireg,ik,:) = ros%evc_multi(ireg,ik,:)*std_chl(ireg)
    enddo
    do ik=ros%kmchl+1,ros%kmchl+ros%kmnit
      ros%evc_multi(ireg,ik,:) = ros%evc_multi(ireg,ik,:)*std_n3n(ireg)
    enddo
  enddo
  
  GlobalCount(1) = ros%nreg
  GlobalCount(2) = ros%neof_multi

  stat = nf90mpi_inq_varid(ncid, 'eva', idvar)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_inq_varid eva", stat)
  stat = nfmpi_get_vara_real_all(ncid,idvar,GlobalStart(1:2), GlobalCount(1:2), x2)
  if (stat /= nf90_noerr) call handle_err("nfmpi_get_vara_real_all", stat)
  ros%eva_multi(:,:) = x2(:,:)
  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  ! ros%evc_multi(:,:,:) = 1.
  ! ros%eva_multi(:,:) = 1.
  
  stat = nf90mpi_close(ncid)
  if (stat /= nf90_noerr) call handle_err("nf90mpi_close", stat)

  DEALLOCATE(x3, x2, std_chl, std_n3n)
  
end subroutine rdeofs_multi


