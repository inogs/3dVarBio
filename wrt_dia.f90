subroutine wrt_dia
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2018 Anna Teruzzi, OGS, Trieste                         !
  !                                                                          !
  !    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                          !
  !                                                                          !
  !    3DVarBio is  free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    3DVarBio is  distributed in the hope that it will be useful,           !
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
  ! Write outputs and diagnostics                                        !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  use pnetcdf
  use filenames
  use mpi_str
  
  implicit none
  
  INTEGER(i4)  :: l,i,j,k
  CHARACTER    :: fgrd
  integer status
  integer            :: ncid,xid,yid,depid,idchl
  integer            :: idvip,idmsk,eofid
  integer(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km

  real(r4), allocatable, dimension(:,:,:) :: Dump_chl

  ALLOCATE ( Dump_chl(grd%im,grd%jm,grd%km) ) ; Dump_chl  = 0.0
  
  ! ---
  ! Innovations
  if(MyId .eq. 0) &
     write(drv%dia,*) 'writes to corrections.dat !!!!!!!!!!!!!!!!!!!!!!!!!'     

  do k=1,grd%km
     do j=1,grd%jm
        do i=1,grd%im
           if (drv%argo_obs .eq. 1) then
              if (grd%msk(i,j,k) .eq. 0) then
                 Dump_chl(i,j,k) = -1.
              else
                 Dump_chl(i,j,k) = REAL(grd%chl(i,j,k), 4)
              endif
           else
              Dump_chl(i,j,k) = REAL(grd%chl(i,j,k), 4 )
           endif
        enddo
     enddo
  enddo
  
  status = nf90mpi_create(Var3DCommunicator, trim(CORR_FILE), NF90_CLOBBER, &
       MPI_INFO_NULL, ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_create ', status)

  global_im = GlobalRow
  global_jm = GlobalCol
  global_km = grd%km
  
  status = nf90mpi_def_dim(ncid,'depth'    ,global_km, depid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim depth ', status)
  status = nf90mpi_def_dim(ncid,'latitude' ,global_jm ,yid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', status)
  status = nf90mpi_def_dim(ncid,'longitude',global_im ,xid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', status)
  
  status = nf90mpi_def_var(ncid,'chl', nf90_float, (/xid,yid,depid/), idchl )
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', status)
  status = nf90mpi_enddef(ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', status)
  
  status = nf90mpi_put_var_all(ncid,idchl,Dump_chl,MyStart,MyCount)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all', status)
  status = nf90mpi_close(ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_close', status)
  
  DEALLOCATE(Dump_chl)

end subroutine wrt_dia
