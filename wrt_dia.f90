subroutine wrt_dia
  
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
  ! Write outputs and diagnostics                                        !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
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
  use bio_str
  
  implicit none
  
  INTEGER(i4)  :: l,i,j,k
  CHARACTER    :: fgrd
  integer status
  integer            :: ncid,xid,yid,depid,idchl,idn3n,idn1p,ido2o
  integer            :: idvip,idmsk,eofid
  integer(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km

  real(r4), allocatable, dimension(:,:,:) :: DumpMatrix

  ALLOCATE ( DumpMatrix(grd%im,grd%jm,grd%km) ) ; DumpMatrix  = 0.0
  
  ! ---
  ! Innovations
  if(MyId .eq. 0) &
     write(drv%dia,*) 'writes to corrections.dat !!!!!!!!!!!!!!!!!!!!!!!!!'     

  
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
  
  if(drv%chl_assim .eq. 1) then
    status = nf90mpi_def_var(ncid,'chl', nf90_float, (/xid,yid,depid/), idchl )
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var chl', status)
    status = nf90mpi_put_att(ncid,idchl   , 'missing_value',1.e+20)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', status)
  endif
  if(drv%nut .eq. 1 .and. bio%n3n .eq. 1) then
    status = nf90mpi_def_var(ncid,'n3n', nf90_float, (/xid,yid,depid/), idn3n )
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var n3n', status)
    status = nf90mpi_put_att(ncid,idn3n   , 'missing_value',1.e+20)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', status)
  endif
  if(drv%nut .eq. 1 .and. bio%n3n .eq. 1 .and. bio%updateN1p .eq. 1) then
    status = nf90mpi_def_var(ncid,'n1p', nf90_float, (/xid,yid,depid/), idn1p )
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var n1p', status)
    status = nf90mpi_put_att(ncid,idn1p   , 'missing_value',1.e+20)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', status)
  endif
  if(drv%nut .eq. 1 .and. bio%o2o .eq. 1) then
    status = nf90mpi_def_var(ncid,'o2o', nf90_float, (/xid,yid,depid/), ido2o )
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var o2o', status)
    status = nf90mpi_put_att(ncid,ido2o   , 'missing_value',1.e+20)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_att', status)
  endif
  
  status = nf90mpi_enddef(ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', status)

  if(drv%chl_assim .eq. 1) then
    do k=1,grd%km
      do j=1,grd%jm
          do i=1,grd%im
            if(grd%msk(i,j,k) .eq. 1) then
              DumpMatrix(i,j,k) = REAL(grd%chl(i,j,k), 4 )
            else
              DumpMatrix(i,j,k) = 1.e20
            endif            
          enddo
      enddo
    enddo
    status = nf90mpi_put_var_all(ncid,idchl,DumpMatrix,MyStart,MyCount)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all chl', status)
  endif


  if(drv%nut .eq. 1 .and. bio%n3n .eq. 1) then
    do k=1,grd%km
      do j=1,grd%jm
          do i=1,grd%im
            if(grd%msk(i,j,k) .eq. 1) then
              DumpMatrix(i,j,k) = REAL(grd%n3n(i,j,k), 4 )
            else
              DumpMatrix(i,j,k) = 1.e20
            endif
          enddo
      enddo
    enddo
    status = nf90mpi_put_var_all(ncid,idn3n,DumpMatrix,MyStart,MyCount)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all n3n', status)
  endif

  if(drv%nut .eq. 1 .and. bio%n3n .eq. 1 .and. bio%updateN1p .eq. 1) then
    do k=1,grd%km
      do j=1,grd%jm
          do i=1,grd%im
            if(grd%msk(i,j,k) .eq. 1) then
              DumpMatrix(i,j,k) = REAL(grd%n3n(i,j,k)*bio%covn3n_n1p(i,j,k), 4 )
            else
              DumpMatrix(i,j,k) = 1.e20
            endif
          enddo
      enddo
    enddo
    status = nf90mpi_put_var_all(ncid,idn1p,DumpMatrix,MyStart,MyCount)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all n1p', status)
  endif

  if(drv%nut .eq. 1 .and. bio%o2o .eq. 1) then
    do k=1,grd%km
      do j=1,grd%jm
          do i=1,grd%im
            if (drv%argo_obs .eq. 1) then
                if (grd%msk(i,j,k) .eq. 1) then
                  DumpMatrix(i,j,k) = REAL(grd%o2o(i,j,k), 4)
                else
                  DumpMatrix(i,j,k) = 1.e20
                endif
            else
                DumpMatrix(i,j,k) = REAL(grd%o2o(i,j,k), 4 )
            endif
          enddo
      enddo
    enddo
    status = nf90mpi_put_var_all(ncid,ido2o,DumpMatrix,MyStart,MyCount)
    if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all o2o', status)
  endif

  status = nf90mpi_close(ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_close', status)
  
  DEALLOCATE(DumpMatrix)

end subroutine wrt_dia
