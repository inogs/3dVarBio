subroutine parallel_wrt_dia
  
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
  
  implicit none
  
  INTEGER(i4)  :: l,i,j,k
  ! REAL         :: Dump_chl(grd%im,grd%jm,grd%km)
  ! REAL         :: Dump_vip(grd%im,grd%jm,ros%neof)
  CHARACTER    :: fgrd
  integer status
  integer            :: ncid,xid,yid,depid,idchl
  integer            :: idvip,idmsk,eofid
  integer(kind=MPI_OFFSET_KIND) :: global_im, global_jm, global_km
  ! integer(kind=MPI_OFFSET_KIND) :: MyStart(3), MyCount(3)
  
  ! ---
  ! Innovations
  if(MyRank .eq. 0) &
     write(drv%dia,*) 'writes to corrections.dat !!!!!!!!!!!!!!!!!!!!!!!!!'     
  print*,""
  print*, "MyRank", MyRank, "WITHIN PARALLEL_WRITE_DIA subroutine"
  print*,""
  write(fgrd,'(i1)')drv%ktr
  
  do l=1,grd%nchl
     do k=1,grd%km
        do j=1,grd%jm
           do i=1,grd%im
              if (drv%argo .eq. 1) then
                 if (grd%msk(i,j,k) .eq. 0) then
                    Dump_chl(i,j,k) = -1.
                 else
                    Dump_chl(i,j,k) = REAL(grd%chl(i,j,k,l), 4)
                 endif
              else
                 Dump_chl(i,j,k) = REAL(grd%chl(i,j,k,l), 4 )
              endif
           enddo
        enddo
     enddo
  enddo
  
  do j=1,grd%jm
     do i=1,grd%im
        Dump_msk(i,j) = real(grd%msk(i,j,1),4);
     enddo
  enddo
  
! #ifdef _USE_MPI
!   if(MyRank .eq. 0) then
!      status = nf90_create("corr0.nc", NF90_CLOBBER, ncid)
!   else
!      status = nf90_create("corr1.nc", NF90_CLOBBER, ncid)
!   end if
! #else
!   status = nf90_create(trim(CORR_FILE), NF90_CLOBBER, ncid)
! #endif
  
  status = nf90mpi_create(MPI_COMM_WORLD, trim(CORR_FILE), NF90_CLOBBER, &
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
  ! status = nf90_sync(ncid)
  status = nf90mpi_close(ncid)
  if (status .ne. NF90_NOERR ) call handle_err('nf90mpi_close', status)
  
  
!   do k=1, ros%neof
!      do j=1, grd%jm
!         do i=1, grd%im
!            Dump_vip(i,j,k) = real(grd%ro(i,j,k),  4)
!         enddo
!      enddo
!   enddo
  
!   status = nf90_create(trim(EIV_FILE), NF90_CLOBBER, ncid) ! Eigenvalues
!   status = nf90_def_dim(ncid,'neof'     ,ros%neof ,eofid)
!   status = nf90_def_dim(ncid,'latitude' ,grd%jm  ,yid)
!   status = nf90_def_dim(ncid,'longitude',grd%im  ,xid)
!   status = nf90_def_var(ncid,'msk' , nf90_float, (/xid,yid/)      , idmsk )
!   status = nf90_def_var(ncid,'vip' , nf90_float, (/xid,yid,eofid/), idvip )
!   status = nf90_enddef(ncid)
!   status = nf90_put_var(ncid,idmsk, Dump_msk )
!   status = nf90_put_var(ncid,idvip, Dump_vip )
!   status = nf90_sync(ncid)
!   status = nf90_close(ncid)
  
  
!   ! ---
!   ! Observations
  
!   !  open(215,file=drv%flag//drv%date//'obs_'//fgrd//'.dat',form='unformatted')

! #ifdef _USE_MPI
!   if(MyRank .eq. 0) then
! #endif

!      open(215,file=trim(OBS_FILE),form='unformatted')
  
!      write(215) chl%no
     
!      close (215)
!      write(*,*)'nchl ',chl%no

! #ifdef _USE_MPI
!   endif
! #endif

end subroutine parallel_wrt_dia
