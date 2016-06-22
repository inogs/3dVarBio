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
 use netcdf
 use filenames
 implicit none

 INTEGER(i4)  :: l,i,j,k
! REAL         :: Dump_chl(grd%im,grd%jm,grd%km)
! REAL         :: Dump_vip(grd%im,grd%jm,ros%neof)
 CHARACTER    :: fgrd
 integer status
 integer            :: ncid,xid,yid,depid,idchl
 integer            :: idvip,idmsk,eofid



! ---
! Innovations

#ifdef __FISICA
   REAL         :: x2(grd%im,grd%jm)
    x2(:,:) = grd%eta(:,:)
    write(101) x2
   do k=1,grd%km
    x2(:,:) = grd%tem(:,:,k)
    write(101) x2
   enddo
   do k=1,grd%km
    x2(:,:) = grd%sal(:,:,k)
    write(101) x2
   enddo
   do k=1,grd%km
    x2(:,:) = grd%uvl(:,:,k)
    write(101) x2
   enddo
   do k=1,grd%km
    x2(:,:) = grd%vvl(:,:,k)
    write(101) x2
   enddo
#endif

   write(drv%dia,*) 'writes to corrections.dat !!!!!!!!!!!!!!!!!!!!!!!!!'

   write(fgrd,'(i1)')drv%ktr

  if(drv%biol.eq.1) then
   do l=1,grd%nchl
   do k=1,grd%km
   do j=1,grd%jm
   do i=1,grd%im
    Dump_chl(i,j,k) = REAL(grd%chl(i,j,k,l), 4 )
   enddo
   enddo
   enddo
   enddo

  endif
   do j=1,grd%jm
   do i=1,grd%im
      Dump_msk(i,j) = real(grd%msk(i,j,1),4);
   enddo
   enddo


status = nf90_create(trim(CORR_FILE), NF90_CLOBBER, ncid)
status = nf90_def_dim(ncid,'depth'    ,grd%km, depid)
status = nf90_def_dim(ncid,'latitude' ,grd%jm ,yid)
status = nf90_def_dim(ncid,'longitude',grd%im ,xid)

status = nf90_def_var(ncid,'chl', nf90_float, (/xid,yid,depid/), idchl )
status = nf90_enddef(ncid)

status = nf90_put_var(ncid,idchl,Dump_chl)
status = nf90_sync(ncid)
status = nf90_close(ncid)


do k=1, ros%neof
do j=1, grd%jm
do i=1, grd%im
   Dump_vip(i,j,k) = real(grd%ro(i,j,k),  4)
enddo
enddo
enddo

status = nf90_create(trim(EIV_FILE), NF90_CLOBBER, ncid) ! Eigenvalues
status = nf90_def_dim(ncid,'neof'     ,ros%neof ,eofid)
status = nf90_def_dim(ncid,'latitude' ,grd%jm  ,yid)
status = nf90_def_dim(ncid,'longitude',grd%im  ,xid)
status = nf90_def_var(ncid,'msk' , nf90_float, (/xid,yid/)      , idmsk )
status = nf90_def_var(ncid,'vip' , nf90_float, (/xid,yid,eofid/), idvip )
status = nf90_enddef(ncid)
status = nf90_put_var(ncid,idmsk, Dump_msk )
status = nf90_put_var(ncid,idvip, Dump_vip )
status = nf90_sync(ncid)
status = nf90_close(ncid)


! ---
! Observations

!  open(215,file=drv%flag//drv%date//'obs_'//fgrd//'.dat',form='unformatted')
  open(215,file=trim(OBS_FILE),form='unformatted')

  write(215) chl%no
  
#ifdef __FISICA
   if(chl%no.ne.0) write (215)                                  &
        gvl%flg(1:gvl%no)                                       &
       ,gvl%dpt(1:gvl%no)                                       &
       ,gvl%err(1:gvl%no), gvl%res(1:gvl%no)                    &
       ,gvl%inc(1:gvl%no)
#endif
  close (215)
  write(*,*)'nchl ',chl%no

end subroutine wrt_dia
