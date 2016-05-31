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
 use bmd_str
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

   write(215) sla%no
   if(sla%no.gt.0)  write(215)                                 &
        sla%ino(1:sla%no), sla%flg(1:sla%no)                   &
       ,sla%lon(1:sla%no), sla%lat(1:sla%no), sla%tim(1:sla%no)&
       ,sla%val(1:sla%no), sla%bac(1:sla%no)                   &
       ,sla%err(1:sla%no), sla%res(1:sla%no)                   &
       ,sla%bia(1:sla%no), sla%inc(1:sla%no)                   &
       ,sla%b_a(1:sla%no)                                      &
       ,sla%dpt(1:sla%no)
   write(215) arg%no
   if(arg%no.ne.0) write (215)                                  &
        arg%ino(1:arg%no), arg%flg(1:arg%no), arg%par(1:arg%no) &
       ,arg%lon(1:arg%no), arg%lat(1:arg%no)                    &
       ,arg%dpt(1:arg%no), arg%tim(1:arg%no)                    &
       ,arg%val(1:arg%no), arg%bac(1:arg%no)                    &
       ,arg%err(1:arg%no), arg%res(1:arg%no)                    &
       ,arg%bia(1:arg%no), arg%inc(1:arg%no)                    &
       ,arg%b_a(1:arg%no)
   write(215) xbt%no
   if(xbt%no.ne.0) write (215)                                  &
        xbt%ino(1:xbt%no), xbt%flg(1:xbt%no), xbt%par(1:xbt%no) &
       ,xbt%lon(1:xbt%no), xbt%lat(1:xbt%no)                    &
       ,xbt%dpt(1:xbt%no), xbt%tim(1:xbt%no)                    &
       ,xbt%val(1:xbt%no), xbt%bac(1:xbt%no)                    &
       ,xbt%err(1:xbt%no), xbt%res(1:xbt%no)                    &
       ,xbt%bia(1:xbt%no), xbt%inc(1:xbt%no)                    &
       ,xbt%b_a(1:xbt%no)
   write(215) gld%no
   if(gld%no.ne.0) write (215)                                  &
        gld%ino(1:gld%no), gld%flg(1:gld%no), gld%par(1:gld%no) &
       ,gld%lon(1:gld%no), gld%lat(1:gld%no)                    &
       ,gld%dpt(1:gld%no), gld%tim(1:gld%no)                    &
       ,gld%val(1:gld%no), gld%bac(1:gld%no)                    &
       ,gld%err(1:gld%no), gld%res(1:gld%no)                    &
       ,gld%bia(1:gld%no), gld%inc(1:gld%no)                    &
       ,gld%b_a(1:gld%no)
   write(215) tra%no
   if(tra%no.ne.0) write (215)                                  &
        tra%dpt                                                 &
       ,tra%ino(1:tra%no), tra%flg(1:tra%no)                    &
       ,tra%loi(1:tra%no), tra%lai(1:tra%no)                    &
       ,tra%lof(1:tra%no), tra%laf(1:tra%no)                    &
       ,tra%lob(tra%nt+1,1:tra%no), tra%lab(tra%nt+1,1:tra%no)  &
       ,tra%rex(1:tra%no), tra%inx(1:tra%no)                    &
       ,tra%rey(1:tra%no), tra%iny(1:tra%no)                    &
       ,tra%loa(1:tra%no), tra%laa(1:tra%no)                    &
       ,tra%erx(1:tra%no), tra%ery(1:tra%no)
   write(215) trd%no
   if(trd%no.ne.0) write (215)                                  &
        trd%dpt                                                 &
       ,trd%ino(1:trd%no), trd%flg(1:trd%no)                    &
       ,trd%loi(1:trd%no), trd%lai(1:trd%no)                    &
       ,trd%lof(1:trd%no), trd%laf(1:trd%no)                    &
       ,trd%lob(trd%nt+1,1:trd%no), trd%lab(trd%nt+1,1:trd%no)  &
       ,trd%rex(1:trd%no), trd%inx(1:trd%no)                    &
       ,trd%rey(1:trd%no), trd%iny(1:trd%no)                    &
       ,trd%loa(1:trd%no), trd%laa(1:trd%no)
   write(215) vdr%no
   if(vdr%no.ne.0) write (215)                                  &
        vdr%ino(1:vdr%no), vdr%flg(1:vdr%no), vdr%par(1:vdr%no) &
       ,vdr%lon(1:vdr%no), vdr%lat(1:vdr%no)                    &
       ,vdr%dpt(1:vdr%no), vdr%tim(1:vdr%no)                    &
       ,vdr%tms(1:vdr%no), vdr%tme(1:vdr%no)                    &
       ,vdr%val(1:vdr%no), vdr%bac(1:vdr%no)                    &
       ,vdr%err(1:vdr%no), vdr%res(1:vdr%no)                    &
       ,vdr%bia(1:vdr%no), vdr%inc(1:vdr%no)                    &
       ,vdr%b_a(1:vdr%no)
   write(215) gvl%no
   if(gvl%no.ne.0) write (215)                                  &
        gvl%ino(1:gvl%no), gvl%flg(1:gvl%no), gvl%par(1:gvl%no) &
       ,gvl%lon(1:gvl%no), gvl%lat(1:gvl%no)                    &
       ,gvl%dpt(1:gvl%no), gvl%tim(1:gvl%no)                    &
       ,gvl%tms(1:gvl%no), gvl%tme(1:gvl%no)                    &
       ,gvl%val(1:gvl%no), gvl%bac(1:gvl%no)                    &
       ,gvl%err(1:gvl%no), gvl%res(1:gvl%no)                    &
       ,gvl%bia(1:gvl%no), gvl%inc(1:gvl%no)                    &
       ,gvl%b_a(1:gvl%no)
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
