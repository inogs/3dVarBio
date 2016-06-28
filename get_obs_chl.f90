subroutine get_obs_chl
  
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
  ! Load Chlorophyll observations                                                !
  !                                                                      !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use netcdf
  use filenames
  
  implicit none
  
  INTEGER(i4)   ::  j,k, kk
  INTEGER(i4)   ::  i
  REAL(r8)      ::  zbo, zbn
  REAL(r4), ALLOCATABLE      ::  chl_mis(:,:),chl_err(:,:)
  INTEGER(i4)   ::  stat, ncid, idvar
  
  chl%no = 0
  chl%nc = 0
  
  
  stat = nf90_open(trim(MISFIT_FILE), NF90_NOWRITE, ncid)
  
  if(stat.ne.0)then
     chl%no = 0
     return
  endif
  
  chl%no = grd%im*grd%jm
  chl%max_val = 10.0
  
  ! ---
  ! Level corresponding to the minimum depth and maximum light propagation
  chl%kdp=grd%km
  ! do k=grd%km, 1, -1
  !    if(grd%dep(k).ge.chl%dep) chl%kdp = k
  ! enddo
  do k=1, grd%km
     if(grd%dep(k).gt.chl%dep) then
        chl%kdp = k
        exit
     endif
  enddo
  
  
  ! ---
  ! Allocate memory for observations 
  
  
  ALLOCATE ( chl_mis(grd%im,grd%jm) ) ; chl_mis = huge(chl_mis(1,1))
  ALLOCATE ( chl_err(grd%im,grd%jm) ) ; chl_err = huge(chl_err(1,1))
  ALLOCATE ( chl%flg(chl%no)) ; chl%flg = huge(chl%flg(1))
  ALLOCATE ( chl%flc(chl%no)) ; chl%flc = huge(chl%flc(1))
  ALLOCATE ( chl%inc(chl%no)) ; chl%inc = huge(chl%inc(1))
  ALLOCATE ( chl%err(chl%no)) ; chl%err = huge(chl%err(1))
  ALLOCATE ( chl%res(chl%no)) ; chl%res = huge(chl%res(1))
  ALLOCATE ( chl%ib(chl%no))  ; chl%ib  = huge(chl%ib(1))
  ALLOCATE ( chl%jb(chl%no))  ; chl%jb  = huge(chl%jb(1))
  ALLOCATE ( chl%pb(chl%no))  ; chl%pb  = huge(chl%pb(1))
  ALLOCATE ( chl%qb(chl%no))  ; chl%qb  = huge(chl%qb(1))
  ALLOCATE ( chl%pq1(chl%no)) ;chl%pq1 = huge(chl%pq1(1))
  ALLOCATE ( chl%pq2(chl%no)) ; chl%pq2 = huge(chl%pq2(1))
  ALLOCATE ( chl%pq3(chl%no)) ; chl%pq3 = huge(chl%pq3(1))
  ALLOCATE ( chl%pq4(chl%no)) ; chl%pq4 = huge(chl%pq4(1))
  ALLOCATE ( chl%dzr(grd%km,chl%no)) ; chl%dzr=huge(chl%dzr(1,1))
  
  
  
  
  stat = nf90_inq_varid (ncid, 'misfchl', idvar)
  if(stat .ne. 0) then
     call f_exit_message(10, 'var misfchl not found')
  endif
  
  stat = nf90_get_var (ncid, idvar, chl_mis)
  if(stat .ne. 0) then
     call f_exit_message(10, 'cannot read var misfchl')
  endif
  
  stat = nf90_inq_varid (ncid, 'errchl', idvar)
  if(stat .ne. 0) then
     call f_exit_message(10, 'var errchl not found')
  endif
  
  stat = nf90_get_var (ncid, idvar, chl_err)
  if(stat .ne. 0) then
     call f_exit_message(10, 'cannot read var errchl')
  endif
  
  do k=1,chl%no
     j = (k-1)/grd%im + 1
     i = k - (j-1)*grd%im
     chl%res(k) = chl_mis(i,j)
     chl%err(k) = chl_err(i,j)
  enddo
  
  
  DEALLOCATE( chl_mis )
  DEALLOCATE( chl_err )
  
  stat = nf90_close (ncid)
  
  !   chl%err(:) =  0.3
  
  ! ---
  ! Initialise quality flag, do residual check, compute vertical integration parameters and count good observations
  chl%nc = 0
  do k=1,chl%no
     j = (k-1)/grd%im + 1
     i = k - (j-1)*grd%im
     if(grd%msk(i,j,chl%kdp).eq.1. )then
        chl%flg(k) = 1
        if(abs(chl%res(k)).gt.chl%max_val) then
           ! residual check
           chl%flg(k) = 0
        else
           ! compute vertical integration parameters
           zbn = grd%dep(1)*2.0
           chl%dzr(1,k) = zbn
           zbo = zbn
           chl%dzr(:,k) = chl%dzr(:,k) / zbo

           ! Update chl%nc variable
           chl%nc = chl%nc + 1
        end if
     else
        chl%flg(k) = 0
     endif
  enddo
  
  print*,'Good chl observations: ',chl%nc
  chl%flc(:) = chl%flg(:)

end subroutine get_obs_chl

subroutine int_par_chl
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Get interpolation parameters for a grid                              !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  
  implicit none
  
  integer(i4)   ::  k
  integer(i4)   ::  i1, kk, j1
  real(r8)      ::  p1, q1
  real(r8)      ::  div_x, div_y
  
  write(drv%dia,*) 'Number of CHL observations:  >>>>>>>>>>>>>',chl%nc
  
  if(chl%nc.gt.0) then
     
     
     ! ---
     ! Interpolation parameters
     do kk = 1,chl%no
        j1 = (kk-1)/grd%im + 1
        i1 = kk - (j1-1)*grd%im
        q1 = 0.0
        p1 = 0.0
        chl%ib(kk) = i1
        chl%jb(kk) = j1
        chl%pb(kk) = p1
        chl%qb(kk) = q1
     enddo
     
     
     ! ---
     ! Horizontal interpolation parameters for each masked grid
     do k = 1,chl%no
        if(chl%flc(k) .eq. 1) then
           
           i1=chl%ib(k)
           p1=chl%pb(k)
           j1=chl%jb(k)
           q1=chl%qb(k)
           
           div_y =  (1.-q1) * max(grd%msk(i1,j1  ,1),grd%msk(i1+1,j1  ,1))     &
                +    q1  * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))
           div_x =  (1.-p1) * grd%msk(i1  ,j1,1) + p1 * grd%msk(i1+1,j1,1)
           chl%pq1(k) = grd%msk(i1,j1,1)                                      &
                * max(grd%msk(i1,j1,1),grd%msk(i1+1,j1,1))             &
                * (1.-p1) * (1.-q1)                                   &
                /( div_x * div_y + 1.e-16 )
           chl%pq2(k) = grd%msk(i1+1,j1,1)                                    &
                * max(grd%msk(i1,j1,1),grd%msk(i1+1,j1,1))             &
                *     p1  * (1.-q1)                                    &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%msk(i1  ,j1+1,1) + p1 * grd%msk(i1+1,j1+1,1)
           chl%pq3(k) = grd%msk(i1,j1+1,1)                                    &
                * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))         &
                * (1.-p1) *     q1                                     &
                /( div_x * div_y + 1.e-16 )
           chl%pq4(k) = grd%msk(i1+1,j1+1,1)                                  &
                * max(grd%msk(i1,j1+1,1),grd%msk(i1+1,j1+1,1))         &
                *     p1  *     q1                                     &
                /( div_x * div_y + 1.e-16 )
           
        endif
     enddo    
  endif
  
end subroutine int_par_chl
