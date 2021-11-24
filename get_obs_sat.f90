subroutine get_obs_sat
  
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
  ! Version 1: A. Teruzzi 2018                                           !
  ! 
  ! Load Chlorophyll observations                                        !
  !                                                                      !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use mpi_str
  use pnetcdf
  use da_params
  use filenames
  
  implicit none
  
  INTEGER(i4)   ::  j,k, kk
  INTEGER(i4)   ::  i
  REAL(r8)      ::  zbo, zbn
  REAL(r4), ALLOCATABLE      ::  chl_mis(:,:),chl_err(:,:)
  REAL(i4), ALLOCATABLE      ::  flagblk(:,:)
  INTEGER(i4)   ::  stat, ncid, idvar, VarId, ierr
  INTEGER(i4)   ::  xid, yid, idF, ii
  INTEGER(KIND=MPI_OFFSET_KIND)   :: global_im, global_jm
  INTEGER(KIND=MPI_OFFSET_KIND)   ::  MyCountTwod(2), MyStartTwod(2)
  CHARACTER(LEN=45)   :: flagFile
  CHARACTER(LEN=15)   :: FlagVarName
  

  global_im = GlobalRow
  global_jm = GlobalCol
  do ii=1,2
    MyCountTwod(ii) = MyCount(ii)
    MyStartTwod(ii) = MyStart(ii)
  enddo



  sat%no = 0
  sat%nc = 0

  stat = nf90mpi_open(Var3DCommunicator, trim(MISFIT_FILE), NF90_NOWRITE, MPI_INFO_NULL, ncid)
  if (stat .ne. NF90_NOERR ) call handle_err('nf90mpi_open '//trim(MISFIT_FILE), stat)
  
  if(stat.ne.0)then
     sat%no = 0
     return
  endif
  
  sat%no = grd%im*grd%jm
  sat%max_val = 10.0
  
  ! ---
  ! Level corresponding to the minimum depth and maximum light propagation
  sat%kdp=grd%km
  ! do k=grd%km, 1, -1
  !    if(grd%dep(k).ge.sat%dep) sat%kdp = k
  ! enddo
  do k=1, grd%km
     if(grd%dep(k).gt.sat%dep) then
        sat%kdp = k
        exit
     endif
  enddo
  
  
  ! ---
  ! Allocate memory for observations   
  
  ALLOCATE ( chl_mis(grd%im,grd%jm) ) ; chl_mis = huge(chl_mis(1,1))
  ALLOCATE ( chl_err(grd%im,grd%jm) ) ; chl_err = huge(chl_err(1,1))
  ALLOCATE ( flagblk(grd%im,grd%jm) )
  ALLOCATE ( sat%flg(sat%no)) ; sat%flg = huge(sat%flg(1))
  ALLOCATE ( sat%flc(sat%no)) ; sat%flc = huge(sat%flc(1))
  ALLOCATE ( sat%inc(sat%no)) ; sat%inc = huge(sat%inc(1))
  ALLOCATE ( sat%err(sat%no)) ; sat%err = huge(sat%err(1))
  ALLOCATE ( sat%res(sat%no)) ; sat%res = huge(sat%res(1))
  ALLOCATE ( sat%ib(sat%no))  ; sat%ib  = huge(sat%ib(1))
  ALLOCATE ( sat%jb(sat%no))  ; sat%jb  = huge(sat%jb(1))
  ALLOCATE ( sat%pb(sat%no))  ; sat%pb  = huge(sat%pb(1))
  ALLOCATE ( sat%qb(sat%no))  ; sat%qb  = huge(sat%qb(1))
  ALLOCATE ( sat%pq1(sat%no)) ;sat%pq1 = huge(sat%pq1(1))
  ALLOCATE ( sat%pq2(sat%no)) ; sat%pq2 = huge(sat%pq2(1))
  ALLOCATE ( sat%pq3(sat%no)) ; sat%pq3 = huge(sat%pq3(1))
  ALLOCATE ( sat%pq4(sat%no)) ; sat%pq4 = huge(sat%pq4(1))
  ALLOCATE ( sat%dzr(grd%km,sat%no)) ; sat%dzr=huge(sat%dzr(1,1))
  

  stat = nf90mpi_inq_varid (ncid, 'misfchl', VarId)
  if (stat .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', stat)
  stat = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, chl_mis)
  if (stat .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all', stat)

  stat = nf90mpi_inq_varid (ncid, 'errchl', VarId)
  if (stat .ne. NF90_NOERR ) call handle_err('nf90mpi_inq_varid', stat)
  stat = nfmpi_get_vara_real_all (ncid, VarId, MyStart, MyCount, chl_err)
  if (stat .ne. NF90_NOERR ) call handle_err('nfmpi_get_vara_real_all', stat)
  
  do k=1,sat%no
     j = (k-1)/grd%im + 1
     i = k - (j-1)*grd%im
     sat%res(k) = chl_mis(i,j)
     sat%err(k) = chl_err(i,j)
  enddo

  
  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  ! sat%res(:) = 0.
  ! sat%err(:) = 1.

  DEALLOCATE( chl_mis )
  DEALLOCATE( chl_err )
  
  stat = nf90mpi_close (ncid)
  if (stat .ne. NF90_NOERR ) call handle_err('nf90mpi_close', stat)
  
  !   sat%err(:) =  0.3
  
  ! ---
  ! Initialise quality flag, do residual check, compute vertical integration parameters and count good observations
  flagblk(:,:) = 0 !blacklisting flag 
  sat%nc = 0
  do k=1,sat%no
     j = (k-1)/grd%im + 1
     i = k - (j-1)*grd%im
     if(grd%msk(i,j,sat%kdp).eq.1. )then
        sat%flg(k) = 1
        if(abs(sat%res(k)).gt.sat%max_val) then
           ! residual check
           sat%flg(k) = 0
           if(abs(sat%res(k)).lt.100) then
             flagblk(i,j) = 1
           endif
        else
           ! compute vertical integration parameters
           zbn = grd%dep(1)*2.0
           sat%dzr(1,k) = zbn
           zbo = zbn
           sat%dzr(:,k) = sat%dzr(:,k) / zbo

           ! Update sat%nc variable
           sat%nc = sat%nc + 1
        end if
     else
        sat%flg(k) = 0
     endif
  enddo

  call MPI_Allreduce(sat%nc, sat%nc_global, 1, MPI_INT, MPI_SUM, Var3DCommunicator, stat)

  if(MyId .eq. 0) then
     print*,'Good chl observations: ',sat%nc_global
     print*,'Saving flag misfit'
  endif

  ! Saving flag misfit sat
  flagFile = 'DA__FREQ_1/flagsat.'//ShortDate//'.nc'
  ierr = nf90mpi_create(Var3DCommunicator, trim(flagFile), NF90_CLOBBER, MPI_INFO_NULL,ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('flagFile ', ierr)

  ierr = nf90mpi_def_dim(ncid,'x',global_im ,xid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim longitude ', ierr)
  ierr = nf90mpi_def_dim(ncid,'y' ,global_jm ,yid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_dim latitude ', ierr)

  FlagVarName='flag_lim_misf'

  ierr = nf90mpi_def_var(ncid, FlagVarName, nf90_int, (/xid,yid/), idF )
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_def_var', ierr)

  ierr = nf90mpi_enddef(ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_enddef', ierr)


  ierr = nf90mpi_put_var_all(ncid,idF,flagblk,MyStartTwod,MyCountTwod)
  if (ierr .ne. NF90_NOERR ) call handle_err('nf90mpi_put_var_all ', ierr)

  ierr = nf90mpi_close(ncid)
  if (ierr .ne. NF90_NOERR ) call handle_err('LimCorrfile ', ierr)



  DEALLOCATE ( flagblk )



  sat%flc(:) = sat%flg(:)

end subroutine get_obs_sat

subroutine int_par_chl
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Get interpolation parameters for a grid                              !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use mpi_str
  
  implicit none
  
  integer(i4)   ::  k
  integer(i4)   ::  i1, kk, j1
  real(r8)      ::  p1, q1
  real(r8)      ::  div_x, div_y

  if(MyId .eq. 0) &
       write(drv%dia,*) 'Number of CHL observations:  >>>>>>>>>>>>>',sat%nc_global

  if(sat%nc.gt.0) then
     
     
     ! ---
     ! Interpolation parameters
     do kk = 1,sat%no
        j1 = (kk-1)/grd%im + 1
        i1 = kk - (j1-1)*grd%im
        q1 = 0.0
        p1 = 0.0
        sat%ib(kk) = i1
        sat%jb(kk) = j1
        sat%pb(kk) = p1
        sat%qb(kk) = q1
     enddo
     
     
     ! ---
     ! Horizontal interpolation parameters for each masked grid
     do k = 1,sat%no
        if(sat%flc(k) .eq. 1) then
           
           i1=sat%ib(k)
           p1=sat%pb(k)
           j1=sat%jb(k)
           q1=sat%qb(k)
           
           div_y =  (1.-q1) * max(grd%global_msk(i1+GlobalRowOffset,j1  ,1),grd%global_msk(i1+GlobalRowOffset+1,j1  ,1))     &
                +    q1  * max(grd%global_msk(i1+GlobalRowOffset,j1+1,1),grd%global_msk(i1+GlobalRowOffset+1,j1+1,1))
           div_x =  (1.-p1) * grd%global_msk(i1+GlobalRowOffset  ,j1,1) + p1 * grd%global_msk(i1+GlobalRowOffset+1,j1,1)
           sat%pq1(k) = grd%global_msk(i1+GlobalRowOffset,j1,1)                                      &
                * max(grd%global_msk(i1+GlobalRowOffset,j1,1),grd%global_msk(i1+GlobalRowOffset+1,j1,1))             &
                * (1.-p1) * (1.-q1)                                   &
                /( div_x * div_y + 1.e-16 )
           sat%pq2(k) = grd%global_msk(i1+GlobalRowOffset+1,j1,1)                                    &
                * max(grd%global_msk(i1+GlobalRowOffset,j1,1),grd%global_msk(i1+GlobalRowOffset+1,j1,1))             &
                *     p1  * (1.-q1)                                    &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%global_msk(i1+GlobalRowOffset  ,j1+1,1) + p1 * grd%global_msk(i1+GlobalRowOffset+1,j1+1,1)
           sat%pq3(k) = grd%global_msk(i1+GlobalRowOffset,j1+1,1)                                    &
                * max(grd%global_msk(i1+GlobalRowOffset,j1+1,1),grd%global_msk(i1+GlobalRowOffset+1,j1+1,1))         &
                * (1.-p1) *     q1                                     &
                /( div_x * div_y + 1.e-16 )
           sat%pq4(k) = grd%global_msk(i1+GlobalRowOffset+1,j1+1,1)                                  &
                * max(grd%global_msk(i1+GlobalRowOffset,j1+1,1),grd%global_msk(i1+GlobalRowOffset+1,j1+1,1))         &
                *     p1  *     q1                                     &
                /( div_x * div_y + 1.e-16 )
           
        endif
     enddo    
  endif
  
  DEALLOCATE ( sat%pb, sat%qb)
  DEALLOCATE ( sat%flg)


end subroutine int_par_chl
