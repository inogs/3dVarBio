subroutine obs_arg_ad
  
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
  !                                                                          !
  !    This file is part of OceanVar.                                        !
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
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
  !                                                                          !
  !---------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply observational operator for ARGO floats (adjoint)               !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use eof_str
  use obs_str
  use mpi_str
  use filenames
  use drv_str
  use bio_str

  implicit none
  
  INTEGER(i4)   ::  i, j, k, kk, condc, condn, my_km
  REAL(r8), DIMENSION(grd%jm,grd%km)  :: slicevar
  REAL(8) :: obsg

  my_km = grd%km
  if(drv%multiv.eq.1) &
    my_km = ros%kmchl
  
  condc = 0
  condn = 0
  if ((drv%chl_assim.eq.1 )  .or. (drv%multiv.eq.1)) then
    condc = 1
    call EXTEND_2D( grd%chl_ad, my_km, ChlExtended_3d )
  endif
  if ((drv%nut.eq.1 .and. bio%N3n.eq.1 )  .or. (drv%multiv.eq.1)) then
    call EXTEND_2D( grd%n3n_ad, grd%km, N3nExtended_3d )
    condn = 1
  endif
  if (drv%nut.eq.1 .and. bio%O2o.eq.1 ) &
    call EXTEND_2D( grd%O2o_ad, grd%km, O2oExtended_3d )


  do kk = 1,arg%no

    i=arg%ib(kk)
    j=arg%jb(kk)
    k=arg%kb(kk)

    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.0 .and. condc.eq.1)then

          obs%k = obs%k + 1
          ! if(drv%multiv .eq. 1) then
          !   obsg = obs%gra(obs%k)*arg%std(kk)
          ! else
          obsg = obs%gra(obs%k)
          ! end if

          ChlExtended_3d(i  ,j  ,k  ) = ChlExtended_3d(i  ,j  ,k  ) + arg%pq1(kk) * obsg
          ChlExtended_3d(i+1,j  ,k  ) = ChlExtended_3d(i+1,j  ,k  ) + arg%pq2(kk) * obsg
          ChlExtended_3d(i  ,j+1,k  ) = ChlExtended_3d(i  ,j+1,k  ) + arg%pq3(kk) * obsg
          ChlExtended_3d(i+1,j+1,k  ) = ChlExtended_3d(i+1,j+1,k  ) + arg%pq4(kk) * obsg
          ChlExtended_3d(i  ,j  ,k+1) = ChlExtended_3d(i  ,j  ,k+1) + arg%pq5(kk) * obsg
          ChlExtended_3d(i+1,j  ,k+1) = ChlExtended_3d(i+1,j  ,k+1) + arg%pq6(kk) * obsg
          ChlExtended_3d(i  ,j+1,k+1) = ChlExtended_3d(i  ,j+1,k+1) + arg%pq7(kk) * obsg
          ChlExtended_3d(i+1,j+1,k+1) = ChlExtended_3d(i+1,j+1,k+1) + arg%pq8(kk) * obsg
    endif
    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.1 .and. condn.eq.1) then

          obs%k = obs%k + 1
          ! if(drv%multiv .eq. 1) then
          !   obsg = obs%gra(obs%k)*arg%std(kk)
          ! else
          obsg = obs%gra(obs%k)
          ! end if

          N3nExtended_3d(i  ,j  ,k  ) = N3nExtended_3d(i  ,j  ,k  ) + arg%pq1(kk) * obsg
          N3nExtended_3d(i+1,j  ,k  ) = N3nExtended_3d(i+1,j  ,k  ) + arg%pq2(kk) * obsg
          N3nExtended_3d(i  ,j+1,k  ) = N3nExtended_3d(i  ,j+1,k  ) + arg%pq3(kk) * obsg
          N3nExtended_3d(i+1,j+1,k  ) = N3nExtended_3d(i+1,j+1,k  ) + arg%pq4(kk) * obsg
          N3nExtended_3d(i  ,j  ,k+1) = N3nExtended_3d(i  ,j  ,k+1) + arg%pq5(kk) * obsg
          N3nExtended_3d(i+1,j  ,k+1) = N3nExtended_3d(i+1,j  ,k+1) + arg%pq6(kk) * obsg
          N3nExtended_3d(i  ,j+1,k+1) = N3nExtended_3d(i  ,j+1,k+1) + arg%pq7(kk) * obsg
          N3nExtended_3d(i+1,j+1,k+1) = N3nExtended_3d(i+1,j+1,k+1) + arg%pq8(kk) * obsg


    endif
    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.2 .and. drv%nut.eq.1 .and. bio%o2o.eq.1) then

          obs%k = obs%k + 1

          O2oExtended_3d(i  ,j  ,k  ) = O2oExtended_3d(i  ,j  ,k  ) + arg%pq1(kk) * obs%gra(obs%k)
          O2oExtended_3d(i+1,j  ,k  ) = O2oExtended_3d(i+1,j  ,k  ) + arg%pq2(kk) * obs%gra(obs%k)
          O2oExtended_3d(i  ,j+1,k  ) = O2oExtended_3d(i  ,j+1,k  ) + arg%pq3(kk) * obs%gra(obs%k)
          O2oExtended_3d(i+1,j+1,k  ) = O2oExtended_3d(i+1,j+1,k  ) + arg%pq4(kk) * obs%gra(obs%k)
          O2oExtended_3d(i  ,j  ,k+1) = O2oExtended_3d(i  ,j  ,k+1) + arg%pq5(kk) * obs%gra(obs%k)
          O2oExtended_3d(i+1,j  ,k+1) = O2oExtended_3d(i+1,j  ,k+1) + arg%pq6(kk) * obs%gra(obs%k)
          O2oExtended_3d(i  ,j+1,k+1) = O2oExtended_3d(i  ,j+1,k+1) + arg%pq7(kk) * obs%gra(obs%k)
          O2oExtended_3d(i+1,j+1,k+1) = O2oExtended_3d(i+1,j+1,k+1) + arg%pq8(kk) * obs%gra(obs%k)

    endif

  enddo

!  we apply contribution in grd%variable_ad

  if((drv%chl_assim.eq.1) .or. (drv%multiv.eq.1))  then
    slicevar(:,1:my_km) = grd%chl_ad(1,:,1:my_km)
    call ADD_PREVCORE_CONTRIB(ChlExtended_3d, my_km, grd%chl_ad, slicevar(:,1:my_km))
    ! call ADD_PREVCORE_CONTRIB(ChlExtended_3d, my_km, grd%chl_ad, grd%chl_ad(1,:,:))
  endif
  if((bio%N3n.eq.1) .or. (drv%multiv.eq.1))  then
    slicevar(:,1:grd%km) = grd%n3n_ad(1,:,:)
    call ADD_PREVCORE_CONTRIB(N3nExtended_3d, grd%km, grd%n3n_ad, slicevar)
    ! call ADD_PREVCORE_CONTRIB(N3nExtended_3d,  grd%km, grd%N3n_ad, grd%n3n_ad(1,:,:))
  endif
  if (bio%O2o.eq.1 ) then
    slicevar(:,1:grd%km) = grd%o2o_ad(1,:,:)
    call ADD_PREVCORE_CONTRIB(O2oExtended_3d, grd%km, grd%o2o_ad, slicevar)
    ! call ADD_PREVCORE_CONTRIB(O2oExtended_3d,  grd%km, grd%O2o_ad, grd%o2o_ad(1,:,:))
  endif



end subroutine obs_arg_ad
