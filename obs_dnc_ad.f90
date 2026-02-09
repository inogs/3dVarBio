subroutine obs_dnc_ad
  
  
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
  ! Apply observational operator for density     (adjoint)               !
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
  ! if(drv%multiv.eq.1) &
  !   my_km = ros%kmchl
  
  ! condc = 0
  ! condn = 0
  ! if ((drv%chl_assim.eq.1 )  .or. (drv%multiv.eq.1)) then
  !   condc = 1
  !   call EXTEND_2D( grd%chl_ad, my_km, ChlExtended_3d )
  ! endif
  ! if ((drv%nut.eq.1 .and. bio%N3n.eq.1 )  .or. (drv%multiv.eq.1)) then
  !   call EXTEND_2D( grd%n3n_ad, grd%km, N3nExtended_3d )
  !   condn = 1
  ! endif
  ! if (drv%nut.eq.1 .and. bio%O2o.eq.1 ) &
  !   call EXTEND_2D( grd%O2o_ad, grd%km, O2oExtended_3d )
  call EXTEND_2D( grd%dnc_ad, grd%km, DncExtended_3d )


  do kk = 1,dnc%no

    i=dnc%ib(kk)
    j=dnc%jb(kk)
    k=dnc%kb(kk)

    if(dnc%flc(kk).eq.1)then

          obs%k = obs%k + 1
          obsg = obs%gra(obs%k)

          DncExtended_3d(i  ,j  ,k  ) = DncExtended_3d(i  ,j  ,k  ) + dnc%pq1(kk) * obsg
          DncExtended_3d(i+1,j  ,k  ) = DncExtended_3d(i+1,j  ,k  ) + dnc%pq2(kk) * obsg
          DncExtended_3d(i  ,j+1,k  ) = DncExtended_3d(i  ,j+1,k  ) + dnc%pq3(kk) * obsg
          DncExtended_3d(i+1,j+1,k  ) = DncExtended_3d(i+1,j+1,k  ) + dnc%pq4(kk) * obsg
          DncExtended_3d(i  ,j  ,k+1) = DncExtended_3d(i  ,j  ,k+1) + dnc%pq5(kk) * obsg
          DncExtended_3d(i+1,j  ,k+1) = DncExtended_3d(i+1,j  ,k+1) + dnc%pq6(kk) * obsg
          DncExtended_3d(i  ,j+1,k+1) = DncExtended_3d(i  ,j+1,k+1) + dnc%pq7(kk) * obsg
          DncExtended_3d(i+1,j+1,k+1) = DncExtended_3d(i+1,j+1,k+1) + dnc%pq8(kk) * obsg
    endif

  enddo


  slicevar(:,1:my_km) = grd%dnc_ad(1,:,1:my_km)
  call ADD_PREVCORE_CONTRIB(DncExtended_3d, my_km, grd%dnc_ad, slicevar(:,1:my_km))



end subroutine obs_dnc_ad
