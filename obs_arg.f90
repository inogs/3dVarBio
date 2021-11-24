subroutine obs_arg
  
  
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
  ! Apply observational operator for ARGO floats                         !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  use eof_str
  use obs_str
  use mpi_str
  use drv_str
  use bio_str
  
  implicit none

  INTEGER(i4)   ::  i, j, k, kk, condc, condn, my_km
  
  my_km = grd%km
  if(drv%multiv.eq.1) &
    my_km = ros%kmchl

  condc = 0
  condn = 0
  if ((drv%chl_assim.eq.1 ) .or. (drv%multiv.eq.1)) then
    condc = 1
    call EXTEND_2D( grd%chl, my_km, ChlExtended_3d )
  endif
  if ((drv%nut.eq.1 .and. bio%N3n.eq.1 ) .or. (drv%multiv.eq.1)) then
    condn = 1
    call EXTEND_2D( grd%n3n, grd%km, N3nExtended_3d )
  endif
  if (bio%O2o.eq.1 ) &
    call EXTEND_2D( grd%O2o, grd%km, O2oExtended_3d )



  do kk = 1,arg%no

      i=arg%ib(kk)
      j=arg%jb(kk)
      k=arg%kb(kk)

    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.0 .and. condc.eq.1) then

        arg%inc(kk) = &
          arg%pq1(kk) * ChlExtended_3d(i  ,j  ,k) +       &
          arg%pq2(kk) * ChlExtended_3d(i+1,j  ,k  ) +       &
          arg%pq3(kk) * ChlExtended_3d(i  ,j+1,k  ) +       &
          arg%pq4(kk) * ChlExtended_3d(i+1,j+1,k  ) +       &
          arg%pq5(kk) * ChlExtended_3d(i  ,j  ,k+1) +       &
          arg%pq6(kk) * ChlExtended_3d(i+1,j  ,k+1) +       &
          arg%pq7(kk) * ChlExtended_3d(i  ,j+1,k+1) +       &
          arg%pq8(kk) * ChlExtended_3d(i+1,j+1,k+1)
        ! if(drv%multiv .eq. 1) &
        !   arg%inc(kk) = arg%inc(kk) * arg%std(kk)
    endif
    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.1 .and. condn.eq.1) then


        arg%inc(kk) = &
          arg%pq1(kk) * N3nExtended_3d(i  ,j  ,k) +       &
          arg%pq2(kk) * N3nExtended_3d(i+1,j  ,k  ) +       &
          arg%pq3(kk) * N3nExtended_3d(i  ,j+1,k  ) +       &
          arg%pq4(kk) * N3nExtended_3d(i+1,j+1,k  ) +       &
          arg%pq5(kk) * N3nExtended_3d(i  ,j  ,k+1) +       &
          arg%pq6(kk) * N3nExtended_3d(i+1,j  ,k+1) +       &
          arg%pq7(kk) * N3nExtended_3d(i  ,j+1,k+1) +       &
          arg%pq8(kk) * N3nExtended_3d(i+1,j+1,k+1)
        ! if(drv%multiv .eq. 1) &
        !   arg%inc(kk) = arg%inc(kk) * arg%std(kk)
    endif
    if(arg%flc(kk).eq.1 .and. arg%par(kk).eq.2 .and. drv%nut.eq.1 .and. bio%o2o.eq.1) then

        arg%inc(kk) = &
          arg%pq1(kk) * O2oExtended_3d(i  ,j  ,k) +       &
          arg%pq2(kk) * O2oExtended_3d(i+1,j  ,k  ) +       &
          arg%pq3(kk) * O2oExtended_3d(i  ,j+1,k  ) +       &
          arg%pq4(kk) * O2oExtended_3d(i+1,j+1,k  ) +       &
          arg%pq5(kk) * O2oExtended_3d(i  ,j  ,k+1) +       &
          arg%pq6(kk) * O2oExtended_3d(i+1,j  ,k+1) +       &
          arg%pq7(kk) * O2oExtended_3d(i  ,j+1,k+1) +       &
          arg%pq8(kk) * O2oExtended_3d(i+1,j+1,k+1)

     endif


  enddo


end subroutine obs_arg
