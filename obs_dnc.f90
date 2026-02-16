subroutine obs_dnc
  
  
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
  use dnc_str
  use mpi_str
  use drv_str
  use bio_str
  
  implicit none

  INTEGER(i4)   ::  i, j, k, kk, condc, condn!, my_km
  
  ! my_km = grd%km
  ! if(drv%multiv.eq.1) &
  !   my_km = ros%kmchl

  if(MyId .eq. 0) then
     print*, "in obs_dnc grd%dnc max and sum", maxval(grd%dnc), sum(grd%dnc)
     write(drv%dia,*) "in obs_dnc grd%dnc max and sum", maxval(grd%dnc), sum(grd%dnc)
  endif
  ! condc = 0
  ! condn = 0
  ! if ((drv%chl_assim.eq.1 ) .or. (drv%multiv.eq.1)) then
  !   condc = 1
  !   call EXTEND_2D( grd%chl, my_km, ChlExtended_3d )
  ! endif
  ! if ((drv%nut.eq.1 .and. drv%dnc.eq.1 ) .or. (drv%multiv.eq.1)) then
    ! condn = 1
  call EXTEND_2D( grd%dnc, grd%km, DncExtended_3d )
  ! endif
  ! if (bio%O2o.eq.1 ) &
  !   call EXTEND_2D( grd%O2o, grd%km, O2oExtended_3d )

  if(MyId .eq. 0) then
     print*, "in obs_dnc DncExtended max and sum", maxval(DncExtended_3d), sum(DncExtended_3d)
     write(drv%dia,*) "in obs_dnc DncExtended max and sum", maxval(DncExtended_3d), sum(DncExtended_3d)
  endif

  if(MyId .eq. 0) then
     print*, "in obs_dnc pq1 dot", dot_product(dnc%pq1,dnc%pq1), maxval(dnc%pq1), sum(dnc%pq1)
     write(drv%dia,*) "in obs_dnc pq1 dot", dot_product(dnc%pq1,dnc%pq1), maxval(dnc%pq1), sum(dnc%pq1)

     print*, "in obs_dnc pq2 dot", dot_product(dnc%pq2,dnc%pq2), maxval(dnc%pq2), sum(dnc%pq2)
     write(drv%dia,*) "in obs_dnc pq2 dot", dot_product(dnc%pq2,dnc%pq2), maxval(dnc%pq2), sum(dnc%pq2)

     print*, "in obs_dnc pq3 dot", dot_product(dnc%pq3,dnc%pq3), maxval(dnc%pq3), sum(dnc%pq3)
     write(drv%dia,*) "in obs_dnc pq3 dot", dot_product(dnc%pq3,dnc%pq3), maxval(dnc%pq3), sum(dnc%pq3)

     print*, "in obs_dnc pq4 dot", dot_product(dnc%pq4,dnc%pq4), maxval(dnc%pq4), sum(dnc%pq4)
     write(drv%dia,*) "in obs_dnc pq4 dot", dot_product(dnc%pq4,dnc%pq4), maxval(dnc%pq4), sum(dnc%pq4)

     print*, "in obs_dnc pq5 dot", dot_product(dnc%pq5,dnc%pq5), maxval(dnc%pq5), sum(dnc%pq5)
     write(drv%dia,*) "in obs_dnc pq5 dot", dot_product(dnc%pq5,dnc%pq5), maxval(dnc%pq5), sum(dnc%pq5)

     print*, "in obs_dnc pq6 dot", dot_product(dnc%pq6,dnc%pq6), maxval(dnc%pq6), sum(dnc%pq6)
     write(drv%dia,*) "in obs_dnc pq6 dot", dot_product(dnc%pq6,dnc%pq6), maxval(dnc%pq6), sum(dnc%pq6)

     print*, "in obs_dnc pq7 dot", dot_product(dnc%pq7,dnc%pq7), maxval(dnc%pq7), sum(dnc%pq7)
     write(drv%dia,*) "in obs_dnc pq7 dot", dot_product(dnc%pq7,dnc%pq7), maxval(dnc%pq7), sum(dnc%pq7)

     print*, "in obs_dnc pq8 dot", dot_product(dnc%pq8,dnc%pq8), maxval(dnc%pq8), sum(dnc%pq8)
     write(drv%dia,*) "in obs_dnc pq8 dot", dot_product(dnc%pq8,dnc%pq8), maxval(dnc%pq8), sum(dnc%pq8)
  endif



  do kk = 1,dnc%no

      i=dnc%ib(kk)
      j=dnc%jb(kk)
      k=dnc%kb(kk)

      if(dnc%flc(kk).eq.1) then
      if(MyId .eq. 0) then
         print*, " dnc%flc(kk)", dnc%flc(kk)
         write(drv%dia,*) " dnc%flc(kk) in obs_dnc ", dnc%flc(kk)
      endif


        dnc%inc(kk) = &
          dnc%pq1(kk) * DncExtended_3d(i  ,j  ,k) +       &
          dnc%pq2(kk) * DncExtended_3d(i+1,j  ,k  ) +       &
          dnc%pq3(kk) * DncExtended_3d(i  ,j+1,k  ) +       &
          dnc%pq4(kk) * DncExtended_3d(i+1,j+1,k  ) +       &
          dnc%pq5(kk) * DncExtended_3d(i  ,j  ,k+1) +       &
          dnc%pq6(kk) * DncExtended_3d(i+1,j  ,k+1) +       &
          dnc%pq7(kk) * DncExtended_3d(i  ,j+1,k+1) +       &
          dnc%pq8(kk) * DncExtended_3d(i+1,j+1,k+1)
    endif
  enddo

  if(MyId .eq. 0) then
     print*, "in obs_dnc dnc%inc dotp, max and sum", dot_product(dnc%inc,dnc%inc), maxval(dnc%inc), sum(dnc%inc)
     write(drv%dia,*) "in obs_dnc dnc%inc dotp, max and sum", dot_product(dnc%inc,dnc%inc), maxval(dnc%inc), sum(dnc%inc)
  endif

end subroutine obs_dnc
