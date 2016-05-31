subroutine min_cfn


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
! Minimise the cost function                                           !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str

 implicit none

 INTEGER(i4)     :: cntr, k
!anna --- Stopping criteria on number of minimization steps
! INTEGER(i4)     :: nstp
!anna ---
 REAL(r8)        :: maxpg

  cntr = 0

!anna
!  nstp=0
!anna



  do while(ctl%task(1:5).eq.'NEW_X' .or. ctl%task(1:2).eq.'FG' .or. ctl%task(1:5).eq.'START')
!anna
!  do while((ctl%task(1:5).eq.'NEW_X' .or. ctl%task(1:2).eq.'FG' .or. ctl%task(1:5).eq.'START').and.nstp.le.300)
!anna


    call setulb(ctl%n, ctl%m, ctl%x_c, ctl%l_c, ctl%u_c, ctl%nbd, ctl%f_c, ctl%g_c,    &
                ctl%factr, ctl%pgtol, ctl%ws, ctl%wy, ctl%sy, ctl%ss, ctl%wt,          &
                ctl%wn, ctl%snd, ctl%z_c, ctl%r_c, ctl%d_c, ctl%t_c, ctl%wa,           &
                ctl%iwa,                                     &
                ctl%task, ctl%iprint,  ctl%csave, ctl%lsave, ctl%isave, ctl%dsave)
!anna
!                ctl%task, ctl%iprint,  ctl%csave, ctl%lsave, ctl%isave, ctl%dsave,nstp)


      if (ctl%task(1:2) .eq. 'FG') then

! Calculate the cost function and its gradient
          call costf

! Modify the stopping criteria
         if(cntr.eq.0 .and. ctl%pgper.ne.0.0 .and. drv%ktr.eq.1 )then
            maxpg = 0.0
           do k=1,ctl%n
            maxpg = max(maxpg,abs(ctl%g_c(k)))
           enddo
            ctl%pgtol = maxpg * ctl%pgper
            cntr = 1
         endif

      endif

  enddo !while

end subroutine min_cfn
