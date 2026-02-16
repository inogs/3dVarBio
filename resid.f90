subroutine resid
  
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
  ! Calculate analysis - observation                                     !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------

  
  use set_knd
  use obs_str
  use drv_str
  use dnc_str

  use mpi_str
  
  implicit none
  
  INTEGER(i4)   :: i, k
  
  k = 0
  
  ! ---
  ! ARGO observations
  if (drv%argo_obs .eq. 1) then
     do i=1,arg%no
        if(arg%flc(i).eq.1)then
           k = k + 1
           obs%inc(k) = arg%inc(i)
           obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
        endif
     enddo
  endif

  
  ! ---
  ! Observations of satellite chlorophyll
  if(drv%sat_obs .eq. 1) then
    do i=1,sat%no
     if(sat%flc(i).eq.1)then
        k = k + 1
        obs%inc(k) = sat%inc(i) 
        obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
     endif
    enddo
  endif
  

  if(MyId .eq. 0) then
     print*, "in resid dnc%inc dotp, sum, max", dot_product(dnc%inc,dnc%inc), sum(dnc%inc), maxval(dnc%inc)
     write(drv%dia,*) "in resid dnc%inc dotp, sum, max", dot_product(dnc%inc,dnc%inc), sum(dnc%inc), maxval(dnc%inc)
  endif
  ! ---
  ! Density increments
  if(drv%dnc .eq. 1) then
   do i=1,dnc%no
      if(dnc%flc(i).eq.1)then
         k = k + 1
         obs%inc(k) = dnc%inc(i) 
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
      endif
   enddo
  endif

  if(MyId .eq. 0) then
      print*, "in resid obs%inc dotp, sum, max", dot_product(obs%inc,obs%inc), sum(obs%inc), maxval(obs%inc)
      write(drv%dia,*) "in resid obs%inc dotp, sum, max", dot_product(obs%inc,obs%inc), sum(obs%inc), maxval(obs%inc)
  endif

  if(MyId .eq. 0) then
      print*, "in resid obs%res dotp, sum, max", dot_product(obs%res,obs%res), sum(obs%res), maxval(obs%res)
      write(drv%dia,*) "in resid obs%res dotp, sum, max", dot_product(obs%res,obs%res), sum(obs%res), maxval(obs%res)
  endif

  if(MyId .eq. 0) then
      print*, "in resid obs%err dotp, sum, max", dot_product(obs%err,obs%err), sum(obs%err), maxval(obs%err)
      write(drv%dia,*) "in resid obs%err dotp, sum, max", dot_product(obs%err,obs%err), sum(obs%err), maxval(obs%err)
  endif

  if(MyId .eq. 0) then
      print*, "in resid dnc%flc dotp, sum, max", dot_product(dnc%flc,dnc%flc), sum(dnc%flc), maxval(dnc%flc)
      write(drv%dia,*) "in resid dnc%flc dotp, sum, max", dot_product(dnc%flc,dnc%flc), sum(dnc%flc), maxval(dnc%flc)
  endif


end subroutine resid
