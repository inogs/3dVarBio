subroutine obs_vec
  
  
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
  ! Create the observational vector                                      !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use drv_str
  use obs_str
  use mpi_str
  
  implicit none
  
  INTEGER(i4)    ::  k, i
  
  ! -------
  ! Define observational vector
  
  obs%no = sat%nc + arg%nc

  if(MyId .eq. 0) &
       write(drv%dia,*) ' Total number of observations: ', obs%no

  ALLOCATE ( obs%inc(obs%no)) ; obs%inc = huge(obs%inc(1))
  ALLOCATE ( obs%amo(obs%no)) ; obs%amo = huge(obs%amo(1))
  ALLOCATE ( obs%res(obs%no)) ; obs%res = huge(obs%res(1))
  ALLOCATE ( obs%err(obs%no)) ; obs%err = huge(obs%err(1))
  ALLOCATE ( obs%gra(obs%no)) ; obs%gra = huge(obs%gra(1))
  
  
 k=0

 if (drv%argo_obs .eq. 1) then
    ! ARGO observations
    do i=1,arg%no
       if(arg%flc(i).eq.1)then
          k=k+1
          obs%res(k) = arg%res(i)
          obs%err(k) = arg%err(i)
       endif
    enddo

    DEALLOCATE(arg%res, arg%err)
    
 endif
 
 ! Observations of satellite chlorophyll
 if(drv%sat_obs .eq. 1) then
  do i=1,sat%no
    if(sat%flc(i).eq.1)then
       k=k+1
       obs%res(k) = sat%res(i)
       obs%err(k) = sat%err(i)
    endif
  enddo
  
  DEALLOCATE(sat%res, sat%err)

 endif
 
 
end subroutine obs_vec
