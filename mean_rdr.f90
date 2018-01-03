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
      
      REAL(r8) FUNCTION mean_rad(k,rad_xy)
      use rcfl
      use grd_str
      use mpi_str

      IMPLICIT NONE
      INTEGER(i4):: j,i,count
      INTEGER(i4),INTENT(IN):: k
      REAL(r8),INTENT(IN):: rad_xy(GlobalRow,GlobalCol,1)
      REAL(r8):: meanxy
     
      ! the Mean Radius is computed 
      ! just on the owned part of the domain
      count=0
      meanxy=0
      do j=1,GlobalCol ! grd%jm
        do i=1,GlobalRow ! grd%im
          if(rad_xy(i,j,1)>0.0001) then
            count=count+1
            meanxy=meanxy+rad_xy(i,j,1)
          end if
        enddo
      enddo
      mean_rad=meanxy/count
      return     
      END FUNCTION mean_rad
