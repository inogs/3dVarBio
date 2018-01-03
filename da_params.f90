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

MODULE DA_PARAMS

  implicit none

  PUBLIC

    character (LEN=17)             :: DA_DATE          != '20130102-120000'
    character (LEN=15)             :: ShortDate        != '20130102-120000'
    integer                        :: jpk_200          != 26
    integer                        :: NBioVar          ! number of biological variables
    CHARACTER(LEN=3), allocatable  :: DA_VarList(:)    ! name of DA biological variables
    double precision               :: DA_JulianDate    ! julian date

  CONTAINS

  SUBROUTINE SET_DA_PARAMS

    DA_DATE = '20130101-12:00:00'
    ShortDate = DA_DATE(1:11)//DA_DATE(13:14)//DA_DATE(16:17)
    jpk_200 = 26
    NBioVar = 17

    allocate(DA_VarList(NBioVar))

    ! DA_VarList init
    ! It must be consistent with NBioVar value

    DA_VarList( 1)='P1l'
    DA_VarList( 2)='P2l'
    DA_VarList( 3)='P3l'
    DA_VarList( 4)='P4l'

    DA_VarList( 5)='P1c'
    DA_VarList( 6)='P2c'
    DA_VarList( 7)='P3c'
    DA_VarList( 8)='P4c'

    DA_VarList( 9)='P1n'
    DA_VarList(10)='P2n'
    DA_VarList(11)='P3n'
    DA_VarList(12)='P4n'

    DA_VarList(13)='P1p'
    DA_VarList(14)='P2p'
    DA_VarList(15)='P3p'
    DA_VarList(16)='P4p'

    DA_VarList(17)='P1s'


  END SUBROUTINE SET_DA_PARAMS

  SUBROUTINE CLEAN_DA_PARAMS

    DEALLOCATE(DA_VarList)

  END SUBROUTINE CLEAN_DA_PARAMS

END MODULE DA_PARAMS

