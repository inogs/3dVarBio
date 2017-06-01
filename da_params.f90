MODULE DA_PARAMS

  implicit none

  PUBLIC

    character (LEN=15)             :: DA_DATE          != '20130102-120000'
    integer                        :: jpk_200          != 26
    integer                        :: NBioVar          ! number of biological variables
    CHARACTER(LEN=3), allocatable  :: DA_VarList(:)    ! name of DA biological variables
    integer                        :: DA_JulianDate    ! julian date

  CONTAINS

  SUBROUTINE SET_DA_PARAMS

    DA_DATE = '20130101-120000'
    jpk_200 = 26
    NBioVar = 17

    allocate(DA_VarList(NBioVar))

    ! DA_VarList init
    ! It must be consistent with NBioVar value

    DA_VarList( 1)='P1l'
    DA_VarList( 2)='P2l'
    DA_VarList( 3)='P3l'
    DA_VarList( 4)='P4l'

    DA_VarList( 5)='P1n'
    DA_VarList( 6)='P2n'
    DA_VarList( 7)='P3n'
    DA_VarList( 8)='P4n'

    DA_VarList( 9)='P1c'
    DA_VarList(10)='P2c'
    DA_VarList(11)='P3c'
    DA_VarList(12)='P4c'

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