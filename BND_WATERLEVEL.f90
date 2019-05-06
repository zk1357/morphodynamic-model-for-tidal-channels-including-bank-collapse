#include 'define.h' 
SUBROUTINE BND_WATERLEVEL(T0,HH2)

    USE HYDRO

    IMPLICIT NONE
    INTEGER   , PARAMETER                :: NNB = SELECTED_REAL_KIND(8)
    REAL(NNB) , DIMENSION(:,:) , POINTER :: HH2
    REAL(NNB)                            :: T0
    INTEGER                              :: I ,J

    DO I=M_START,M_END
        HH2(I , N_START)=AMP*SIN(T0/43200.D0*2*3.1415926D0) + DEP(I,N_START)
    END DO

END SUBROUTINE BND_WATERLEVEL