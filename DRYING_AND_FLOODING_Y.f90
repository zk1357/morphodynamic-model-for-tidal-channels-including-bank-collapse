#include 'define.h' 
SUBROUTINE DRYING_AND_FLOODING_Y(PP2,HH2,QQ2)
    
    USE HYDRO
    USE INTERFACE_DEFINITIONS
    
    IMPLICIT NONE
    INTEGER , PARAMETER        :: NND=SELECTED_REAL_KIND(8)
    REAL(NND) , PARAMETER      :: DF_LIMIT=0.5D0
    REAL(NND) , DIMENSION(:,:) :: PP2, HH2, QQ2
    REAL(NND)                  :: BDEPTH, TDEPTH, WATERLEVEL
    INTEGER :: I, J

    DO I= M_START+1 , M_END-1 !X-SWEEP流速点判断
DAF0:   DO J= N_START , N_END
            IF (NDD(I,J)==-1) CYCLE DAF0

            IF (QQ2(I,J)>0.D0) THEN !UP-WINDING
                TDEPTH=HH2(I,J)
            ELSE IF (QQ2(I,J)<0.D0) THEN
                TDEPTH=HH2(I,J+1)
            ELSE IF (QQ2(I,J)==0.D0) THEN
                WATERLEVEL = MAX( (HH2(I,J)-DEP(I,J)), (HH2(I,J+1)-DEP(I,J+1)) )
                BDEPTH = MIN( DEP(I,J), DEP(I,J+1) )
                TDEPTH = WATERLEVEL+BDEPTH                 
            END IF

            IF (VDDX(I,J)== 1) THEN
                IF(TDEPTH<=0.5D0*FLOODINGDEPTHX) THEN
                   VDDX(I,J) =-1
                   QQ2(I,J)  =0.D0
                END IF
            END IF

        END DO DAF0
    END DO

    DO I= N_START+1 , N_END-1 !Y-SWEEP流速点判断
DAF:    DO J= M_START+1 , M_END-1
            IF (NDD(J,I)==-1) CYCLE DAF

            IF (PP2(J,I)>0.D0) THEN  !UP-WINDINF
                TDEPTH=HH2(J,I)
            ELSE IF (PP2(J,I)<0.D0) THEN
                TDEPTH=HH2(J+1,I)
            ELSE IF (PP2(J,I)==0.D0) THEN 
                WATERLEVEL = MAX( (HH2(J,I)-DEP(J,I)), (HH2(J+1,I)-DEP(J+1,I)) )
                BDEPTH = MIN( DEP(J,I), DEP(J+1,I) )
                TDEPTH = WATERLEVEL+BDEPTH
            END IF

            IF (VDDY(J,I)==-1) THEN
                IF (TDEPTH>FLOODINGDEPTHY)       VDDY(J,I) = 1
            ELSE IF (VDDY(J,I)==1) THEN
                IF (TDEPTH<0.5D0*FLOODINGDEPTHY) VDDY(J,I) =-1
            END IF
        END DO DAF
    END DO


    DO I= N_START+1 , N_END-1 !水位点干湿判断
DAF1:   DO J= M_START+1 , M_END-1
            IF(NDD(J,I)==-1) CYCLE DAF1
            IF(VDDY(J,I)== 1  .OR. VDDY(J-1,I)== 1  .OR. VDDX(J,I)== 1  .OR. VDDX(J,I-1)== 1) NDD(J,I)=0
            IF(VDDY(J,I)==-1 .AND. VDDY(J-1,I)==-1 .AND. VDDX(J,I)==-1 .AND. VDDX(J,I-1)==-1) NDD(J,I)=1
        END DO DAF1
    END DO

    RETURN
END SUBROUTINE DRYING_AND_FLOODING_Y