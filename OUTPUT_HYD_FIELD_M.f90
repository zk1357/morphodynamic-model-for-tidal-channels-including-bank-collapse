#include 'define.h' 
!MATLAB流速场信息作图
SUBROUTINE OUTPUT_HYD_FIELD_M
    USE TIME_DIN
#ifdef PRE_HYD        
    USE HYDRO
#endif

#ifdef PRE_AD
    USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
    USE BANK_EROSION
    USE EXCHANGE
    USE SED_EXCHANGE
#endif
    
    IMPLICIT NONE
    
    INTEGER            :: I , J , K, P
    REAL*8             :: U , V , UV, TB, CHEZY
    CHARACTER*40       :: FILENAME
    
    FILENAME=''
    WRITE(FILENAME(1:4) , '(I4.4)') GLOBTIME%GLOBYEAR
    WRITE(FILENAME(5:6) , '(I2.2)') GLOBTIME%GLOBMON
    WRITE(FILENAME(7:8) , '(I2.2)') GLOBTIME%GLOBDAY
    WRITE(FILENAME(9:10) , '(I2.2)') GLOBTIME%GLOBHOUR
    WRITE(FILENAME(11:12) , '(I2.2)') GLOBTIME%GLOBMIN
    WRITE(FILENAME(13:14) , '(I2.2)') GLOBTIME%GLOBSEC
    
    FILENAME(15:22)='HYDM.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME
    
    OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    
    K=1    
    !DO I = M_START , M_END
    !    DO J = N_START , N_END-1 ,1
    !        WRITE(600 , REC = K ) P1(I,J)
    !        K=K+1   
    !    END DO
    !END DO   
    !
    !DO I = M_START , M_END
    !    DO J = N_START , N_END-1 ,1
    !        WRITE(600 , REC = K ) Q1(I,J)
    !        K=K+1   
    !    END DO
    !END DO    
    !
    !DO I = M_START , M_END
    !    DO J = N_START , N_END-1 ,1
    !        WRITE(600 , REC = K ) H1(I,J)
    !        K=K+1   
    !    END DO
    !END DO 
    
    DO I = M_START+1 , M_END-1
        DO J = N_START+1 , N_END-1 ,1
            CALL UV_CACULATION(I,J,U,V,UV)           
            WRITE(600 , REC = K ) U
            K=K+1   
        END DO
    END DO 
    
    DO I = M_START+1 , M_END-1
        DO J = N_START+1 , N_END-1 ,1
            CALL UV_CACULATION(I,J,U,V,UV)           
            WRITE(600 , REC = K ) v
            K=K+1   
        END DO
    END DO     

    DO I = M_START+1 , M_END-1
        DO J = N_START+1 , N_END-1 ,1          
            WRITE(600 , REC = K ) h1(i,j)
            K=K+1   
        END DO
    END DO 
    
    DO I = M_START+1 , M_END-1
        DO J = N_START+1 , N_END-1 ,1
            !计算流速
            CALL UV_CACULATION(I,J,U,V,UV)    
#ifdef PRE_CHEZY 
            !谢才公式
            CHEZY=CHE
#endif

#ifdef PRE_MANNING
            !MANNING公式
            CHEZY=MANNING*H1(I,J)**(1.D0/6.D0)
#endif               
            TB=1000*9.8D0/CHEZY**2*UV*U          
            WRITE(600 , REC = K ) TB
            K=K+1   
        END DO
    END DO    
    
    
    CLOSE(600)
    
END SUBROUTINE OUTPUT_HYD_FIELD_M 