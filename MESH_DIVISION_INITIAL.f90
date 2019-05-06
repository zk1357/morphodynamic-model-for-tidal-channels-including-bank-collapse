#include 'define.h'  
#ifdef PRE_BANK_EROSION  

SUBROUTINE MESH_DIVISION_INITIAL
            
#ifdef PRE_HYD        
    USE HYDRO
#endif

#ifdef PRE_AD
    USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
    USE BANK_EROSION
    USE EXCHANGE
    USE DEF_EXCHANGE
    USE SED_EXCHANGE
#endif
    
    IMPLICIT NONE


    INTEGER   , PARAMETER                    :: NFI=SELECTED_REAL_KIND(8)
    INTEGER                                  :: I, J
    REAL(NFI)                                :: X1, X2, X3, X4, Y1, Y2, Y3, Y4
    
LOOP0:  DO I = 1 , NUMBER_BV
            
            EXCHANGE1%EXCHANGE_BANK_SIZE(I,:,:)=0.D0  
                
            DO J = 1 , EXCHANGE1%ELEMENT
                X1=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,1),1) 
                X2=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,2),1) 
                X3=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,3),1) 
                Y1=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,1),2) 
                Y2=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,2),2) 
                Y3=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,3),2)
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,10)=0.5D0*ABS(Y3*(X2-X1)+Y1*(X3-X2)+Y2*(X1-X3))
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,1)=EXCHANGE1%MESH(J,1) 
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,2)=EXCHANGE1%MESH(J,2) 
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,3)=EXCHANGE1%MESH(J,3)
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,4)=X1 
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,5)=Y1
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,6)=X2
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,7)=Y2
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,8)=X3 
                EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,9)=Y3
            END DO            

        END DO LOOP0
 
        RETURN
        
    END SUBROUTINE MESH_DIVISION_INITIAL  
    
#endif

