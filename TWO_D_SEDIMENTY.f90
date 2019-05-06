#include 'define.h'
#ifdef PRE_AD
SUBROUTINE TWO_D_SEDIMENTY
    
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
    REAL*8                              :: DTT , TEMP , QD , TB , QE , &
                                        &  A , B , CX , CY , U , V , UV  , &
                                        &  DXX , DYY , HL , HR, CHEZY
    INTEGER                             :: I , J , K , JC
    
    DTT=0.5D0*DT

LOOP0:  DO J = N_START+1 , N_END-1
            AMT3=0.D0
            AMT3(2,:)=1.D0
LOOP1:      DO I = M_START , M_END
                !固壁边界
                IF ( NDD(I,J)==-1 .OR. NDD(I,J)==1 ) THEN
                    AMT3(4,I)=C1(I,J)
                    CYCLE LOOP1 
                END IF

                !水深小于临界值，不算
                IF (H1(I,J)<= FLOODINGDEPTHY) THEN
                    AMT3(4,J)=C1(I,J)
                    CYCLE LOOP1 
                END IF                 
                
                !开边界，梯度为0
                IF ( NDD(I,J)==4 ) THEN
                    IF ( I==M_START ) THEN
                        AMT3(3,I)=-1 
                    ELSE 
                        AMT3(1,I)=-1
                    END IF
                    CYCLE LOOP1
                END IF

                CALL UV_CACULATION(I,J,U,V,UV)
                CX= 0.5D0*U*DTT/DX1(I,J)
                CY= 0.5D0*V*DTT/DY1(I,J)
                A = SIGN(1.D0 , U)
                B = SIGN(1.D0 , V)
#ifdef PRE_CHEZY                
                !谢才公式
                CHEZY=CHE 
#endif

#ifdef PRE_MANNING
                !MANNING公式
                CHEZY=MANNING*H1(I,J)**(1.D0/6.D0)
#endif
                TB= FLOWDENSITY*9.8D0/CHEZY**2*UV**2
                QE= MAX( 0.D0 , ME*(TB/TE-1.D0) )

                !扩散项
                !DXX
                IF ( ABS(NDD(I,J+1))==1 .OR. NDD(I,J+1)==2) THEN
                    DXX= KX*DTT/DX1(I,J)**2*( -C1(I,J)+C1(I,J-1) )
                ELSE IF ( ABS(NDD(I,J-1))==1 .OR. NDD(I,J-1)==2) THEN
                    DXX= KX*DTT/DX1(I,J)**2*( -C1(I,J)+C1(I,J+1) )
                ELSE
                    DXX= KX*DTT/DX1(I,J)**2*( C1(I,J+1)-2*C1(I,J)+C1(I,J-1) )
                END IF
                
                !DYY
                IF ( ABS(NDD(I+1,J))==1 .OR. NDD(I+1,J)==2) THEN
                    DYY= KY*DTT/DY1(I,J)**2*( -C1(I,J)+C1(I-1,J) )
                ELSE IF ( ABS(NDD(I-1,J))==1 .OR. NDD(I-1,J)==2) THEN
                    DYY= KY*DTT/DY1(I,J)**2*( -C1(I,J)+C1(I+1,J) )
                ELSE
                    DYY= KY*DTT/DY1(I,J)**2*( C1(I+1,J)-2*C1(I,J)+C1(I-1,J) )
                END IF 
                
                
                
!                AMT(1,J)= -(1.D0+B)*CY*H2(I-1,J)
!                AMT(2,J)=  (1.D0+2.D0*CY*B)*H2(I,J)+WS*DTT*(1.D0-TB/TD)
!                AMT(3,J)=  (1.D0-B)*CY*H2(I+1,J)
!                AMT(4,J)=  C1(I,J)*H1(I,J)+QE*DTT-CX*( (1.D0-A)*C1(I,J+1)*H1(I,J+1)+2.D0*A*C1(I,J)*H1(I,J)-(1.D0+A)*C1(I,J-1)*H1(I,J-1) )
!               可能Y向流速过小,导致AMT(2,J)太小,计算时无法做除法,因此提前除以H2(I,J),上面注释部分为初始离散
                AMT3(1,I)= -(1.D0+B)*CY
                AMT3(2,I)=  (1.D0+2.D0*CY*B)+WS*DTT*MAX(0.D0,(1.D0-TB/TD))/H1(I,J)
                AMT3(3,I)=  (1.D0-B)*CY
                AMT3(4,I)=  C1(I,J)+QE*DTT/H1(I,J)-CX*( (1.D0-A)*C1(I,J+1)+2.D0*A*C1(I,J)-(1.D0+A)*C1(I,J-1) )+DXX+DYY
               
                
!LOOP2:          DO K = 1 , NUMBER_BV
!                    IF ( I==BVP(K,1) .AND. J==BVP(K,2)  ) THEN
!                        AMT3(4,J)=AMT3(4,J)+S_BANK_EROSION(K)/H2(I,J)*DTT/DY1(1,1)
!                        EXIT LOOP2
!                    END IF
!                END DO LOOP2

            END DO LOOP1


            JC=M_END-M_START+1
!            WRITE(*,*)I1,I2,JC
!            PAUSE
!    OPEN(1,FILE='AMT3.TXT')
!    DO K = M_START ,   M_END
!        WRITE(1,'(4E15.4)') (AMT3(I,K),I=1,4)
!    END DO
!    CLOSE(1)
!    PAUSE  
            CALL SOLVE3(JC)
!    OPEN(1,FILE='AMT3.TXT')
!    DO K = ML ,   MU
!        WRITE(1,'(4E15.4)') (AMT3(J,K),J=1,4)
!    END DO
!    CLOSE(1)
!    PAUSE  
            DO I = M_START+1 , M_END-1
                C1(I,J)=AMT3(4,I)
                !IF ( C1(I,J) .LT. 0.D0 .OR.  C1(I,J)>=C_LIMIT) C1(I,J)=0.D0
            END DO
#ifdef PRE_MOR
            !DO I = M_START+1 , M_END-1
            !    !U = P2(I,J)/(H2(I,J)+H2(I,J+1))+P2(I,J-1)/(H2(I,J)+H2(I,J-1))
            !    !V = Q2(I,J)/(H2(I,J)+H2(I+1,J))+Q2(I-1,J)/(H2(I,J)+H2(I-1,J))
            !    IF ( NDD(I,J)==-1 ) CYCLE
            !    CALL UV_CACULATION(I,J,U,V,UV)
            !    TB= FLOWDENSITY*CD*UV**2
            !    QE= MAX( 0.D0 , ME*(TB/TE-1.D0) )
            !    QD = C1(I,J)*WS*MAX( 0.D0 , (1.D0-TB/TD) )
            !    !DEP(I,J)=DEP(I,J)+MOR_FACTOR*(QE-QD)*DTT/SANDDENSITY
            !    DEP(I,J)=DEP(I,J)+MOR_FACTOR*(QE-QD)*DT*0.001
            !END DO
#endif


        END DO LOOP0
                
        RETURN
    END SUBROUTINE TWO_D_SEDIMENTY
#endif