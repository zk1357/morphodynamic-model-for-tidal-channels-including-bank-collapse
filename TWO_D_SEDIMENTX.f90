#include 'define.h'
#ifdef PRE_AD

SUBROUTINE TWO_D_SEDIMENTX  !写守恒形式的,加水深
    
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
    REAL*8                              :: DTT , TEMP , QD , TB , TB_BANK , QE , QE_BANK , &
                                        &  A , B , CX , CY , U , V  , UV , CHEZY, &
                                        &  DXX , DYY , HL , HR
    INTEGER                             :: I , J , K , P , JC
    
    DTT=0.5D0*DT

LOOP0:  DO I = M_START+1 , M_END-1  
            AMT3=0.D0
            AMT3(2,:)=1.D0
LOOP1:      DO J = N_START , N_END
                !固壁边界
                IF ( NDD(I,J)==-1 .OR. NDD(I,J)==1 ) THEN
                    AMT3(4,J)=C1(I,J)
                    CYCLE LOOP1 
                END IF
                
                !水深小于临界值，不算
                IF (H1(I,J)<= FLOODINGDEPTHX) THEN
                    AMT3(4,J)=C1(I,J)
                    CYCLE LOOP1 
                END IF                    
                
                !开边界，梯度为0,尚未调试好
                IF ( NDD(I,J)==4 ) THEN
                    IF ( J==N_START ) THEN
                        AMT3(3,J)=-1 
                    ELSE 
                        AMT3(1,J)=-1
                    END IF
                    CYCLE LOOP1
                END IF
                
                !水位边界
                IF ( NDD(I,J)==2) THEN 
                    IF (P1(I,1)>=0.D0) THEN
                        AMT3(4,J)=C_BOUNDARY
                        !AMT3(3,J)=-1
                    ELSE
                         AMT3(3,J)=-1
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
                   
                AMT3(1,J)= -(1.D0+A)*CX
                AMT3(2,J)=  (1.D0+2.D0*CX*A)+WS*DTT*MAX(0.D0,(1.D0-TB/TD))/H1(I,J)
                AMT3(3,J)=  (1.D0-A)*CX
                AMT3(4,J)=  C1(I,J)+QE*DTT/H1(I,J)-CY*( (1.D0-B)*C1(I+1,J)+2.D0*B*C1(I,J)-(1.D0+B)*C1(I-1,J) )+DXX+DYY
                IF ( NDD(I,J-1)==-1 ) AMT3(1,J)=0.D0
                IF ( NDD(I,J+1)==-1 ) AMT3(3,J)=0.D0
                
#ifdef PRE_BANK_EROSION
LOOP2:          DO K = 1 , NUMBER_BV   !加上的含沙量不应该算在底床冲刷里面
                    IF ( I==BVP(K,1) .AND. J==BVP(K,2)  ) THEN
                        DO P = 1 , 5
                            AMT3(4,J-3+P)=AMT3(4,J-3+P)+S_BANK_EROSION(K)
                        END DO
                        
                        !IF (GLOBTIME%GLOBSTEP%DAYSTEP/10==GLOBTIME%GLOBSTEP%DAYSTEP/10.0) THEN
                        !    QE_OUT( 3*K-1, GLOBTIME%GLOBSTEP%DAYSTEP/10+1 ) = & 
                        !&   QE_OUT( 3*K-1, GLOBTIME%GLOBSTEP%DAYSTEP/10+1 )+ QE*DTT/H1(I,J)*2
                        !END IF
                        
                        EXIT LOOP2
                    END IF
                END DO LOOP2
#endif
            END DO LOOP1


            JC=N_END-N_START+1
!            WRITE(*,*)I1,I2,JC
!            PAUSE
    !OPEN(1,FILE='AMT.TXT')
    !DO K = 1 ,  N_END 
    !    WRITE(1,'(4E15.4)') (AMT(J,K),J=1,4)
    !END DO
    !CLOSE(1)
    !PAUSE  
            CALL SOLVE3(JC)
    !OPEN(1,FILE='AMT.TXT')
    !DO K = 1 ,  N_END
    !    WRITE(1,'(4E15.4)') (AMT(J,K),J=1,4)
    !END DO
    !CLOSE(1)
    !PAUSE  
            DO J = N_START , N_END
                C1(I,J)=AMT3(4,J)
                !IF ( C1(I,J) .LT. 0.D0 .OR.  C1(I,J)>=C_LIMIT) C1(I,J)=0.D0
            END DO


        END DO LOOP0
        
        !DO J = N_START , N_END
        !    C1(M_START,J)=C1(M_START+1,J)
        !    C1(M_END,J)=C1(M_END-1,J)
        !END DO
        
!        OPEN(100,FILE='C.TXT')
!        WRITE(100,'(100E12.4)')( C1(2,J) , J=1,100)  
!        CLOSE(100)
!        PAUSE
             
        RETURN
END SUBROUTINE TWO_D_SEDIMENTX

subroutine solve3(JC)
    USE HYDRO
    IMPLICIT NONE
    INTEGER                             :: I , J , JC
    real*8                              :: a
    !REAL*8                              :: AMT(4,MAX(M_END,N_END)+2)
!        WRITE(*,*)JC
!        PAUSE
      do  i=2,jc
        a=amt3(2,i-1)
        do  j=2,4
          amt3(j,i-1)=amt3(j,i-1)/a
        end do
        a=-amt3(1,i)
        amt3(2,i)=amt3(2,i)+amt3(3,i-1)*a
        amt3(4,i)=amt3(4,i)+amt3(4,i-1)*a
      end do
      a=amt3(2,jc)
      do  j=2,4
        amt3(j,jc)=amt3(j,jc)/a
      end do

      a=amt3(3,jc-1)
      do  i=jc-1,1,-1
        amt3(4,i)=amt3(4,i)-amt3(3,i)*amt3(4,i+1)
      end do

       RETURN
    end subroutine solve3
#endif