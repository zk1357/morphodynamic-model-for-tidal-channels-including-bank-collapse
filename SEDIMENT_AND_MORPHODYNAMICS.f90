SUBROUTINE SEDIMENT_AND_MORPHODYNAMICS
!    USE INTERFACE_DEFINITIONS
    USE ARRAYSIZE
    USE HYDRO
    USE DEPTH
    USE H_PARAMETERS
    USE Morphodynamics
    IMPLICIT NONE
    INTEGER   , PARAMETER                :: SAM=SELECTED_REAL_KIND(8)
    INTEGER                              :: I , J , STEP
!    REAL(SAM) , DIMENSION(:,:) , POINTER :: PP1 , QQ1 , HH1 , CC1 , CC2 , DXX , DYY 
    REAL(SAM)                            :: DTDX , DTDY , TEMP
    REAL(SAM)                            :: CX , CY , TX , TY
    REAL(SAM)                            :: A1 , A2 , A3 , A4 , A5 , &
                                          & B1 , B2 , B3 , B4 , B5 , &
                                          & QE , QD , TB , A , B , C  , D  ,  D1 ,   &
                                          & TXW , TXE , TYN , TYS , UL , UR , VU , VD

!        TX   = KX*DT/DX1(1,1)**2
!        TY   = KY*DT/DY1(1,1)**2
!        DTDX = DT/DX1(1,1)
!        DTDY = DT/DY1(1,1)
        
!        IF ( P1(10,2)>=0.D0 ) C1(M_START+1:M_END-1,N_START:N_START+1)=0.5D0
!        IF ( P1(10,2)< 0.D0 ) C1(M_START+1:M_END-1,N_START:N_START+1)=0.D0
        DO I = M_START+1 , M_END-1
            IF (P1(I,2)>=0.D0) THEN
                C1(I,1:2)=0.5D0
            ELSE
                C1(I,1:2)=0.D0
            END IF
        END DO

SAM0:   DO I = M_START+1 , M_END-1
           DO J = N_START+2 , N_END-1  
!---------------------------------------------------------------------------------------------------------------
                UR=2.D0*P1(I,J  )/(H1(I,J)+H1(I,J+1))
                UL=2.D0*P1(I,J-1)/(H1(I,J)+H1(I,J-1))
                VU=2.D0*Q1(I,J  )/(H1(I,J)+H1(I+1,J))
                VD=2.D0*Q1(I-1,J)/(H1(I,J)+H1(I-1,J))
                A =DT/DX1(1,1)*UR
                A1=DT/DY1(1,1)*VU
                B =0.5D0*( C1(I,J+1)+C1(I,J) )
                B1=0.5D0*( C1(I+1,J)+C1(I,J) )
                IF (UR>0.D0) THEN
                    C=-0.125D0*( C1(I,J+1)-2.D0*C1(I,J  )+C1(I,J-1) )
                    TXE=A*(B+C)
                ELSE IF (UR==0.D0) THEN
                    TXE=0.D0
                ELSE
                    C=-0.125D0*( C1(I,J+2)-2.D0*C1(I,J+1)+C1(I,J  ) )
                    TXE=A*(B+C)
                END IF
                
                IF (VU>0.D0) THEN
                    C=-0.125D0*( C1(I+1,J)-2.D0*C1(I,J  )+C1(I-1,J) )
                    TYN=A1*(B1+C)
                ELSE IF (VU==0.D0) THEN
                    TYN=0.D0
                ELSE
                    C=-0.125D0*( C1(I+2,J)-2.D0*C1(I+1,J)+C1(I,J  ) )
                    TYN=A1*(B1+C)
                END IF

                A =DT/DX1(1,1)*UL
                A1=DT/DY1(1,1)*VD
                B =0.5D0*( C1(I,J-1)+C1(I,J) )
                B1=0.5D0*( C1(I-1,J)+C1(I,J) )
                IF (UL>0.D0) THEN
                    C=-0.125D0*( C1(I,J  )-2.D0*C1(I,J-1)+C1(I,J-2) )
                    TXW=A*(B+C)
                ELSE IF (UL==0.D0) THEN
                    TXW=0.D0
                ELSE
                    C=-0.125D0*( C1(I,J+1)-2.D0*C1(I,J  )+C1(I,J-1) )
                    TXW=A*(B+C)
                END IF

                IF (VD>0.D0) THEN
                    C=-0.125D0*( C1(I  ,J)-2.D0*C1(I-1,J  )+C1(I-2,J) )
                    TYS=A1*(B1+C)
                ELSE IF (VD==0.D0) THEN
                    TYS=0.D0
                ELSE
                    C=-0.125D0*( C1(I+1,J)-2.D0*C1(I  ,J  )+C1(I-1,J ) )
                    TYS=A1*(B1+C)
                END IF


                D =KX*DT/DX1(1,1)**2*( C1(I,J+1)-C1(I,J)-( C1(I,J)-C1(I,J-1) ) )
                D1=KY*DT/DY1(1,1)**2*( C1(I+1,J)-C1(I,J)-( C1(I,J)-C1(I-1,J) ) )
                TB = FLOWDENSITY*CD*( ( P1(I,J  )/(H1(I,J)+H1(I,J+1))+P1(I,J-1)/(H1(I,J)+H1(I,J-1)) )**2 )
!                QE = ME*(TB/TE-1)
                QD = C1(I,J)*WS*(1-TB/TD)
                C2(I,J)=C1(I,J)-TXE+TXW-TYN+TYS+D+D1-QD



!                A1 = ( CX**2/6.D0-0.5D0*CX+1.D0/3.D0+TX )*CX
!                A2 = ( -CX**2/3.D0+0.5D0*CX+5.D0/6.D0-0.5D0*CX*CY-0.5D0*CY**2+0.5D0*CY-2.D0*TX-2.D0*TY )*CX
!                A3 = ( -1.D0/6.D0+CX**2/6.D0+TX )*CX
!                A4 = ( -0.5D0*CY+0.5D0*CY**2+TY )*CX
!                A5 = ( 0.5D0*CX*CY+TY )*CX
!                B1 = ( CY**2/6.D0-0.5D0*CY+1.D0/3.D0+TY )*CY
!                B2 = ( -CY**2/3.D0+0.5D0*CY+5.D0/6.D0-0.5D0*CX*CY-0.5D0*CX**2+0.5D0*CX-2.D0*TY-2.D0*TX )*CY
!                B3 = ( -1.D0/6.D0+CY**2/6.D0+TY )*CY
!                B4 = ( -0.5D0*CX+0.5D0*CX**2+TX )*CY
!                B5 = ( 0.5D0*CX*CY+TX )*CY
!                CX =2.D0*P1(I,J  )/( H1(I,J)+H1(I,J+1) )*DTDX
!                IF (CX==0.D0) THEN
!                    TXE=0.D0
!                ELSE IF (CX>0.D0) THEN
!                    TXE = A1*C1(I,J+1)+A2*C1(I,J  )+A3*C1(I,J-1)+A4*C1(I+1,J  )+A5*C1(I-1,J  )-TX*C1(I,J+1)+TX*C1(I,J  )
!                ELSE 
!                    TXE = A1*C1(I,J  )+A2*C1(I,J+1)+A3*C1(I,J+2)+A4*C1(I+1,J+1)+A5*C1(I-1,J+1)-TX*C1(I,J+1)+TX*C1(I,J  )
!                END IF
!
!                IF (CY==0.D0) THEN
!                    TYN=0.D0
!                ELSE
!                    TYN = B1*C1(I+1,J)+B2*C1(I  ,J)+B3*C1(I-1,J)+B4*C1(I  ,J+1)+B5*C1(I  ,J-1)-TY*C1(I+1,J)+TY*C1(I  ,J)
!                END IF

!                CX = 2.D0*P1(I,J-1)/( H1(I,J)+H1(I,J-1) )*DTDX
!                CY = 2.D0*Q1(I-1,J)/( H1(I,J)+H1(I-1,J) )*DTDY
!                
!                A1 = ( CX**2/6.D0-0.5D0*CX+1.D0/3.D0+TX )*CX
!                A2 = ( -CX**2/3.D0+0.5D0*CX+5.D0/6.D0-0.5D0*CX*CY-0.5D0*CY**2+0.5D0*CY-2.D0*TX-2.D0*TY )*CX
!                A3 = ( -1.D0/6.D0+CX**2/6.D0+TX )*CX
!                A4 = ( -0.5D0*CY+0.5D0*CY**2+TY )*CX
!                A5 = ( 0.5D0*CX*CY+TY )*CX
!                B1 = ( CY**2/6.D0-0.5D0*CY+1.D0/3.D0+TY )*CY
!                B2 = ( -CY**2/3.D0+0.5D0*CY+5.D0/6.D0-0.5D0*CX*CY-0.5D0*CX**2+0.5D0*CX-2.D0*TY-2.D0*TX )*CY
!                B3 = ( -1.D0/6.D0+CY**2/6.D0+TY )*CY
!                B4 = ( -0.5D0*CX+0.5D0*CX**2+TX )*CY
!                B5 = ( 0.5D0*CX*CY+TX )*CY
!                CX = 2.D0*P1(I,J-1)/( H1(I,J)+H1(I,J-1) )*DTDX
!                IF (CX==0.D0) THEN
!                    TXW=0.D0
!                ELSE IF (CX>0.D0) THEN
!                    TXW = A1*C1(I,J  )+A2*C1(I,J-1)+A3*C1(I,J-2)+A4*C1(I+1,J-1)+A5*C1(I-1,J-1)-TX*C1(I,J  )+TX*C1(I,J-1)
!                ELSE
!                    TXW = A1*C1(I,J-1)+A2*C1(I,J  )+A3*C1(I,J+1)+A4*C1(I+1,J  )+A5*C1(I-1,J  )-TX*C1(I,J  )+TX*C1(I,J-1)
!                END IF
!                IF (CY==0.D0) THEN
!                    TYS=0.D0
!                ELSE
!                    TYS = B1*C1(I  ,J)+B2*C1(I-1,J)+B3*C1(I-2,J)+B4*C1(I-1,J+1)+B5*C1(I-1,J-1)-TY*C1(I  ,J)+TY*C1(I-1,J)
!                END IF

!                TB = FLOWDENSITY*CD*( ( PP1(I,J  )/(HH1(I,J)+HH1(I,J+1))+PP1(I,J-1)/(HH1(I,J)+HH1(I,J-1)) )**2 + & 
!                &                     ( QQ1(I,J  )/(HH1(I,J)+HH1(I+1,J))+QQ1(I-1,J)/(HH1(I,J)+HH1(I-1,J)) )**2 )
!                QE = ME*(TB/TE-1)
!                QD = CC1(I,J)*WS*(1-TB/TD)
!                C2(I,J)=C1(I,J)+TXW-TXE+TYS-TYN

!---------------------------------------------------------------------------------------------------------------

            END DO
        END DO SAM0

SAM1:   DO I = M_START+1 , M_END-1
           DO J = N_START+2 , N_END-1  
                TEMP=C1(I,J)
                C1(I,J)=C2(I,J)   
                C2(I,J)=TEMP
            END DO
        END DO SAM1
        !写泥沙源函数，即冲淤公式，以及床面变形方程
        !检查上述公式是否正确
        !边界条件书写，固壁前的处理，规划起算点位置
OPEN (1,FILE='含沙量.TXT')
    WRITE(1,'(<N_END>F12.4)')((C1(I,J),J=N_START,N_END),I=M_START,M_END)
CLOSE(1)


END SUBROUTINE SEDIMENT_AND_MORPHODYNAMICS