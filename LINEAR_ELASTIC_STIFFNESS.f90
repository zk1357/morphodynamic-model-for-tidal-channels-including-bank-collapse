#include 'define.h' 
#ifdef PRE_BANK_COLLAPSE   
SUBROUTINE LINEAR_ELASTIC_STIFFNESS(NUMBER_I,NUMBER_J)
    
#ifdef PRE_BANK_EROSION
        USE BANK_EROSION
        !USE EXCHANGE
        USE DEF_EXCHANGE
        !USE SED_EXCHANGE
#endif

#ifdef PRE_BANK_COLLAPSE
        USE BANK_COLLAPSE
#endif
        IMPLICIT NONE
        
        INTERFACE
            SUBROUTINE pardiso_sym_f90( N , NNZ , A , JA , IA , B , X )
                IMPLICIT NONE
                INTEGER , PARAMETER                        :: PSF=SELECTED_REAL_KIND(8)
                INTEGER                                    :: N , NNZ
                REAL(PSF) , ALLOCATABLE , DIMENSION(:)     :: A , B , X
                INTEGER   , ALLOCATABLE , DIMENSION(:)     :: JA, IA
            END SUBROUTINE pardiso_sym_f90
        END INTERFACE
        
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: NODE_ST     !记录单元劲度矩阵在整体劲度矩阵的位置
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: B           !应变矩阵
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: BT          !转置
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: D           !弹性矩阵
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: H
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KE
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: KT
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: STRESS 
        REAL*8, ALLOCATABLE, DIMENSION(:)     :: A_LE        !一维数组存储非零值
        REAL*8, ALLOCATABLE, DIMENSION(:)     :: F, WY       !节点力与位移
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: GRID , MESH ,EXPORT
        INTEGER, ALLOCATABLE, DIMENSION(:)    :: B_LE , C_LE !一维数值存储行列
        INTEGER, ALLOCATABLE, DIMENSION(:)    :: TEMP_C
        INTEGER, ALLOCATABLE, DIMENSION(:)    :: DIG_NUMBER
        INTEGER                               :: N , NNZ , NUMBER_I , NUMBER_J
        INTEGER                               :: I , J , K
        REAL*8                                :: TEMP
        CHARACTER*40                          :: FILENAME , FILENAME_PATH

        
        !FILENAME=''
        !WRITE(FILENAME(1:4),'(I4.4)') NUMBER_I
        !FILENAME(5:9)='F.RES'
        !FILENAME_PATH='MATLAB_PLOT\'//FILENAME
        !OPEN (UNIT=100 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2*2)
        !DO I = 1 , EXCHANGE1%NODE
        !    WRITE(100 , REC=I )BANKCOLLAPSE%F(NUMBER_I,2*I-1),BANKCOLLAPSE%F(NUMBER_I,2*I)
        !END DO
        !CLOSE(100)   
        !!WRITE(*,*)1 
        !PAUSE
        
        !WRITE(*,'(12F6.2)') (BANKCOLLAPSE%F(NUMBER_I,I),I=1,12)
        !WRITE(*,'(12I6)') (BANKCOLLAPSE%RESTRICTION(I),I=1,12)
        !PAUSE         
        
        ! CHECK CHANGE OF ECTERNAL LOAD
        IF (NUMBER_J .NE. -1) THEN 
        
            !求单元、整体劲度矩阵
        
            ALLOCATE( B(3,6) )
            ALLOCATE( BT(6,3) )
            ALLOCATE( D(3,3) )
            ALLOCATE( NODE_ST(6,EXCHANGE1%ELEMENT) )
            ALLOCATE( H(3,6,EXCHANGE1%ELEMENT) )
            ALLOCATE( KE(6,6,EXCHANGE1%ELEMENT) )
            ALLOCATE( KT(2*EXCHANGE1%NODE,2*EXCHANGE1%NODE) )
            ALLOCATE( F(2*EXCHANGE1%NODE) )
            ALLOCATE( WY(2*EXCHANGE1%NODE) )
        
            B      =0.D0
            BT     =0.D0
            NODE_ST=0.D0
            KT     =0.D0
        
            !求解每个单元弹性矩阵
            DO I = 1 , EXCHANGE1%ELEMENT 
                !开挖单元以空气填充
                IF ( EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,11) .NE. 0 ) &
                &    BANKCOLLAPSE%MODULUS(NUMBER_I,I)=0.D0 
                D=0.D0
                D(1,1)=BANKCOLLAPSE%MODULUS(NUMBER_I,I)/(1.D0-BANKCOLLAPSE%MU**2)
                D(2,2)=BANKCOLLAPSE%MODULUS(NUMBER_I,I)/(1.D0-BANKCOLLAPSE%MU**2)
                D(1,2)=BANKCOLLAPSE%MODULUS(NUMBER_I,I)/(1.D0-BANKCOLLAPSE%MU**2)*BANKCOLLAPSE%MU
                D(2,1)=BANKCOLLAPSE%MODULUS(NUMBER_I,I)/(1.D0-BANKCOLLAPSE%MU**2)*BANKCOLLAPSE%MU
                D(3,3)=BANKCOLLAPSE%MODULUS(NUMBER_I,I)/(1.D0+BANKCOLLAPSE%MU)*0.5D0
            
                B=0.D0
                B(1,1)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,7)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,9)     
                B(1,3)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,9)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,5)     
                B(1,5)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,5)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,7)
                B(2,2)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,8)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,6)     
                B(2,4)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,4)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,8)     
                B(2,6)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,6)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,4)
                B(3,1)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,8)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,6)     
                B(3,2)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,7)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,9)     
                B(3,3)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,4)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,8)
                B(3,4)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,9)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,5)     
                B(3,5)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,6)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,4)     
                B(3,6)=EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,5)-EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,7)

                NODE_ST(1,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,1)-1 
                NODE_ST(2,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,1)   
                NODE_ST(3,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,2)-1
                NODE_ST(4,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,2)   
                NODE_ST(5,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,3)-1 
                NODE_ST(6,I)=2.D0*EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,3)

                BT=TRANSPOSE(B)

                DO J = 1 , 3
                    DO K = 1 , 6
                        H(J,K,I)=D(J,1)*B(1,K)+D(J,2)*B(2,K)+D(J,3)*B(3,K)
                    END DO
                END DO
                
                DO J = 1 , 6
                    DO K = 1 , 6
                        KE(J,K,I)=0.25D0*BANKCOLLAPSE%T_LE*(BT(J,1)*H(1,K,I)+BT(J,2)*H(2,K,I)+BT(J,3)*H(3,K,I))/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)
                    END DO
                END DO
                
            END DO

            DO I = 1 , EXCHANGE1%ELEMENT
                DO J = 1 , 6
                    DO K = 1 , 6
                        KT(NODE_ST(J,I),NODE_ST(K,I))=KT(NODE_ST(J,I),NODE_ST(K,I))+KE(J,K,I)
                    END DO
                END DO
            END DO
            
            
            !WRITE(*,'(12F6.2)') ((KT(I,J)*4,I=1,12),J=1,12)
            !PAUSE            
            
            DO I = 1 , 2*EXCHANGE1%NODE
                IF (BANKCOLLAPSE%RESTRICTION(I)==0) THEN
                    BANKCOLLAPSE%F(NUMBER_I,I)=0.D0
                    DO J = 1 , 2*EXCHANGE1%NODE
                        IF (I .NE. J) THEN
                            KT(I,J)=0.D0
                            KT(J,I)=0.D0
                        ELSE
                            KT(I,J)=1.D0
                        END IF
                    END DO
                END IF
            END DO
        
            DO I = 1, 2*EXCHANGE1%NODE
                IF (KT(I,I)==0.D0) THEN
                    BANKCOLLAPSE%F(NUMBER_I,I)=0.D0
                    KT(I,I)=1.D0
                END IF
            END DO

            !WRITE(*,'(12F6.2)') ((KT(I,J)*4,I=1,12),J=1,12)
            !PAUSE 
            
            !WRITE(*,'(12F6.2)') (BANKCOLLAPSE%F(NUMBER_I,I),I=1,12)
            !PAUSE  
            
            
            !以下为用一维数组存储大型稀疏矩阵
            K=0                 !找出非零元素个数   
            DO I = 1 , 2*EXCHANGE1%NODE
                DO J = I , 2*EXCHANGE1%NODE
                    IF (KT(I,J) .NE. 0.D0) THEN
                        K=K+1
                    END IF
                END DO
            END DO        
        
            ALLOCATE(A_LE(K),B_LE(K),TEMP_C(K),C_LE(2*EXCHANGE1%NODE+1))
            NNZ=K
            A_LE     =0.D0
            B_LE     =0
            B_LE     =0
            TEMP_C   =0
        
            K=1
            DO I = 1 , 2*EXCHANGE1%NODE   !找出非零元素
                DO J = I , 2*EXCHANGE1%NODE !对称矩阵只需存储上三角矩阵
                    IF (KT(I,J) .NE. 0.D0) THEN
                        A_LE(K)  =KT(I,J)
                        B_LE(K)  =J
                        TEMP_C(K)=I
                        K=K+1
                    END IF
                END DO
            END DO  
        
            K      =2
            C_LE(1)=1
            DO I = 2 , NNZ
                IF (TEMP_C(I-1) .NE. TEMP_C(I)) THEN
                    C_LE(K)=I         
                    K=K+1
                END IF
            END DO
            C_LE(2*EXCHANGE1%NODE+1)=NNZ+1

            N=2*EXCHANGE1%NODE    
            DO I = 1, 2*EXCHANGE1%NODE
                !WY(I) = 0.D0
                F(I)  = BANKCOLLAPSE%F(NUMBER_I,I)
            END DO
            
            !WRITE(*,'(12F6.2)') (F(I),I=1,12)
            !WRITE(*,'(12F6.2)') (BANKCOLLAPSE%F(NUMBER_I,I),I=1,12)
            !PAUSE 
            
            
            !FILENAME=''
            !WRITE(FILENAME(1:4),'(I4.4)') NUMBER_I
            !FILENAME(5:9)='F.RES'
            !FILENAME_PATH='MATLAB_PLOT\'//FILENAME
            !OPEN (UNIT=100 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2*2)
            !DO I = 1 , EXCHANGE1%NODE
            !    WRITE(100 , REC=I )BANKCOLLAPSE%F(NUMBER_I,2*I-1),BANKCOLLAPSE%F(NUMBER_I,2*I)
            !END DO
            !CLOSE(100)
            
            
            !矩阵求解
            CALL pardiso_sym_f90( N , NNZ , A_LE , B_LE , C_LE , F , WY )

            !WRITE(*,'(12F6.2)') (WY(I),I=1,N)
            !PAUSE 
            
            !WRITE(FILENAME(1:4),'(I4.4)') NUMBER_I
            !FILENAME(5:10)='WY.RES'
            !FILENAME_PATH='MATLAB_PLOT\'//FILENAME
            !OPEN (UNIT=100 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2*4)
            !DO I = 1 , EXCHANGE1%NODE
            !    WRITE(100 , REC=I ) WY(2*I-1), WY(2*I), F(2*I-1), F(2*I)
            !END DO
            !CLOSE(100)  
            
            !DO I=1,EXCHANGE1%ELEMENT
            !    WRITE(*,'(6F6.2)') ((H(J,K,I),K=1,6),J=1,3)
            !    WRITE(*,'(3F6.2)')            
            !END DO
            !WRITE(*,'(11F6.2)')((EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,J),J=1,11),I=1,EXCHANGE1%ELEMENT )
            !PAUSE 
            
            !求各单元应力
            DO I = 1 , EXCHANGE1%ELEMENT 
                
                DO J = 1 , 3
                    !单位为KPA
                    BANKCOLLAPSE%TSTRESS(NUMBER_I,I,J)=BANKCOLLAPSE%TSTRESS(NUMBER_I,I,J)                                          &
                &  +5.D-1*H(J,1,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,1)*2-1)/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)  &          
                &  +5.D-1*H(J,2,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,1)*2  )/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)  &
                &  +5.D-1*H(J,3,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,2)*2-1)/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)  &
                &  +5.D-1*H(J,4,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,2)*2  )/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)  &
                &  +5.D-1*H(J,5,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,3)*2-1)/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10)  &
                &  +5.D-1*H(J,6,I)*WY(EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,3)*2  )/EXCHANGE1%EXCHANGE_BANK_SIZE(NUMBER_I,I,10) 
                END DO

                IF (ISNAN(BANKCOLLAPSE%TSTRESS(NUMBER_I,I,3))) THEN
                    WRITE(*,*)NUMBER_I,NUMBER_J , I
                END IF
                
            END DO
            
            !求各单元累积位移
            DO I = 1, 2*EXCHANGE1%NODE
                BANKCOLLAPSE%WY(NUMBER_I,I) = BANKCOLLAPSE%WY(NUMBER_I,I) + WY(I)
            END DO
            
            ! CHECK FAILURE
            CALL FAILURE_CRITERION(NUMBER_I)            
            
            ! OUTPUT STRESS FIELD
            !CALL OUTPUT_STRESS_M(NUMBER_I)
            
            ! NUMBER OF LOAD +1 
			EXCHANGE1%NUMBER_LOAD(NUMBER_I) = EXCHANGE1%NUMBER_LOAD(NUMBER_I)+1
            
            IF ( ALLOCATED( B ) )      DEALLOCATE( B )
            IF ( ALLOCATED( BT ) )      DEALLOCATE( BT )
            IF ( ALLOCATED( D ) )      DEALLOCATE( D )
            IF ( ALLOCATED( NODE_ST ) )      DEALLOCATE( NODE_ST )
            IF ( ALLOCATED( KE ) )      DEALLOCATE( KE )
            IF ( ALLOCATED( KT ) )      DEALLOCATE( KT )
            IF ( ALLOCATED( TEMP_C ) )      DEALLOCATE( TEMP_C )
            IF ( ALLOCATED( H ) )      DEALLOCATE( H )
            !位移wy、F是否需要移除，以及总位移是否需要重新计算
            IF ( ALLOCATED( WY ) )      DEALLOCATE( WY )
            IF ( ALLOCATED( F ) )      DEALLOCATE( F )
            
        END IF

        RETURN

    END SUBROUTINE LINEAR_ELASTIC_STIFFNESS
#endif