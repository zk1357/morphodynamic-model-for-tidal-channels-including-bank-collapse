#include 'define.h' 
    
! �����½��ڻ�����    
SUBROUTINE OUTPUT(BEFORE_STEP)
        USE TIME_DIN
        
#ifdef PRE_HYD        
        USE HYDRO
#endif

#ifdef PRE_AD
        USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
        USE BANK_EROSION
        USE DEF_EXCHANGE
        USE SED_EXCHANGE
#endif

#ifdef PRE_BANK_COLLAPSE
        USE BANK_COLLAPSE
#endif
    
        IMPLICIT NONE
        REAL*8               :: TEMP, SHIELD, X
        INTEGER              :: I, J, K, BEFORE_STEP
        CHARACTER*40         :: FILENAME, FILENAME1
        
        ! EVERY TIME STEP 
        IF ( GLOBTIME%GLOBSTEP%SECSTEP >  BEFORE_STEP  ) THEN               
            !CALCULATING Q~A            
            !IF ( GLOBTIME%GLOBSTEP%DAYSTEP/10== &
            !&    GLOBTIME%GLOBSTEP%DAYSTEP/10.0 ) THEN
            !    !���û�ж�����ɳģ�飬���MOR_FACTOR_BANK�ĵ�
            !    CALL CACULATE_Q_A(GLOBTIME%GLOBSTEP%DAYSTEP/10+1)    
            !END IF
                
            !CALL CACULATE_Q_A(1)    

            !MATLABˮ����Ϣ��ͼ�����Գ���ʱʹ��
            !IF ( GLOBTIME%GLOBSTEP%SECSTEP/100==GLOBTIME%GLOBSTEP%SECSTEP/100.0) CALL OUTPUT_HYD_FIELD_M
            
#ifdef PRE_AD  

            !IF ( GLOBTIME%GLOBSTEP%SECSTEP/100==GLOBTIME%GLOBSTEP%SECSTEP/100.0) CALL OUTPUT_SSC_FIELD_M
            !MATLAB���㺬ɳ����Ϣ��ͼ
            !IF ( GLOBTIME%GLOBSTEP%SECSTEP/100==GLOBTIME%GLOBSTEP%SECSTEP/100.0) THEN
            !    DO J = 150 , 390 , 10
            !        DO I = 47 , 55
            !            CALL OUTPUT_SED_M(I,J,GLOBTIME%GLOBSTEP%SECSTEP/100,0)
            !        END DO
            !    END DO
            !END IF
                
            !DO J = 150 , 230 , 10
            !    CALL OUTPUT_SED_M(50,J,GLOBTIME%GLOBSTEP-BEFORE_STEP,0)
            !END DO
            !DO I = 51 , 54 
            !    CALL OUTPUT_SED_M(I,220,GLOBTIME%GLOBSTEP-BEFORE_STEP,0)
            !END DO
            !
            !CALL OUTPUT_SED_M(52,0,GLOBTIME%GLOBSTEP-BEFORE_STEP,1)
#endif
             
            !EVERY HOUR                
            IF ( GLOBTIME%GLOBMIN==0 .AND. GLOBTIME%GLOBSEC==0 ) THEN  
                WRITE(*,'(2X, 1A, 1I4, 1A, 2X, 1I2, 1A)')'��'  , GLOBTIME%GLOBSTEP%DAYSTEP , '��' , GLOBTIME%GLOBHOUR , 'Сʱ'
                
                !!BVP OUTPUT
                !CALL OUTPUT_BVP
                
                !DO I = 20,400,50
                !    CALL OUTPUT_HYD_M(51,I,GLOBTIME%GLOBSTEP%HOURSTEP,0)
                !END DO
                   
                !MATLAB������ͼ
            !     IF (DAY_STEP/5==DAY_STEP/5.0) &
            !&    
                !CALL OUTPUT_HYD_FIELD_M                                      

#ifdef PRE_AD                    
                !MATLAB��ɳ������ͼ
                !CALL OUTPUT_SSC_FIELD_M

                !DO I = 1, NUMBER_BV/2
                !    CALL OUTPUT_SED_M(I,BVP(I,1),BVP(I,2),GLOBTIME%GLOBSTEP%HOURSTEP,0)
                !END DO
                !DO I = 300,340,5
                !    CALL OUTPUT_SED_M(51,I,GLOBTIME%GLOBSTEP%HOURSTEP,0)
                !END DO
                !DO I = 50,250,50
                !    CALL OUTPUT_SED_M(51,I,GLOBTIME%GLOBSTEP%HOURSTEP,0)
                !END DO

#endif

#ifdef PRE_BANK_EROSION   
                ! OUTPUT BANK PROFILE
                !CALL OUTPUT_POSITION
                
                !��¼�����߱����ٹ���
                !CALL OUTPUT_HYD_BVP_M                
#endif
                    
#ifdef PRE_BANK_COLLAPSE

                !���̮�����롢ˮ����ʴ�����������
                !CALL OUTPUT_COLLAPSE_M
                    
                !!���ÿ������ÿ��COLLAPSE EROSION 
                !CALL OUTPUT_CONTRIBUTION_M                   
#endif                                                  
                      
                !EVERY DAY
                IF ( GLOBTIME%GLOBHOUR==0) THEN 
                        
                    !OUTPUT Q~A RELATIONSHIP
                    !CALL OUTPUT_Q_A  
                        
                    !OUTPUT accumulated QE AND BANK HEGIHT
                    !CALL OUTPUT_QE
                        
                    !ƽ�����Σ��Գƶ���
                    !CALL BED_SY_SMOOTH
                    
                    !MATLAB��ò�ݱ���ͼ
                    CALL OUTPUT_MOR_M
                    
                    !BVP OUTPUT
                    !CALL OUTPUT_BVP
                    
                    !CALL OUTPUT_HYD_FIELD_M 
#ifdef PRE_AD                      
                    !CALL OUTPUT_SSC_FIELD_M
#endif  

#ifdef PRE_BANK_EROSION                    

                    !���ÿ������ÿ��COLLAPSE EROSION 
                    CALL OUTPUT_CONTRIBUTION_M
                        
                    !���bvp����
                    !CALL OUTPUT_BVPN
#endif  

#ifdef PRE_BANK_COLLAPSE
                    ! OUTPUT CBC
                    CALL OUTPUT_BANKEROSION_VOLUME
#endif 
                    !WRITE(*,*)'��',GLOBTIME%GLOBSTEP%DAYSTEP,'��'
                END IF
            END IF

        END IF 
        

        

        RETURN
END SUBROUTINE OUTPUT