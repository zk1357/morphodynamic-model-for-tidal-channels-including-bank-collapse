#include 'define.h'

MODULE TIME_DIN
    TYPE TIMESTEP_TYPE
        INTEGER :: SECSTEP, HOURSTEP, DAYSTEP, BEFSTEP     ! definition of second, minute, hour, BEFORE
    END TYPE

    TYPE TIME_TYPE
        TYPE(TIMESTEP_TYPE) :: GLOBSTEP
        INTEGER             :: GLOBMIN, GLOBHOUR, &           ! definition of minute, hour
        &                      GLOBHALF,                   &  ! definition of half-day
        &                      GLOBDAY, GLOBMON, GLOBYEAR     ! definition of day, month, year
        REAL                :: GLOBSEC                        ! definition of second
        END TYPE
        TYPE(TIME_TYPE)     :: GLOBTIME
END MODULE    
    
#ifdef PRE_HYD    
MODULE HYDRO 
    IMPLICIT NONE
    INTEGER,PARAMETER                           :: NP1=SELECTED_REAL_KIND(8)
    !������������ʼ����    
    INTEGER,PARAMETER                           :: M_START=1 , M_END=41 , N_START=1 , N_END=276       
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: P1,P2
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: Q1,Q2 
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: H1,H2 
    !X,Y�������񲽳�    
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: DX1,DY1 
    !ʱ�䲽��
    REAL(NP1),PARAMETER                         :: DT=6.D0
    !�¶ȱ߽�,��λ���,�߽�ˮ��,��ʪ�ж��ٽ�ˮ��
    REAL(NP1),PARAMETER                         :: SLOPE=2.D-3 , SLOPE1=0.D-3 , AMP=1.D0 , DEPOUTSEA=3.D0 , FLOODINGDEPTHX=0.1D0 , FLOODINGDEPTHY=0.1D0
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: DEP
    !��Ħ�輰��ճϵ��,Ĭ��ʹ�ó���л�Ź�ʽ
    REAL(NP1),PARAMETER                         :: MANNING=38 , CHE=38, CHE1=10, EDDY_X=2.0D0 , EDDY_Y=1.D0 !���ݼ��㣬��ֵС��1.9����
    !ˮ��,���ٽڵ�״̬
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: NDD
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: VDDX   !VELOCITY-NDD
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: VDDY
    REAL  ,   ALLOCATABLE,TARGET,DIMENSION(:,:) :: IDN
    REAL  ,   ALLOCATABLE,TARGET,DIMENSION(:,:) :: JDN   
    !Q~A��ϵ
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: QA_Q , QA_A
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: QA_Q_P , QA_A_P !��������Q~A��ϵ
    !����λ��,���ڸ���ֵ    
    INTEGER                                     :: TCI_START , TCI_END  , TCJ_START , TCJ_END
    !���శ��BVP�ܸ���
    INTEGER   ,PARAMETER                        :: NUMBER_BV=100                      
    !�������    
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:)   :: BANK_H
    !��������
    REAL(NP1),PARAMETER                         :: WIDTH_CREEK_INI=20.D0
    !�������    
    REAL(NP1)                                   :: AMT(6,3001)        
    REAL(NP1)                                   :: AMT3(4,3001)
END MODULE HYDRO 
#endif

#ifdef PRE_AD
MODULE Morphodynamics
    IMPLICIT NONE
    INTEGER,PARAMETER                                 :: NP2=SELECTED_REAL_KIND(8)
    !��ɳ��    
    REAL(NP2),ALLOCATABLE,TARGET,DIMENSION(:,:)       :: C1 , C2 
    !��������
    INTEGER                                           :: MOR_FACTOR
    !X,Y�������ɳ��ɢϵ��
    REAL(NP2),PARAMETER                               :: KX=10.D0  , KY=2.D0    
    !��ɳ����                                           MEȡֵΪDELFT3dĬ��ֵ��WSΪĬ��ֵһ�롢��ɳ�ܶ�ΪĬ��ֵ
    REAL(NP2),PARAMETER                               :: ME=5.D-5 , TE=0.5D0 , TD=1000.D0 , WS=5.D-4 , &
                                                       &  CD=0.002 , FLOWDENSITY=1000 , SANDDENSITY=2650 , &
                                                       &  C_BOUNDARY=0.D0 , C_LIMIT=1000.5D0 , PORO=0.43
END MODULE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
MODULE BANK_EROSION
    IMPLICIT NONE
    INTEGER   ,PARAMETER                     :: NP6=SELECTED_REAL_KIND(8)   
    INTEGER   ,PARAMETER                     :: COLLAPSE_START=1 , COLLAPSE_END=39 !һ���̮����ʼ��
    INTEGER   ,PARAMETER                     :: NUMBER_OB=21                       !���ڼ����ܸ���,������Ŵ�������,����20�ȷ�
    INTEGER                                  :: MOR_FACTOR_BANK
    REAL(NP6) ,PARAMETER                     :: WIDTH_BANK=10.D0
    REAL(NP6) ,PARAMETER                     :: COL_PER_BED=0.1D0                   !COLLAPSE BECOME BED IN PERCENTAGE
    REAL(NP6)                                :: BANK_START_VELOCITY                     !��������
    REAL(NP6) ,PARAMETER                     :: ME_BANK=2.D-4 , CD_BANK=0.002 , TE_BANK=0.062D0 , &
                                              & FLOWDENSITY_BANK=1000 , SANDDENSITY_BANK=1940 , PORO_BANK=0.
    REAL(NP6) ,PARAMETER                     :: BANK_EROSION_LIMIT=0.1D0 , BANK_COLLAPSE_LIMIT=0.3D0
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: POSITION 
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: QE_OUT !��¼�����ڸ��������ʴͨ��
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: BVP_VELOCITY !��¼BVP���ٹ���,��������ˮ����Ӧ��
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) :: BVP , BVP_IN        !��¼�����ڣ��������Ե����ٵ�λ��,�ֱ�ΪI,J
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: BVP_JUDGE  !�ж��Ƿ���°��ڱ߱ڵ㣬�������������һ�� 
END MODULE BANK_EROSION

MODULE EXCHANGE
    IMPLICIT NONE
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: JUDGE_BANK,& !�ж��Ƿ����Ϊ��ʼ��������
    &                                           JUDGE_TIME   !���䳱ʱ���ж� 0�ǳ� -1�䳱
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) :: BF           !��¼�ƻ��������ٷֱ�,1-erosion���� 2-CO���� 3-cycʶ�� 4-e/c�ж�
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: NUMBER_LOAD  !����ÿһ��������˵,Ϊ�Ӻ��ش�����ֱ���ƻ����°���
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: NUMBER_BANKCYCLE !��¼�ǵڼ�������,����������+1
    !ÿ�鰶����Ϣ
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_RETREATS  !ÿ�鰶�ں��˾���
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_COLLAPSE  !ÿ�鰶��̮������
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_FLOW      !ÿ�鰶����ʴ����,��������̮������ֵΪˮ����ʴ�ܾ��룬��Ϊ���ڸ����ж���Mesh�ӳ���
    !ÿ�찶����Ϣ
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LENGTH_EROSION_DAY  !ÿ��BANK EROSION����
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LENGTH_COLLAPSE_DAY  !ÿ��COLLAPSE    ����   
    !���ڳ�ʼ�������������,Ӧ�Ա߱ڸ߶ȼ���
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: ZOOMIN           
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LHL , LHLB , &   !LHL��¼���������ˮλ�߶�
    &                                           STATION_BANK , & !����ˮ����ˢ�϶�İ���,����2mû��,���°���,0����,1��
    &                                           LENGTH_COLLAPSE_BEFORE , & !��¼ǰһ������仯ʱ��̮�����룬�Ӷ�����ǰһ������ĸ߳�����ֵ
    &                                           RETREATS_TOTAL         , & !��¼ÿ��CS������ܺ��˾���
    &                                           UPDATE_C                   !�ۼ�ÿ���������ǰ��̮�������
    INTEGER                                  :: JUDGE_TEMP       !��ʱ�ж��Ƿ񻻰���
END MODULE EXCHANGE

MODULE SED_EXCHANGE
    IMPLICIT NONE
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: S_BANK_EROSION !���ڼ�¼���ڳ�ˢ����ɳ,ÿһ��ʱ�䲽������,Ϊ��ɳ��
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: S_BANK_COLLAPSE!���ڼ�¼̮��ת��Ϊ�½��ڻ�����ɳ,Ϊ��ɳ���
ENDMODULE SED_EXCHANGE
    
MODULE DEF_EXCHANGE
    IMPLICIT NONE
    TYPE EXCHANGE_TYPE
        !��¼ģ���а��ڳߴ���Ϣ
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: EXCHANGE_BANK_SIZE
        !��¼ģ�Ͱ��������Ԥ������ߴ�����ű���
        REAL*8, ALLOCATABLE, DIMENSION(:)    :: ZOOMIN 
        !�����Ƴ��ж�
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DIG_NUMBER         
        !��¼����ڵ����ꡢ��Ԫ��Ӧ�ڵ��
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: GRID, MESH 
        !��¼̮����ɳ���
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_VOLUME 
        !��¼ˮ����ʴ��ɳ���,��������ʴ�߱��ڻ�����ɳ
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: FLOWEROSION_VOLUME
        !��¼�����ڻ���ɳ���
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: SHIELD_VOLUME        
        !��¼�ƻ�ԭ��1��̮����2������ˮ����ˢ
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_DEF
        !��¼�����ƻ�ʱ��
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_TIME 
        !��¼�ۼƵ�ˮ����ʴ������̮�����,ÿ�����; ��һ����flow;�ڶ���collapse;������collapse shield 
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: BANKEROSION_VOLUME 
        !�ж��Ƿ��������,��̮������
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: JUDGE_BANK
        !�ж��Ƿ����FLOW EROSION
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: FLOWEROSION_CHECK        
        !���ڼӺ��ش�������0�μ��ؽ��Ϊ��ʼӦ����
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NUMBER_LOAD  
        !��¼�ǵڼ�������,����������+1
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NUMBER_BANKCYCLE 
        !����ˮ����ˢ�϶�İ���,����һ������û��,���°���,0����,1��
        REAL*8, ALLOCATABLE, DIMENSION(:)    :: FLOW_EXCESS_JUDGE 
        !��¼����Ԫ���ڵ����
        INTEGER                              :: ELEMENT, NODE, &
        ! �ж��Ƿ���Ҫ�������Ӻ���
        &                                       JUDGE_TEMP
        !Ԥ��İ��ڳߴ�
        REAL*8                               :: GRID_BANK_H, GRID_BANK_W  
        
    END TYPE
    
    TYPE(EXCHANGE_TYPE) :: EXCHANGE1
    
END MODULE     
       
#endif

#ifdef PRE_BANK_COLLAPSE
MODULE BANK_COLLAPSE !�ļ���Ŵ�100��ʼ
    IMPLICIT NONE
    
    TYPE BANK_COLLAPSE_TYPE
        !��¼�����ڡ������鵥Ԫ����ģ��
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: MODULUS
        !��¼�����ڡ������鵥ԪӦ��ֵ��1~7�ֱ�Ϊ��������Ӧ����һ����Ӧ������С��Ӧ�����нǡ��ƻ�����
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TSTRESS
        !��¼�����ڡ������鵥Ԫ�ƻ����, 0�ȶ���1���ƻ���2���ƻ�
        INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: STATION   
        !��¼����������ڵ��������λ��ֵ
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: F, WY   
        !��¼�������Բ��������ɱȡ��ھ�������Ħ����
        REAL*8                                :: MU, COH, FRI
        !��¼����ĵ�Ԫ����
        REAL*8                                :: T_LE
        !��¼����ģ����ʼֵ
        REAL*8                                :: MODULUS_INITIAL  
        !������ɳת��Ϊ�״��ٷֱ�
        REAL*8                                :: COLL_TO_BED         
        !��¼̮������߽�����
        INTEGER, ALLOCATABLE, DIMENSION(:)  :: RESTRICTION  
        !��¼�ڵ���Χ��Ԫ��Ϣ
        INTEGER, ALLOCATABLE, DIMENSION(:,:):: MATCH        
    END TYPE
    
    TYPE(BANK_COLLAPSE_TYPE) :: BANKCOLLAPSE
        
END MODULE

#endif

MODULE INTERFACE_DEFINITIONS
    INTERFACE
        ! X-MOMENTUM SIMULATION MODULE
        SUBROUTINE CONTINUITY_MOMENTUM_X(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DXX,DYY,I1,I2,IP,UPS,DWS)
            IMPLICIT NONE
            INTEGER   , PARAMETER                :: NNX = SELECTED_REAL_KIND(8)
            REAL(NNX) , DIMENSION(:,:) , POINTER :: PP1 , PP2 , HH1 , HH2 , QQ1 , QQ2 , DXX , DYY
            REAL(NNX)                            :: T0 
            INTEGER                              :: I1 , I2 , IP
            LOGICAL                              :: UPS , DWS
        END SUBROUTINE CONTINUITY_MOMENTUM_X
        ! Y-MOMENTUM SIMUTION MODULE
        SUBROUTINE CONTINUITY_MOMENTUM_Y(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DXX,DYY,I1,I2,IP,UPS,DWS)
            IMPLICIT NONE
            INTEGER , PARAMETER                  :: NNY = SELECTED_REAL_KIND(8)
            REAL(NNY) , DIMENSION(:,:) , POINTER :: PP1 , PP2 , HH1 , HH2 , QQ1 , QQ2 , DXX , DYY
            REAL(NNY)                            :: T0 
            INTEGER                              :: I1 , I2 , IP
            LOGICAL                              :: UPS , DWS
        END SUBROUTINE CONTINUITY_MOMENTUM_Y
        !BOUNDARY CONDITION
        SUBROUTINE BND_WATERLEVEL(T0,HH2)
            IMPLICIT NONE
            INTEGER   , PARAMETER                :: NNB = SELECTED_REAL_KIND(8)
            REAL(NNB) , DIMENSION(:,:) , POINTER :: HH2
            REAL(NNB)                            :: T0
        END SUBROUTINE BND_WATERLEVEL
        !DRYING AND FLOODING
        SUBROUTINE DRYING_AND_FLOODING_X(PP2,HH2,QQ2)
            IMPLICIT NONE
            INTEGER   , PARAMETER                :: NNB = SELECTED_REAL_KIND(8)
            REAL(NNB) , DIMENSION(:,:) , POINTER :: PP2, HH2, QQ2
        END SUBROUTINE DRYING_AND_FLOODING_X
        SUBROUTINE DRYING_AND_FLOODING_Y(PP2,HH2,QQ2)
            IMPLICIT NONE
            INTEGER   , PARAMETER                :: NNB = SELECTED_REAL_KIND(8)
            REAL(NNB) , DIMENSION(:,:) , POINTER :: PP2, HH2, QQ2
        END SUBROUTINE DRYING_AND_FLOODING_Y
    END INTERFACE
END MODULE INTERFACE_DEFINITIONS