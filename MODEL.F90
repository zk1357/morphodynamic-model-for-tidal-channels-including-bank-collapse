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
    !计算域网格起始点编号    
    INTEGER,PARAMETER                           :: M_START=1 , M_END=41 , N_START=1 , N_END=276       
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: P1,P2
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: Q1,Q2 
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: H1,H2 
    !X,Y方向网格步长    
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: DX1,DY1 
    !时间步长
    REAL(NP1),PARAMETER                         :: DT=6.D0
    !坡度边界,潮位振幅,边界水深,干湿判断临界水深
    REAL(NP1),PARAMETER                         :: SLOPE=2.D-3 , SLOPE1=0.D-3 , AMP=1.D0 , DEPOUTSEA=3.D0 , FLOODINGDEPTHX=0.1D0 , FLOODINGDEPTHY=0.1D0
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: DEP
    !底摩阻及涡粘系数,默认使用常数谢才公式
    REAL(NP1),PARAMETER                         :: MANNING=38 , CHE=38, CHE1=10, EDDY_X=2.0D0 , EDDY_Y=1.D0 !根据计算，该值小于1.9可行
    !水深,流速节点状态
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: NDD
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: VDDX   !VELOCITY-NDD
    INTEGER  ,ALLOCATABLE,TARGET,DIMENSION(:,:) :: VDDY
    REAL  ,   ALLOCATABLE,TARGET,DIMENSION(:,:) :: IDN
    REAL  ,   ALLOCATABLE,TARGET,DIMENSION(:,:) :: JDN   
    !Q~A关系
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:) :: QA_Q , QA_A
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: QA_Q_P , QA_A_P !求各个点的Q~A关系
    !潮沟位置,便于赋初值    
    INTEGER                                     :: TCI_START , TCI_END  , TCJ_START , TCJ_END
    !两侧岸壁BVP总个数
    INTEGER   ,PARAMETER                        :: NUMBER_BV=100                      
    !潮沟深度    
    REAL(NP1),ALLOCATABLE,TARGET,DIMENSION(:)   :: BANK_H
    !潮沟宽度
    REAL(NP1),PARAMETER                         :: WIDTH_CREEK_INI=20.D0
    !矩阵求解    
    REAL(NP1)                                   :: AMT(6,3001)        
    REAL(NP1)                                   :: AMT3(4,3001)
END MODULE HYDRO 
#endif

#ifdef PRE_AD
MODULE Morphodynamics
    IMPLICIT NONE
    INTEGER,PARAMETER                                 :: NP2=SELECTED_REAL_KIND(8)
    !含沙量    
    REAL(NP2),ALLOCATABLE,TARGET,DIMENSION(:,:)       :: C1 , C2 
    !加速因子
    INTEGER                                           :: MOR_FACTOR
    !X,Y方向的泥沙扩散系数
    REAL(NP2),PARAMETER                               :: KX=10.D0  , KY=2.D0    
    !泥沙参数                                           ME取值为DELFT3d默认值、WS为默认值一半、泥沙密度为默认值
    REAL(NP2),PARAMETER                               :: ME=5.D-5 , TE=0.5D0 , TD=1000.D0 , WS=5.D-4 , &
                                                       &  CD=0.002 , FLOWDENSITY=1000 , SANDDENSITY=2650 , &
                                                       &  C_BOUNDARY=0.D0 , C_LIMIT=1000.5D0 , PORO=0.43
END MODULE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
MODULE BANK_EROSION
    IMPLICIT NONE
    INTEGER   ,PARAMETER                     :: NP6=SELECTED_REAL_KIND(8)   
    INTEGER   ,PARAMETER                     :: COLLAPSE_START=1 , COLLAPSE_END=39 !一侧的坍塌起始点
    INTEGER   ,PARAMETER                     :: NUMBER_OB=21                       !岸壁监测点总个数,监测点序号从下往上,岸壁20等分
    INTEGER                                  :: MOR_FACTOR_BANK
    REAL(NP6) ,PARAMETER                     :: WIDTH_BANK=10.D0
    REAL(NP6) ,PARAMETER                     :: COL_PER_BED=0.1D0                   !COLLAPSE BECOME BED IN PERCENTAGE
    REAL(NP6)                                :: BANK_START_VELOCITY                     !启动流速
    REAL(NP6) ,PARAMETER                     :: ME_BANK=2.D-4 , CD_BANK=0.002 , TE_BANK=0.062D0 , &
                                              & FLOWDENSITY_BANK=1000 , SANDDENSITY_BANK=1940 , PORO_BANK=0.
    REAL(NP6) ,PARAMETER                     :: BANK_EROSION_LIMIT=0.1D0 , BANK_COLLAPSE_LIMIT=0.3D0
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: POSITION 
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: QE_OUT !记录潮沟内各断面的侵蚀通量
    REAL(NP6) , ALLOCATABLE , DIMENSION(:,:) :: BVP_VELOCITY !记录BVP流速过程,用来计算水流切应力
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) :: BVP , BVP_IN        !记录潮沟内，潮沟壁旁的流速点位置,分别为I,J
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: BVP_JUDGE  !判断是否更新岸壁边壁点，即岸壁网格后退一格 
END MODULE BANK_EROSION

MODULE EXCHANGE
    IMPLICIT NONE
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: JUDGE_BANK,& !判断是否更换为初始岸壁网格
    &                                           JUDGE_TIME   !涨落潮时期判断 0涨潮 -1落潮
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) :: BF           !记录破坏次数及百分比,1-erosion次数 2-CO次数 3-cyc识别 4-e/c判断
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: NUMBER_LOAD  !对于每一个岸壁来说,为加荷载次数，直至破坏换新岸壁
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: NUMBER_BANKCYCLE !记录是第几个岸壁,即更换次数+1
    !每块岸壁信息
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_RETREATS  !每块岸壁后退距离
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_COLLAPSE  !每块岸壁坍塌距离
    REAL*8    , ALLOCATABLE , DIMENSION(:,:) :: LENGTH_FLOW      !每块岸壁侵蚀距离,若不计算坍塌，此值为水流侵蚀总距离，因为岸壁更换判断在Mesh子程序
    !每天岸壁信息
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LENGTH_EROSION_DAY  !每天BANK EROSION距离
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LENGTH_COLLAPSE_DAY  !每天COLLAPSE    距离   
    !用于初始岸壁网格的缩放,应对边壁高度减少
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: ZOOMIN           
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: LHL , LHLB , &   !LHL记录各断面最低水位高度
    &                                           STATION_BANK , & !用于水流冲刷较多的岸壁,冲完2m没塌,换新岸壁,0不换,1换
    &                                           LENGTH_COLLAPSE_BEFORE , & !记录前一次网格变化时的坍塌距离，从而计算前一个网格的高程增加值
    &                                           RETREATS_TOTAL         , & !记录每个CS断面的总后退距离
    &                                           UPDATE_C                   !累计每次网格更新前的坍塌距离合
    INTEGER                                  :: JUDGE_TEMP       !临时判断是否换岸壁
END MODULE EXCHANGE

MODULE SED_EXCHANGE
    IMPLICIT NONE
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: S_BANK_EROSION !用于记录岸壁冲刷的泥沙,每一步时间步长带入,为含沙量
    REAL*8    , ALLOCATABLE , DIMENSION(:)   :: S_BANK_COLLAPSE!用于记录坍塌转化为坡脚掩护的泥沙,为泥沙体积
ENDMODULE SED_EXCHANGE
    
MODULE DEF_EXCHANGE
    IMPLICIT NONE
    TYPE EXCHANGE_TYPE
        !记录模型中岸壁尺寸信息
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: EXCHANGE_BANK_SIZE
        !记录模型岸壁相对于预设网格尺寸的缩放比例
        REAL*8, ALLOCATABLE, DIMENSION(:)    :: ZOOMIN 
        !网格移除判断
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DIG_NUMBER         
        !记录网格节点坐标、单元对应节点号
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: GRID, MESH 
        !记录坍塌泥沙体积
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_VOLUME 
        !记录水流侵蚀泥沙体积,不包含侵蚀边壁掩护的泥沙
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: FLOWEROSION_VOLUME
        !记录岸壁掩护泥沙体积
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: SHIELD_VOLUME        
        !记录破坏原因：1：坍塌；2：过量水流冲刷
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_DEF
        !记录岸壁破坏时间
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: COLLAPSE_TIME 
        !记录累计的水流侵蚀、岸壁坍塌体积,每天输出; 第一列是flow;第二列collapse;第三列collapse shield 
        REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: BANKEROSION_VOLUME 
        !判断是否更换岸壁,由坍塌引起
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: JUDGE_BANK
        !判断是否继续FLOW EROSION
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: FLOWEROSION_CHECK        
        !岸壁加荷载次数，第0次加载结果为初始应力场
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NUMBER_LOAD  
        !记录是第几个岸壁,即更换次数+1
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NUMBER_BANKCYCLE 
        !用于水流冲刷较多的岸壁,冲完一定距离没塌,换新岸壁,0不换,1换
        REAL*8, ALLOCATABLE, DIMENSION(:)    :: FLOW_EXCESS_JUDGE 
        !记录网格单元、节点个数
        INTEGER                              :: ELEMENT, NODE, &
        ! 判断是否需要继续叠加荷载
        &                                       JUDGE_TEMP
        !预设的岸壁尺寸
        REAL*8                               :: GRID_BANK_H, GRID_BANK_W  
        
    END TYPE
    
    TYPE(EXCHANGE_TYPE) :: EXCHANGE1
    
END MODULE     
       
#endif

#ifdef PRE_BANK_COLLAPSE
MODULE BANK_COLLAPSE !文件序号从100开始
    IMPLICIT NONE
    
    TYPE BANK_COLLAPSE_TYPE
        !记录各岸壁、各土块单元弹性模量
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: MODULUS
        !记录各岸壁、各土块单元应力值：1~7分别为：两个正应力、一个剪应力、大、小主应力、夹角、破坏类型
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TSTRESS
        !记录各岸壁、各土块单元破坏与否, 0稳定、1拉破话、2剪破坏
        INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: STATION   
        !记录各岸壁网格节点的外力、位移值
        REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: F, WY   
        !记录土体特性参数：泊松比、内聚力、内摩擦角
        REAL*8                                :: MU, COH, FRI
        !记录假设的单元体厚度
        REAL*8                                :: T_LE
        !记录弹性模量初始值
        REAL*8                                :: MODULUS_INITIAL  
        !塌落泥沙转化为底床百分比
        REAL*8                                :: COLL_TO_BED         
        !记录坍塌土体边界条件
        INTEGER, ALLOCATABLE, DIMENSION(:)  :: RESTRICTION  
        !记录节点周围单元信息
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