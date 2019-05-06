#include 'define.h'
    
#ifdef PRE_BANK_COLLAPSE 
MODULE LINEAR_ELASTIC !文件序号从100开始
    IMPLICIT NONE
    INTEGER , PARAMETER ::NLE=SELECTED_REAL_KIND(8)
    REAL(NLE) , PARAMETER                    ::MODULUS_INITIAL=5.D6 !初始弹性模量
    REAL(NLE) , PARAMETER                    ::T_LE=1.D0 !厚度
    REAL(NLE) , ALLOCATABLE , DIMENSION(:,:) ::MEA !记录网格单元的信息:坐标、节点、面积等等
    REAL(NLE) , ALLOCATABLE , DIMENSION(:)   ::MODULUS
    REAL(NLE) , ALLOCATABLE , DIMENSION(:,:) ::TSTRESS
    REAL(NLE) , ALLOCATABLE , DIMENSION(:)   ::F
    REAL(NLE) , ALLOCATABLE , DIMENSION(:)   ::WY
    INTEGER   , ALLOCATABLE , DIMENSION(:)   ::RESTRICTION    
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) ::MATCH
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: STATION  !记录单元状态
    INTEGER                                  ::NODE , ELEMENT
    REAL(NLE)                                ::MU , COH , FRI
   
    CONTAINS
    
    SUBROUTINE LINEAR_ELASTIC_STRESS
        !求应力的主程序，包含各个子程序
        USE BANK_EROSION
        USE EXCHANGE
        USE DEF_EXCHANGE
        USE HYDRO
        INTEGER  :: I , J , K , NUMBER_J
        CHARACTER*40 :: FILENAME , FILENAME_PATH
             
        ALLOCATE( MODULUS(ELEMENT) )
        ALLOCATE( STATION(ELEMENT) ) !0稳定、1拉破话、2剪破坏
        ALLOCATE( F(2*NODE)        )

LOOP0:  DO I = 1 , NUMBER_BV/2
            
            !给固定断面计算坍塌
            !IF ( I < COLLAPSE_START                                    .OR. &
            !&    I > COLLAPSE_END .AND. I < COLLAPSE_START+NUMBER_BV/2 .OR. &
            !&    I > COLLAPSE_END+NUMBER_BV/2   ) CYCLE
            !漫滩不计算坍塌
            !IF ( LHL(I)>BANK_H(I) .AND. LHLB(I)>BANK_H(I) ) THEN
            !    IF ( JUDGE_BANK(I)==1 ) THEN
            !        !BF(I,3)=1
            !    END IF
            !    CYCLE LOOP0
            !END IF
            FILENAME(1:7)='MODULUS'
            WRITE(FILENAME(8:11),'(I4)')I
            FILENAME(12:15)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2)
            DO J = 1 , ELEMENT
                READ( 101 , REC = J ) MODULUS(J)
            END DO
            CLOSE(101)
            
            FILENAME(1:7)='STATION'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2)
            DO J = 1 , ELEMENT
                READ( 101 , REC = J ) STATION(J)
            END DO
            CLOSE(101)

            CALL LINEAR_ELASTIC_EXCAVATION( I , NUMBER_J )
            
            !WRITE(*,*)I,NUMBER_J
            
 LOOP101:  DO J = 1 , NUMBER_J
                IF (JUDGE_TEMP==1) THEN
                    JUDGE_TEMP=0
                    EXIT LOOP101
                END IF
                !WRITE(*,*)'I=',I,'J=',J,NUMBER_J
                CALL LINEAR_ELASTIC_STIFFNESS( I , J )
            END DO LOOP101
            JUDGE_TEMP=0
            
            FILENAME(1:7)='MODULUS'
            WRITE(FILENAME(8:11),'(I4)')I
            FILENAME(12:15)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2)
            DO J = 1 , ELEMENT
                WRITE( 101 , REC = J ) MODULUS(J)
            END DO            
            CLOSE(101)

            FILENAME(1:7)='STATION'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2)
            DO J = 1 , ELEMENT
                WRITE( 101 , REC = J ) STATION(J)
            END DO
            CLOSE(101)

!            DEALLOCATE( MODULUS )
!            DEALLOCATE( F )

!            PAUSE
        END DO LOOP0

        DEALLOCATE( MODULUS )
        DEALLOCATE( STATION )
        DEALLOCATE( F )
!        DEALLOCATE( TSTRESS )

!        PAUSE
    END SUBROUTINE LINEAR_ELASTIC_STRESS

    SUBROUTINE LINEAR_ELASTIC_GRID
        !读取网格以及网格剖分，外来可以加入网格自动剖分
        USE HYDRO
        USE EXCHANGE
        USE BANK_EROSION
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:) ::GRID , MESH !记录网格节点信息
        INTEGER :: I , J , K , MEA_NUMBER
        REAL(NLE) :: X1 , X2 , X3 , Y1 , Y2 , Y3
        CHARACTER*40 :: FILENAME , FILENAME_PATH
        OPEN (UNIT=100 , FILE='初始网格节点个数、单元个数.TXT' ,STATUS='OLD')
        OPEN (UNIT=103 , FILE='初始网格节点坐标.TXT'         ,STATUS='OLD')
        OPEN (UNIT=101 , FILE='初始网格单元节点号.TXT'       ,STATUS='OLD')
        OPEN (UNIT=102 , FILE='土体参数.TXT'               ,STATUS='OLD')

        READ (100,*) NODE , ELEMENT
        READ (102,*) MU , COH , FRI
        ALLOCATE(MEA(ELEMENT,11),GRID(NODE,2),MESH(ELEMENT,3),MATCH(NODE,10),RESTRICTION(2*NODE),F(2*NODE))
        READ (103,*) ((GRID(I,J),J=1,2),I=1,NODE)
        READ (101,*) ((MESH(I,J),J=1,3),I=1,ELEMENT)
        MEA=0.D0
!-----------------------------------------------------------------------------记录各单元信息
        DO I = 1 , ELEMENT
            X1=GRID(MESH(I,1),1) ; X2=GRID(MESH(I,2),1) ; X3=GRID(MESH(I,3),1) 
            Y1=GRID(MESH(I,1),2) ; Y2=GRID(MESH(I,2),2) ; Y3=GRID(MESH(I,3),2)
            MEA(I,10)=0.5D0*ABS(Y3*(X2-X1)+Y1*(X3-X2)+Y2*(X1-X3))
            MEA(I,1)=MESH(I,1) ; MEA(I,2)=MESH(I,2) ; MEA(I,3)=MESH(I,3)
            MEA(I,4)=GRID(MESH(I,1),1) ; MEA(I,5)=GRID(MESH(I,1),2)
            MEA(I,6)=GRID(MESH(I,2),1) ; MEA(I,7)=GRID(MESH(I,2),2)
            MEA(I,8)=GRID(MESH(I,3),1) ; MEA(I,9)=GRID(MESH(I,3),2)
        END DO
!        FILENAME(1:3)='MEA'
!        WRITE(FILENAME(4:5),'(I2)') MEA_NUMBER
!        FILENAME(6:9)='.DAT'
!        OPEN (UNIT=101 , FILE=FILENAME ,STATUS='OLD')
!
!        READ(101,'(F12.4)')((MEA(I,J),J=1,11),I=1,ELEMENT)
!        READ(102,)
!-----------------------------------------------------------------------------记录节点周围单元
        MATCH=0
        K=0
        DO I = 1 , NODE
            K=1
            DO J = 1 , ELEMENT
                IF (I==MEA(J,1) .OR. I==MEA(J,2) .OR. I==MEA(J,3)) THEN
                    MATCH(I,K)=J    
                    K=K+1
                END IF
            END DO
        END DO
!-----------------------------------------------------------------------------边界约束及外力
        RESTRICTION=1
        DO I = 1 , NODE
            IF (GRID(I,1)==0 ) THEN
                RESTRICTION(2*I-1)=0
!                RESTRICTION(2*I  )=0
            END IF
            IF (GRID(I,2)==0) THEN
                RESTRICTION(2*I-1)=0
                RESTRICTION(2*I  )=0
            END IF
        END DO

        F=0.D0
        DO I = 1 , NODE
            DO J = 1 , 10
                IF (MATCH(I,J)==0) CYCLE
                F(2*I)=F(2*I)-MEA(MATCH(I,J),10)*3.13D0*T_LE !单位为KN
            END DO
        END DO
        
        X1=0.D0
        OPEN (104,FILE='自重节点力.DAT')
        DO I = 1 , NODE
            X1=X1+F(2*I)
            WRITE(104,'(2F12.4)') F(2*I-1) , F(2*I)
        END DO
        WRITE(104,*)X1
        CLOSE(104)
        
        ALLOCATE( MODULUS(ELEMENT) )
        MODULUS=MODULUS_INITIAL
        DO I = 1 , NUMBER_BV   !给MODULUS赋初值
            FILENAME(1:7)='MODULUS'
            WRITE(FILENAME(8:11),'(I4)')I
            FILENAME(12:15)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=104 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2)
            DO J=1,ELEMENT
                WRITE( 104 , REC = J) MODULUS(J)
            END DO
            CLOSE(104)   
        END DO        

        ALLOCATE( STATION(ELEMENT) )
        STATION=0
        DO I = 1 , NUMBER_BV   !给STATION赋初值 0稳定、1拉破话、2剪破坏
            FILENAME(1:7)='STATION'
            WRITE(FILENAME(8:11),'(I4)')I
            FILENAME(12:15)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=104 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='REPLACE' , ACCESS='DIRECT' , RECL=2)
            DO J=1,ELEMENT
                WRITE( 104 , REC = J) STATION(J)
            END DO
            CLOSE(104)   
        END DO  
        DEALLOCATE( STATION )

!        PAUSE


!        F(2)=-10.D0   !单位为KN
        DEALLOCATE(GRID,MESH)
        CLOSE(100)
        CLOSE(101)
        CLOSE(103)
        CLOSE(102)

        !OPEN(100 , FILE= 'GRID.RES' , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        !DO I = 1 , NODE
        !    WRITE (100 , REC=I) GRID(I,1)
        !END DO
        !CLOSE(100)
        !
        !OPEN(600 , FILE= 'MESH.RES' , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*3)
        !DO I = 1 , ELEMENT
        !    WRITE (600 , REC=I) MESH(I,1), MESH(I,2), MESH(I,3)
        !END DO
        !CLOSE(600)       
             
        RETURN
    END SUBROUTINE LINEAR_ELASTIC_GRID


    SUBROUTINE LINEAR_ELASTIC_EXCAVATION(NUMBER_I,NUMBER_J)
        USE EXCHANGE
        USE HYDRO
        
        IMPLICIT NONE
        INTEGER , PARAMETER                        :: LEE=SELECTED_REAL_KIND(8)
        REAL(LEE) , ALLOCATABLE , DIMENSION(:,:)   :: B , BT , GRID
        INTEGER   , ALLOCATABLE , DIMENSION(:)     :: DIG_NODE !开挖面上的节点
        INTEGER                                    :: I , J , K , DIG_NUMBER , NUMBER_I , NUMBER_J
        LOGICAL                                    :: TEMP1 , TEMP2
        REAL(LEE)                                  :: TEMP
        CHARACTER*40                               :: FILENAME , FILENAME1 , FILENAME_PATH

        ALLOCATE( B (3,6) )
        ALLOCATE( BT(6,3) )
!        ALLOCATE( F(2*NODE) )
!        ALLOCATE( MEA(ELEMENT,11) )
        ALLOCATE( TSTRESS(ELEMENT,3))

        FILENAME(1:3)='MEA'
        WRITE(FILENAME(4:7),'(I4)') NUMBER_I
        FILENAME(8:11)='.RES'
        FILENAME_PATH='TEMP_FILE\'//FILENAME
        OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2*11)
!        OPEN (UNIT=102 , FILE='MEA.DAT' )

        DO I=1,ELEMENT
            READ( 101 , REC = I) (MEA(I,J),J=1,11)
        END DO
        CLOSE(101)

!        WRITE(*,*)1
!        PAUSE

        IF ( JUDGE_BANK(NUMBER_I)==1 ) THEN !变量进行初始化
!            JUDGE_BANK(NUMBER_I)=0
            STATION=0
            MODULUS=MODULUS_INITIAL
            OPEN (UNIT=102 , FILE='初始应力分布.RES' , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2*3)
            DO I = 1 , ELEMENT
                READ( 102 , REC=I )(TSTRESS(I,J),J=1,3)
            END DO
            CLOSE(102)
        ELSE
            FILENAME1(1:8)='总应力场'
            WRITE(FILENAME1(9:12),'(I4)') NUMBER_I
            FILENAME1(13:16)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME1
            OPEN (UNIT=102 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2*3)
            DO I = 1 , ELEMENT
                READ(102 , REC=I )(TSTRESS(I,J),J=1,3)
            END DO
            CLOSE(102)
        END IF

!        PAUSE
        F=0.D0
        DO I = 1 , ELEMENT       !加开挖荷载    
            IF (MEA(I,11) .NE. 1.D0) CYCLE
!            WRITE(*,*)I,NODE,MEA(I,1),MEA(I,2),MEA(I,3)
            B =0.D0
            BT=0.D0
            B(1,1)=MEA(I,7)-MEA(I,9)     ; B(1,3)=MEA(I,9)-MEA(I,5)     ; B(1,5)=MEA(I,5)-MEA(I,7)
            B(2,2)=MEA(I,8)-MEA(I,6)     ; B(2,4)=MEA(I,4)-MEA(I,8)     ; B(2,6)=MEA(I,6)-MEA(I,4)
            B(3,1)=MEA(I,8)-MEA(I,6)     ; B(3,2)=MEA(I,7)-MEA(I,9)     ; B(3,3)=MEA(I,4)-MEA(I,8)
            B(3,4)=MEA(I,9)-MEA(I,5)     ; B(3,5)=MEA(I,6)-MEA(I,4)     ; B(3,6)=MEA(I,5)-MEA(I,7)
            BT=TRANSPOSE(B)

            F(2*INT(MEA(I,1))-1)=0.5D0*T_LE*( BT(1,1)*TSTRESS(I,1) + BT(1,2)*TSTRESS(I,2) + BT(1,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,1))-1)
            F(2*INT(MEA(I,1))  )=0.5D0*T_LE*( BT(2,1)*TSTRESS(I,1) + BT(2,2)*TSTRESS(I,2) + BT(2,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,1))  )
            
            F(2*INT(MEA(I,2))-1)=0.5D0*T_LE*( BT(3,1)*TSTRESS(I,1) + BT(3,2)*TSTRESS(I,2) + BT(3,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,2))-1)
            F(2*INT(MEA(I,2))  )=0.5D0*T_LE*( BT(4,1)*TSTRESS(I,1) + BT(4,2)*TSTRESS(I,2) + BT(4,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,2))  )       

            F(2*INT(MEA(I,3))-1)=0.5D0*T_LE*( BT(5,1)*TSTRESS(I,1) + BT(5,2)*TSTRESS(I,2) + BT(5,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,3))-1)
            F(2*INT(MEA(I,3))  )=0.5D0*T_LE*( BT(6,1)*TSTRESS(I,1) + BT(6,2)*TSTRESS(I,2) + BT(6,3)*TSTRESS(I,3) ) + F(2*INT(MEA(I,3))  )

            !是否需加土荷载？？？
            F(2*INT(MEA(I,1))  )=F(2*INT(MEA(I,1))  )-MEA(I,10)*5227*1.D-3*T_LE  !单位为KN
            F(2*INT(MEA(I,2))  )=F(2*INT(MEA(I,2))  )-MEA(I,10)*5227*1.D-3*T_LE  !单位为KN
            F(2*INT(MEA(I,3))  )=F(2*INT(MEA(I,3))  )-MEA(I,10)*5227*1.D-3*T_LE  !单位为KN


        END DO
       
!        TEMP=0.d0
!        OPEN (UNIT=101 , FILE='开挖后节点力.DAT')
!        DO I = 1 , NODE
!            TEMP=TEMP+DMIN1(F(2*I),0.D0)
!            IF (F(2*I-1)==0.D0 .OR. F(2*I)==0.D0) CYCLE
!            WRITE(101,'(2F12.4)') F(2*I-1) , F(2*I)
!        END DO
!        CLOSE(101)
!        WRITE(*,*)1 , TEMP
!        PAUSE

!---------------------------
        ALLOCATE( GRID(NODE,2) )
        OPEN (UNIT=101 , FILE='初始网格节点坐标.TXT'         ,STATUS='OLD')
        READ(101,*)((GRID(I,J),J=1,2),I=1,NODE)
        CLOSE(101)
!        ALLOCATE( DIG_NODE(NODE) )
!        FILENAME(1:3)='DIG'
!        WRITE(FILENAME(4:5),'(I2)') NUMBER_I
!        FILENAME(6:9)='.DAT'
!        WRITE(*,*)FILENAME
!        OPEN (UNIT=40 , FILE=FILENAME ,STATUS='OLD')
!        READ(40,'(I6)') (DIG_NODE(K),K=1,NODE)
!        CLOSE(40)
        
!        OPEN(1,FILE='TEST.DAT')
!!        ALLOCATE( NODE_DIG(NODE) )
!        DO I=1,NODE
!            IF (DIG_NODE(I)==0)CYCLE
!            WRITE(1,'(4F12.2)')GRID(I,1) , GRID(I,2) , F(2*I-1) , F(2*I)
!        END DO
!        CLOSE(1)
!        DEALLOCATE( DIG_NODE )
!        PAUSE
!---------------------------



!        PAUSE
        ALLOCATE( DIG_NODE(100) )
        DIG_NODE  =0
        DIG_NUMBER=0
LOOP1:  DO I = 1 , NODE
            TEMP1=.FALSE. !未开挖
            TEMP2=.FALSE. !刚开挖
LOOP2:      DO J = 1 , 10
                IF (MATCH(I,J)==0) CYCLE
                IF (MEA(MATCH(I,J),11)==0.D0) TEMP1=.TRUE.
                IF (MEA(MATCH(I,J),11)==1.D0) TEMP2=.TRUE.
                IF (TEMP1 .AND. TEMP2) THEN
                   DIG_NUMBER=DIG_NUMBER+1 
                   DIG_NODE(DIG_NUMBER)=I
                   EXIT LOOP2
                END IF
            END DO LOOP2
        END DO LOOP1
!        WRITE(*,*)DIG_NUMBER
!        DO I = 1 , DIG_NUMBER
!            WRITE(*,*)F(2*DIG_NODE(I)-1),F(2*DIG_NODE(I))
!        END DO
!        PAUSE
        DO I = 1 , NODE  !只施加开挖面上的荷载
            K=0
            DO J = 1 , DIG_NUMBER
                IF (DIG_NODE(J)==I) K=1    
            END DO
            IF (K==0) F(2*I-1)=0.D0
            IF (K==0) F(2*I  )=0.D0
        END DO  
          
!        DO I = 1 , DIG_NUMBER
!            WRITE(*,*)DIG_NODE(I) , F(2*DIG_NODE(I)-1) , F(2*DIG_NODE(I))
!        END DO
!        PAUSE

!        ALLOCATE( GRID(NODE,2) )
!        OPEN (UNIT=101 , FILE='初始网格节点坐标.TXT'         ,STATUS='OLD')
!        READ(101,*)((GRID(I,J),J=1,2),I=1,NODE)
!        CLOSE(101)
        
        !!饱和土不考虑静水压力的变化
        !IF ( LHL(NUMBER_I)>LHLB(NUMBER_I) .AND. LHL(NUMBER_I)<=BANK_H(NUMBER_I) ) THEN   !涨潮过程
        !!IF ( LHL(NUMBER_I)>LHLB(NUMBER_I)  ) THEN   !涨潮过程    
        !    DO I = 1 , DIG_NUMBER  !考虑水体静压力
        !        IF ( GRID(DIG_NODE(I),2) < LHL(NUMBER_I) ) THEN
        !            IF ( GRID(DIG_NODE(I),2) > LHLB(NUMBER_I) ) THEN
        !                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)-T_LE*490.D0*(LHL(NUMBER_I)-GRID(DIG_NODE(I),2))*1.D-3 !化为KN
        !            ELSE
        !                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)-T_LE*490.D0*(LHL(NUMBER_I)-LHLB(NUMBER_I)     )*1.D-3
        !            END IF
        !        END IF
        !    END DO
        !ELSE IF (LHL(NUMBER_I)<LHLB(NUMBER_I) .AND. LHLB(NUMBER_I)<=BANK_H(NUMBER_I) ) THEN !落潮过程
        !!ELSE IF (LHL(NUMBER_I)<LHLB(NUMBER_I)  ) THEN !落潮过程
        !     DO I = 1 , DIG_NUMBER  !考虑水体静压力
        !        IF ( GRID(DIG_NODE(I),2) < LHLB(NUMBER_I) ) THEN
        !            IF ( GRID(DIG_NODE(I),2) > LHL(NUMBER_I) ) THEN
        !                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)+T_LE*490.D0*(LHLB(NUMBER_I)-GRID(DIG_NODE(I),2))*1.D-3
        !            ELSE
        !                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)+T_LE*490.D0*(LHLB(NUMBER_I)-LHL (NUMBER_I)     )*1.D-3
        !            END IF
        !        END IF 
        !     END DO          
        !END IF

!        IF ( LHLB(NUMBER_I)<=2.D0 .AND. LHL(NUMBER_I)>=2.D0 ) THEN !涨潮漫滩，加浮力
!!            DO I = 1 , NODE
!!                DO J = 1 , 10
!!                    IF (MATCH(I,J)==0) CYCLE
!!                    F(2*I)=F(2*I)+MEA(MATCH(I,J),10)*3345*1.D-3*T_LE  !单位为KN,海水密度1024千克每立方米
!!                END DO
!!            END DO
!            DO I = 1 , DIG_NUMBER  !漫滩去水体静压力
!                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)+T_LE*490.D0*(LHLB(NUMBER_I)-GRID(DIG_NODE(I),2))*1.D-3
!            END DO 
!        ELSE IF ( LHLB(NUMBER_I)>=2.D0 .AND. LHL(NUMBER_I)<=2.D0  ) THEN !落潮过程
!!            DO I = 1 , NODE    
!!                DO J = 1 , 10
!!                    IF (MATCH(I,J)==0) CYCLE
!!                    F(2*I)=F(2*I)-MEA(MATCH(I,J),10)*3345*1.D-3*T_LE  !单位为KN,海水密度1024千克每立方米
!!                END DO
!!            END DO
!            DO I = 1 , DIG_NUMBER  !退滩加水体静压力
!                F(2*DIG_NODE(I)-1)=F(2*DIG_NODE(I)-1)-T_LE*490.D0*(LHL(NUMBER_I)-GRID(DIG_NODE(I),2))*1.D-3
!            END DO 
!        END IF                
        
!        DO I = 1 , NODE  !加浮力
!            IF ( GRID(I,2)>LHLB(NUMBER_I) .AND. GRID(I,2)<LHL(NUMBER_I) ) THEN !涨潮过程
!                DO J = 1 , 10
!                    IF (MATCH(I,J)==0) CYCLE
!                    F(2*I)=F(2*I)+MEA(MATCH(I,J),10)*3345*1.D-3*T_LE  !单位为KN,海水密度1024千克每立方米
!                END DO
!            ELSE IF ( GRID(I,2)<LHLB(NUMBER_I) .AND. GRID(I,2)>LHL(NUMBER_I) ) THEN !落潮过程
!                DO J = 1 , 10
!                    IF (MATCH(I,J)==0) CYCLE
!                    F(2*I)=F(2*I)-MEA(MATCH(I,J),10)*3345*1.D-3*T_LE  !单位为KN,海水密度1024千克每立方米
!                END DO
!            END IF                
!        END DO

!        WRITE(*,*)LHL(NUMBER_I)
!        PAUSE

!        PAUSE
        TEMP=0.D0
        DO I = 1 , DIG_NUMBER
            IF ( DABS(F(2*DIG_NODE(I)-1)) > DABS(TEMP) ) TEMP=F(2*DIG_NODE(I)-1)
            IF ( DABS(F(2*DIG_NODE(I)  )) > DABS(TEMP) ) TEMP=F(2*DIG_NODE(I)  )
        END DO
!        WRITE(*,*)NUMBER_I,TEMP
!        PAUSE
        NUMBER_J=CEILING(ABS(TEMP)/2)   !荷载分多少份加载
!        WRITE(*,*)NUMBER_J
!        PAUSE
        IF (NUMBER_J==0) NUMBER_J=1 !荷载小于10KN,不划分
!        WRITE(*,*)NUMBER_J
!        PAUSE
!        NUMBER_J=1
        F=F/NUMBER_J
!        WRITE(*,*)NUMBER_I
!        DO I=1,DIG_NUMBER
!            WRITE(*,*)F(DIG_NODE(I))
!        END DO
!        PAUSE
!        DO I = 1,2*NODE
!        IF (F(I)==0.D0)CYCLE
!        WRITE(*,*)F(I)
!        END DO
!        WRITE(*,*)NUMBER_J 
!        PAUSE
        FILENAME(1:3)='MEA'        !修改MEA(I,11),1=>2
        WRITE(FILENAME(4:7),'(I4)') NUMBER_I
        FILENAME(8:11)='.RES'
        FILENAME_PATH='TEMP_FILE\'//FILENAME
!        OPEN (UNIT=101 , FILE=FILENAME ,STATUS='REPLACE')
        OPEN (UNIT=101 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*11)
        DO I = 1 , ELEMENT
            IF (MEA(I,11)==1.D0) MEA(I,11)=2.D0
        END DO
        DO I = 1 , ELEMENT
            WRITE(101 , REC=I ) (MEA(I,J),J=1,11)
        END DO
        CLOSE(101)

!        WRITE(*,*)1
!        PAUSE

        DEALLOCATE( B )
        DEALLOCATE( BT )      
        DEALLOCATE( GRID )
        DEALLOCATE( DIG_NODE )
        DEALLOCATE( TSTRESS  )
        RETURN
    END SUBROUTINE LINEAR_ELASTIC_EXCAVATION

    SUBROUTINE LINEAR_ELASTIC_STIFFNESS(NUMBER_I,NUMBER_J)
        !求单元、整体劲度矩阵
        IMPLICIT NONE
        INTERFACE
            SUBROUTINE pardiso_sym_f90( N , NNZ , A , JA , IA , B , X )
                IMPLICIT NONE
                INTEGER , PARAMETER                        :: PSF=SELECTED_REAL_KIND(8)
                INTEGER                                    :: N , NNZ
                REAL(PSF) , ALLOCATABLE , DIMENSION(:)     :: A , B , X
                INTEGER   , ALLOCATABLE , DIMENSION(:)     :: JA, IA
            END SUBROUTINE pardiso_sym_f90

            SUBROUTINE FAILURE_CRITERION(STRESS , NUMBER_I , NUMBER_J)
                IMPLICIT NONE
                INTEGER , PARAMETER                        :: FC=SELECTED_REAL_KIND(8)
                REAL(FC) , ALLOCATABLE , DIMENSION(:,:)    :: STRESS
                INTEGER                                    :: NUMBER_I , NUMBER_J
            END SUBROUTINE FAILURE_CRITERION
   
!            SUBROUTINE OUTPLOT(STRESS , MAXSTRESS , MINSTRESS , NUMBER_I , NUMBER_J)
!                IMPLICIT NONE
!                INTEGER , PARAMETER                        :: OT=SELECTED_REAL_KIND(8)
!                REAL(OT) , ALLOCATABLE , DIMENSION(:,:)    :: STRESS
!                REAL(OT) , ALLOCATABLE , DIMENSION(:)      :: MAXSTRESS , MINSTRESS
!                INTEGER                                    :: NUMBER_I , NUMBER_J
!            END SUBROUTINE OUTPLOT
        END INTERFACE

        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: NODE_ST     !记录单元劲度矩阵在整体劲度矩阵的位置
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: B           !应变矩阵
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: BT          !转置
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: D           !弹性矩阵
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:,:) :: H
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:,:) :: KE
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: KT
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: STRESS 
        REAL(NLE) , ALLOCATABLE , DIMENSION(:)     :: A_LE        !一维数组存储非零值
        REAL(NLE) , ALLOCATABLE , DIMENSION(:,:)   :: GRID , MESH ,EXPORT
        INTEGER   , ALLOCATABLE , DIMENSION(:)     :: B_LE , C_LE !一维数值存储行列
        INTEGER   , ALLOCATABLE , DIMENSION(:)     :: TEMP_C
        INTEGER   , ALLOCATABLE , DIMENSION(:)     :: DIG_NUMBER
        INTEGER                                    :: N , NNZ , NUMBER_I , NUMBER_J
        INTEGER                                    :: I , J , K
        REAL(NLE)                                  :: TEMP
        CHARACTER*40                               :: FILENAME , FILENAME_PATH
!        OPEN (UNIT=100 , FILE='材料参数.TXT' ,STATUS='OLD')
!        
!        READ(100,*) MU
        ALLOCATE(B(3,6),BT(6,3),D(3,3),NODE_ST(6,ELEMENT),H(3,6,ELEMENT),KE(6,6,ELEMENT),KT(2*NODE,2*NODE))
!        MODULUS=5.D6   !后面需替换
        D      =0.D0 
        B      =0.D0
        BT     =0.D0
        NODE_ST=0.D0
        KT     =0.D0
        
        DO I = 1 , ELEMENT !求解每个单元弹性矩阵
            IF (MEA(I,11) .NE. 0) MODULUS(I)=0.D0 !开挖单元以空气填充
            D=0.D0
            D(1,1)=MODULUS(I)/(1.D0-MU**2)
            D(2,2)=MODULUS(I)/(1.D0-MU**2)
            D(1,2)=MODULUS(I)/(1.D0-MU**2)*MU
            D(2,1)=MODULUS(I)/(1.D0-MU**2)*MU
            D(3,3)=MODULUS(I)/(1.D0+MU)*0.5D0
            
            B=0.D0
            B(1,1)=MEA(I,7)-MEA(I,9)     ; B(1,3)=MEA(I,9)-MEA(I,5)     ; B(1,5)=MEA(I,5)-MEA(I,7)
            B(2,2)=MEA(I,8)-MEA(I,6)     ; B(2,4)=MEA(I,4)-MEA(I,8)     ; B(2,6)=MEA(I,6)-MEA(I,4)
            B(3,1)=MEA(I,8)-MEA(I,6)     ; B(3,2)=MEA(I,7)-MEA(I,9)     ; B(3,3)=MEA(I,4)-MEA(I,8)
            B(3,4)=MEA(I,9)-MEA(I,5)     ; B(3,5)=MEA(I,6)-MEA(I,4)     ; B(3,6)=MEA(I,5)-MEA(I,7)

            NODE_ST(1,I)=2.D0*MEA(I,1)-1 ; NODE_ST(2,I)=2.D0*MEA(I,1)   ; NODE_ST(3,I)=2.D0*MEA(I,2)-1
            NODE_ST(4,I)=2.D0*MEA(I,2)   ; NODE_ST(5,I)=2.D0*MEA(I,3)-1 ; NODE_ST(6,I)=2.D0*MEA(I,3)

            BT=TRANSPOSE(B)

            DO J = 1 , 3
                DO K = 1 , 6
                    H(J,K,I)=D(J,1)*B(1,K)+D(J,2)*B(2,K)+D(J,3)*B(3,K)
                END DO
            END DO

            DO J = 1 , 6
                DO K = 1 , 6
                    KE(J,K,I)=0.25D0*T_LE*(BT(J,1)*H(1,K,I)+BT(J,2)*H(2,K,I)+BT(J,3)*H(3,K,I))/MEA(I,10)
                END DO
            END DO
        END DO

        DO I = 1 , ELEMENT
            DO J = 1 , 6
                DO K = 1 , 6
                    KT(NODE_ST(J,I),NODE_ST(K,I))=KT(NODE_ST(J,I),NODE_ST(K,I))+KE(J,K,I)
                END DO
            END DO
        END DO
!!        ----------
!        OPEN(10,FILE='1.TXT')
!        WRITE(10,'(<2*NODE>F10.2)')((KT(I,J),J=1,2*NODE),I=1,2*NODE)
!        CLOSE(10)
!        PAUSE
!!        ----------
        DO I = 1 , 2*NODE
            IF (RESTRICTION(I)==0) THEN
                F(I)=0.D0
                DO J = 1 , 2*NODE
                    IF (I .NE. J) THEN
                        KT(I,J)=0.D0
                        KT(J,I)=0.D0
                    ELSE
                        KT(I,J)=1.D0
                    END IF
                END DO
            END IF
        END DO
        
        IF (NUMBER_I .NE. 0) THEN  !NUMBER_I==0表示在BEFORE里调用,此时网格为初始网格
            ALLOCATE( DIG_NUMBER( NODE ) )
            DIG_NUMBER=0
            FILENAME(1:3)='DIG'
            WRITE(FILENAME(4:7),'(I4)') NUMBER_I
            FILENAME(8:11)='.RES'
            FILENAME_PATH='TEMP_FILE\'//FILENAME
            OPEN (UNIT=100 , FILE=FILENAME_PATH , FORM='UNFORMATTED' , STATUS='OLD' , ACCESS='DIRECT' , RECL=2)
!            K=0
!            DO WHILE ( .NOT. EOF(100))
!                K=K+1
!                READ (100,*) DIG_NUMBER(K)
!!                WRITE (*,*) DIG_NUMBER(K)
!            END DO
            DO I = 1 , NODE
                READ (100 , REC=I ) DIG_NUMBER(I)
!                IF (DIG_NUMBER(I).NE.0) WRITE(*,*)DIG_NUMBER(I)
            END DO
            CLOSE(100)

!            DO I= 1 , NODE   !将挖去节点的对角线设为1，避免NAN
!                IF (DIG_NUMBER(I)==0) CYCLE
!!                WRITE(*,*)DIG_NUMBER(I)
!!                KT(2*I-1,:    )=0.D0
!!                KT(:    ,2*I-1)=0.D0
!!                KT(2*I  ,:    )=0.D0
!!                KT(:    ,2*I  )=0.D0
!                KT(2*I-1,2*I-1)=1.D0
!                KT(2*I  ,2*I  )=1.D0
!!                F(2*I)         =0.D0
!!                F(2*I-1)       =0.D0  
!            END DO
!            PAUSE
!            DO I = 1 , NODE  !将周围单元全破坏节点的对角线设为1，避免NAN,后面修改，用STATION判断
!                K=0
!LOOP1:          DO J = 1 , 10
!
!                    IF (MATCH(I,J)==0) CYCLE LOOP1
!                    IF ( STATION(MATCH(I,J))==0 ) THEN
!                        K=1
!                        EXIT LOOP1
!                    END IF
!                END DO LOOP1
!                IF (K==0) THEN
!                    KT(2*I-1,:    )=0.D0
!                    KT(:    ,2*I-1)=0.D0
!                    KT(2*I  ,:    )=0.D0
!                    KT(:    ,2*I  )=0.D0                    
!                    KT(2*I-1,2*I-1)=1.D0
!                    KT(2*I  ,2*I  )=1.D0 
!                END IF
!            END DO 
        END IF
        
        DO I = 1, 2*NODE
            IF (KT(I,I)==0.D0) THEN
                F(I)=0.D0
                KT(I,I)=1.D0
            END IF
        END DO

!        DO I = 1 , 2*NODE
!            IF (KT(I,I)==0.D0) WRITE(*,*) I
!        END DO
!        PAUSE      

!        ----------
!        OPEN(10,FILE='1.TXT')
!!        WRITE(10,'(<2*NODE>F10.2)')((KT(I,J),J=1,2*NODE),I=1,2*NODE)
!        WRITE(10,'(2F12.4)')(F(I),I=1,2*NODE)
!        CLOSE(10)
!        PAUSE
!        ----------
        DEALLOCATE(B,BT,D,NODE_ST,KE)

!------------------------------------------------------------------以下为用一维数组存储大型稀疏矩阵
        
        K=0                 !找出非零元素个数   
        DO I = 1 , 2*NODE
            DO J = I , 2*NODE
                IF (KT(I,J) .NE. 0.D0) THEN
                    K=K+1
                END IF
            END DO
        END DO        
        
        ALLOCATE(A_LE(K),B_LE(K),TEMP_C(K),C_LE(2*NODE+1))
        NNZ=K
        A_LE     =0.D0
        B_LE     =0
        B_LE     =0
        TEMP_C   =0
        
        K=1
        DO I = 1 , 2*NODE   !找出非零元素
            DO J = I , 2*NODE !对称矩阵只需存储上三角矩阵
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
        C_LE(2*NODE+1)=NNZ+1

        DEALLOCATE(KT,TEMP_C)

        N=2*NODE
        ALLOCATE(WY(N))

        CALL pardiso_sym_f90( N , NNZ , A_LE , B_LE , C_LE , F , WY )!矩阵求解

!------------------------------------------------------------------求各单元应力
        ALLOCATE(STRESS(ELEMENT,3))
        STRESS=0.D0
        DO I = 1 , ELEMENT
            DO J = 1 , 3
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,1,I)*WY(MEA(I,1)*2-1)/MEA(I,10)  !单位为KPA
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,2,I)*WY(MEA(I,1)*2  )/MEA(I,10)
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,3,I)*WY(MEA(I,2)*2-1)/MEA(I,10)
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,4,I)*WY(MEA(I,2)*2  )/MEA(I,10)
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,5,I)*WY(MEA(I,3)*2-1)/MEA(I,10)
                STRESS(I,J)=STRESS(I,J)+5.D-1*H(J,6,I)*WY(MEA(I,3)*2  )/MEA(I,10)
            END DO
            IF (ISNAN(STRESS(I,3))) THEN
                WRITE(*,*)NUMBER_I,NUMBER_J , I
                !PAUSE
            END IF
        END DO
        
        !写总应力
        !WRITE(*,*)5
        
        IF (NUMBER_I==0) THEN
            OPEN(UNIT=101 , FILE= '初始应力分布.RES' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*3)
            DO I = 1 , ELEMENT
                WRITE(101 , REC=I)(STRESS(I,J),J=1,3)
            END DO
            CLOSE(101)
            DEALLOCATE(MODULUS)
! --------------------------------------
            OPEN(UNIT=101 , FILE= '初始应力分布作图.DAT')
            OPEN (UNIT=103 , FILE='初始网格节点坐标.TXT'         ,STATUS='OLD')
            OPEN (UNIT=102 , FILE='初始网格单元节点号.TXT'       ,STATUS='OLD')
            ALLOCATE(EXPORT(NODE,3))
            ALLOCATE(GRID(NODE,2))
            ALLOCATE(MESH(ELEMENT,3))
            READ (103,*) ((GRID(I,J),J=1,2),I=1,NODE)
            READ (102,*) ((MESH(I,J),J=1,3),I=1,ELEMENT)
            EXPORT=0.D0
            DO I=1,NODE
                TEMP=0.D0
                DO J=1,10
                    IF (MATCH(I,J)==0)CYCLE
                    EXPORT(I,1)=EXPORT(I,1)+STRESS(MATCH(I,J),1)*MEA(MATCH(I,J),10)
                    EXPORT(I,2)=EXPORT(I,2)+STRESS(MATCH(I,J),2)*MEA(MATCH(I,J),10)
                    EXPORT(I,3)=EXPORT(I,3)+STRESS(MATCH(I,J),3)*MEA(MATCH(I,J),10)
                    TEMP=TEMP+MEA(MATCH(I,J),10)
                END DO
!                WRITE(*,*)TEMP
                EXPORT(I,1)=EXPORT(I,1)/TEMP
                EXPORT(I,2)=EXPORT(I,2)/TEMP
                EXPORT(I,3)=EXPORT(I,3)/TEMP
            END DO
            WRITE(101,*)" title= ""triangle"" "  
            WRITE(101,*)" variables=x,y,z1,z2,z3,z4,z5"
            WRITE(101,*)"ZONE N=",NODE,", E=",ELEMENT,", F=FEPOINT, ET=triangle"
            DO I=1,NODE
                WRITE(101,'(7E12.4)')GRID(I,1),GRID(I,2),(EXPORT(I,J),J=1,3),WY(2*I-1),WY(2*I)
            END DO
            WRITE (101,*) ((MESH(I,J),J=1,3),I=1,ELEMENT)
            CLOSE(101)
            CLOSE(102)
            CLOSE(103)
            DEALLOCATE(EXPORT,GRID,MESH)
!            PAUSE
!---------------------------------------
        ELSE
            !WRITE(*,*)5
            !CALL FAILURE_CRITERION(STRESS , NUMBER_I , NUMBER_J)
!            DEALLOCATE( TSTRESS )
        END IF

!        DEALLOCATE( MEA )
!        DEALLOCATE( MATCH )
!        DEALLOCATE( RESTRICTION )
        IF (NUMBER_I==0) DEALLOCATE( F )
        DEALLOCATE( H )
!        DEALLOCATE( MODULUS )
!        DEALLOCATE( A_LE )
!        DEALLOCATE( B_LE )
!        DEALLOCATE( C_LE )
        DEALLOCATE( WY )
        DEALLOCATE( STRESS )
!        CALL OUTPLOT(STRESS , MAXSTRESS , MINSTRESS , NUMBER_I , NUMBER_J)
        


!!        ----------
!        OPEN(11,FILE='3.TXT')
!        WRITE(11,'(3F12.6)')((STRESS(I,J),J=1,3),I=1,ELEMENT)
!        CLOSE(11)
!        PAUSE
!!        ---------- 
        RETURN
!       IF (NUMBER_I/=0) DEALLOCATE(TSTRESS)
    END SUBROUTINE LINEAR_ELASTIC_STIFFNESS

END MODULE LINEAR_ELASTIC
#endif 