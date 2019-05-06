#include 'define.h'
SUBROUTINE CONTINUITY_MOMENTUM_X(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DXX,DYY,I1,I2,IP,UPS,DWS)

USE HYDRO

implicit none
integer,parameter::nnx=selected_real_kind(8)
real(nnx),dimension(:,:),pointer ::pp1,pp2,qq1,qq2,hh1,hh2
real(nnx),dimension(:,:),pointer :: dxx,dyy
real(nnx)::t0
integer::i1,i2,ip
logical::ups,dws,RE_CAL

real(nnx) :: fr , al
real(nnx) :: adve , advw
real(nnx) :: vn , vs , pn , ps , un , us , um , hn , hs , hm , cadl , cadr , cdiffl , cdiffr
real(nnx) :: resist , chezy  , qav , hav , grv
real(nnx) :: diffx , diffy, eddy, temp_u 
real(nnx) , pointer :: dx , dx2 , dy , dy2
integer   :: ml , mu , nl , nu
integer   :: i , j , k , jc , jjc , j0
real(nnx) :: gra
character(len=80)::err_msg
integer::status

ml=lbound(pp1,2)
mu=ubound(pp1,2)
nl=lbound(pp1,1)
nu=ubound(pp1,1)
gra=9.8d0
dx=>dxx(1,1)
dy=>dyy(1,1)

loop1:do i= i1 , i2 , ip
 100  CONTINUE   
      RE_CAL=.FALSE.  
      !continuity equation
      jc=-1
      loop2:do j= ml , mu , 1
                jc=jc+2
                if(ndd(i,j)==-1 .OR. ndd(i,j)== 1) then !陆地
                    !IF (HH1(I,J)<0.1) HH1(I,J)=0.1D0
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = hh1(i,j)                
                    cycle loop2
                end if
                if(ndd(i,j)== 2) then !水位边界
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = hh2(i,j)
                    cycle loop2
                end if
                !if(ndd(i,j)== 4) then !NEUMAN边界
                !    if (j==ml) then
                !        amt(1,jc) = 0.d0
                !        amt(2,jc) = 0.d0
                !        amt(3,jc) = 1.d0
                !        amt(4,jc) = 0.d0
                !        amt(5,jc) =-1.d0
                !        amt(6,jc) = 0.d0
                !    else if (j==mu) then
                !        amt(1,jc) =-1.d0
                !        amt(2,jc) = 0.d0
                !        amt(3,jc) = 1.d0
                !        amt(4,jc) = 0.d0
                !        amt(5,jc) = 0.d0
                !        amt(6,jc) = 0.d0                        
                !    end if
                !    cycle loop2
                !end if                
                amt(1,jc) =  0.d0
                amt(2,jc) = -dt*dy
                amt(3,jc) = +4.d0*dx*dy
                amt(4,jc) = +dt*dy
                amt(5,jc) = +0.d0
                amt(6,jc) = +4.d0*dx*dy*hh1(i,j) &
                &            -dt*dy*(pp1(i,j)-pp1(i,j-1)) &
                &            -dt*dx*(qq2(i,j)-qq2(i-1,j)+qq1(i,j)-qq1(i-1,j))
            end do loop2
            !momentum equation
            jc=0
      loop3:do j=ml , mu , 1
                jc=jc+2
                !close boundary condition
                if(ndd(i,j)==-1 .OR. VDDX(I,J)==-1) then
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = 0.d0
                    cycle loop3
                end if
                if(ndd(i,j+1)==-1 ) then
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = 0.d0
                    cycle loop3
                end if                
                !waterlevel boundary condition   
                if( ndd(i,j)== 2 .or. ndd(i,j+1)==2 ) then
                    grv=gra*0.5d0*(hh1(i,j)+hh1(i,j+1))
                    amt(1,jc) =   0.d0
                    amt(2,jc) = - grv*dt*dy
                    amt(3,jc) =   dx*dy
                    amt(4,jc) =   grv*dt*dy
                    amt(5,jc) =   0.d0
                    amt(6,jc) =   pp1(i,j)*dx*dy &
                    &           + grv*dt*dy*(dep(i,j+1)-dep(i,j))  
                    cycle loop3
                end if
               
                 !froude number
                 fr=abs(pp1(i,j))/sqrt(0.125d0*gra*(hh1(i,j)+hh1(i,j+1))**3) 
                 if(fr>=1.d0) then
                   al=1.d0
                 else if(fr<1.d0 .and. fr>0.25d0) then
                    al=(fr-0.25d0)*4.d0/3.d0 
                 else
                    al=0.d0
                 end if
                 if(pp1(i,j)<0.d0) al=-al
                 
                 !advection momentum
                 adve=0.25d0*((1.d0-al)*pp1(i,j+1)+(1.d0+al)*pp1(i,j))  &
                   & /((1.d0-dabs(al))*hh1(i,j+1)+dmax1(al,0.d0)*hh1(i,j)-dmin1(al,0.d0)*hh1(i,j+2))  !
                 advw=0.25d0*((1d0-al)*pp1(i,j)+(1d0+al)*pp1(i,j-1)) &
                   & /((1.d0-dabs(al))*hh1(i,j)+dmax1(al,0.d0)*hh1(i,j-1)-dmin1(al,0.d0)*hh1(i,j+1))  
                 
                 !cross momentum and cross diffusion
                 if(ndd(i+1,j)==0 .and. ndd(i+1,j+1)==0) then
                    vn = 2.d0*(qq2(i,j)+qq2(i,j+1))/(hh1(i,j)+hh1(i+1,j)+hh1(i,j+1)+hh1(i+1,j+1))
                    pn = pp1(i+1,j)
                    un = 2.d0*pn/(hh1(i+1,j)+hh1(i+1,j+1))
                    hn = (hh1(i,j)+hh1(i+1,j)+hh1(i,j+1)+hh1(i+1,j+1))/4.d0
                 else if(ndd(i+1,j)/=0 .and. ndd(i+1,j+1)==0) then
                    vn = 2.d0*(qq2(i,j)+qq2(i,j+1))/(hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i+1,j+1))
                    pn = pp1(i,j)
                    un = 2.d0*pn/(hh1(i,j)+hh1(i,j+1))
                    hn = (hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i+1,j+1))/4.d0
                 else if(ndd(i+1,j)==0 .and. ndd(i+1,j+1)/=0) then
                    vn = 2.d0*(qq2(i,j)+qq2(i,j+1))/(hh1(i,j)+hh1(i+1,j)+hh1(i,j+1)+hh1(i,j+1))
                    pn = pp1(i,j)
                    un = 2.d0*pn/(hh1(i,j)+hh1(i,j+1))
                    hn = (hh1(i,j)+hh1(i+1,j)+hh1(i,j+1)+hh1(i,j+1))/4.d0                    
                 else
                    vn = 0.d0
                    pn = pp1(i,j)
                    un = 2.d0*pn/(hh1(i,j)+hh1(i,j+1))
                    hn = (hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i,j+1))/4.d0                    
                 end if
                 
                 if(ndd(i-1,j)==0 .and. ndd(i-1,j+1)==0) then
                    vs = 2.d0*(qq2(i-1,j)+qq2(i-1,j+1))/(hh1(i,j)+hh1(i-1,j)+hh1(i,j+1)+hh1(i-1,j+1))
                    ps = pp1(i-1,j)
                    us = 2.d0*ps/(hh1(i-1,j)+hh1(i-1,j+1))
                    hs = (hh1(i,j)+hh1(i-1,j)+hh1(i,j+1)+hh1(i-1,j+1))/4.d0                    
                 else if(ndd(i-1,j)/=0 .and. ndd(i-1,j+1)==0) then
                    vs = 2.d0*(qq2(i-1,j)+qq2(i-1,j+1))/(hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i-1,j+1))
                    ps = pp1(i,j)
                    us = 2.d0*ps/(hh1(i,j)+hh1(i,j+1))
                    hs = (hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i-1,j+1))/4.d0                     
                 else if(ndd(i-1,j)==0 .and. ndd(i-1,j+1)/=0) then
                    vs = 2.d0*(qq2(i-1,j)+qq2(i-1,j+1))/(hh1(i,j)+hh1(i-1,j)+hh1(i,j+1)+hh1(i,j+1))
                    ps = pp1(i,j)
                    us = 2.d0*ps/(hh1(i,j)+hh1(i,j+1))
                    hs = (hh1(i,j)+hh1(i-1,j)+hh1(i,j+1)+hh1(i,j+1))/4.d0                     
                 else
                    vs = 0.d0
                    ps = pp1(i,j)
                    us = 2.d0*ps/(hh1(i,j)+hh1(i,j+1))
                    hs = (hh1(i,j)+hh1(i,j)+hh1(i,j+1)+hh1(i,j+1))/4.d0                     
                 end if
                 
                 !diffy=0.5d0*(0.5d0*vn+0.5d0*vs)**2*dt+eddy_y
                 diffy=0.5d0*(0.5d0*vn+0.5d0*vs)**2*dt
                 
                 if(dws) then
                    if(ndd(i+1,j)==0 .and. ndd(i+1,j+1)==0) then
                        cadl   = -0.5d0*vs
                        cadr   =  0.5d0*(-vn*pp2(i+1,j)+vs*ps-vn*pp1(i,j))
                        cdiffl =  diffy
                        cdiffr =  diffy*(pp2(i+1,j)-pp1(i,j)+ps) 
                    else
                        cadl   =  0.5d0*(vn-vs)
                        cadr   =  0.5d0*(vs*ps-vn*pp1(i,j))
                        cdiffl =  0.d0
                        cdiffr =  diffy*(-pp1(i,j)+ps)  
                    end if
                    
                  else if(ups) then
                    if(ndd(i-1,j)==0 .and. ndd(i-1,j+1)==0) then
                        cadl  = 0.5d0*vn
                        cadr  = 0.5d0*(-vn*pn+vs*pp2(i-1,j)+vs*pp1(i,j))
                        cdiffl= diffy
                        cdiffr= diffy*(pn-pp1(i,j)+pp2(i-1,j))     
                    else
                        cadl  = 0.5d0*(vn-vs)
                        cadr  = 0.5d0*(-vn*pn+vs*pp1(i,j))
                        cdiffl= 0.d0
                        cdiffr= diffy*(pn-pp1(i,j))  
                    end if
                  end if   
                  
                  !eddy viscosity
                  hm = (hh1(i,j)+hh1(i,j+1))/2.d0 
                  um = pp1(i,j)/hm
                  !hn = max(hh1(i,j),hh1(i+1,j),hh1(i,j+1),hh1(i+1,j+1))
                  !hs = max(hh1(i,j),hh1(i-1,j),hh1(i,j+1),hh1(i-1,j+1))
                  hn = hm
                  hs = hm
                  eddy = eddy_y*(hn*(un-um)/dy-hs*(um-us)/dy)/dy + & 
                &        eddy_x*(hm/dx*( 2*pp1(i,j+1)/(hh1(i,j+1)+hh1(i,j+2))-2*pp1(i,j  )/(hh1(i,j)+hh1(i,j+1)) )- &
                &                hm/dx*( 2*pp1(i,j  )/(hh1(i,j  )+hh1(i,j+1))-2*pp1(i,j-1)/(hh1(i,j)+hh1(i,j-1)) ))/dx
                  
                  !gravity
                  grv=gra*0.5d0*(hh1(i,j)+hh1(i,j+1))
                  
                  !resistance term
                  qav=0.125d0*(qq1(i,j)+qq1(i-1,j)+qq1(i,j+1)+qq1(i-1,j+1) &
                  &           +qq2(i,j)+qq2(i-1,j)+qq2(i,j+1)+qq2(i-1,j+1))  
                  if(pp1(i,j)>=0.d0) then
                    hav=hh1(i,j)
                  else if(pp1(i,j)<0.d0) then
                    hav=hh1(i,j+1)
                  end if
#ifdef PRE_CHEZY 
                  !常系数谢才公式
                  resist=gra/CHE**2*sqrt(pp1(i,j)**2+qav**2)/hav**2.D0
                  ! bianjie jia zao
                  !if (i<=m_start+5 .or. i>=m_end-5) then
                  !    resist=gra/CHE1**2*sqrt(pp1(i,j)**2+qav**2)/hav**2.D0
                  !end if
#endif

#ifdef PRE_MANNING                   
                  !曼宁公式
                  resist=gra/manning**2*sqrt(pp1(i,j)**2+qav**2)/hav**(7.d0/3.d0)
                  !if (i<=m_start+5 .or. i>=m_end-5) then
                  !    resist=gra/che1**2*sqrt(pp1(i,j)**2+qav**2)/hav**(7.d0/3.d0)
                  !end if
#endif                     

                  !x-diffusion term
                  !diffx=0.5d0*(2.d0*pp1(i,j)/(hh1(i,j+1)+hh1(i,j)))**2*dt+eddy_x 
                  diffx=0.5d0*(2.d0*pp1(i,j)/(hh1(i,j+1)+hh1(i,j)))**2*dt
                  
                  amt(1,jc) = - (1.d0+al)*advw*dt*dy &
                            & - diffx*dt*dy/dx
                  
                  amt(2,jc) = - grv*dt*dy
                  
                  amt(3,jc) = + dx*dy &
                            & + ((1.d0+al)*adve-(1.d0-al)*advw)*dt*dy &
                            & + cadl*dt*dx                            &
                            & + resist*dt*dy*dx                       &
                            & + 2.d0*diffx*dt*dy/dx                   &
                            & + cdiffl*dt*dx/dy 
                                            
                  amt(4,jc) = + grv*dt*dy
                  
                  amt(5,jc) = + (1d0-al)*adve*dt*dy                   &
                            & - diffx*dt*dy/dx
                  
                  amt(6,jc) = + dx*dy*pp1(i,j)                        &
                            & + cadr*dt*dx                            &
                            & + cdiffr*dt*dx/dy                       &
                            & + grv*dt*dy*(dep(i,j+1)-dep(i,j))       &       
                            & + eddy*dt*dx*dy
            end do loop3
            
            !jc=jc+1 !最后一行为连续性方程   
            
            call solve5(jc)   
            
            j0 = 0
            do j= 2 , jc , 2
                j0 = j0 + 1
                pp2(i,j0)=amt(6,j)
            end do
            j0 = 0
            do j= 1 , jc , 2
                j0 = j0 + 1 
                hh2(i,j0)=amt(6,j)
            end do
!----------------------------------------------------------------------------------
            DO J=N_START+1,N_END-1
                IF (NDD(I,J)==-1) CYCLE
                IF (HH2(I,J)<0.04D0) THEN
                    !temp_u=2*pp2(i,j)/(hh2(i,j)+hh2(i,j+1))
                    NDD(I,J)=1
                    VDDX(I,J)=-1
                    VDDX(I,J-1)=-1
                    VDDY(I,J)=-1
                    VDDY(I-1,J)=-1
                    !QQ2(I,J)=0.D0
                    !QQ2(I-1,J)=0.D0
                    RE_CAL=.TRUE.
                    !write(*,*)i,j,temp_u,1, hh2(i,j), hh1(i,j)
                    !pause
                END IF
            END DO
            IF (RE_CAL) GOTO 100
!----------------------------------------------------------------------------------
    end do loop1     

    ! neuman边界赋值
    !do j=ml, mu
    !    hh2(m_start,j)=hh2(m_start+1,j)
    !    hh2(m_end,j)=hh2(m_end-1,j)
    !    pp2(m_start,j)=pp2(m_start+1,j)
    !    pp2(m_end,j)=pp2(m_end-1,j)
    !end do

    end subroutine CONTINUITY_MOMENTUM_X                    
    
SUBROUTINE UV_CACULATION(I,J,U,V,UV)
    USE TIME_DIN
    USE HYDRO
    IMPLICIT NONE
    INTEGER  ::  I , J 
    REAL*8   ::  U , V , UV
    
    
    
    IF (GLOBTIME%GLOBSTEP%SECSTEP/2==GLOBTIME%GLOBSTEP%SECSTEP/2.0) THEN
        
        !CULCULATING
        IF (ABS(NDD(I,J))==1) THEN
            U=0.D0
            V=0.D0
            UV=0.D0
        ELSE 
            U=(P1(I,J)+P1(I,J-1))/H1(I,J)/2.D0
            V=(Q1(I,J)+Q1(I-1,J))/H1(I,J)/2.D0
            UV=SQRT(U**2+V**2)
        END IF        
        
    ELSE       
        
        !CULCULATING
        IF (ABS(NDD(I,J))==1) THEN
            U=0.D0
            V=0.D0
            UV=0.D0
        ELSE 
            U=(P2(I,J)+P2(I,J-1))/H1(I,J)/2.D0
            V=(Q2(I,J)+Q2(I-1,J))/H1(I,J)/2.D0
            UV=SQRT(U**2+V**2)
        END IF  
    
    END IF  
    
END SUBROUTINE