#include 'define.h'
SUBROUTINE CONTINUITY_MOMENTUM_y(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DXX,DYY,I1,I2,IP,UPS,DWS)

use HYDRO

implicit none
integer,parameter::nny=selected_real_kind(8)
real(nny),dimension(:,:),pointer ::pp1,pp2,qq1,qq2,hh1,hh2
real(nny),dimension(:,:),pointer :: dxx,dyy
real(nny)::t0
integer::i1,i2,ip
logical::ups,dws,RE_CAL

real(nny) :: fr , al
real(nny) :: adve , advw
real(nny) :: vn , vs , pn , ps , un , us , um , hn , hs , hm , cadl , cadr , cdiffl , cdiffr
real(nny) :: resist , chezy  , qav , hav , grv
real(nny) :: diffx , diffy , eddy, temp_u
real(nny) , pointer :: dx , dx2 , dy , dy2
integer   :: ml , mu , nl , nu
integer   :: i , j , k , jc , jjc , j0
real(nny) :: gra
character(len=80)::err_msg
integer::status

ml=lbound(pp1,1)
mu=ubound(pp1,1)
nl=lbound(pp1,2)
nu=ubound(pp1,2)
gra=9.8d0
!eddy_x=0.d0
!eddy_y=0.d0
!manning=32.d0
dx=>dxx(1,1)
dy=>dyy(1,1)
 !       open(1,file='ndd.txt')
 !write(1,'(101i4)')((VddY(i,j),i=1,101),j=1,401)
 !close(1)
 !pause
loop1:do i= i1 , i2 , ip
 101  CONTINUE
      RE_CAL=.FALSE.     
      !continuity equation
      jc=-1
      loop2:do j= ml , mu , 1
                jc=jc+2
                if(ndd(j,i)==-1 .OR. ndd(j,i)== 1) then !陆地
                    !IF (HH1(J,I)<0.01) HH1(J,I)=0.01D0
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = hh1(j,i)
                    cycle loop2
                end if
                if(ndd(j,i)== 2) then !水位边界
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = hh2(j,i)
                    cycle loop2
                end if
                if(ndd(j,i)== 4) then !NEUMAN边界
                    if (j==ml) then
                        amt(1,jc) = 0.d0
                        amt(2,jc) = 0.d0
                        amt(3,jc) = 1.d0
                        amt(4,jc) = 0.d0
                        amt(5,jc) =-1.d0
                        amt(6,jc) = 0.d0
                    else if (j==mu) then
                        amt(1,jc) =-1.d0
                        amt(2,jc) = 0.d0
                        amt(3,jc) = 1.d0
                        amt(4,jc) = 0.d0
                        amt(5,jc) = 0.d0
                        amt(6,jc) = 0.d0
                    end if
                    cycle loop2
                end if                
                amt(1,jc) =  0.d0
                amt(2,jc) = -dt*dy
                amt(3,jc) =  4.d0*dx*dy
                amt(4,jc) =  dt*dy
                amt(5,jc) =  0.d0
                amt(6,jc) =  4.d0*dx*dy*hh1(j,i) &
                &            -dt*dy*(pp1(j,i)-pp1(j-1,i)) &
                &            -dt*dx*(qq2(j,i)-qq2(j,i-1)+qq1(j,i)-qq1(j,i-1))
            end do loop2
            !momentum equation
            jc=0
      loop3:do j=ml , mu , 1
                jc=jc+2
                !close boundary condition
                if(ndd(j,i)==-1 .OR. VDDY(J,I)==-1) then
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = 0.d0
                    cycle loop3
                end if
                if(ndd(j+1,i)==-1) then
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = 0.d0
                    cycle loop3
                end if                
                !waterlevel boundary condition   
                if(ndd(j,i)== 2 .or. ndd(j+1,i)==2 ) then
                    grv=gra*0.5d0*(hh1(j,i)+hh1(j+1,i))
                    amt(1,jc) =   0.d0
                    amt(2,jc) = - grv*dt*dy
                    amt(3,jc) =   dx*dy
                    amt(4,jc) =   grv*dt*dy
                    amt(5,jc) =   0.d0
                    amt(6,jc) =   pp1(j,i)*dx*dy &
                    &           + grv*dt*dy*(dep(j+1,i)-dep(j,i))  
                    cycle loop3
                end if
                if(ndd(j,i)== 4 .and. j==ml) then !水位边界,目前不可用，没调试好
                    amt(1,jc) = 0.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) =-1.d0
                    amt(6,jc) = 0.d0
                    cycle loop3
                end if
                
                if(ndd(j,i)== 4 .and. j==mu) then !水位边界
                    amt(1,jc) =-1.d0
                    amt(2,jc) = 0.d0
                    amt(3,jc) = 1.d0
                    amt(4,jc) = 0.d0
                    amt(5,jc) = 0.d0
                    amt(6,jc) = 0.d0  
                    cycle loop3
                end if   
                 
                 !froude number
                 fr=abs(pp1(j,i))/sqrt(0.125d0*gra*(hh1(j,i)+hh1(j+1,i))**3) 
                 if(fr>=1.d0) then
                   al=1.d0
                 else if(fr<1.d0 .and. fr>0.25d0) then
                    al=(fr-0.25d0)*4.d0/3.d0 
                 else
                    al=0.d0
                 end if
                 if(pp1(j,i)<0.d0) al=-al
                 
                 !advection momentum
                 adve=0.25d0*((1.d0-al)*pp1(j+1,i)+(1.d0+al)*pp1(j,i))  &
                   & /((1.d0-dabs(al))*hh1(j+1,i)+dmax1(al,0.d0)*hh1(j,i)-dmin1(al,0.d0)*hh1(j+2,i))  !
                 advw=0.25d0*((1d0-al)*pp1(j,i)+(1d0+al)*pp1(j-1,i)) &
                   & /((1.d0-dabs(al))*hh1(j,i)+dmax1(al,0.d0)*hh1(j-1,i)-dmin1(al,0.d0)*hh1(j+1,i))  
                 
                 !cross momentum and cross diffusion
                 if(ndd(j,i+1)==0 .and. ndd(j+1,i+1)==0) then
                    vn = 2.d0*(qq2(j,i)+qq2(j+1,i))/(hh1(j,i)+hh1(j,i+1)+hh1(j+1,i)+hh1(j+1,i+1))
                    pn = pp1(j,i+1)
                    un = 2.d0*pn/(hh1(j,i+1)+hh1(j+1,i+1))
                    hn = (hh1(j,i)+hh1(j,i+1)+hh1(j+1,i)+hh1(j+1,i+1))/4.d0
                 else if(ndd(j,i+1)/=0 .and. ndd(j+1,i+1)==0) then
                    vn = 2.d0*(qq2(j,i)+qq2(j+1,i))/(hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i+1))
                    pn = pp1(j,i)
                    un = 2.d0*pn/(hh1(j,i)+hh1(j+1,i))
                    hn = (hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i+1))/4.d0                    
                 else if(ndd(j,i+1)==0 .and. ndd(j+1,i+1)/=0) then
                    vn = 2.d0*(qq2(j,i)+qq2(j+1,i))/(hh1(j,i)+hh1(j,i+1)+hh1(j+1,i)+hh1(j+1,i))
                    pn = pp1(j,i)
                    un = 2.d0*pn/(hh1(j,i)+hh1(j+1,i))
                    hn = (hh1(j,i)+hh1(j,i+1)+hh1(j+1,i)+hh1(j+1,i))/4.d0                     
                 else
                    vn = 0.d0
                    pn = pp1(j,i)
                    un = 2.d0*pn/(hh1(j,i)+hh1(j+1,i))
                    hn = (hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i))/4.d0                    
                 end if
                 
                 if(ndd(j,i-1)==0 .and. ndd(j+1,i-1)==0) then
                    vs = 2.d0*(qq2(j,i-1)+qq2(j+1,i-1))/(hh1(j,i)+hh1(j,i-1)+hh1(j+1,i)+hh1(j+1,i-1))
                    ps = pp1(j,i-1)
                    us = 2.d0*ps/(hh1(j,i-1)+hh1(j+1,i-1))
                    hs = (hh1(j,i)+hh1(j,i-1)+hh1(j+1,i)+hh1(j+1,i-1))/4.d0                     
                 else if(ndd(j,i-1)/=0 .and. ndd(j+1,i-1)==0) then
                    vs = 2.d0*(qq2(j,i-1)+qq2(j+1,i-1))/(hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i-1))
                    ps = pp1(j,i)
                    us = 2.d0*ps/(hh1(j,i)+hh1(j+1,i))
                    hs = (hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i-1))/4.d0                     
                 else if(ndd(j,i-1)==0 .and. ndd(j+1,i-1)/=0) then
                    vs = 2.d0*(qq2(j,i-1)+qq2(j+1,i-1))/(hh1(j,i)+hh1(j,i-1)+hh1(j+1,i)+hh1(j+1,i))
                    ps = pp1(j,i)
                    us = 2.d0*ps/(hh1(j,i)+hh1(j+1,i))
                    hs = (hh1(j,i)+hh1(j,i-1)+hh1(j+1,i)+hh1(j+1,i))/4.d0                     
                 else
                    vs = 0.d0
                    ps = pp1(j,i)
                    us = 2.d0*ps/(hh1(j,i)+hh1(j+1,i))
                    hs = (hh1(j,i)+hh1(j,i)+hh1(j+1,i)+hh1(j+1,i))/4.d0                     
                 end if
                 
                 !diffy=0.5d0*(0.5d0*vn+0.5d0*vs)**2*dt+eddy_y
                 diffy=0.5d0*(0.5d0*vn+0.5d0*vs)**2*dt
                 
                 if(dws) then
                    if(ndd(j,i+1)==0 .and. ndd(j+1,i+1)==0) then
                        cadl   = -0.5d0*vs
                        cadr   =  0.5d0*(-vn*pp2(j,i+1)+vs*ps-vn*pp1(j,i))
                        cdiffl =  diffy
                        cdiffr =  diffy*(pp2(j,i+1)-pp1(j,i)+ps) 
                    else
                        cadl   =  0.5d0*(vn-vs)
                        cadr   =  0.5d0*(vs*ps-vn*pp1(j,i))
                        cdiffl =  0.d0
                        cdiffr =  diffy*(-pp1(j,i)+ps)  
                    end if
                  else if(ups) then
                    if(ndd(j,i-1)==0 .and. ndd(j+1,i-1)==0) then
                        cadl  = 0.5d0*vn
                        cadr  = 0.5d0*(-vn*pn+vs*pp2(j,i-1)+vs*pp1(j,i))
                        cdiffl= diffy
                        cdiffr= diffy*(pn-pp1(j,i)+pp2(j,i-1))     
                    else
                        cadl  = 0.5d0*(vn-vs)
                        cadr  = 0.5d0*(-vn*pn+vs*pp1(j,i))
                        cdiffl= 0.d0
                        cdiffr= diffy*(pn-pp1(j,i))  
                    end if
                  end if   
                  
                  !eddy viscosity
                  hm = (hh1(j,i)+hh1(j+1,i))/2.d0 
                  um = pp1(j,i)/hm
                  !hn = max(hh1(j,i),hh1(j,i+1),hh1(j+1,i),hh1(j+1,i+1))
                  !hs = max(hh1(j,i),hh1(j,i-1),hh1(j+1,i),hh1(j+1,i-1))
                  hn = hm
                  hs = hm
                  eddy = eddy_y*(hn*(un-um)/dy-hs*(um-us)/dy)/dy + & 
                &        eddy_x*(hm/dx*( 2*pp1(j+1,i)/(hh1(j+1,i)+hh1(j+2,i))-2*pp1(j,i  )/(hh1(j,i)+hh1(j+1,i)) )- &
                &                hm/dx*( 2*pp1(j,i  )/(hh1(j,i  )+hh1(j+1,i))-2*pp1(j-1,i)/(hh1(j,i)+hh1(j-1,i)) ))/dx 
                  
                  !gravity
                  grv=gra*0.5d0*(hh1(j,i)+hh1(j+1,i))
                  
                  !resistance term
                  qav=0.125d0*(qq1(j,i)+qq1(j,i-1)+qq1(j+1,i)+qq1(j+1,i-1) &
                  &           +qq2(j,i)+qq2(j,i-1)+qq2(j+1,i)+qq2(j+1,i-1))  
                  if(pp1(j,i)>=0.d0) then
                    hav=hh1(j,i)
                  else if(pp1(j,i)<0.d0) then
                    hav=hh1(j+1,i)
                  end if
#ifdef PRE_CHEZY                   
                  !常系数谢才公式
                  resist=gra/CHE**2*sqrt(pp1(j,i)**2+qav**2)/hav**2.D0
                  ! bianjie jia zao
                  if (j<=m_start+5 .or. j>=m_end-5) then
                      resist=gra/CHE1**2*sqrt(pp1(j,i)**2+qav**2)/hav**2.D0
                  end if                  
#endif

#ifdef PRE_MANNING                   
                  !曼宁公式                  
                  !resist=gra/manning**2*sqrt(pp1(j,i)**2+qav**2)/hav**(7.d0/3.d0)  
                  !if (j<=m_start+5 .or. j>=m_end-5) then
                  !    resist=gra/che1**2*sqrt(pp1(j,i)**2+qav**2)/hav**(7.d0/3.d0)
                  !end if                  
#endif
                  
                  
                  !x-diffusion term
                  !diffx=0.5d0*(2.d0*pp1(j,i)/(hh1(j+1,i)+hh1(j,i)))**2*dt+eddy_x 
                  diffx=0.5d0*(2.d0*pp1(j,i)/(hh1(j+1,i)+hh1(j,i)))**2*dt
                  
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
                  
                  amt(6,jc) = + dx*dy*pp1(j,i)                        &
                            & + cadr*dt*dx                            &
                            & + cdiffr*dt*dx/dy                       &
                            & + grv*dt*dy*(dep(j+1,i)-dep(j,i))       &      
                            & + eddy*dt*dx*dy
                      
            end do loop3
            
            !jc=jc+1 !最后一行为连续性方程  
      
            !OPEN(1,FILE='AMT5.TXT')
            !DO K = 1 ,  jc 
            !    WRITE(1,'(6E15.4)') (AMT(J,K),J=1,6)
            !END DO
            !CLOSE(1)
            !PAUSE      
            
            call solve5(jc)   
            
            !OPEN(1,FILE='AMT5.TXT')
            !DO K = 1 ,  jc 
            !    WRITE(1,'(6E15.4)') (AMT(J,K),J=1,6)
            !END DO
            !CLOSE(1)
            !PAUSE
            
            j0 = 0
            do j= 2 , jc , 2
                j0 = j0 + 1
                pp2(j0,i)=amt(6,j)
            end do
            j0 = 0
            do j= 1 , jc , 2
                j0 = j0 + 1 
                hh2(j0,i)=amt(6,j)
            end do
!--------------------------------------------------------------------------------
            DO J=M_START+1,M_END-1
                IF (NDD(J,I)==-1) CYCLE
                IF (HH2(J,I)<=0.04D0) THEN
                    !temp_u=2*pp2(j,i)/(hh2(j,i)+hh2(j+1,i))
                    !if (temp_u>1.5) then
                        NDD(J,I)=1
                        VDDY(J,I)=-1
                        VDDY(J-1,I)=-1
                        VDDX(J,I)=-1
                        VDDX(J,I-1)=-1
                        RE_CAL=.TRUE.
                        !write(*,*)i,j,temp_u,2
                    !end if
                END IF
            END DO
            IF (RE_CAL) GOTO 101
!----------------------------------------------------------------------------------

      end do loop1      
end subroutine             