#include 'define.h' 
    subroutine solve5(jc)
      use HYDRO
      implicit none
      integer :: i , j , jc
      DOUBLE PRECISION :: a

      do  i=3,jc
        a=amt(3,i-2)
        do  j=3,6
          amt(j,i-2)=amt(j,i-2)/a
        end do
        a=-amt(2,i-1)
        amt(3,i-1)=amt(3,i-1)+amt(4,i-2)*a
        amt(4,i-1)=amt(4,i-1)+amt(5,i-2)*a
        amt(6,i-1)=amt(6,i-1)+amt(6,i-2)*a
        a=-amt(1,i)
        amt(2,i)=amt(2,i)+amt(4,i-2)*a
        amt(3,i)=amt(3,i)+amt(5,i-2)*a
        amt(6,i)=amt(6,i)+amt(6,i-2)*a
      end do
      a=amt(3,jc-1)
      do  i=3,6
        amt(i,jc-1)=amt(i,jc-1)/a
      end do
      a=-amt(2,jc)
      amt(3,jc)=amt(3,jc)+a*amt(4,jc-1)
      amt(6,jc)=(amt(6,jc)+a*amt(6,jc-1))/amt(3,jc)
      amt(6,jc-1)=amt(6,jc-1)-amt(4,jc-1)*amt(6,jc)
      do  i=jc-2,1,-1
        amt(6,i)=amt(6,i)-amt(4,i)*amt(6,i+1)-amt(5,i)*amt(6,i+2)          
      end do

      return
      end