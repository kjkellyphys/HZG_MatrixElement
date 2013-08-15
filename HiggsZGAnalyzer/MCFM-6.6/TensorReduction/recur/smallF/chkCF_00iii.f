      subroutine chkCF_00iii(i1,i2,i3,f,Gr,Shat5,N0)
C---  Expression for rearrangement of extension of Eq. 5.70
C---  Checks C00iii
C---  Small terms of order f(i)*Cijkl,Gr(i,j)*Cijklm
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=2)
      double precision f(np),Gr(np,np),den
      double complex Shat5(np,z4max,-2:0),diff
       
      do ep=-2,0
      if     ((i1 .eq. i2) .and. (i1 .eq. i3)) then
        den=8d0
	k=i1
      elseif (i1 .eq. i2) then
        den=6d0
	k=i1
      elseif (i1 .eq. i3) then
        den=6d0
	k=i1
      elseif (i2 .eq. i3) then
        den=6d0
	k=i2
      else
        den=4d0
	k=i1
      endif      
      diff=Cv(czziii(z3(i1,i2,i3))+N0,ep)*den-
     . (Shat5(k,z4(k,i1,i2,i3),ep)
     . -f(k)*Cv(ciiii(z4(k,i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Cv(ciiiii(z5(1,k,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiii(z5(2,k,i1,i2,i3))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCF_00iii',i1,i2,i3,diff
     
      enddo
      
      return
      end
