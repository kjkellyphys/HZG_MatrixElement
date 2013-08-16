      subroutine chkY_00i(k,l,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for DD Eq. 5.56b
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),diff

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     . 2*Gtwiddle(k,l)*Dv(dzzi(i1)+N0,ep)
     . -(-2*Gtwiddle(k,i1)*Dv(dzzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat3(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat3(2,z2(l,i1),ep)
     . +Gtwiddle(k,3)*Shat3(3,z2(l,i1),ep)
     . +Xtwiddle(k,0)*Dv(dii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diii(z3(k,l,i1))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chk5_56b',k,l,i1,diff

      enddo

      return
      end
  


