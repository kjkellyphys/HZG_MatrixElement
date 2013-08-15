      subroutine chkC5_45(j,i1,i2,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      implicit none
      include 'constants.f'  
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,j1,j2,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),Shat3s(np,z2max,-2:0),diff
      double complex res
       
      do ep=-2,0
      do n=1,np
      Shat3s(n,z2(i1,i2),ep)=Shat3(n,z2(i1,i2),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czzi(i2),ep)
     .        +delta(n,i2)*Cv(N0+czzi(i1),ep))
      enddo 
      res=         
     . -(
     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)
     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)
     . -DetGr*Cv(ciii(z3(j,i1,i2))+N0,ep))/Xtwiddle0(j)
      diff=
     .  +Xtwiddle0(j)*Cv(cii(z2(i1,i2))+N0,ep)
     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)
     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)
     . -DetGr*Cv(ciii(z3(j,i1,i2))+N0,ep)
      if (abs(diff) .gt. weenumber) write(6,*) 'chkC5_45',j,i1,i2,diff
     , ,res
c      write(6,*) 'chkC5_45',j,i1,i2,diff,res
      enddo

      return
      end
