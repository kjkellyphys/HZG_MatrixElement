      subroutine chkC_0000i(k,l,i1,DetGr,f,Gtwiddle,Gtt, 
     . Shat4zz,Shat5zzzz,S0000i,Shat5zz,N0) 
      implicit none 
      include 'constants.f'   
      include 'Cnames.f'   
      include 'Cv.f'   
      include 'Carraydef.f'   
      include 'Carrays.f'   
      include 'weenumber.f'   
      integer ep,N0,k,l,n,m,i1,np 
      parameter(np=2) 
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np), 
     . f(np) 
      double complex Shat4zz(np,z1max,-2:0),Shat5zzzz(np,-2:0), 
     . S0000i(np,-2:0),Shat5zz(np,z2max,-2:0),bit,pole,diff
        
      do ep=-2,0
      bit=czip 
      do n=1,np 
      do m=1,np 
      bit=bit 
     . +Gtt(k,n,l,m)*(f(n)*Shat4zz(m,i1,ep) 
     . +2*(delta(n,i1)*Shat5zzzz(m,ep)) 
     . -f(n)*f(m)*Cv(czzi(i1)+N0,ep) 
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Cv(cc0000+N0,ep))
      enddo 
      enddo 
      pole=czip 
      if (ep .gt. -2)  
     . pole=-4*Gtwiddle(k,l)*Cv(czzzzi(i1)+N0,ep-1) 
     
      diff=14d0*Gtwiddle(k,l)*Cv(czzzzi(i1)+N0,ep)
     . +pole 
     . +DetGr*Cv(czziii(z3(k,l,i1))+N0,ep) 
     . -Gtwiddle(k,l)*S0000i(i1,ep) 
     . -Gtwiddle(1,l)*Shat5zz(1,z2(k,i1),ep) 
     . -Gtwiddle(2,l)*Shat5zz(2,z2(k,i1),ep) 
     . +Gtwiddle(k,l) 
     . *(Shat5zz(1,z2(1,i1),ep) 
     .  +Shat5zz(2,z2(2,i1),ep)) 
     . +bit
 
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))  
     . write(6,*) 'chkG_0000i',k,l,i1,diff 
     
      enddo 
      
      return
      end
