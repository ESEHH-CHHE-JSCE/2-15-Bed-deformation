      module allo1
      real(8),dimension(:),allocatable :: B1,h1,Zb1,Se1,Ts1,Qs1
     &,Zb2,Qup,WLd
      end module
      use allo1
      OPEN(9901,FILE='Results.dat',STATUS='UNKNOWN')
      OPEN(9902,FILE='Conditions.dat',STATUS='UNKNOWN')
      OPEN(9903,FILE='Hydrograph.dat',STATUS='UNKNOWN')
      OPEN(9904,FILE='WLevel.dat',STATUS='UNKNOWN')
      read(9902,*)Nx
      read(9902,*)Nqw
      read(9902,*)g
      read(9902,*)Dt
      read(9902,*)Dx
      read(9902,*)Bu
      read(9902,*)Bm
      read(9902,*)Bd
      read(9902,*)Seu
      read(9902,*)Sem
      read(9902,*)Sed
      read(9902,*)Nxu
      read(9902,*)Nxm
      read(9902,*)Dm
      read(9902,*)Rs
      read(9902,*)Porosity
      read(9902,*)Cm
      read(9902,*)Tend
      read(9902,*)Tqup
      read(9902,*)LSTEP
      read(9902,*)TdataS
      read(9902,*)TscreenS
      Sws=Rs-1.
      T=0.
      Tdata=0.
      Tscreen=0.

      allocate (B1(Nx),h1(Nx),Zb1(Nx),Se1(Nx),Ts1(Nx),Qs1(Nx)
     &,Zb2(Nx),Qup(Nqw),WLd(Nqw))

      do 1099 I=1,Nxu
      B1(I)=Bu
 1099 continue

      do 1100 I=Nxu+1,Nxm
      B1(I)=Bm
 1100 continue

      do 1101 I=Nxm+1,Nx
      B1(I)=Bd
 1101 continue

      do 1105 I=1,Nx
      if(I.le.Nxu)then
      Se1(I)=Seu
      elseif(I.le.Nxm)then
      Se1(I)=Sem
      else
      Se1(I)=Sed
      endif
 1105 continue

*---- Hydrograph ----*
      READ(9903,*)
      do L=1, LSTEP
      READ(9903,*)Dmmy1,Qup(L)
      enddo
      Q=Qup(1)

*---- WElevation ----*
      READ(9904,*)
      do L=1, LSTEP
      READ(9904,*)Dmmy1,WLd(L)
      enddo
      WLdown=WLd(1)

      h1(Nx)=WLdown
      do 1102 I=1,Nx-1
      h1(I)=((cm*Q)/(B1(I)*Se1(I)**0.5))**(3./5.)
 1102 continue

      Zb1(Nx)=0.
      do 1104 I=Nx-1,1,-1
      Zb1(I)=Zb1(I+1)+Se1(I)*Dx
 1104 continue

      do 1109 I=1,Nx
      Zb2(I)=Zb1(I)
 1109 continue

 1113 Continue

      I=T/Tqup+1
      Q=(Qup(I+1)-Qup(I))/Tqup*(T-Tqup*(I-1))+Qup(I)
      WLdown=(WLd(I+1)-WLd(I))/Tqup*(T-Tqup*(I-1))+WLd(I)
      h1(Nx)=WLdown-Zb1(Nx)

      do 1110 I=Nx-1,1,-1

 1112 Down=Zb1(I+1)+h1(I+1)+Q**2./(2.*g*B1(I+1)**2.*h1(I+1)**2.)
     &+(Q**2.*cm**2.*Dx)/(2.*B1(I+1)**2.*h1(I+1)**(10./3.))
      Up=Zb1(I)+h1(I)+Q**2./(2.*g*B1(I)**2.*h1(I)**2.)
     &-(Q**2.*cm**2.*Dx)/(2.*B1(I)**2.*h1(I)**(10./3.))

      if(Up.lt.Down)then
      if(J1.eq.2)go to 1111
      h1(I)=h1(I)+0.001
      J1=1
      else
      if(J1.eq.1)go to 1111
      h1(I)=h1(I)-0.001
      J1=2
      endif
      go to 1112

 1111 J1=0

 1110 continue

      if(T.le.Dt*0.5)go to 2001

      do 1106 I=1,Nx
      if(I.le.Nx-1)then
      Se2=(Q**2.*cm**2.)/(2.*B1(I+1)**2.*h1(I+1)**(10./3.))
     &+(Q**2.*cm**2.)/(2.*B1(I)**2.*h1(I)**(10./3.))
      else
      Se2=(Q**2.*cm**2.)/(B1(I)**2.*h1(I)**(10./3.))
      endif
      Us=(g*h1(I)*Se2)**0.5
      Ts1(I)=Us**2./(Sws*g*Dm)
 1106 continue

      Tsc=0.05
      do 1107 I=1,Nx
      U_flux=Q/(B1(I)*h1(I))
      Rough=Dm*(1.+2.*Ts1(I))
      If(h1(I).gt.Rough)then
      Use=U_flux/(6.+1./0.4*log(h1(I)/Rough))
      Else
      Use=U_flux/6.
      Endif
      Tse=Use**2./(Sws*g*Dm)
      if(Tsc.le.Ts1(I))then
      Qs1(I)=17.*(Sws*g*Dm**3.)**0.5*Tse**1.5*(1.-Tsc/Ts1(I))
     &*(1.-(Tsc/Ts1(I))**0.5)*B1(I)
      else
      Qs1(I)=0.
      endif
 1107 continue

      do 1108 I=2,Nx
      Zb1(I)=Zb1(I)-Dt/(1.-Porosity)*(Qs1(I)-Qs1(I-1))/(Dx*B1(I))
 1108 continue
      Zb1(1)=Zb1(2)+Se1(1)*Dx

 2001 continue

      Tdata1=TdataS-Dt*0.5
      if((Tdata.ge.Tdata1).or.(T.le.Dt*0.5))then
      write(9901,*)'T=',T
      write(9901,1001)'I','B1(m)','h1(m)','Z1(m)','Zb1(m)','Se1','Ts1'
     &,'Qs1(m3/s)','Zb2(m)'
      do 1103 I=1,Nx
      write(9901,1002)I,B1(I),h1(I),Zb1(I)+h1(I),Zb1(I),Se1(I),Ts1(I)
     &,Qs1(I),Zb2(I)
 1103 continue
      Tdata=0.
      endif

      Tscreen1=TscreenS-Dt*0.5
      if((Tscreen.ge.Tscreen1).or.(T.le.0.001))then
      write(*,*)'T=',T,'Q=',Q
      Tscreen=0.
      endif

      T=T+Dt
      Tdata=Tdata+Dt
      Tscreen=Tscreen+Dt

      if(T.le.Tend)go to 1113

 1001 format(A5,8A9)
 1002 format(I5,8F9.4)

      deallocate (B1,h1,Zb1,Se1,Ts1,Qs1
     &,Zb2,Qup,WLd)

      stop
      end
