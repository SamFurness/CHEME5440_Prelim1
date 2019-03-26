function X=compute(LX1,LX2,LX3,eX,RXT,Gj,KX,tauX,kdX,mu,WI1,W11,W12,W13,W22,W23,W33,eL,RLT,KL,tauL,kdL)

%Begin by defining all parameters given:
%TRANSCRIPTION
LX1=LX1;
LX2=LX2;
LX3=LX3;
LX=1000;
eX=eX;
kEX1=60*eX/LX1;%hr^-1
kEX2=60*eX/LX2;
kEX3=60*eX/LX3;
GDW=(10^9)/((6.022*(10^23))*(0.3*2.8*(10^-13)));%Converts molecules/cell into
    %nmol/gDW
RXT=RXT*GDW;
Gj=Gj*GDW;
KX=KX;%nmol/gDW
tauX=tauX;
kdX=kdX;%min^-1
mu=mu;%min^-1
%Now defining the weights and parameters for a  Moon/Voigt formulation:
WI1=WI1;
W11=W11;
W12=W12;
W13=W13;
W22=W22;
W23=W23;
W33=W33;

%TRANSLATION
LL1=LX1/3;
LL2=LX2/3;
LL3=LX3/3;
LL=333;
eL=eL;
kEL1=60*eL/LL1;
kEL2=60*eL/LL2;
kEL3=60*eL/LL3;%hr^-1
RLT=RLT*GDW;%nmol/gDW
KL=KL;%nmol/gDW
tauL=tauL;
kdL=kdL;%min^-1

%------------------------------------------------------

%Now we will define reaction rates for TRANSCRIPTION:
rX1=kEX1*RXT*Gj/(KX*tauX+Gj*(tauX+1));
rX2=kEX2*RXT*Gj/(KX*tauX+Gj*(tauX+1));
rX3=kEX3*RXT*Gj/(KX*tauX+Gj*(tauX+1));

%Set-up differential equations matrices:
A=[-(mu+kdX),0,0,0,0,0;
   0,-(mu+kdX),0,0,0,0;
   0,0,-(mu+kdX),0,0,0;
   0,0,0,-(mu+kdL),0,0;
   0,0,0,0,-(mu+kdL),0;
   0,0,0,0,0,-(mu+kdL)];
S=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0;
   0,0,0,1,0,0;
   0,0,0,0,1,0;
   0,0,0,0,0,1];

%Define a time vector:
T = 1;%in hours
tsim=linspace(0,300,((300/T)+1));
x=zeros(6,length(tsim));%preallocate x vector

%Set initial conditions
x(:,1)=[0;0;0;0;0;0];%Initially no mRNA or protein
m1_0=x(1,1);
m2_0=x(2,1);
m3_0=x(3,1);
p1_0=x(4,1);
p2_0=x(5,1);
%p3_0=x(6,1);
I_0=0;
%Initial values for r. Uses function bind_fxn for fI values
u1_0=(W11+(WI1*bind_fxn(I_0,1,1)))/(1+W11+(WI1*bind_fxn(I_0,1,1)));
u2_0=(W22+(W12*bind_fxn(p1_0,1,2)))/(1+W22+(W12*bind_fxn(p1_0,1,2)));
u3_0=(W33+(W13*bind_fxn(p1_0,1,3)))/(1+W33+(W13*bind_fxn(p1_0,1,3))+(W23*bind_fxn(p2_0,2,3)));
rL1_0=kEL1*RLT*(m1_0/(KL*tauL+m1_0*(tauL+1)));
rL2_0=kEL2*RLT*(m2_0/(KL*tauL+m2_0*(tauL+1)));
rL3_0=kEL3*RLT*(m3_0/(KL*tauL+m3_0*(tauL+1)));
TX1_0=rX1*u1_0;
TX2_0=rX2*u2_0;
TX3_0=rX3*u3_0;
TL4_0=rL1_0;%assume at kinetic limit
TL5_0=rL2_0;
TL6_0=rL3_0;

%Set up r vector initial conditions
r0=[TX1_0;TX2_0;TX3_0;TL4_0;TL5_0;TL6_0];
r=zeros(6,length(tsim));
r(1,1)=r0(1);
r(2,1)=r0(2);
r(3,1)=r0(3);
r(4,1)=r0(4);
r(5,1)=r0(5);
r(6,1)=r0(6);

%Performing the Discretization
A_hat=expm(A*T);
A_inverse=inv(A);
A_sub=A_hat-eye(6);
S_sub=A\A_sub;
S_hat=S_sub*S;

%Solve Discretization
for i=1:(length(tsim)-1)
   Ix=0;
   %Get the current r and x vectors
   xk=zeros(6,1);
   rk=zeros(6,1);
   xk(1)=x(1,i);
   xk(2)=x(2,i);
   xk(3)=x(3,i);
   xk(4)=x(4,i);
   xk(5)=x(5,i);
   xk(6)=x(6,i);
  
   rk(1)=r(1,i);
   rk(2)=r(2,i);
   rk(3)=r(3,i);
   rk(4)=r(4,i);
   rk(5)=r(5,i);
   rk(6)=r(6,i);
   
   %Calculate the new entry:
   xkplus1=A_hat*xk+S_hat*rk;
   
   %Update x with new values
   x(1,i+1)=xkplus1(1);
   x(2,i+1)=xkplus1(2);
   x(3,i+1)=xkplus1(3);
   x(4,i+1)=xkplus1(4);
   x(5,i+1)=xkplus1(5);
   x(6,i+1)=xkplus1(6);
   
   %Update r 
   m1=xkplus1(1);
   m2=xkplus1(2);
   m3=xkplus1(3);
   p1=xkplus1(4);
   p2=xkplus1(5);
   p3=xkplus1(6);
   u1=(W11+(WI1*bind_fxn(Ix,1,1)))/(1+W11+(WI1*bind_fxn(Ix,1,1)));
   u2=(W22+(W12*bind_fxn(p1,1,2)))/(1+W22+(W12*bind_fxn(p1,1,2)));
   u3=(W33+(W13*bind_fxn(p1,1,3)))/(1+W33+(W13*bind_fxn(p1,1,3))+(W23*bind_fxn(p2,2,3)));
   rL1=kEL1*RLT*(m1/(KL*tauL+m1*(tauL+1)));
   rL2=kEL2*RLT*(m2/(KL*tauL+m2*(tauL+1)));
   rL3=kEL3*RLT*(m3/(KL*tauL+m3*(tauL+1)));
   TX1=rX1*u1;
   TX2=rX2*u2;
   TX3=rX3*u3;
   TL4=rL1;
   TL5=rL2;
   TL6=rL3;
%now r
   r(1,i+1)=TX1;
   r(2,i+1)=TX2;
   r(3,i+1)=TX3;
   r(4,i+1)=TL4;
   r(5,i+1)=TL5;
   r(6,i+1)=TL6;
end

%Once at steady state, add more time to add the inducer
%New time simulated over
tSim=linspace(tsim(end),tsim(end)+300,((300)/T+1));
%Set initial conditions
xI=zeros(6,length(tSim));
xI(:,1)=x(:,length(tsim));%Initially values from no inducer
m1_0=xI(1,1);
m2_0=xI(2,1);
m3_0=xI(3,1);
p1_0=xI(4,1);
p2_0=xI(5,1);
p3_0=xI(6,1);
I_0=0;
%Initial values for r. Uses function bind_fxn for fI values
u1_0=(W11+(WI1*bind_fxn(I_0,1,1)))/(1+W11+(WI1*bind_fxn(I_0,1,1)));
u2_0=(W22+(W12*bind_fxn(p1_0,1,2)))/(1+W22+(W12*bind_fxn(p1_0,1,2)));
u3_0=(W33+(W13*bind_fxn(p1_0,1,3)))/(1+W33+(W13*bind_fxn(p1_0,1,3))+(W23*bind_fxn(p2_0,2,3)));
rL1_0=kEL1*RLT*(m1_0/(KL*tauL+m1_0*(tauL+1)));
rL2_0=kEL2*RLT*(m2_0/(KL*tauL+m2_0*(tauL+1)));
rL3_0=kEL3*RLT*(m3_0/(KL*tauL+m3_0*(tauL+1)));
TX1_0=rX1*u1_0;
TX2_0=rX2*u2_0;
TX3_0=rX3*u3_0;
TL4_0=rL1_0;%assume at kinetic limit
TL5_0=rL2_0;
TL6_0=rL3_0;
%Set up r vector initial conditions
r0=[TX1_0;TX2_0;TX3_0;TL4_0;TL5_0;TL6_0];
rI=zeros(6,length(tSim));
rI(1,1)=r0(1);
rI(2,1)=r0(2);
rI(3,1)=r0(3);
rI(4,1)=r0(4);
rI(5,1)=r0(5);
rI(6,1)=r0(6);

for i=1:(length(tSim)-1)
   Ix=10;
   %Get the current r and x vectors
   xk=zeros(6,1);
   rk=zeros(6,1);
   xk(1)=xI(1,i);
   xk(2)=xI(2,i);
   xk(3)=xI(3,i);
   xk(4)=xI(4,i);
   xk(5)=xI(5,i);
   xk(6)=xI(6,i);
  
   rk(1)=rI(1,i);
   rk(2)=rI(2,i);
   rk(3)=rI(3,i);
   rk(4)=rI(4,i);
   rk(5)=rI(5,i);
   rk(6)=rI(6,i);
   
   %Calculate the new entry:
   xkplus1=A_hat*xk+S_hat*rk;
   
   %Update x with new values
   xI(1,i+1)=xkplus1(1);
   xI(2,i+1)=xkplus1(2);
   xI(3,i+1)=xkplus1(3);
   xI(4,i+1)=xkplus1(4);
   xI(5,i+1)=xkplus1(5);
   xI(6,i+1)=xkplus1(6);
   
   %Update r 
   m1=xkplus1(1);
   m2=xkplus1(2);
   m3=xkplus1(3);
   p1=xkplus1(4);
   p2=xkplus1(5);
   p3=xkplus1(6);
   u1=(W11+(WI1*bind_fxn(Ix,1,1)))/(1+W11+(WI1*bind_fxn(Ix,1,1)));
   u2=(W22+(W12*bind_fxn(p1,1,2)))/(1+W22+(W12*bind_fxn(p1,1,2)));
   u3=(W33+(W13*bind_fxn(p1,1,3)))/(1+W33+(W13*bind_fxn(p1,1,3))+(W23*bind_fxn(p2,2,3)));
   rL1=kEL1*RLT*(m1/(KL*tauL+m1*(tauL+1)));
   rL2=kEL2*RLT*(m2/(KL*tauL+m2*(tauL+1)));
   rL3=kEL3*RLT*(m3/(KL*tauL+m3*(tauL+1)));
   TX1=rX1*u1;
   TX2=rX2*u2;
   TX3=rX3*u3;
   TL4=rL1;
   TL5=rL2;
   TL6=rL3;
%now r
   rI(1,i+1)=TX1;
   rI(2,i+1)=TX2;
   rI(3,i+1)=TX3;
   rI(4,i+1)=TL4;
   rI(5,i+1)=TL5;
   rI(6,i+1)=TL6;
end

%Concatenate the values to give the full array of species as output
X=[x xI(:,2:end)];
R=[r rI(:,2:end)];
T=[tsim tSim(2:end)];
    

end

