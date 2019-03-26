%Sam Furness
%CHEME5440 Prelim1 Question 2 Part B,C
%3/22/19

%This code calls the function compute.m that runs the whole simulation at
    %every time step for given parameters. Creates s array that gives the
    %sensitivity of each parameter. This s array is integrated to create a 
    %time averaged sensitivity array contained every paramter
    %and this undergoes SVD to yield the U array.
%--------------------------------------------------
%Begin my calculating sensitivity for each phase for each parameter
    %with central difference:
    
%Length Gene 1
LX1_low=1200*0.95;
LX1=1200;
LX1_high=1200*1.05;
X_LX1low=compute(LX1_low,2400,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX1=compute(LX1,2400,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX1high=compute(LX1_high,2400,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient 
%Phase 1
for i=280:300
    for j=1:6
    s_LX1(j,i)=(LX1/X_LX1(j,i))*(X_LX1high(j,i)-(X_LX1low(j,i)))/(0.1*LX1);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_LX1(j,i)=(LX1/X_LX1(j,i))*(X_LX1high(j,i)-(X_LX1low(j,i)))/(0.1*LX1);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_LX1(j,i)=(LX1/X_LX1(j,i))*(X_LX1high(j,i)-(X_LX1low(j,i)))/(0.1*LX1);
    end
end

%Length Gene 2
LX2_low=2400*0.95;
LX2=2400;
LX2_high=2400*1.05;
X_LX2low=compute(1200,LX2_low,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX2=compute(1200,LX2,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX2high=compute(1200,LX2_high,600,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_LX2(j,i)=(LX2/X_LX2(j,i))*(X_LX2high(j,i)-(X_LX2low(j,i)))/(0.1*LX2);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_LX2(j,i)=(LX2/X_LX2(j,i))*(X_LX2high(j,i)-(X_LX2low(j,i)))/(0.1*LX2);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_LX2(j,i)=(LX2/X_LX2(j,i))*(X_LX2high(j,i)-(X_LX2low(j,i)))/(0.1*LX2);
    end
end

%Length Gene 3
LX3=600;
LX3_low=LX3*0.95;
LX3_high=LX3*1.05;
X_LX3low=compute(1200,2400,LX3_low,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX3=compute(1200,2400,LX3,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_LX3high=compute(1200,2400,LX3_high,60,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_LX3(j,i)=(LX3/X_LX3(j,i))*(X_LX3high(j,i)-(X_LX3low(j,i)))/(0.1*LX3);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_LX3(j,i)=(LX3/X_LX3(j,i))*(X_LX3high(j,i)-(X_LX3low(j,i)))/(0.1*LX3);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_LX3(j,i)=(LX3/X_LX3(j,i))*(X_LX3high(j,i)-(X_LX3low(j,i)))/(0.1*LX3);
    end
end

%Elongation rate
eX=60;
eX_low=eX*0.95;
eX_high=eX*1.05;
X_eXlow=compute(1200,2400,600,eX_low,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_eX=compute(1200,2400,600,eX,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_eXhigh=compute(1200,2400,600,eX_high,1150,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_eX(j,i)=(eX/X_eX(j,i))*(X_eXhigh(j,i)-(X_eXlow(j,i)))/(0.1*eX);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_eX(j,i)=(eX/X_eX(j,i))*(X_eXhigh(j,i)-(X_eXlow(j,i)))/(0.1*eX);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_eX(j,i)=(eX/X_eX(j,i))*(X_eXhigh(j,i)-(X_eXlow(j,i)))/(0.1*eX);
    end
end

%RNAP Concentration
RXT=1150;
RXT_low=RXT*0.95;
RXT_high=RXT*1.05;
X_RXTlow=compute(1200,2400,600,60,RXT_low,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_RXT=compute(1200,2400,600,60,RXT,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_RXThigh=compute(1200,2400,600,60,RXT_high,200,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_RXT(j,i)=(RXT/X_RXT(j,i))*(X_RXThigh(j,i)-(X_RXTlow(j,i)))/(0.1*RXT);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_RXT(j,i)=(RXT/X_RXT(j,i))*(X_RXThigh(j,i)-(X_RXTlow(j,i)))/(0.1*RXT);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_RXT(j,i)=(RXT/X_RXT(j,i))*(X_RXThigh(j,i)-(X_RXTlow(j,i)))/(0.1*RXT);
    end
end

%Gene Concentration
Gj=200;
Gj_low=Gj*0.95;
Gj_high=Gj*1.05;
X_Gjlow=compute(1200,2400,600,60,1150,Gj_low,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_Gj=compute(1200,2400,600,60,1150,Gj,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_Gjhigh=compute(1200,2400,600,60,1150,Gj_high,.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_Gj(j,i)=(Gj/X_Gj(j,i))*(X_Gjhigh(j,i)-(X_Gjlow(j,i)))/(0.1*Gj);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_Gj(j,i)=(Gj/X_Gj(j,i))*(X_Gjhigh(j,i)-(X_Gjlow(j,i)))/(0.1*Gj);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_Gj(j,i)=(Gj/X_Gj(j,i))*(X_Gjhigh(j,i)-(X_Gjlow(j,i)))/(0.1*Gj);
    end
end

%Transcription Saturation Constant
KX=0.24;
KX_low=KX*0.95;
KX_high=KX*1.05;
X_KXlow=compute(1200,2400,600,60,1150,200,KX_low,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_KX=compute(1200,2400,600,60,1150,200,KX,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_KXhigh=compute(1200,2400,600,60,1150,200,KX_high,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_KX(j,i)=(KX/X_KX(j,i))*(X_KXhigh(j,i)-(X_KXlow(j,i)))/(0.1*KX);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_KX(j,i)=(KX/X_KX(j,i))*(X_KXhigh(j,i)-(X_KXlow(j,i)))/(0.1*KX);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_KX(j,i)=(KX/X_KX(j,i))*(X_KXhigh(j,i)-(X_KXlow(j,i)))/(0.1*KX);
    end
end


%Transcription time constant
tauX=2.7;
tauX_low=tauX*0.95;
tauX_high=tauX*1.05;
X_tauXlow=compute(1200,2400,600,60,1150,200,0.24,tauX_low,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_tauX=compute(1200,2400,600,60,1150,200,0.24,tauX,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_tauXhigh=compute(1200,2400,600,60,1150,200,0.24,tauX_high,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_tauX(j,i)=(tauX/X_tauX(j,i))*(X_tauXhigh(j,i)-(X_tauXlow(j,i)))/(0.1*tauX);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_tauX(j,i)=(tauX/X_tauX(j,i))*(X_tauXhigh(j,i)-(X_tauXlow(j,i)))/(0.1*tauX);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_tauX(j,i)=(tauX/X_tauX(j,i))*(X_tauXhigh(j,i)-(X_tauXlow(j,i)))/(0.1*tauX);
    end
end

%mRNA degradation rate
kdX=log(2)/2.1;
kdX_low=kdX*0.95;
kdX_high=kdX*1.05;
X_kdXlow=compute(1200,2400,600,60,1150,200,0.24,2.7,kdX_low,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_kdX=compute(1200,2400,600,60,1150,200,0.24,2.7,kdX,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_kdXhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,kdX_high,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_kdX(j,i)=(kdX/X_kdX(j,i))*(X_kdXhigh(j,i)-(X_kdXlow(j,i)))/(0.1*kdX);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_kdX(j,i)=(kdX/X_kdX(j,i))*(X_kdXhigh(j,i)-(X_kdXlow(j,i)))/(0.1*kdX);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_kdX(j,i)=(kdX/X_kdX(j,i))*(X_kdXhigh(j,i)-(X_kdXlow(j,i)))/(0.1*kdX);
    end
end

%Growth rate
mu=log(2)/40;
mu_low=mu*0.95;
mu_high=mu*1.05;
X_mulow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,mu_low,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_mu=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,mu,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_muhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,mu_high,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_mu(j,i)=(mu/X_mu(j,i))*(X_muhigh(j,i)-(X_mulow(j,i)))/(0.1*mu);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_mu(j,i)=(mu/X_mu(j,i))*(X_muhigh(j,i)-(X_mulow(j,i)))/(0.1*mu);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_mu(j,i)=(mu/X_mu(j,i))*(X_muhigh(j,i)-(X_mulow(j,i)))/(0.1*mu);
    end
end

%WI1
WI1=100;
WI1_low=WI1*0.95;
WI1_high=WI1*1.05;
X_WI1low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,WI1_low,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_WI1=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,WI1,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_WI1high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,WI1_high,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_WI1(j,i)=(WI1/X_WI1(j,i))*(X_WI1high(j,i)-(X_WI1low(j,i)))/(0.1*WI1);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_WI1(j,i)=(WI1/X_WI1(j,i))*(X_WI1high(j,i)-(X_WI1low(j,i)))/(0.1*WI1);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_WI1(j,i)=(WI1/X_WI1(j,i))*(X_WI1high(j,i)-(X_WI1low(j,i)))/(0.1*WI1);
    end
end

%W11
W11=1e-10;
W11_low=W11*0.95;
W11_high=W11*1.05;
X_W11low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,W11_low,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W11=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,W11,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W11high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,W11_high,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W11(j,i)=(W11/X_W11(j,i))*(X_W11high(j,i)-(X_W11low(j,i)))/(0.1*W11);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W11(j,i)=(W11/X_W11(j,i))*(X_W11high(j,i)-(X_W11low(j,i)))/(0.1*W11);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W11(j,i)=(W11/X_W11(j,i))*(X_W11high(j,i)-(X_W11low(j,i)))/(0.1*W11);
    end
end

%W12
W12=100;
W12_low=W12*0.95;
W12_high=W12*1.05;
X_W12low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,W12_low,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W12=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,W12,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W12high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,W12_high,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W12(j,i)=(W12/X_W12(j,i))*(X_W12high(j,i)-(X_W12low(j,i)))/(0.1*W12);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W12(j,i)=(W12/X_W12(j,i))*(X_W12high(j,i)-(X_W12low(j,i)))/(0.1*W12);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W12(j,i)=(W12/X_W12(j,i))*(X_W12high(j,i)-(X_W12low(j,i)))/(0.1*W12);
    end
end

%W13
W13=0.5;
W13_low=W13*0.95;
W13_high=W13*1.05;
X_W13low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,W13_low,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W13=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,W13,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W13high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,W13_high,1e-10,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W13(j,i)=(W13/X_W13(j,i))*(X_W13high(j,i)-(X_W13low(j,i)))/(0.1*W13);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W13(j,i)=(W13/X_W13(j,i))*(X_W13high(j,i)-(X_W13low(j,i)))/(0.1*W13);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W13(j,i)=(W13/X_W13(j,i))*(X_W13high(j,i)-(X_W13low(j,i)))/(0.1*W13);
    end
end

%W22
W22=1e-10;
W22_low=W22*0.95;
W22_high=W22*1.05;
X_W22low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,W22_low,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W22=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,W22,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W22high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,W22_high,500000,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W22(j,i)=(W22/X_W22(j,i))*(X_W22high(j,i)-(X_W22low(j,i)))/(0.1*W22);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W22(j,i)=(W22/X_W22(j,i))*(X_W22high(j,i)-(X_W22low(j,i)))/(0.1*W22);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W22(j,i)=(W22/X_W22(j,i))*(X_W22high(j,i)-(X_W22low(j,i)))/(0.1*W22);
    end
end

%W23
W23=500000;
W23_low=W23*0.95;
W23_high=W23*1.05;
X_W23low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,W23_low,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W23=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,W23,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W23high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,W23_high,1e-10,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W23(j,i)=(W23/X_W23(j,i))*(X_W23high(j,i)-(X_W23low(j,i)))/(0.1*W23);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W23(j,i)=(W23/X_W23(j,i))*(X_W23high(j,i)-(X_W23low(j,i)))/(0.1*W23);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W23(j,i)=(W23/X_W23(j,i))*(X_W23high(j,i)-(X_W23low(j,i)))/(0.1*W23);
    end
end

%W33
W33=1e-10;
W33_low=W33*0.95;
W33_high=W33*1.05;
X_W33low=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,W33_low,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W33=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,W33,16.5,45000,454.64,0.8,log(20)/(24*60));
X_W33high=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,W33_high,16.5,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_W33(j,i)=(W33/X_W33(j,i))*(X_W33high(j,i)-(X_W33low(j,i)))/(0.1*W33);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_W33(j,i)=(W33/X_W33(j,i))*(X_W33high(j,i)-(X_W33low(j,i)))/(0.1*W33);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_W33(j,i)=(W33/X_W33(j,i))*(X_W33high(j,i)-(X_W33low(j,i)))/(0.1*W33);
    end
end

%eL
eL=16.5;
eL_low=eL*0.95;
eL_high=eL*1.05;
X_eLlow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,eL_low,45000,454.64,0.8,log(20)/(24*60));
X_eL=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,eL,45000,454.64,0.8,log(20)/(24*60));
X_eLhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,eL_high,45000,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_eL(j,i)=(eL/X_eL(j,i))*(X_eLhigh(j,i)-(X_eLlow(j,i)))/(0.1*eL);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_eL(j,i)=(eL/X_eL(j,i))*(X_eLhigh(j,i)-(X_eLlow(j,i)))/(0.1*eL);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_eL(j,i)=(eL/X_eL(j,i))*(X_eLhigh(j,i)-(X_eLlow(j,i)))/(0.1*eL);
    end
end

%RLT
RLT=45000;
RLT_low=RLT*0.95;
RLT_high=RLT*1.05;
X_RLTlow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,RLT_low,454.64,0.8,log(20)/(24*60));
X_RLT=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,RLT,454.64,0.8,log(20)/(24*60));
X_RLThigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,RLT_high,454.64,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_RLT(j,i)=(RLT/X_RLT(j,i))*(X_RLThigh(j,i)-(X_RLTlow(j,i)))/(0.1*RLT);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_RLT(j,i)=(RLT/X_RLT(j,i))*(X_RLThigh(j,i)-(X_RLTlow(j,i)))/(0.1*RLT);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_RLT(j,i)=(RLT/X_RLT(j,i))*(X_RLThigh(j,i)-(X_RLTlow(j,i)))/(0.1*RLT);
    end
end

%KL
KL=454.64;
KL_low=KL*0.95;
KL_high=KL*1.05;
X_KLlow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,KL_low,0.8,log(20)/(24*60));
X_KL=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,KL,0.8,log(20)/(24*60));
X_KLhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,KL_high,0.8,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_KL(j,i)=(KL/X_KL(j,i))*(X_KLhigh(j,i)-(X_KLlow(j,i)))/(0.1*KL);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_KL(j,i)=(KL/X_KL(j,i))*(X_KLhigh(j,i)-(X_KLlow(j,i)))/(0.1*KL);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_KL(j,i)=(KL/X_KL(j,i))*(X_KLhigh(j,i)-(X_KLlow(j,i)))/(0.1*KL);
    end
end

%tauL
tauL=0.8;
tauL_low=tauL*0.95;
tauL_high=tauL*1.05;
X_tauLlow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,tauL_low,log(20)/(24*60));
X_tauL=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,tauL,log(20)/(24*60));
X_tauLhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,tauL_high,log(20)/(24*60));
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_tauL(j,i)=(tauL/X_tauL(j,i))*(X_tauLhigh(j,i)-(X_tauLlow(j,i)))/(0.1*tauL);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_tauL(j,i)=(tauL/X_tauL(j,i))*(X_tauLhigh(j,i)-(X_tauLlow(j,i)))/(0.1*tauL);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_tauL(j,i)=(tauL/X_tauL(j,i))*(X_tauLhigh(j,i)-(X_tauLlow(j,i)))/(0.1*tauL);
    end
end

%kdL
kdL=log(2)/(24*60);
kdL_low=kdL*0.95;
kdL_high=kdL*1.05;
X_kdLlow=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,kdL_low);
X_kdL=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,kdL);
X_kdLhigh=compute(1200,2400,600,60,1150,200,0.24,2.7,log(2)/2.1,log(2)/40,100,1e-10,100,0.5,1e-10,500000,1e-10,16.5,45000,454.64,0.8,kdL_high);
T=linspace(0,600,601);
%Compute the sensitivity coefficient
%Phase 1
for i=280:300
    for j=1:6
    s_kdL(j,i)=(kdL/X_kdL(j,i))*(X_kdLhigh(j,i)-(X_kdLlow(j,i)))/(0.1*kdL);
    end
end
%Phase 2 Early
for i=300:320
    for j=1:6
    s_kdL(j,i)=(kdL/X_kdL(j,i))*(X_kdLhigh(j,i)-(X_kdLlow(j,i)))/(0.1*kdL);
    end
end
%Phase 2 Late
for i=580:600
    for j=1:6
    s_kdL(j,i)=(kdL/X_kdL(j,i))*(X_kdLhigh(j,i)-(X_kdLlow(j,i)))/(0.1*kdL);
    end
end

%----------------------------------------------------------

%Now to calculate the time averaged sensitivity array we will integrate
    %the sensitivity coefficient arrays using the trapezoidal rule and
    %divide by the time windows (19 minutes). Note the time step is
    %1 minute so the width of each trapezoid is 1.
    
%____PHASE1______
%Calculating time averaged sensitivity vector for LX1
for i=280:299
    for j=1:6
    sum_LX1(j,i) = (s_LX1(j,i)+s_LX1(j,(i+1)))/2;
    end
end
S_LX1(1,1)=sum(sum_LX1(1,280:299))/19;
S_LX1(2,1)=sum(sum_LX1(2,280:299))/19;
S_LX1(3,1)=sum(sum_LX1(3,280:299))/19;
S_LX1(4,1)=sum(sum_LX1(4,280:299))/19;
S_LX1(5,1)=sum(sum_LX1(5,280:299))/19;
S_LX1(6,1)=sum(sum_LX1(6,280:299))/19;

%LX2
for i=280:299
    for j=1:6
    sum_LX2(j,i) = (s_LX2(j,i)+s_LX2(j,(i+1)))/2;
    end
end
S_LX2(1,1)=sum(sum_LX2(1,280:299))/19;
S_LX2(2,1)=sum(sum_LX2(2,280:299))/19;
S_LX2(3,1)=sum(sum_LX2(3,280:299))/19;
S_LX2(4,1)=sum(sum_LX2(4,280:299))/19;
S_LX2(5,1)=sum(sum_LX2(5,280:299))/19;
S_LX2(6,1)=sum(sum_LX2(6,280:299))/19;

%LX3
for i=280:299
    for j=1:6
    sum_LX3(j,i) = (s_LX3(j,i)+s_LX3(j,(i+1)))/2;
    end
end
S_LX3(1,1)=sum(sum_LX3(1,280:299))/19;
S_LX3(2,1)=sum(sum_LX3(2,280:299))/19;
S_LX3(3,1)=sum(sum_LX3(3,280:299))/19;
S_LX3(4,1)=sum(sum_LX3(4,280:299))/19;
S_LX3(5,1)=sum(sum_LX3(5,280:299))/19;
S_LX3(6,1)=sum(sum_LX3(6,280:299))/19;

%eX
for i=280:299
    for j=1:6
    sum_eX(j,i) = (s_eX(j,i)+s_eX(j,(i+1)))/2;
    end
end
S_eX(1,1)=sum(sum_eX(1,280:299))/19;
S_eX(2,1)=sum(sum_eX(2,280:299))/19;
S_eX(3,1)=sum(sum_eX(3,280:299))/19;
S_eX(4,1)=sum(sum_eX(4,280:299))/19;
S_eX(5,1)=sum(sum_eX(5,280:299))/19;
S_eX(6,1)=sum(sum_eX(6,280:299))/19;

%RXT
for i=280:299
    for j=1:6
    sum_RXT(j,i) = (s_RXT(j,i)+s_RXT(j,(i+1)))/2;
    end
end
S_RXT(1,1)=sum(sum_RXT(1,280:299))/19;
S_RXT(2,1)=sum(sum_RXT(2,280:299))/19;
S_RXT(3,1)=sum(sum_RXT(3,280:299))/19;
S_RXT(4,1)=sum(sum_RXT(4,280:299))/19;
S_RXT(5,1)=sum(sum_RXT(5,280:299))/19;
S_RXT(6,1)=sum(sum_RXT(6,280:299))/19;

%Gj
for i=280:299
    for j=1:6
    sum_Gj(j,i) = (s_Gj(j,i)+s_Gj(j,(i+1)))/2;
    end
end
S_Gj(1,1)=sum(sum_Gj(1,280:299))/19;
S_Gj(2,1)=sum(sum_Gj(2,280:299))/19;
S_Gj(3,1)=sum(sum_Gj(3,280:299))/19;
S_Gj(4,1)=sum(sum_Gj(4,280:299))/19;
S_Gj(5,1)=sum(sum_Gj(5,280:299))/19;
S_Gj(6,1)=sum(sum_Gj(6,280:299))/19;

%KX
for i=280:299
    for j=1:6
    sum_KX(j,i) = (s_KX(j,i)+s_KX(j,(i+1)))/2;
    end
end
S_KX(1,1)=sum(sum_KX(1,280:299))/19;
S_KX(2,1)=sum(sum_KX(2,280:299))/19;
S_KX(3,1)=sum(sum_KX(3,280:299))/19;
S_KX(4,1)=sum(sum_KX(4,280:299))/19;
S_KX(5,1)=sum(sum_KX(5,280:299))/19;
S_KX(6,1)=sum(sum_KX(6,280:299))/19;

%tauX
for i=280:299
    for j=1:6
    sum_tauX(j,i) = (s_tauX(j,i)+s_tauX(j,(i+1)))/2;
    end
end
S_tauX(1,1)=sum(sum_tauX(1,280:299))/19;
S_tauX(2,1)=sum(sum_tauX(2,280:299))/19;
S_tauX(3,1)=sum(sum_tauX(3,280:299))/19;
S_tauX(4,1)=sum(sum_tauX(4,280:299))/19;
S_tauX(5,1)=sum(sum_tauX(5,280:299))/19;
S_tauX(6,1)=sum(sum_tauX(6,280:299))/19;

%kdX
for i=280:299
    for j=1:6
    sum_kdX(j,i) = (s_kdX(j,i)+s_kdX(j,(i+1)))/2;
    end
end
S_kdX(1,1)=sum(sum_kdX(1,280:299))/19;
S_kdX(2,1)=sum(sum_kdX(2,280:299))/19;
S_kdX(3,1)=sum(sum_kdX(3,280:299))/19;
S_kdX(4,1)=sum(sum_kdX(4,280:299))/19;
S_kdX(5,1)=sum(sum_kdX(5,280:299))/19;
S_kdX(6,1)=sum(sum_kdX(6,280:299))/19;

%mu
for i=280:299
    for j=1:6
    sum_mu(j,i) = (s_mu(j,i)+s_mu(j,(i+1)))/2;
    end
end
S_mu(1,1)=sum(sum_mu(1,280:299))/19;
S_mu(2,1)=sum(sum_mu(2,280:299))/19;
S_mu(3,1)=sum(sum_mu(3,280:299))/19;
S_mu(4,1)=sum(sum_mu(4,280:299))/19;
S_mu(5,1)=sum(sum_mu(5,280:299))/19;
S_mu(6,1)=sum(sum_mu(6,280:299))/19;

%WI1
for i=280:299
    for j=1:6
    sum_WI1(j,i) = (s_WI1(j,i)+s_WI1(j,(i+1)))/2;
    end
end
S_WI1(1,1)=sum(sum_WI1(1,280:299))/19;
S_WI1(2,1)=sum(sum_WI1(2,280:299))/19;
S_WI1(3,1)=sum(sum_WI1(3,280:299))/19;
S_WI1(4,1)=sum(sum_WI1(4,280:299))/19;
S_WI1(5,1)=sum(sum_WI1(5,280:299))/19;
S_WI1(6,1)=sum(sum_WI1(6,280:299))/19;

%W11
for i=280:299
    for j=1:6
    sum_W11(j,i) = (s_W11(j,i)+s_W11(j,(i+1)))/2;
    end
end
S_W11(1,1)=sum(sum_W11(1,280:299))/19;
S_W11(2,1)=sum(sum_W11(2,280:299))/19;
S_W11(3,1)=sum(sum_W11(3,280:299))/19;
S_W11(4,1)=sum(sum_W11(4,280:299))/19;
S_W11(5,1)=sum(sum_W11(5,280:299))/19;
S_W11(6,1)=sum(sum_W11(6,280:299))/19;

%W12
for i=280:299
    for j=1:6
    sum_W12(j,i) = (s_W12(j,i)+s_W12(j,(i+1)))/2;
    end
end
S_W12(1,1)=sum(sum_W12(1,280:299))/19;
S_W12(2,1)=sum(sum_W12(2,280:299))/19;
S_W12(3,1)=sum(sum_W12(3,280:299))/19;
S_W12(4,1)=sum(sum_W12(4,280:299))/19;
S_W12(5,1)=sum(sum_W12(5,280:299))/19;
S_W12(6,1)=sum(sum_W12(6,280:299))/19;

%W13
for i=280:299
    for j=1:6
    sum_W13(j,i) = (s_W13(j,i)+s_W13(j,(i+1)))/2;
    end
end
S_W13(1,1)=sum(sum_W13(1,280:299))/19;
S_W13(2,1)=sum(sum_W13(2,280:299))/19;
S_W13(3,1)=sum(sum_W13(3,280:299))/19;
S_W13(4,1)=sum(sum_W13(4,280:299))/19;
S_W13(5,1)=sum(sum_W13(5,280:299))/19;
S_W13(6,1)=sum(sum_W13(6,280:299))/19;

%W22
for i=280:299
    for j=1:6
    sum_W22(j,i) = (s_W22(j,i)+s_W22(j,(i+1)))/2;
    end
end
S_W22(1,1)=sum(sum_W22(1,280:299))/19;
S_W22(2,1)=sum(sum_W22(2,280:299))/19;
S_W22(3,1)=sum(sum_W22(3,280:299))/19;
S_W22(4,1)=sum(sum_W22(4,280:299))/19;
S_W22(5,1)=sum(sum_W22(5,280:299))/19;
S_W22(6,1)=sum(sum_W22(6,280:299))/19;

%W23
for i=280:299
    for j=1:6
    sum_W23(j,i) = (s_W23(j,i)+s_W23(j,(i+1)))/2;
    end
end
S_W23(1,1)=sum(sum_W23(1,280:299))/19;
S_W23(2,1)=sum(sum_W23(2,280:299))/19;
S_W23(3,1)=sum(sum_W23(3,280:299))/19;
S_W23(4,1)=sum(sum_W23(4,280:299))/19;
S_W23(5,1)=sum(sum_W23(5,280:299))/19;
S_W23(6,1)=sum(sum_W23(6,280:299))/19;

%W33
for i=280:299
    for j=1:6
    sum_W33(j,i) = (s_W33(j,i)+s_W33(j,(i+1)))/2;
    end
end
S_W33(1,1)=sum(sum_W33(1,280:299))/19;
S_W33(2,1)=sum(sum_W33(2,280:299))/19;
S_W33(3,1)=sum(sum_W33(3,280:299))/19;
S_W33(4,1)=sum(sum_W33(4,280:299))/19;
S_W33(5,1)=sum(sum_W33(5,280:299))/19;
S_W33(6,1)=sum(sum_W33(6,280:299))/19;

%eL
for i=280:299
    for j=1:6
    sum_eL(j,i) = (s_eL(j,i)+s_eL(j,(i+1)))/2;
    end
end
S_eL(1,1)=sum(sum_eL(1,280:299))/19;
S_eL(2,1)=sum(sum_eL(2,280:299))/19;
S_eL(3,1)=sum(sum_eL(3,280:299))/19;
S_eL(4,1)=sum(sum_eL(4,280:299))/19;
S_eL(5,1)=sum(sum_eL(5,280:299))/19;
S_eL(6,1)=sum(sum_eL(6,280:299))/19;

%RLT
for i=280:299
    for j=1:6
    sum_RLT(j,i) = (s_RLT(j,i)+s_RLT(j,(i+1)))/2;
    end
end
S_RLT(1,1)=sum(sum_RLT(1,280:299))/19;
S_RLT(2,1)=sum(sum_RLT(2,280:299))/19;
S_RLT(3,1)=sum(sum_RLT(3,280:299))/19;
S_RLT(4,1)=sum(sum_RLT(4,280:299))/19;
S_RLT(5,1)=sum(sum_RLT(5,280:299))/19;
S_RLT(6,1)=sum(sum_RLT(6,280:299))/19;

%KL
for i=280:299
    for j=1:6
    sum_KL(j,i) = (s_KL(j,i)+s_KL(j,(i+1)))/2;
    end
end
S_KL(1,1)=sum(sum_KL(1,280:299))/19;
S_KL(2,1)=sum(sum_KL(2,280:299))/19;
S_KL(3,1)=sum(sum_KL(3,280:299))/19;
S_KL(4,1)=sum(sum_KL(4,280:299))/19;
S_KL(5,1)=sum(sum_KL(5,280:299))/19;
S_KL(6,1)=sum(sum_KL(6,280:299))/19;

%tauL
for i=280:299
    for j=1:6
    sum_tauL(j,i) = (s_tauL(j,i)+s_tauL(j,(i+1)))/2;
    end
end
S_tauL(1,1)=sum(sum_tauL(1,280:299))/19;
S_tauL(2,1)=sum(sum_tauL(2,280:299))/19;
S_tauL(3,1)=sum(sum_tauL(3,280:299))/19;
S_tauL(4,1)=sum(sum_tauL(4,280:299))/19;
S_tauL(5,1)=sum(sum_tauL(5,280:299))/19;
S_tauL(6,1)=sum(sum_tauL(6,280:299))/19;

%kdL
for i=280:299
    for j=1:6
    sum_kdL(j,i) = (s_kdL(j,i)+s_kdL(j,(i+1)))/2;
    end
end
S_kdL(1,1)=sum(sum_kdL(1,280:299))/19;
S_kdL(2,1)=sum(sum_kdL(2,280:299))/19;
S_kdL(3,1)=sum(sum_kdL(3,280:299))/19;
S_kdL(4,1)=sum(sum_kdL(4,280:299))/19;
S_kdL(5,1)=sum(sum_kdL(5,280:299))/19;
S_kdL(6,1)=sum(sum_kdL(6,280:299))/19;

%_____PHASE2___________
%Calculating for LX1
for i=300:319
    for j=1:6
    sum_LX1(j,i) = (s_LX1(j,i)+s_LX1(j,(i+1)))/2;
    end
end
S2_LX1(1,1)=sum(sum_LX1(1,300:319))/19;
S2_LX1(2,1)=sum(sum_LX1(2,300:319))/19;
S2_LX1(3,1)=sum(sum_LX1(3,300:319))/19;
S2_LX1(4,1)=sum(sum_LX1(4,300:319))/19;
S2_LX1(5,1)=sum(sum_LX1(5,300:319))/19;
S2_LX1(6,1)=sum(sum_LX1(6,300:319))/19;

%LX2
for i=300:319
    for j=1:6
    sum_LX2(j,i) = (s_LX2(j,i)+s_LX2(j,(i+1)))/2;
    end
end
S2_LX2(1,1)=sum(sum_LX2(1,300:319))/19;
S2_LX2(2,1)=sum(sum_LX2(2,300:319))/19;
S2_LX2(3,1)=sum(sum_LX2(3,300:319))/19;
S2_LX2(4,1)=sum(sum_LX2(4,300:319))/19;
S2_LX2(5,1)=sum(sum_LX2(5,300:319))/19;
S2_LX2(6,1)=sum(sum_LX2(6,300:319))/19;

%LX3
for i=300:319
    for j=1:6
    sum_LX3(j,i) = (s_LX3(j,i)+s_LX3(j,(i+1)))/2;
    end
end
S2_LX3(1,1)=sum(sum_LX3(1,300:319))/19;
S2_LX3(2,1)=sum(sum_LX3(2,300:319))/19;
S2_LX3(3,1)=sum(sum_LX3(3,300:319))/19;
S2_LX3(4,1)=sum(sum_LX3(4,300:319))/19;
S2_LX3(5,1)=sum(sum_LX3(5,300:319))/19;
S2_LX3(6,1)=sum(sum_LX3(6,300:319))/19;

%eX
for i=300:319
    for j=1:6
    sum_eX(j,i) = (s_eX(j,i)+s_eX(j,(i+1)))/2;
    end
end
S2_eX(1,1)=sum(sum_eX(1,300:319))/19;
S2_eX(2,1)=sum(sum_eX(2,300:319))/19;
S2_eX(3,1)=sum(sum_eX(3,300:319))/19;
S2_eX(4,1)=sum(sum_eX(4,300:319))/19;
S2_eX(5,1)=sum(sum_eX(5,300:319))/19;
S2_eX(6,1)=sum(sum_eX(6,300:319))/19;

%RXT
for i=300:319
    for j=1:6
    sum_RXT(j,i) = (s_RXT(j,i)+s_RXT(j,(i+1)))/2;
    end
end
S2_RXT(1,1)=sum(sum_RXT(1,300:319))/19;
S2_RXT(2,1)=sum(sum_RXT(2,300:319))/19;
S2_RXT(3,1)=sum(sum_RXT(3,300:319))/19;
S2_RXT(4,1)=sum(sum_RXT(4,300:319))/19;
S2_RXT(5,1)=sum(sum_RXT(5,300:319))/19;
S2_RXT(6,1)=sum(sum_RXT(6,300:319))/19;

%Gj
for i=300:319
    for j=1:6
    sum_Gj(j,i) = (s_Gj(j,i)+s_Gj(j,(i+1)))/2;
    end
end
S2_Gj(1,1)=sum(sum_Gj(1,300:319))/19;
S2_Gj(2,1)=sum(sum_Gj(2,300:319))/19;
S2_Gj(3,1)=sum(sum_Gj(3,300:319))/19;
S2_Gj(4,1)=sum(sum_Gj(4,300:319))/19;
S2_Gj(5,1)=sum(sum_Gj(5,300:319))/19;
S2_Gj(6,1)=sum(sum_Gj(6,300:319))/19;

%KX
for i=300:319
    for j=1:6
    sum_KX(j,i) = (s_KX(j,i)+s_KX(j,(i+1)))/2;
    end
end
S2_KX(1,1)=sum(sum_KX(1,300:319))/19;
S2_KX(2,1)=sum(sum_KX(2,300:319))/19;
S2_KX(3,1)=sum(sum_KX(3,300:319))/19;
S2_KX(4,1)=sum(sum_KX(4,300:319))/19;
S2_KX(5,1)=sum(sum_KX(5,300:319))/19;
S2_KX(6,1)=sum(sum_KX(6,300:319))/19;

%tauX
for i=300:319
    for j=1:6
    sum_tauX(j,i) = (s_tauX(j,i)+s_tauX(j,(i+1)))/2;
    end
end
S2_tauX(1,1)=sum(sum_tauX(1,300:319))/19;
S2_tauX(2,1)=sum(sum_tauX(2,300:319))/19;
S2_tauX(3,1)=sum(sum_tauX(3,300:319))/19;
S2_tauX(4,1)=sum(sum_tauX(4,300:319))/19;
S2_tauX(5,1)=sum(sum_tauX(5,300:319))/19;
S2_tauX(6,1)=sum(sum_tauX(6,300:319))/19;

%kdX
for i=300:319
    for j=1:6
    sum_kdX(j,i) = (s_kdX(j,i)+s_kdX(j,(i+1)))/2;
    end
end
S2_kdX(1,1)=sum(sum_kdX(1,300:319))/19;
S2_kdX(2,1)=sum(sum_kdX(2,300:319))/19;
S2_kdX(3,1)=sum(sum_kdX(3,300:319))/19;
S2_kdX(4,1)=sum(sum_kdX(4,300:319))/19;
S2_kdX(5,1)=sum(sum_kdX(5,300:319))/19;
S2_kdX(6,1)=sum(sum_kdX(6,300:319))/19;

%mu
for i=300:319
    for j=1:6
    sum_mu(j,i) = (s_mu(j,i)+s_mu(j,(i+1)))/2;
    end
end
S2_mu(1,1)=sum(sum_mu(1,300:319))/19;
S2_mu(2,1)=sum(sum_mu(2,300:319))/19;
S2_mu(3,1)=sum(sum_mu(3,300:319))/19;
S2_mu(4,1)=sum(sum_mu(4,300:319))/19;
S2_mu(5,1)=sum(sum_mu(5,300:319))/19;
S2_mu(6,1)=sum(sum_mu(6,300:319))/19;

%WI1
for i=300:319
    for j=1:6
    sum_WI1(j,i) = (s_WI1(j,i)+s_WI1(j,(i+1)))/2;
    end
end
S2_WI1(1,1)=sum(sum_WI1(1,300:319))/19;
S2_WI1(2,1)=sum(sum_WI1(2,300:319))/19;
S2_WI1(3,1)=sum(sum_WI1(3,300:319))/19;
S2_WI1(4,1)=sum(sum_WI1(4,300:319))/19;
S2_WI1(5,1)=sum(sum_WI1(5,300:319))/19;
S2_WI1(6,1)=sum(sum_WI1(6,300:319))/19;

%W11
for i=300:319
    for j=1:6
    sum_W11(j,i) = (s_W11(j,i)+s_W11(j,(i+1)))/2;
    end
end
S2_W11(1,1)=sum(sum_W11(1,300:319))/19;
S2_W11(2,1)=sum(sum_W11(2,300:319))/19;
S2_W11(3,1)=sum(sum_W11(3,300:319))/19;
S2_W11(4,1)=sum(sum_W11(4,300:319))/19;
S2_W11(5,1)=sum(sum_W11(5,300:319))/19;
S2_W11(6,1)=sum(sum_W11(6,300:319))/19;

%W12
for i=300:319
    for j=1:6
    sum_W12(j,i) = (s_W12(j,i)+s_W12(j,(i+1)))/2;
    end
end
S2_W12(1,1)=sum(sum_W12(1,300:319))/19;
S2_W12(2,1)=sum(sum_W12(2,300:319))/19;
S2_W12(3,1)=sum(sum_W12(3,300:319))/19;
S2_W12(4,1)=sum(sum_W12(4,300:319))/19;
S2_W12(5,1)=sum(sum_W12(5,300:319))/19;
S2_W12(6,1)=sum(sum_W12(6,300:319))/19;

%W13
for i=300:319
    for j=1:6
    sum_W13(j,i) = (s_W13(j,i)+s_W13(j,(i+1)))/2;
    end
end
S2_W13(1,1)=sum(sum_W13(1,300:319))/19;
S2_W13(2,1)=sum(sum_W13(2,300:319))/19;
S2_W13(3,1)=sum(sum_W13(3,300:319))/19;
S2_W13(4,1)=sum(sum_W13(4,300:319))/19;
S2_W13(5,1)=sum(sum_W13(5,300:319))/19;
S2_W13(6,1)=sum(sum_W13(6,300:319))/19;

%W22
for i=300:319
    for j=1:6
    sum_W22(j,i) = (s_W22(j,i)+s_W22(j,(i+1)))/2;
    end
end
S2_W22(1,1)=sum(sum_W22(1,300:319))/19;
S2_W22(2,1)=sum(sum_W22(2,300:319))/19;
S2_W22(3,1)=sum(sum_W22(3,300:319))/19;
S2_W22(4,1)=sum(sum_W22(4,300:319))/19;
S2_W22(5,1)=sum(sum_W22(5,300:319))/19;
S2_W22(6,1)=sum(sum_W22(6,300:319))/19;

%W23
for i=300:319
    for j=1:6
    sum_W23(j,i) = (s_W23(j,i)+s_W23(j,(i+1)))/2;
    end
end
S2_W23(1,1)=sum(sum_W23(1,300:319))/19;
S2_W23(2,1)=sum(sum_W23(2,300:319))/19;
S2_W23(3,1)=sum(sum_W23(3,300:319))/19;
S2_W23(4,1)=sum(sum_W23(4,300:319))/19;
S2_W23(5,1)=sum(sum_W23(5,300:319))/19;
S2_W23(6,1)=sum(sum_W23(6,300:319))/19;

%W33
for i=300:319
    for j=1:6
    sum_W33(j,i) = (s_W33(j,i)+s_W33(j,(i+1)))/2;
    end
end
S2_W33(1,1)=sum(sum_W33(1,300:319))/19;
S2_W33(2,1)=sum(sum_W33(2,300:319))/19;
S2_W33(3,1)=sum(sum_W33(3,300:319))/19;
S2_W33(4,1)=sum(sum_W33(4,300:319))/19;
S2_W33(5,1)=sum(sum_W33(5,300:319))/19;
S2_W33(6,1)=sum(sum_W33(6,300:319))/19;

%eL
for i=300:319
    for j=1:6
    sum_eL(j,i) = (s_eL(j,i)+s_eL(j,(i+1)))/2;
    end
end
S2_eL(1,1)=sum(sum_eL(1,300:319))/19;
S2_eL(2,1)=sum(sum_eL(2,300:319))/19;
S2_eL(3,1)=sum(sum_eL(3,300:319))/19;
S2_eL(4,1)=sum(sum_eL(4,300:319))/19;
S2_eL(5,1)=sum(sum_eL(5,300:319))/19;
S2_eL(6,1)=sum(sum_eL(6,300:319))/19;

%RLT
for i=300:319
    for j=1:6
    sum_RLT(j,i) = (s_RLT(j,i)+s_RLT(j,(i+1)))/2;
    end
end
S2_RLT(1,1)=sum(sum_RLT(1,300:319))/19;
S2_RLT(2,1)=sum(sum_RLT(2,300:319))/19;
S2_RLT(3,1)=sum(sum_RLT(3,300:319))/19;
S2_RLT(4,1)=sum(sum_RLT(4,300:319))/19;
S2_RLT(5,1)=sum(sum_RLT(5,300:319))/19;
S2_RLT(6,1)=sum(sum_RLT(6,300:319))/19;

%KL
for i=300:319
    for j=1:6
    sum_KL(j,i) = (s_KL(j,i)+s_KL(j,(i+1)))/2;
    end
end
S2_KL(1,1)=sum(sum_KL(1,300:319))/19;
S2_KL(2,1)=sum(sum_KL(2,300:319))/19;
S2_KL(3,1)=sum(sum_KL(3,300:319))/19;
S2_KL(4,1)=sum(sum_KL(4,300:319))/19;
S2_KL(5,1)=sum(sum_KL(5,300:319))/19;
S2_KL(6,1)=sum(sum_KL(6,300:319))/19;

%tauL
for i=300:319
    for j=1:6
    sum_tauL(j,i) = (s_tauL(j,i)+s_tauL(j,(i+1)))/2;
    end
end
S2_tauL(1,1)=sum(sum_tauL(1,300:319))/19;
S2_tauL(2,1)=sum(sum_tauL(2,300:319))/19;
S2_tauL(3,1)=sum(sum_tauL(3,300:319))/19;
S2_tauL(4,1)=sum(sum_tauL(4,300:319))/19;
S2_tauL(5,1)=sum(sum_tauL(5,300:319))/19;
S2_tauL(6,1)=sum(sum_tauL(6,300:319))/19;

%kdL
for i=300:319
    for j=1:6
    sum_kdL(j,i) = (s_kdL(j,i)+s_kdL(j,(i+1)))/2;
    end
end
S2_kdL(1,1)=sum(sum_kdL(1,300:319))/19;
S2_kdL(2,1)=sum(sum_kdL(2,300:319))/19;
S2_kdL(3,1)=sum(sum_kdL(3,300:319))/19;
S2_kdL(4,1)=sum(sum_kdL(4,300:319))/19;
S2_kdL(5,1)=sum(sum_kdL(5,300:319))/19;
S2_kdL(6,1)=sum(sum_kdL(6,300:319))/19;

%_____PHASE3____________
%Calculating for LX1
for i=580:599
    for j=1:6
    sum_LX1(j,i) = (s_LX1(j,i)+s_LX1(j,(i+1)))/2;
    end
end
S3_LX1(1,1)=sum(sum_LX1(1,580:599))/19;
S3_LX1(2,1)=sum(sum_LX1(2,580:599))/19;
S3_LX1(3,1)=sum(sum_LX1(3,580:599))/19;
S3_LX1(4,1)=sum(sum_LX1(4,580:599))/19;
S3_LX1(5,1)=sum(sum_LX1(5,580:599))/19;
S3_LX1(6,1)=sum(sum_LX1(6,580:599))/19;

%LX2
for i=580:599
    for j=1:6
    sum_LX2(j,i) = (s_LX2(j,i)+s_LX2(j,(i+1)))/2;
    end
end
S3_LX2(1,1)=sum(sum_LX2(1,580:599))/19;
S3_LX2(2,1)=sum(sum_LX2(2,580:599))/19;
S3_LX2(3,1)=sum(sum_LX2(3,580:599))/19;
S3_LX2(4,1)=sum(sum_LX2(4,580:599))/19;
S3_LX2(5,1)=sum(sum_LX2(5,580:599))/19;
S3_LX2(6,1)=sum(sum_LX2(6,580:599))/19;

%LX3
for i=580:599
    for j=1:6
    sum_LX3(j,i) = (s_LX3(j,i)+s_LX3(j,(i+1)))/2;
    end
end
S3_LX3(1,1)=sum(sum_LX3(1,580:599))/19;
S3_LX3(2,1)=sum(sum_LX3(2,580:599))/19;
S3_LX3(3,1)=sum(sum_LX3(3,580:599))/19;
S3_LX3(4,1)=sum(sum_LX3(4,580:599))/19;
S3_LX3(5,1)=sum(sum_LX3(5,580:599))/19;
S3_LX3(6,1)=sum(sum_LX3(6,580:599))/19;

%eX
for i=580:599
    for j=1:6
    sum_eX(j,i) = (s_eX(j,i)+s_eX(j,(i+1)))/2;
    end
end
S3_eX(1,1)=sum(sum_eX(1,580:599))/19;
S3_eX(2,1)=sum(sum_eX(2,580:599))/19;
S3_eX(3,1)=sum(sum_eX(3,580:599))/19;
S3_eX(4,1)=sum(sum_eX(4,580:599))/19;
S3_eX(5,1)=sum(sum_eX(5,580:599))/19;
S3_eX(6,1)=sum(sum_eX(6,580:599))/19;

%RXT
for i=580:599
    for j=1:6
    sum_RXT(j,i) = (s_RXT(j,i)+s_RXT(j,(i+1)))/2;
    end
end
S3_RXT(1,1)=sum(sum_RXT(1,580:599))/19;
S3_RXT(2,1)=sum(sum_RXT(2,580:599))/19;
S3_RXT(3,1)=sum(sum_RXT(3,580:599))/19;
S3_RXT(4,1)=sum(sum_RXT(4,580:599))/19;
S3_RXT(5,1)=sum(sum_RXT(5,580:599))/19;
S3_RXT(6,1)=sum(sum_RXT(6,580:599))/19;

%Gj
for i=580:599
    for j=1:6
    sum_Gj(j,i) = (s_Gj(j,i)+s_Gj(j,(i+1)))/2;
    end
end
S3_Gj(1,1)=sum(sum_Gj(1,580:599))/19;
S3_Gj(2,1)=sum(sum_Gj(2,580:599))/19;
S3_Gj(3,1)=sum(sum_Gj(3,580:599))/19;
S3_Gj(4,1)=sum(sum_Gj(4,580:599))/19;
S3_Gj(5,1)=sum(sum_Gj(5,580:599))/19;
S3_Gj(6,1)=sum(sum_Gj(6,580:599))/19;

%KX
for i=580:599
    for j=1:6
    sum_KX(j,i) = (s_KX(j,i)+s_KX(j,(i+1)))/2;
    end
end
S3_KX(1,1)=sum(sum_KX(1,580:599))/19;
S3_KX(2,1)=sum(sum_KX(2,580:599))/19;
S3_KX(3,1)=sum(sum_KX(3,580:599))/19;
S3_KX(4,1)=sum(sum_KX(4,580:599))/19;
S3_KX(5,1)=sum(sum_KX(5,580:599))/19;
S3_KX(6,1)=sum(sum_KX(6,580:599))/19;

%tauX
for i=580:599
    for j=1:6
    sum_tauX(j,i) = (s_tauX(j,i)+s_tauX(j,(i+1)))/2;
    end
end
S3_tauX(1,1)=sum(sum_tauX(1,580:599))/19;
S3_tauX(2,1)=sum(sum_tauX(2,580:599))/19;
S3_tauX(3,1)=sum(sum_tauX(3,580:599))/19;
S3_tauX(4,1)=sum(sum_tauX(4,580:599))/19;
S3_tauX(5,1)=sum(sum_tauX(5,580:599))/19;
S3_tauX(6,1)=sum(sum_tauX(6,580:599))/19;

%kdX
for i=580:599
    for j=1:6
    sum_kdX(j,i) = (s_kdX(j,i)+s_kdX(j,(i+1)))/2;
    end
end
S3_kdX(1,1)=sum(sum_kdX(1,580:599))/19;
S3_kdX(2,1)=sum(sum_kdX(2,580:599))/19;
S3_kdX(3,1)=sum(sum_kdX(3,580:599))/19;
S3_kdX(4,1)=sum(sum_kdX(4,580:599))/19;
S3_kdX(5,1)=sum(sum_kdX(5,580:599))/19;
S3_kdX(6,1)=sum(sum_kdX(6,580:599))/19;

%mu
for i=580:599
    for j=1:6
    sum_mu(j,i) = (s_mu(j,i)+s_mu(j,(i+1)))/2;
    end
end
S3_mu(1,1)=sum(sum_mu(1,580:599))/19;
S3_mu(2,1)=sum(sum_mu(2,580:599))/19;
S3_mu(3,1)=sum(sum_mu(3,580:599))/19;
S3_mu(4,1)=sum(sum_mu(4,580:599))/19;
S3_mu(5,1)=sum(sum_mu(5,580:599))/19;
S3_mu(6,1)=sum(sum_mu(6,580:599))/19;

%WI1
for i=580:599
    for j=1:6
    sum_WI1(j,i) = (s_WI1(j,i)+s_WI1(j,(i+1)))/2;
    end
end
S3_WI1(1,1)=sum(sum_WI1(1,580:599))/19;
S3_WI1(2,1)=sum(sum_WI1(2,580:599))/19;
S3_WI1(3,1)=sum(sum_WI1(3,580:599))/19;
S3_WI1(4,1)=sum(sum_WI1(4,580:599))/19;
S3_WI1(5,1)=sum(sum_WI1(5,580:599))/19;
S3_WI1(6,1)=sum(sum_WI1(6,580:599))/19;

%W11
for i=580:599
    for j=1:6
    sum_W11(j,i) = (s_W11(j,i)+s_W11(j,(i+1)))/2;
    end
end
S3_W11(1,1)=sum(sum_W11(1,580:599))/19;
S3_W11(2,1)=sum(sum_W11(2,580:599))/19;
S3_W11(3,1)=sum(sum_W11(3,580:599))/19;
S3_W11(4,1)=sum(sum_W11(4,580:599))/19;
S3_W11(5,1)=sum(sum_W11(5,580:599))/19;
S3_W11(6,1)=sum(sum_W11(6,580:599))/19;

%W12
for i=580:599
    for j=1:6
    sum_W12(j,i) = (s_W12(j,i)+s_W12(j,(i+1)))/2;
    end
end
S3_W12(1,1)=sum(sum_W12(1,580:599))/19;
S3_W12(2,1)=sum(sum_W12(2,580:599))/19;
S3_W12(3,1)=sum(sum_W12(3,580:599))/19;
S3_W12(4,1)=sum(sum_W12(4,580:599))/19;
S3_W12(5,1)=sum(sum_W12(5,580:599))/19;
S3_W12(6,1)=sum(sum_W12(6,580:599))/19;

%W13
for i=580:599
    for j=1:6
    sum_W13(j,i) = (s_W13(j,i)+s_W13(j,(i+1)))/2;
    end
end
S3_W13(1,1)=sum(sum_W13(1,580:599))/19;
S3_W13(2,1)=sum(sum_W13(2,580:599))/19;
S3_W13(3,1)=sum(sum_W13(3,580:599))/19;
S3_W13(4,1)=sum(sum_W13(4,580:599))/19;
S3_W13(5,1)=sum(sum_W13(5,580:599))/19;
S3_W13(6,1)=sum(sum_W13(6,580:599))/19;

%W22
for i=580:599
    for j=1:6
    sum_W22(j,i) = (s_W22(j,i)+s_W22(j,(i+1)))/2;
    end
end
S3_W22(1,1)=sum(sum_W22(1,580:599))/19;
S3_W22(2,1)=sum(sum_W22(2,580:599))/19;
S3_W22(3,1)=sum(sum_W22(3,580:599))/19;
S3_W22(4,1)=sum(sum_W22(4,580:599))/19;
S3_W22(5,1)=sum(sum_W22(5,580:599))/19;
S3_W22(6,1)=sum(sum_W22(6,580:599))/19;

%W23
for i=580:599
    for j=1:6
    sum_W23(j,i) = (s_W23(j,i)+s_W23(j,(i+1)))/2;
    end
end
S3_W23(1,1)=sum(sum_W23(1,580:599))/19;
S3_W23(2,1)=sum(sum_W23(2,580:599))/19;
S3_W23(3,1)=sum(sum_W23(3,580:599))/19;
S3_W23(4,1)=sum(sum_W23(4,580:599))/19;
S3_W23(5,1)=sum(sum_W23(5,580:599))/19;
S3_W23(6,1)=sum(sum_W23(6,580:599))/19;

%W33
for i=580:599
    for j=1:6
    sum_W33(j,i) = (s_W33(j,i)+s_W33(j,(i+1)))/2;
    end
end
S3_W33(1,1)=sum(sum_W33(1,580:599))/19;
S3_W33(2,1)=sum(sum_W33(2,580:599))/19;
S3_W33(3,1)=sum(sum_W33(3,580:599))/19;
S3_W33(4,1)=sum(sum_W33(4,580:599))/19;
S3_W33(5,1)=sum(sum_W33(5,580:599))/19;
S3_W33(6,1)=sum(sum_W33(6,580:599))/19;

%eL
for i=580:599
    for j=1:6
    sum_eL(j,i) = (s_eL(j,i)+s_eL(j,(i+1)))/2;
    end
end
S3_eL(1,1)=sum(sum_eL(1,580:599))/19;
S3_eL(2,1)=sum(sum_eL(2,580:599))/19;
S3_eL(3,1)=sum(sum_eL(3,580:599))/19;
S3_eL(4,1)=sum(sum_eL(4,580:599))/19;
S3_eL(5,1)=sum(sum_eL(5,580:599))/19;
S3_eL(6,1)=sum(sum_eL(6,580:599))/19;

%RLT
for i=580:599
    for j=1:6
    sum_RLT(j,i) = (s_RLT(j,i)+s_RLT(j,(i+1)))/2;
    end
end
S3_RLT(1,1)=sum(sum_RLT(1,580:599))/19;
S3_RLT(2,1)=sum(sum_RLT(2,580:599))/19;
S3_RLT(3,1)=sum(sum_RLT(3,580:599))/19;
S3_RLT(4,1)=sum(sum_RLT(4,580:599))/19;
S3_RLT(5,1)=sum(sum_RLT(5,580:599))/19;
S3_RLT(6,1)=sum(sum_RLT(6,580:599))/19;

%KL
for i=580:599
    for j=1:6
    sum_KL(j,i) = (s_KL(j,i)+s_KL(j,(i+1)))/2;
    end
end
S3_KL(1,1)=sum(sum_KL(1,580:599))/19;
S3_KL(2,1)=sum(sum_KL(2,580:599))/19;
S3_KL(3,1)=sum(sum_KL(3,580:599))/19;
S3_KL(4,1)=sum(sum_KL(4,580:599))/19;
S3_KL(5,1)=sum(sum_KL(5,580:599))/19;
S3_KL(6,1)=sum(sum_KL(6,580:599))/19;

%tauL
for i=580:599
    for j=1:6
    sum_tauL(j,i) = (s_tauL(j,i)+s_tauL(j,(i+1)))/2;
    end
end
S3_tauL(1,1)=sum(sum_tauL(1,580:599))/19;
S3_tauL(2,1)=sum(sum_tauL(2,580:599))/19;
S3_tauL(3,1)=sum(sum_tauL(3,580:599))/19;
S3_tauL(4,1)=sum(sum_tauL(4,580:599))/19;
S3_tauL(5,1)=sum(sum_tauL(5,580:599))/19;
S3_tauL(6,1)=sum(sum_tauL(6,580:599))/19;

%kdL
for i=580:599
    for j=1:6
    sum_kdL(j,i) = (s_kdL(j,i)+s_kdL(j,(i+1)))/2;
    end
end
S3_kdL(1,1)=sum(sum_kdL(1,580:599))/19;
S3_kdL(2,1)=sum(sum_kdL(2,580:599))/19;
S3_kdL(3,1)=sum(sum_kdL(3,580:599))/19;
S3_kdL(4,1)=sum(sum_kdL(4,580:599))/19;
S3_kdL(5,1)=sum(sum_kdL(5,580:599))/19;
S3_kdL(6,1)=sum(sum_kdL(6,580:599))/19;

%---------------------------------------------------

%Now the full sensitivity matrix will be concatentated and analyzed
    %for each phase

sensitivity1=[S_LX1 S_LX2 S_LX3 S_eX S_RXT S_Gj S_KX S_tauX S_kdX S_mu S_WI1 S_W11 S_W12 S_W13 S_W22 S_W23 S_W33 S_eL S_RLT S_KL S_tauL S_kdL];
sensitivity2=[S2_LX1 S2_LX2 S2_LX3 S2_eX S2_RXT S2_Gj S2_KX S2_tauX S2_kdX S2_mu S2_WI1 S2_W11 S2_W12 S2_W13 S2_W22 S2_W23 S2_W33 S2_eL S2_RLT S2_KL S2_tauL S2_kdL];
sensitivity3=[S3_LX1 S3_LX2 S3_LX3 S3_eX S3_RXT S3_Gj S3_KX S3_tauX S3_kdX S3_mu S3_WI1 S3_W11 S3_W12 S3_W13 S3_W22 S3_W23 S3_W33 S3_eL S3_RLT S3_KL S3_tauL S3_kdL];

%Performing the singular value decomposition for each phase:

[U1,S1,V1]=svd(sensitivity1);
[U2,S2,V2]=svd(sensitivity2);
[U3,S3,V3]=svd(sensitivity3);

disp('U1')
disp(U1(:,1))
disp('U2')
disp(U2(:,1))
disp('U3')
disp(U3(:,1))


    
    





