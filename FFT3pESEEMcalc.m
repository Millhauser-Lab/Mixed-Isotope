%% Section 1: Plot 3pESEEM time domain and spectra for one and two equivalent 15N. This result is shown in SI Figure 4B, calculated using Equation S9-S12. 

clc;clear;

%Input parameters
tau=.210;
va=0.3;
vb=2.5;
k=0.5; %k=sin^2(2*eta)=(B*wI/wa/wb)^2, where eta is the half the angle between the two nuclear quantization axes, B is the hyperfine coupling, wI is the nuclear Larmor frequency, and wa and wb are the nuclear frequencies in the alpha and beta electron spin manifolds.
step=0.016; %(change in T)
points=2048;
deadtime=0.012; %T(0)
tdecay=.8; %exponential decay constant for time domain

%calculation
wa=va*2*pi;
wb=vb*2*pi;
T=deadtime:step:points*step-step+deadtime;

for i=1:points
    V1(i)=0.5*(1-k/2*((1-cos(2*pi*vb*tau))*(1-cos(2*pi*va*(tau+T(i)))))) + 0.5*(1-k/2*((1-cos(2*pi*va*tau))*(1-cos(2*pi*vb*(tau+T(i)))))); %one coupled 15N
    V2(i)=0.5*(1-k/2*((1-cos(2*pi*vb*tau))*(1-cos(2*pi*va*(tau+T(i))))))^2 + 0.5*(1-k/2*((1-cos(2*pi*va*tau))*(1-cos(2*pi*vb*(tau+T(i))))))^2; %two coupled 15N
    Va(i)=((1-cos(wb*tau))*(1-cos(wa*(T(i)+tau))))^2; %two nuclei, just the alpha terms
    Vb(i)=((1-cos(wa*tau))*(1-cos(wb*(T(i)+tau))))^2; %two nuclei, just the beta terms
end

V1a=V1;
V1=V1-mean(V1);
V2a=V2;
V2=V2-mean(V2);
decay=exp(-T*tdecay);
V1=V1.*decay;
%V1=V1/max(V1);
V2=V2.*decay;
%V2=V2/max(V2);

y1 = fft(V1);                               % Compute DFT of x
m1 = abs(y1);                               % Magnitud
p1 = unwrap(angle(y1));                     % Phase
f = (0:length(y1)-1)/step/length(y1);        % Frequency vector
y2 = fft(V2);                               % Compute DFT of x
m2 = abs(y2);                               % Magnitud
p2 = unwrap(angle(y2));                     % Phase
Va=Va-mean(Va);
Va=Va.*decay;
Va=Va/max(Va);
ya = fft(Va);                               % Compute DFT of x
ma = abs(ya);                               % Magnitud
pa = unwrap(angle(ya));                     % Phase

subplot(2,1,1)
p=plot(T,V1,T,V2)
p(1).LineWidth=1.3;
p(2).LineWidth=1.3;
legend('one 15N','two 15N',FontSize=12)
xlabel('T (\mus)',FontSize=12)
xlim([0 8])
subplot(2,1,2)
p=plot(f,m1,f,m2)
p(1).LineWidth=1.3;
p(2).LineWidth=1.3;
xlim([0 8])
legend('one 15N','two 15N',FontSize=12)
xlabel('Frequency (MHz)',FontSize=12)

%% Section 2: Check to see if Expanded equation, Eq S11, is correct: V3 = (1-k/2*(1-cos(wb*tau))*(1-cos(wa*(T(i)+tau))))^2... Only considers the first half of the equation dependent on wa, since the equation is identical when dependent on wb.
% Must be run after Section 1
for i=1:points
    VA(i)=(1-k/2*(1-cos(wb*tau))*(1-cos(wa*(T(i)+tau))))^2;
    V3(i)=(1-k*(1-cos(wb*tau))+k*k/4*(2.25-3*cos(wb*tau)+3/4*cos(2*wb*tau))+(k*(1-cos(wb*tau))-k*k/4*(3-4*cos(wb*tau)+cos(2*wb*tau)))*cos(wa*(T(i)+tau))+k*k/16*(3-4*cos(wb*tau)+cos(2*wb*tau))*cos(2*wa*(T(i)+tau)));
end

norm=1; %change the relative intensity of the two 15N ESEEM spectra

V3=V3-mean(V3);
decay=exp(-T*tdecay);
V3=V3.*decay;
V3=V3/max(V3)*norm;
y3 = fft(V3);                               % Compute DFT of x
m3 = abs(y3);                               % Magnitud
p3 = unwrap(angle(y3));                     % Phase
VA=VA-mean(VA);
VA=VA.*decay;
VA=VA/max(VA);
yA = fft(VA);                               % Compute DFT of x
mA = abs(yA);                               % Magnitud
pA = unwrap(angle(yA));                     % Phase

p=plot(f,mA,f,m3,'--')
p(1).LineWidth=1.3;
p(2).LineWidth=1.3;
xlim([0 8])
legend('one 15N','two 15N',FontSize=12)
xlabel('Frequency (MHz)',FontSize=12)

%% Section 3: Visualize relative intensity for va, vb, 2va, and 2vb at various tau values. This is shown in SI Figure 5
%this section works by MANUALLY picking the location of va, vb, 2va, and
%2vb. Check lines 153-156 if you significantly change the va and vb input.
%It's best to run section 4 after this, as this section has a lot of
%outputs that can be overwhelming

%input parameters
va=0.3;
vb=2.5;
k=0.5;
step=0.016;
points=2048;
deadtime=0.012;
tdecay=0.3;
taumin=.0;  %minimum tau to calculate, microseconds
taumax=3;  %maximum tau to calculate, microseconds

%calculation

wa=va*2*pi;
wb=vb*2*pi;
T=deadtime:step:points*step-step+deadtime;
decay=exp(-T*tdecay);


tauint=0.002;
taunumber=(taumax-taumin)/tauint+1;

tau=zeros(1,taunumber);
V2=zeros(points,taunumber);

for i=1:taunumber
    tau(i)=taumin+tauint*(i-1);
end

for j=1:points
    for i=1:taunumber
        V(j,i)=0.5*(1-k/2*((1-cos(wb*tau(i)))*(1-cos(wa*(tau(i)+T(j))))))^2 + 0.5*(1-k/2*((1-cos(wa*tau(i)))*(1-cos(wb*(tau(i)+T(j))))))^2; 
    end
end

V=V-mean(V);

for i=1:points
    V(i,:)=V(i,:)*decay(i);
end

ft=0*V;
ft=abs(fft(V));
fta=ft;

yoffset=01;

for i=1:taunumber
    fta(:,i)=fta(:,i)+yoffset*(i-1);
end

FAtau=ft(12,:);
F2Atau=ft(22,:)-ft(19,:); %Because 2va is on the shoulder of va, we subtract the contribution from va to best the isolated intensity of 2va
FBtau=ft(83,:);
F2Btau=ft(164,:);

f = (0:length(ft)-1)/step/length(ft);        % Frequency vector

subplot(2,2,1)
plot(T,V)
title('All time domain signals',FontSize=12)
xlabel('T (\mus)',FontSize=12)
xlim([0 10])
subplot(2,2,2)
plot(f,ft)
title('All frequency spectra',FontSize=12)
xlabel('Frequency (MHz)',FontSize=12)
xlim([0 8])
subplot(2,2,3)
p=plot(tau,FAtau,'r',tau,F2Atau,'r--',tau,FBtau,'b',tau,F2Btau,'b--')
p(1).LineWidth=1.3;
p(2).LineWidth=1.3;
p(3).LineWidth=1.3;
p(4).LineWidth=1.3;
title('Relative peak intensity as a function of \tau',FontSize=12)
xlabel('\tau (\mus)',FontSize=12)
legend('\nu_\alpha = 0.33 MHz','2\nu_\alpha = 0.64 MHz','\nu_\beta = 2.50 MHz','2\nu_\beta = 5.00 MHz')
subplot(2,2,4)
p=plot(f,ft(:,700))
p(1).LineWidth=1.3;
xlim([0 8])
title('Visualize a specific spectra',FontSize=12)
xlabel('Frequency (MHz)',FontSize=12)

%% Section 4: Replot previous figure for Relative peak intensitys as a function of tau.
%Requires Section 3 to be run. Lines 200 and 201 can be used to input vertical lines at experimental tau
tau1000=1000*tau;

subplot(1,1,1)
p=plot(tau1000,FAtau,'r',tau1000,F2Atau,'r--',tau1000,FBtau,'b',tau1000,F2Btau,'b--')
p(1).LineWidth=1.3;
p(2).LineWidth=1.3;
p(3).LineWidth=1.3;
p(4).LineWidth=1.3;
xlim([0 1500])
title('Relative Peak Intensity as a Function of \tau',FontSize=12)
xlabel('\tau (ns)',FontSize=12)
ylabel('ESEEM intensity')
xline(144,'--','Color','#77AC30') %line one
xline(210,'Color','#77AC30')      %line two
legend('\nu_\alpha = 0.33 MHz','2\nu_\alpha = 0.64 MHz','\nu_\beta = 2.50 MHz','2\nu_\beta = 5.00 MHz','\tau = 144 ns','\tau = 210 ns')

