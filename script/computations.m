clear all
close all
clc

E = 2.06e11;
rho = 7800;

ome_max=2*pi*20;

Iyb = 3.892e-5;
Ab = 3.912e-3;
Iyr = 4.8210e-4;
Ar = 0.01155;

mb=rho*Ab
mr=rho*Ar;
E*Ar
E*Iyr
E*Ab
E*Iyb

Lmaxb=sqrt((pi^2/(1.5*ome_max))*sqrt((E*Iyb/mb)));
Lmaxr=sqrt((pi^2/(1.5*ome_max))*sqrt((E*Iyr/mr)));

%% Q-4 
load("structural_mat_mkr.mat");
%Matrices partition
ndof=60;
ntot=3*22;

MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);
CFF=R(1:ndof,1:ndof);

MFC=M(1:ndof,ndof+1:ntot);
KFC=K(1:ndof,ndof+1:ntot);
CFC=R(1:ndof,ndof+1:ntot);

MCF=M(ndof+1:ntot,1:ndof);
KCF=K(ndof+1:ntot,1:ndof);
CCF=R(ndof+1:ntot,1:ndof);

MCC=M(ndof+1:ntot,ndof+1:ntot);
KCC=K(ndof+1:ntot,ndof+1:ntot);
CCC=R(ndof+1:ntot,ndof+1:ntot);



i=sqrt(-1); %the imaginary operator i has to be defined in the script
vect_f=0:0.01:20; %vector that collects all frequencies
xC0=[1;0;0;1;0;0];
for k=1:length(vect_f)
    omega=vect_f(k)*2*pi; 
    A=-omega^2*MFF+i*omega*CFF+KFF;
    vect_fC0=-(-omega^2*MFC+i*omega*CFC+KFC)*xC0; 
    x0=A\vect_fC0; %computation of the response of free coordinates
    Qc=(-omega^2*MCF+i*omega*CCF+KCF)*x0+(-omega^2*MCC+i*omega*CCC+KCC)*xC0; %constraint
    H1=Qc(1);
    M1=Qc(3);
    mod1(k)=abs(H1);
    phase1(k)=angle(H1);
    mod2(k)=abs(M1);
    phase2(k)=angle(M1);

end


figure;
subplot 211;plot(vect_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('H')
subplot 212;plot(vect_f,phase1);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod2);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('M')
subplot 212;plot(vect_f,phase2);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')

%% Q-5
T = 0.3; % period
A = 2000; % amplitude

t = 0:1e-3:1.5; % time vector
ty=0:1e-3:T;
y = zeros(1,length(ty));
% Create the triangular wave
y = A*(2/T)*abs(mod(ty-T/4,T) - T/2) - A/2;



figure()
plot(ty, y);
xlabel('Time (seconds)');
ylabel('Amplitude');



fftout=fft(y);
N=length(y);
df=1/T;
fmax=(N/2-1)*df;
vett_freq=0:df:fmax;
modf(1)=1/N*abs(fftout(1));
modf(2:N/2)=2/N*abs(fftout(2:N/2));
fasf(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211;bar(vett_freq,modf);
subplot 212;plot(vett_freq,fasf);

dof_yD=idb(11,2); %vertical displacement for point D -->node 11


ome0=2*pi/T; %foundamental frequency
output_yD=zeros(1,length(t));
output_yDdd=zeros(1,length(t));

vect_force=zeros(60,1);
vect_force(dof_yD,1)=1;

for k = 1:N/2
   omega(k)=2*pi*vett_freq(k);
   if omega(k)<=(pi*fmax*2)
    A=-omega(k)^2*MFF+i*omega(k)*CFF+KFF;
    x0 = A\(modf(k)*vect_force);
    yD=x0(dof_yD);
    yDdd=-omega(k)^2*yD;
    output_yD=output_yD+abs(yD)*cos(omega(k)*t+angle(yD));
    output_yDdd=output_yDdd+abs(yDdd)*cos(omega(k)*t+angle(yDdd));
    modyD(k)=abs(yD);
    phaseyD(k)=angle(yD);
   end
end

figure;
plot(t,output_yD,'LineWidth',1)
title("Vertical displacement of point D due to vertical periodic force");
figure;
plot(t,output_yDdd,'LineWidth',1)
title("Vertical acceleration of point D due to vertical periodic force");

figure;
subplot 211;bar(vett_freq,modyD);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('yD')
subplot 212;plot(vett_freq,phaseyD);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')


%% Q-6 Unbalanced rotating mass

dof_xB=idb(12,1); %horizontal displacement for point B -->node 12
dof_yB=idb(12,2); %vertical displacement for point B -->node 12

dof_xC=idb(15,1); %horizontal displacement for point C -->node 12
dof_yC=idb(15,2); %vertical displacement for point C -->node 12

dof_yA=idb(7,2); %vertical displacement for point C -->node 12
dof_xA=idb(7,1); %horizontal displacement for point C -->node 12


epsilon=0.005; %eccentricity of the rotating 
ms=5;

vect_force=zeros(60,1);
vect_force(dof_yA,1)=1; %the real value of the sinusoidal forcing input has to be computed for each value of frequency omega
vect_force(dof_xA,1)=1; %the real value of the sinusoidal forcing input has to be computed for each value of frequency omega

for k=1:length(vect_f)
    omega=vect_f(k)*2*pi; 
    A=-omega^2*MFF+i*omega*CFF+KFF;
    vect_f0=vect_force*ms*epsilon*omega^2;
    x0=A\vect_f0; %computation of the response of free coordinates
    yA=x0(dof_yA);
    xB=x0(dof_xB);
    xBdd=-omega^2*xB;
    yB=x0(dof_yB);
    yBdd=-omega^2*yB;
    xC=x0(dof_xC);
    xCdd=-omega^2*xC;
    yC=x0(dof_yC);
    yCdd=-omega^2*yC;
%     mod1(k)=abs(xB);
%     phase1(k)=angle(xB);
%     mod2(k)=abs(xBdd);
%     phase2(k)=angle(xBdd);
%     mod3(k)=abs(yB);
%     phase3(k)=angle(yB);
%     mod4(k)=abs(yBdd);
%     phase4(k)=angle(yBdd);
    mod1(k)=abs(xC);
    phase1(k)=angle(xC);
    mod2(k)=abs(xCdd);
    phase2(k)=angle(xCdd);
    mod3(k)=abs(yC);
    phase3(k)=angle(yC);
    mod4(k)=abs(yCdd);
    phase4(k)=angle(yCdd);

end
figure;
subplot 211;plot(vect_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('xC')
subplot 212;plot(vect_f,phase1);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod2);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('xCdd')
subplot 212;plot(vect_f,phase2);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')

figure;
subplot 211;plot(vect_f,mod3);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('yC')
subplot 212;plot(vect_f,phase3);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod4);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('yCdd')
subplot 212;plot(vect_f,phase4);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')


