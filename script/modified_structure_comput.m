clear all
close all
clc


% Per risolvere il punto 7 aggiungo una trave nel pillar di destra, tra il nodo centrale del pillar e che 
% si collega a terra ad una distanza di 4 metri dalla base del pillar di
% destra. Ho usato una sezione IPE 270
E = 2.06e11;
rho = 7800;

ome_max=2*pi*20;

Iyb = 5.790e-5;
Ab = 4.595e-3;
Iyr = 9.2080e-4;
Ar = 0.01560;

mb=rho*Ab
mr=rho*Ar;
E*Ab
E*Iyb
%% Q-7

load("modified_structure_mkr.mat");

%Matrices partition
ndof=60;
ntot=3*23;

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
xC0=[1;0;0;1;0;0;0;0;0];
for k=1:length(vect_f)
    omega=vect_f(k)*2*pi; 
    A=-omega^2*MFF+i*omega*CFF+KFF;
    vect_fC0=-(-omega^2*MFC+i*omega*CFC+KFC)*xC0; 
    x0=A\vect_fC0; %computation of the response of free coordinates
    Qc=(-omega^2*MCF+i*omega*CCF+KCF)*x0+(-omega^2*MCC+i*omega*CCC+KCC)*xC0; %constraint
    H1_new=Qc(1);
    M1_new=Qc(3);
    mod1(k)=abs(H1_new);
    phase1(k)=angle(H1_new);
    mod2(k)=abs(M1_new);
    phase2(k)=angle(M1_new);

end


figure;
subplot 211;plot(vect_f,mod1);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('H')
subplot 212;plot(vect_f,phase1);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')
figure;
subplot 211;plot(vect_f,mod2);grid;xlabel('Frequency [Hz]');ylabel('Amplitude');title('M')
subplot 212;plot(vect_f,phase2);grid;xlabel('Frequency [Hz]');ylabel('Phase [rad]')