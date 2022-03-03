clear all;close all
simulationname='E00002noOscIAI0.5vth1600me'
L=2*pi; %unit of de
NT=100000;NTOUT=100;
NG=64*4;%mesh number
DT=L/NG*0.7;   %unit of 1/wpe
PPC=128; %Particle per cell, 128 is good.
N=NG*PPC;%total particle number
WP=1;%weight?
QM=-1;
QMion=1/1600;
V0=0;
VTe=0.025;  %unit of c
VTi=VTe*sqrt(QMion)/4; %unit of c it is actually 5 keV, maybe I should reduced it to 1 keV, like 0.05
%VTi=0;
XP1e=1;

V1=0.5*VTe;
mode=1;
Q=WP^2/(QM*N/L); % Don't know, charge per partical
Qi=WP^2/(QMion*N/L)*QMion;

%rho_back=-Q*N/L; % background density
dx=L/NG;
% initial loading for the 2 Stream instability
xpe=linspace(L/N/2,L-L/N/2,N)';
xpi=xpe;
vpe=VTe*randn(N,1);
vpey=VTe*randn(N,1);
vpi=VTi*randn(N,1);
vpiy=VTi*randn(N,1);

%pm=[1:N]';
%pm=1-2*mod(pm,2);
%vpe=vpe+pm.*V0;
% Perturbation
vpe=vpe+V1;
vpi=vpi;%+V1;
%xpe=xpe+XP1e*(L/N)*cos(2*pi*xpe/L*mode);
%xpi=xpe;

vpiOUT=zeros(N,NT/NTOUT);
vpeOUT=vpiOUT;
rhoiOUT=zeros(NG,NT/NTOUT); %unit of wpe^2
rhoeOUT=rhoiOUT;
Eout=rhoiOUT;

Veedge=linspace(-VTe*10, VTe*10, 100)+V1;
Viedge=linspace(-VTi*10,VTi*10,100);
gedge=0:2:NG;
pexOUT=zeros(length(Veedge)-1,length(gedge)-1,NT/NTOUT);
pixOUT=zeros(length(Viedge)-1,length(gedge)-1,NT/NTOUT);



p=1:N;p=[p p];
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
% Main computational cycle


for it=1:NT

    % update xp
    xpe=xpe+vpe*DT;
    xpi=xpi+vpi*DT;
    % apply bc on the particle positions
    out=(xpe<0); xpe(out)=xpe(out)+L;
    %vpe(out)=VTe*randn(1,1)+V1;
    out=(xpi<0); xpi(out)=xpi(out)+L;
    %vpi(out)=VTi*randn(1,1);


    out=(xpe>=L);xpe(out)=xpe(out)-L;
    %vpe(out)=VTe*randn(1,1)+V1;

    out=(xpi>=L);xpi(out)=xpi(out)-L;
    %vpi(out)=VTi*randn(1,1);


    % projection p->g
    g1=floor(xpe/dx-.5)+1;g=[g1;g1+1]; %mark which cell the particle is in
    g1i=floor(xpi/dx-.5)+1;gi=[g1i;g1i+1]; %mark which cell the particle is in

    fraz1=1-abs(xpe/dx-g1+.5);fraz=[fraz1;1-fraz1]; %position in the cell
    fraz1i=1-abs(xpi/dx-g1i+.5);frazi=[fraz1i;1-fraz1i]; %position in the cell

    % apply bc on the projection
    out=(g<1);g(out)=g(out)+NG;
    out=(gi<1);gi(out)=gi(out)+NG;

    out=(g>NG);g(out)=g(out)-NG;
    out=(gi>NG);gi(out)=gi(out)-NG;

    mat=sparse(p,g,fraz,N,NG);
    mati=sparse(p,gi,frazi,N,NG);
    rhoe=full((Q/dx).*sum(mat))';
    rhoi=full((Qi/dx).*sum(mati))';
    rho=rhoe+rhoi; %rho is normalized charge? normalized to n0?
    % computing fields
    Phi=Poisson\(-rho(1:NG-1)*dx^2);%matrix left division: A * x=B, grad^2 phi=-rho
    Phi=[Phi;0];
    Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
    %Eext=0;
    %Eext=-0.00001+0.0001*sin(-2*pi*it*DT/500+2*pi*0*L*(1:NG)'/NG/500);
    Eext=-0.00002;

    Eg=Eg+Eext;
    %Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx)+0.0001; %modified to add external electric field.
    % projection q->p and update of vp
    B=0;
    %B=0.1;
    vpe=vpe+mat*QM*Eg*DT+vpey*QM*B*DT;
    vpey=vpey-vpe*B*QM*DT;
    vpi=vpi+mati*QMion*Eg*DT+vpiy*QMion*B*DT;
    vpiy=vpiy-vpi*B*QMion*DT;
    
    if ~mod(it,NTOUT)
        vpiOUT(:,it/NTOUT)=vpi;
        vpeOUT(:,it/NTOUT)=vpe;
        rhoiOUT(:,it/NTOUT)=rhoi;
        rhoeOUT(:,it/NTOUT)=rhoe;
        Eout(:,it/NTOUT)=Eg;
        disp(it);
        pexOUT(:,:,it/NTOUT)=histcounts2(vpe,g1,Veedge,gedge);
        pixOUT(:,:,it/NTOUT)=histcounts2(vpi,g1i,Viedge,gedge);
    end
    

   
end
%%
clear Ve Vi;
for i=1:NT/NTOUT
    [Ve(:,i),Veedge]=histcounts(vpeOUT(:,i),Veedge);
    [Vi(:,i),Viedge]=histcounts(vpiOUT(:,i),Viedge);
end
%%
Veaxis=Veedge(1:end-1)+(Veedge(2)-Veedge(1))/2;
Taxis=(1:NT/NTOUT)*DT*NTOUT;
Viaxis=Viedge(1:end-1)+(Viedge(2)-Viedge(1))/2;
figure;imagesc(Taxis, Veaxis,Ve);
figure;imagesc(Taxis,Viaxis,Vi);
figure;imagesc(rhoiOUT);
figure;imagesc(rhoeOUT);
%%
window=64;
tstart=150;
kwindow=window;
kaxis4=2*pi*kwindow/L*(-kwindow/2:kwindow/2)/kwindow;
kaxis=2*pi*NG/L*(-NG/2:NG/2)/NG;
omegaaxis=2*pi/DT/NTOUT*(-window/2:window/2)/window;
rhoeOUT4=(rhoeOUT(1:4:256,:)+rhoeOUT(2:4:256,:)+rhoeOUT(3:4:256,:)+rhoeOUT(4:4:256,:))/4;
figure;imagesc(omegaaxis,kaxis4,(abs(fftshift(fft2(rhoeOUT4(:,tstart:tstart+window))))));caxis([0 100])
xlabel('\omega/\omega_{pe}'), ylabel('k \times d_e')

figure;imagesc(omegaaxis,kaxis,log10(abs(fftshift(fft2(Eout(:,tstart:tstart+window))))));caxis([0.5 2])
xlabel('\omega/\omega_{pe}'), ylabel('k \times d_e')


%%

Ve1000=mean(Ve(:,1000:1100),2);
Ve1=mean(Ve(:,1:100),2);
Ve700=mean(Ve(:,700:800),2);
Ve500=mean(Ve(:,500:600),2);
Ve300=mean(Ve(:,300:300),2);
Ve300=mean(Ve(:,300:400),2);
Ve200=mean(Ve(:,200:300),2);
Ve100=mean(Ve(:,100:200),2);

%%
Ve1=mean(Ve(:,1:100),2);
Ve10000=mean(Ve(:,10000:15000),2);
Ve50000=mean(Ve(:,50000:55000),2);
Ve70000=mean(Ve(:,70000:75000),2);
Ve90000=mean(Ve(:,90000:95000),2);
figure;plot(Veaxis,Ve1);hold on;
plot(Veaxis,Ve10000);
plot(Veaxis,Ve50000);
plot(Veaxis,Ve70000);
plot(Veaxis,Ve90000);
xlabel('v/c');
