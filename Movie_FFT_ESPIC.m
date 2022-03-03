
projectpath='/Users/shuzhang/Workspace/LsrCapCoil20A/TSdata'
movievariable='Rhoe';
filename=[simulationname,movievariable];

window=64;
start=1; lastframe=length(rhoeOUT(1,:))-window;
Eout4=(Eout(1:4:NG,:)+Eout(2:4:NG,:)+Eout(3:4:NG,:)+Eout(4:4:NG,:))/4;

kwindow=size(Eout4,1);
kaxis4=-2*pi*kwindow/L*(-kwindow/2:kwindow/2)/kwindow;% I believe the yaxis is wrong in 2D fft
omegaaxis=2*pi/DT/NTOUT*(-window/2:window/2)/window;
Spectrumtime=zeros(lastframe,window+1);




clear F;
F(lastframe) = struct('cdata',[],'colormap',[]);
for timenumber=start:lastframe
    tstart=timenumber;
    rhoeOUT4=(rhoeOUT(1:4:NG,:)+rhoeOUT(2:4:NG,:)+rhoeOUT(3:4:NG,:)+rhoeOUT(4:4:NG,:))/4;
    Eout4=(Eout(1:4:NG,:)+Eout(2:4:NG,:)+Eout(3:4:NG,:)+Eout(4:4:NG,:))/4;
    PowerSpectrum=abs(fftshift(fft2(rhoeOUT4(:,tstart:tstart+window))));
    Spectrumtime(tstart,:)=PowerSpectrum(10,:)+PowerSpectrum(11,:)+PowerSpectrum(12,:);
    figure(100);imagesc(omegaaxis,kaxis4,log10(PowerSpectrum));
    colorbar;
    axis xy;
    caxis([0 2]);
    xlabel('\omega/\omega_{pe}'), ylabel('k \times d_e');
    title([num2str((timenumber+window/2)*NTOUT*DT),'\omega_{pe}^{-1}']);
  
    if(timenumber==start)
        w = waitforbuttonpress
    end
    F(timenumber)=getframe(gcf);
end
F(1:start-1)=[];
movie(F,1)
%%
myVideo=VideoWriter([projectpath,'/',filename,'.mp4'],'MPEG-4');
myVideo.FrameRate = 60;  % Default 30
myVideo.Quality = 95;    % Default 75
open(myVideo);
writeVideo(myVideo, F);
close(myVideo)

%%
figure;plot(Taxis(window/2:end-window/2-1),mean(Spectrumtime(:,34),2),Taxis(window/2:end-window/2-1),mean(Spectrumtime(:,38:48),2))