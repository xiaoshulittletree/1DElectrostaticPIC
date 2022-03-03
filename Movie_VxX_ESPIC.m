
projectpath='/Users/shuzhang/Workspace/LsrCapCoil20A/TSdata'
movievariable='Vex';
filename=[simulationname,movievariable];

start=1; lastframe=size(pexOUT,3);
%Veaxis=Veedge(1:end-1)+(Veedge(2)-Veedge(1))/2;
xaxis=(1:NG)/NG*L;





clear F;
F(lastframe) = struct('cdata',[],'colormap',[]);
switch movievariable
    case 'Vex'
        for timenumber=start:lastframe
            figure(101);imagesc(xaxis,Veaxis,pexOUT(:,:,timenumber));    colorbar;
            title([num2str(timenumber*NTOUT*DT),'\omega_{pe}^{-1}']);
            xlabel('X (de)'), ylabel('Vx (c)');
            
            if(timenumber==start)
                w = waitforbuttonpress;
            end
            F(timenumber)=getframe(gcf);
        end
    case 'Vix'
        for timenumber=start:lastframe
            figure(101);imagesc(xaxis,Viaxis,pixOUT(:,:,timenumber));    colorbar;
            title([num2str(timenumber*NTOUT*DT),'\omega_{pe}^{-1}']);
            xlabel('X (de)'), ylabel('Vx (c)');
            
            if(timenumber==start)
                w = waitforbuttonpress;
            end
            F(timenumber)=getframe(gcf);
        end
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