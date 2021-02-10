%%
clc
close all
clear all
%monofasico
filepath='E:\dados\rockFlow\us\teste0\';
filename='2020-12-17 18_46_59_,636testeSync_1_5rad.xml';
filename='2020-12-17 18_53_55_,939testeSync1_1_0rad.xml';
%filename='2020-12-17 18_59_16_,188testeSync0_0_5rad.xml';

file=[filepath filename];
%%
channels=1;
d = System.IO.File.GetCreationTime(file);
usObj=Ultrasonic.loadData(file,channels);
% show Ultrasonic Object
usObj
%%
%figure
%imagesc(abs(hilbert(usObj.data(1900:3600,1:100:40000))))

%% data transformation step
% crop
% flip
usObj.workData=(usObj.data(1700:3300,:));
 
%% fill doppler properties

c=1400; %chute
clc
ns=1*floor(2*usObj.cycles*usObj.fs/usObj.fc); %janela espacial de med. de vel.
dec=1; 
ns=32;
nc=32; %janela temporal de medição de velocidade
dif=0; 
ovt=1; ovs=1;
% create doppler object
usdObj=UltrasonicDoppler(usObj,ns,nc,ovs,ovt)


%% velocity estimation
NpRange=[-2 -1 0 1 2];
EAM=usdObj.ExtAutoCorrelation(c,NpRange);

%% convert angle
angle=10
vmax=2.0
vr=EAM.vel;
flowr=vr/sind(angle);
flow=(-flowr);

%filter outliers
flow=medfilt2(flow,[3 3]);

figure
subplot 211
imagesc(abs(hilbert(usObj.workData)))
subplot 212
imagesc(flow)

%% calculate tof to make a bubble mask

clc
% threshould for a valid echo signal reflection
thr=1500;

usData=usObj.workData;
ndata=size(usData,2);
range1=1:200; %second echo
range2=200:1600;%afterwards
%ndata=3000;
R=12.85;
D=25.7;
diamMin=3;


% mean signal (use maybe full or void)
sinal=mean(usData(range1,:),2);
sinalRef=abs(hilbert(sinal));

for i=1:ndata
    %function US= tof3(this,dados,periodo_de_interesse1,periodo_de_interesse2,vel_som, corte,out,sinalRef)
    USDO(i)=usObj.tofMaxPeak(usData(:,i),range1,range2,c, thr,0,sinalRef);    
end

%
for i=1:ndata
    if(USDO(i).valid==1)
        distO(i)=USDO(i).dist;    
        maxValO(i)=USDO(i).max2;
        distSamples(i)=USDO(i).pos2;  
    else
        % neste caso, o filtro NaN vai eliminar o final
        distSamples(i)=range2(end);  
    end
end



%%
% compensate angle
ds=distSamples*cosd(10);
% mov average for smooth
distSamplesM=movmean(ds,11);


do=distO*cosd(10);
distOM=movmean(do,11);

%% see if interface is correct
  

fi=figure
fi.Color='w';
inter=1:1500;
re2=range2(end)
range1(1)=50;
re1=range1(1);
re2=ceil((re2+re1)/2)

%re=ceil(re/2);
hold all

plot(1000*do(inter),'LineWidth',1)
plot(1000*distOM(inter),'LineWidth',1)
legend('original','média móvel 11')

ts=(range1(1)-range1(1):re2-range1(1))/usObj.fs;
y=1000*c*ts/2;
x=inter;
im=imagesc(gca,x,y,(abs(usData(range1(1):re2,inter))));

axis tight
im.AlphaData= .5;
ax  = gca;
%plot(timerange,repmat(mean(FD.mean),1,2),'--b');
%title([filename(1:end-15) ' novo=' num2str(mean(dn*1000)) 'mm velho=' num2str(mean(do*1000)) 'mm' ])
title([filename(36:end) ' novo=' num2str(mean(distOM*1000)) 'mm'])
title(['1.5 radianos/s filme médio=' num2str(mean(distOM*1000)) 'mm'])
xlabel('samples')
ylabel('distance(mm)')
%% if not OK debug

figure
imagesc(flipud(abs(hilbert(usData(:,1:2000)))))
%%
figure
subplot 311
plot(abs(hilbert(usData(:,200))))
subplot 312
plot(abs(hilbert(usData(:,700))))
subplot 313
plot(abs(hilbert(usData(:,1600))))



%% erase bubble positions in slug velocity data

%- flow velocity have windows 
% 1 - samples Ns 
bubblePosUndersampled=ceil(distSamplesM/ns);
% 2 - time Nc
bubblePos=bubblePosUndersampled(1:nc:end);

filtFlowTemp=util.slugBubbleNanFilter(flow,bubblePos);
%cortar começo
filtFlow=filtFlowTemp(5:end,:)

%% show and compare new velocity map
x=linspace(0,10,size(filtFlow,2));
y=linspace(0,25,size(filtFlow,1));

figure
subplot 211
imagesc(x,y,flipud(flow))
h = colorbar;
set(get(h,'label'),'string','velocity(m/s)');
colormap jet
caxis([-0.75 0.75])
xlabel('time(s)')
ylabel('distance(mm)')

subplot 212
imagesc(x,y,flipud(filtFlow))
h = colorbar;
set(get(h,'label'),'string','velocity(m/s)');
colormap jet
caxis([-0.75 0.75])
xlabel('time(s)')
ylabel('distance(mm)')

%title('1.5 rad/s')


%% get image from the video directly
videopath = 'E:\dados\rockFlow\camera\';
videoname = '2020-12-17 17-46-58-1_5_rad.mkv';
videoname = '2020-12-17 17-53-54.mkv';
%videoname= '2020-12-17 17-59-15.mkv';
videoSource=VideoReader([videopath videoname])


%%

delete([videopath videoname(1:end-4) '_new3.mp4'])
myVideo = VideoWriter([videopath videoname(1:end-4) '_new'],'MPEG-4');
myVideo.FrameRate = 7;  %can adjust this, 5 - 10 works well for me
open(myVideo)
%
flowMap=-flipud(filtFlow);
step=1;
vInit=1;
usWindow=ceil(500/32);
vEnd=vInit+usWindow;
%%
y=linspace(0,25,size(filtFlow,1));

fprf=63
figure
x0=500;
y0=500;
width=1200;
height=500
set(gcf,'position',[x0,y0,width,height])
offsetCamera=5;
iterFrame=offsetCamera+ceil(vInit/2);
videoFrame = read(videoSource,iterFrame);
%cam=TTImg.Var1(iterFrame);
%videoFrame=imread([IMAGE.path char(cam)]);  
imgR = videoFrame(180:400,size(videoFrame,2)/2+250:size(videoFrame,2)-350,:); % hor


iterCount=1;
iterSync=1;

%se 30 -> congela 1
%se 60 -> congela 2

fileCell='E:\dados\rockFlow\Cell.jpg';

imCell=imread(fileCell);

fps=30;
revolutionTime=5;
revolutionSamples=fps*revolutionTime;
%velocity samples
revolutionSamples=200;

maxAngle=10;
minAngle=-10;
angle=linspace(minAngle,maxAngle,revolutionSamples)
initVal=5;

iterAngle=1;

%or a specific init
[ d, iterAngle ] = min( abs( angle-initVal ) );
% subindo ou descendo?
iterAngle=154;
reverse=1;


for i=1:600%size(flowMap,2)
     a=round(vInit/fprf,2);
     b=round(((vInit+vEnd)/(2*fprf)),2);
     c=round(vEnd/fprf,2);
     xx=[a b c];
    
    x=linspace(a,c,usWindow);
    as=sprintf('%.2f', a);
    bs=sprintf('%.2f', b);
    cs=sprintf('%.2f', c);
    xl={as, bs, cs};
    subplot(2,5,[1 2 3 4])        
    imagesc(x,y,flowMap(:,vInit:vEnd))
    ax=gca;
    ax.XTickMode = 'manual';
    ax.XTick=xx;
    ax.XTickLabel=xl;
    
    h = colorbar;
    set(get(h,'label'),'string','vel(m/s)');
    colormap jet
    caxis([-0.75 0.75])
    xlabel('time(s)')
    ylabel('distance(mm)')
    subplot(2,5,5)
    
    plot(flipud(mean(flowMap(:,vInit:vEnd),2)),y)
    title(['mean v:' util.roundTxt(mean(mean(flowMap(:,vInit:vEnd),'omitnan'),'omitnan'),2) 'm/s'])
    ylim([0 25])
    xlim([-0.8 0.8])    
    xlabel('vel(m/s)')
    grid on
    
    
     
    subplot(2,5,[6 7 8 9])
   
    %cam=TTImg.Var1(iterFrame);
    %videoFrame=imread([IMAGE.path char(cam)]);  
    %
    hold on
    imshow(imgR)    
    hold on

    %# define points (in matrix coordinates)
    % (x,y)
    p1 = [10,320];
    p2 = [180,350];
    
    plot([p1(2),p2(2)],[p1(1),p2(1)],'Color','r','LineWidth',1)
    
    p1 = [10,360];
    p2 = [170,388];
    
    plot([p1(2),p2(2)],[p1(1),p2(1)],'Color','r','LineWidth',1)
    
    subplot(2,5,10)
    
    imr=imrotate(imCell,angle(iterAngle));
    imshow(imr)
    title(['angle=' util.roundTxt(angle(iterAngle),2) 'º'])
    
    if(reverse==0)
        iterAngle=iterAngle+1
    else
        iterAngle=iterAngle-1
    end
    
    if(angle(iterAngle)==maxAngle)
        reverse=1;
    end
    if(angle(iterAngle)==minAngle)
        reverse=0;
    end    
    
    
       
    


    
    %rectangle('Position',[320,1,70, 100],'Edgecolor', 'r');
    %%
    vInit=vInit+step;
    vEnd=vEnd+step;
    iterCount=iterCount+1;
    %le a imagem a cada dois counts
    if(iterCount>2)
        iterCount=1;
        iterFrame=iterFrame+1
        videoFrame = read(videoSource,iterFrame);
        if(iterSync==30)
            iterFrame=iterFrame-1;
            iterSync=1;        
        end
        imgR = videoFrame(180:400,size(videoFrame,2)/2+250:size(videoFrame,2)-350,:); % hor
    end

    if(vEnd>size(flowMap,2))
        return;
    end
    pause(.01)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
end
close(myVideo)

%imagesc(flowMap(:,vInit:vEnd))

%%