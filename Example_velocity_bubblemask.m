

%% clear all data
clc
close all
clear all

%monofasico
filepath='E:\dados\rockFlow\testes_11_02\us\50oil_50air\175_degrees\';
filename='2021-02-11 18_59_39_,710_50oil_50air_175_degrees_1_rad.xml';
%filename='2021-02-11 18_54_51_,735_50oil_50air_175_degrees_05_rad.xml';
%filename='2021-02-11 19_01_29_,114_50oil_50air_175_degrees_15_rad.xml';
file=[filepath filename];

%% extract from the filename parameters of the plot titles
clc
str=filename;
out=util.getNumbersFromFilenames(filename,'back','rad','degrees');
rad=out{1};
degrees=out{2};
expTitle=[degrees ' degrees and ' rad ' rad/s']



%% load ultrasonic data
channels=1;
d = System.IO.File.GetCreationTime(file);
usObj=Ultrasonic.loadData(file,channels);
% show Ultrasonic Object
usObj


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Analysis
% in this step we evaluate A and B scan of ultrasonic data

%% evaluate the B-mode data
% in this step we evaluate the data for cropping afterwards
figure
imagesc(abs(hilbert(usObj.data(1900:3600,1:20000))))
xlabel('pulses')
ylabel('samples')

%% evaluate A-mode data
figure
hold all
plot(abs(hilbert(usObj.data(1500:end,14000))))
plot(abs(hilbert(usObj.data(1500:end,8000))))
grid on
xlabel('samples')
ylabel('amplitude')
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Selection
% in this step we select what part of data we are interested

usObj.workData=(usObj.data(1700:3300,1:10000));

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doppler Velocity Parameters
% in this step we fill doppler properties and verify its parameters

c=1410; %sound speed in the medium
clc
%ns=1*floor(2*usObj.cycles*usObj.fs/usObj.fc); %janela espacial de med. de vel.
dec=1; %decinmation
ns=32; %spatial window
ns_mm=1000*c*(ns/usObj.fs)/2;
nc=32; %temporal window
nc_msec=1000*nc/usObj.fprf;
v_nc=0.9;
nc_mm=v_nc*nc_msec;

dif=0; 
ovt=1; ovs=1; %override if needed

vmax=c*usObj.fprf/(4*usObj.fc);
pmax=c/(2*usObj.fprf);
angle=10;
disp('Size of data')
disp(size(usObj.workData))
disp(['Sound Speed in the media c: ' num2str(c) ' m/s'])
disp(['Spatial Window Ns: ' num2str(ns_mm) ' mm'])
disp(['Temporal Window Nc: ' num2str(nc_msec) ' msec'])
disp(['Spatial Window Nc considering v_nc=' num2str(v_nc) ' m/s: ' num2str(nc_mm) ' mm'])
disp(['Max Velocity for Auto Correlation: ' num2str(vmax) ' m/s'])
disp(['Transducer Angle: ' num2str(angle) ' degrees'])
disp(['Max Velocity for Auto Correlation angle component: ' num2str(vmax/sind(angle)) ' m/s'])
disp(['Max depth for fprf=' num2str(usObj.fprf) ' pulses/sec: ' num2str(pmax) ' meters'])


%% create doppler object
% Matched filter for signal enhancement
usdObj=UltrasonicDoppler(usObj,ns,nc,ovs,ovt)


%% velocity estimation using EAM
NpRange=[-2 -1 0 1 2];
EAM=usdObj.ExtAutoCorrelation(c,NpRange);

%% post-processing of velocity data
angle=10
vmax=2.0
vr=EAM.vel;
flowr=vr/sind(angle);
flow=(flowr);

%filter outliers
flow=medfilt2(flow,[3 3]);

figure
subplot 211
imagesc(abs(hilbert(usObj.workData)))
subplot 212
imagesc(flow)

% end of velocity code

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bubble Mask Code
% in this step we use time-of-fligth techniques to recognize liquid gas
% interfaces

%% evaluate A-mode data again but on the workdata
figure
subplot 211
%imagesc(flipud(abs(hilbert(usObj.workData))))
imagesc((abs(hilbert(usObj.workData))))
subplot 212
hold all
pos1=2000
pos2=5000
plot(abs(hilbert(usObj.workData(:,pos1))))
plot(abs(hilbert(usObj.workData(:,pos2))))
legend(num2str(pos1),num2str(pos2))
grid on
xlabel('samples')
ylabel('amplitude')

%% select appropriate threshould
clc
usData=usObj.workData;
ndata=size(usData,2);

thr=0.30; 
range1=1:200; %second echo
range2=200:1600;%afterwards

display(['range of the reference peak is: ' num2str(range1(1)) ' to ' num2str(range1(end)) ' samples'])
display(['range of the water gas peak is: ' num2str(range2(1)) ' to ' num2str(range2(end)) ' samples'])
display(['Amplitude Threshold is: ' num2str(thr)])


%%
%ndata=3000;
R=12.85;
D=25.7;
diamMin=3;

for i=1:ndata
    % using max values 
    USDO(i)=usObj.tofMaxPeak(usData(:,i),range1,range2,thr,'max');    
    %if not valid
    if(USDO(i).valid==1)
        distO(i)=c*USDO(i).tt./2;            
        maxValO(i)=USDO(i).max2;
        distSamples(i)=USDO(i).pos2;  
    else
        % if not valid use the NaN filter to eliminate the end
        % range2(end) is the last position (should be the pipe)
        tt=(range2(end)-USDO(i).pos1)./usObj.fs;
        distO(i)=c*tt./2;
        distSamples(i)=range2(end);  
    end
end


% compensate angle
do=distO*cosd(10);
display(['Compensate for angle: ' num2str(angle)])
dmean=movmean(do,11);
%d(find(d>0.024))=0;

%% see if the distance is correct
x=usObj.tofMaxPeak(usData(:,800),range1,range2,thr) 
%% plot to compare B-mode and the distance found
fi=figure
fi.Color='w';
inter=1:10000;

%choose initial B-mode
%should match the peak of the reference echo
re1=50;
%choose end of B-mode
re2=ceil((re2+re1)/2); %center
re2=re2+500;

display(['B-mode for samples:' num2str(re1) ' to ' num2str(re2) ' samples'])

hold all

plot(1000*do(inter),'LineWidth',1)
plot(1000*dmean(inter),'LineWidth',1)
legend('original','mov avg 11')

ts=(range1(1)-range1(1):re2-range1(1))/usObj.fs;
y=1000*c*ts/2;
x=inter;
im=imagesc(gca,x,y,(abs(usData(range1(1):re2,inter))));

axis tight
im.AlphaData= .5;
ax  = gca;
%plot(timerange,repmat(mean(FD.mean),1,2),'--b');


title(['B-mode x Distance - Mean Thickness=' num2str(mean(dmean*1000)) 'mm'])
xlabel('samples')
ylabel('distance(mm)')

%% -------------------
% mov average for smooth

distSamplesM=movmean(distSamples,15);
figure
plot(distSamplesM)

% end of distance processing

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Bubble Mask Code
% in this step undersample the distance vector to match velocity map (as it
% was reduced based on Ns and Nc)

%- flow velocity have windows 
% 1 - samples Ns 
bubblePosUndersampled=ceil(distSamples/ns);
% 2 - time Nc
bubblePos=bubblePosUndersampled(1:nc:end);
distanceVel=dmean(1:nc:end)*1000;

figure;plot(bubblePos)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Bubble Mask Code
% in this step we use a NanFilter for the data after the measured distance 
filtFlowTemp=util.slugBubbleNanFilter(flow,bubblePos);

%cortar começo
filtFlow=filtFlowTemp(1:end,:)

time_t=size(filtFlow,2)*nc/usObj.fprf;
x=linspace(0,time_t,size(filtFlow,2));
y=linspace(0,25,size(filtFlow,1));

figure
subplot 211
imagesc(x,y,flipud(flow))
h = colorbar;
set(get(h,'label'),'string','velocity(m/s)');
colormap jet
%caxis([-0.75 0.75])
xlabel('time(s)')
ylabel('distance(mm)')

 subplot 212
imagesc(x,y,flipud(filtFlow))
h = colorbar;
set(get(h,'label'),'string','velocity(m/s)');
colormap jet
%caxis([-0.75 0.75])
xlabel('time(s)')
ylabel('distance(mm)')

title('7.5 degree - 1 rad/s')
title(expTitle)

%%
figure
imagesc(x,y,flipud(filtFlow))
set(gcf, 'Position',  [100, 200, 800, 200])
h = colorbar;
set(get(h,'label'),'string','velocity(m/s)');
colormap jet

caxis([-0.75 0.75])
xlabel('time(s)')
ylabel('distance(mm)')
title(expTitle)
savefig([filepath filename(1:end-4)])

