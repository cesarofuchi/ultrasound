%%
clc
close all
clear all
%monofasico
filepath='E:\dados\rockFlow\us\teste0\';
filename='2020-12-17 18_46_59_,636testeSync_1_5rad.xml';
%filename='2020-12-17 18_53_55_,939testeSync1_1_0rad.xml';
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
corte=20;
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
    USDO(i)=usObj.tofMaxPeak(usData(:,i),range1,range2,c, corte,0,sinalRef);    
end

%
for i=1:ndata
    distO(i)=USDO(i).dist;    
    maxValO(i)=USDO(i).max2;
    distSamples(i)=USDO(i).pos2;  
end

% compensate angle
do=distO*cosd(10);
dmean=movmean(do,11);
%d(find(d>0.024))=0;

%% see if the distance is correct
x=usObj.tofMaxPeak(usData(:,800),range1,range2,c, corte,0,sinalRef) 
%%
fi=figure
fi.Color='w';
inter=1:20000;
re2=range2(end)
range1(1)=50;
re1=range1(1);
re2=ceil((re2+re1)/2)

%re=ceil(re/2);
hold all

plot(1000*do(inter),'LineWidth',1)
plot(1000*dmean(inter),'LineWidth',1)
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
title([filename(36:end) ' novo=' num2str(mean(dmean*1000)) 'mm'])
title(['1.5 radianos/s filme médio=' num2str(mean(dmean*1000)) 'mm'])
xlabel('samples')
ylabel('distance(mm)')

%% -------------------
% mov average for smooth
distSamplesM=movmean(distSamples,15);
figure
plot(distSamplesM)
%% erase bubble positions in slug velocity data

%- flow velocity have windows 
% 1 - samples Ns 
bubblePosUndersampled=ceil(distSamples/ns);
% 2 - time Nc
bubblePos=bubblePosUndersampled(1:nc:end);

figure;plot(bubblePos)

filtFlowTemp=util.slugBubbleNanFilter(flow,bubblePos);

%cortar começo
filtFlow=filtFlowTemp(5:end,:)

%%
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

title('1.5 rad/s')

%%

