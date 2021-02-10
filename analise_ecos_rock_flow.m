%%
clear all
close all
clc

%
tic

filename='date_2020_10_28_time_19_17_05_,547_min_JL_max_JG_2k.xml';
filename='date_2020_10_28_time_20_06_25_,716_max_JL_min_JG_2k.xml';


%filename='date_2020_10_28_time_19_33_59_,273_min_JL_min_JG_1k.xml'


filepath='E:\dados\anular\ultrassom\'


filepath='E:\dados\anular\ultrassom\maxJLmaxJG\';
filename='date_2020_10_28_time_19_55_20_,829_max_JL_max_JG_1k.xml';

filepath='E:\dados\rockFlow\us\teste0\';
filename='2020-12-17 18_46_59_,636testeSync_1_5rad.xml';

filename='2020-12-17 18_53_55_,939testeSync1_1_0rad.xml';

filename='2020-12-17 18_59_16_,188testeSync0_0_5rad.xml';


file=[filepath filename];


%% read data

channels=1;
d = System.IO.File.GetCreationTime(file);
usObj=Ultrasonic.loadData(file,channels);

% show Ultrasonic Object
usObj
%%
temperatura=23;
c=util.water_c_temp(temperatura)

lambda=1000*c/usObj.fc
%crop the data of this experiment
fprf=usObj.fprf

%%
roi=1700:usObj.samples;
usData=usObj.data(roi,:);

%% uses mat file from mean full pipe file

tuboCheioFile=[filepath 'calib.mat'];
load(tuboCheioFile)
    %%
usDataMeanFullRoi=usDataMeanFull(roi);

wave=104;

usWave=usData(:,wave);
figure
hold all
plot(usWave)
plot(usDataMeanFullRoi)
legend(['wave=' num2str(wave)],'tubo cheio')

figure
hold all
plot(abs(hilbert(usWave)))
plot(abs(hilbert(usDataMeanFullRoi)))
legend(['wave=' num2str(wave)],'tubo cheio')


%%
figure
imagesc(abs(hilbert(usData)))

%% regiao de interesse de procura
% eco parede liquido
range1=148:200;

%parede - final do tubo
%10 mm em metros
dm=24e-3;
samplesRange=usObj.fs*2*dm/c

range2=range1(end):samplesRange;


%repmat(mean(data)',[1 size(data,1)])';
intervalo=1:usObj.samples;
t=(intervalo)*(c/2)/usObj.fs;
t=t*1000;
intervalo_plot=1:1000;
%  figure 
%  plot(usData)
% 
% figure
% plot(abs(usData(:,intervalo_plot)))
%% modo A

interval=range1(1):range2(end);
data=(usData(interval,:));
%datas=detrend(data);
datah=abs(hilbert(usData(interval,:)));
% %%
figure
imagesc(datah)
%%
samples=631-37;
ts=samples/usObj.fs;
dist=1000*c*ts/2

samples=711-37;
ts=samples/usObj.fs;
dtubo=1000*c*ts/2


%%
clc
fc=4e6
cycles=4;
% US=sensor.Ultrasonic
% US.fc=fc;
% US.fs=fs;
% US.cycles=cycles;

corte=20;
clear tt1
clear dist_zeros
ndata=size(usData,2);
range1=148:200; %second echo
range2=200:3301;%afterwards
%ndata=3000;
R=12.85;
D=25.7;
diamMin=3;
out=0.0784;

sinal=mean(usData(range1,:),2);
sinalRef=abs(hilbert(sinal));

%1005 0.003
tic
for i=1:ndata
    %function US= tof3(this,dados,periodo_de_interesse1,periodo_de_interesse2,vel_som, corte,out,sinalRef)
    USDO(i)=usObj.tofMaxPeak(usData(:,i),range1,range2,c, corte,0,sinalRef);    
end

%
for i=1:ndata
    distO(i)=USDO(i).dist;    
    maxValO(i)=USDO(i).max2;
end
% %%
% out=0.0033
% 
% for i=1:ndata
%     USDN(i)=US.tof3(usData(:,i),range1,range2,c, corte,out);    
% end
% 
% 
% for i=1:ndata
%     distn(i)=USDN(i).dist;    
%     maxValn(i)=USDN(i).max2;
% end
% toc
%

%%
do=distO*cosd(10);
%d(find(d>0.024))=0;
dmean=movmean(do,11);

%% dados em samples


for i=1:ndata
    distSamples(i)=USDO(i).pos2;  
    
end
%%
save([file(1:end-4) 'dist'],'distSamples')

%
% dn=distn*25.7/26;
% %d(find(d>0.024))=0;
% dn=movmean(dn,11);
% mean(dn)*1000
%
% figure
%  subplot 211
%  hold all
% % %yyaxis left
%  inter=1:3000;
%  plot(dn(inter)*1000)
%  plot(do(inter)*1000)
% % %plot(linspace(0.0125,0.0125,1500))
%  legend('novo','velho')
%  
% %yyaxis right
% %plot(maxVal(1:1500))
% ylim([0 26])
% grid on
% 
% subplot 212
% r1=50;
% 
% ts=(range1(1)-range1(1):range2(end)-range1(1))/fs;
% 
% y=1000*c*ts/2;
% x=inter;
% imagesc(x,y,(abs(usData(range1(1):range2(end),inter))))
% %imagesc((abs(usData(range1(1)+600:range2(end),inter))))
% set(gca,'YDir','normal') 
%  title([filename(1:end-15) ' novo=' num2str(mean(dn*1000)) 'mm velho=' num2str(mean(do*1000)) 'mm' ])
 
%%
%
fi=figure
fi.Color='w';
inter=1:20000;
re2=range2(end)
range1(1)=196;
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
% %%
% USDN=USDO
% figure
% iter=1830;
% subplot 211
% hold all
% t=1:size(usData,1);
% plot(t,abs(usData(:,iter)))
% plot(range1(1):range1(end),abs(usData(range1(1):range1(end),iter)),'g')
% plot(range2(1):range2(end),abs(usData(range2(1):range2(end),iter)),'r')
% 
% scatter(USDN(iter).pos1,abs(usData(USDN(iter).pos1,iter)))
% scatter(USDN(iter).pos2,abs(usData(USDN(iter).pos2,iter)))
% grid on
% title(num2str(iter))
% 
% subplot 212
% hold all
% iter=2644;
% t=1:size(usData,1);
% plot(t,abs(usData(:,iter)))
% 
% plot(range1(1):range1(end),abs(usData(range1(1):range1(end),iter)),'g')
% plot(range2(1):range2(end),abs(usData(range2(1):range2(end),iter)),'r')
% 
% 
% scatter(USDN(iter).pos1,abs(usData(USDN(iter).pos1,iter)))
% scatter(USDN(iter).pos2,abs(usData(USDN(iter).pos2,iter)))
% grid on
% title(num2str(iter))
% 
% %%
% figure
% iter=1845
% %iter=2644
% t=(10^6)*(1:size(usData,1)-1)%./fs;
% plot(t,abs(hilbert(usData(1:end-1,iter))))
% ylabel('Amplitude (V)')
% xlabel('time(\mus)')
% grid on
% %plot(range1(1):range1(end),abs(usData(range1(1):range1(end),iter)),'g')
% %plot(range2(1):range2(end),abs(usData(range2(1):range2(end),iter)),'r')
% 
% %i=224;
% %X=US.tof3(usData(:,i),range1,range2,c, corte,out);
% 
% %%
% figure
% 
% mEnd=mean(usData(6850:7050,:)');
% max(mEnd)
% plot(mEnd)
% figure
% mIni=mean(usData(3230:3330,:)');
% plot(mIni)
% max(mIni)
% %%
% figure
% m=detrend(usData(1:end-1,:));
% mEnd=mean(m,2);
% max(mEnd)
% plot(mEnd)
% %% ABRIR DADOS DO TUBO CHEIO
% tuboCheioFile='E:\dados\anular\ultrassom\tuboCheio.mat';
% load(tuboCheioFile)
% %%
% 
% 
% 
% meanCalib(3200:3400)=meanCalib(3400:3600);
% %%
% 
% figure
% iter=1845
% %iter=2644 %bom exemplo
% t=(10^6)*(1:size(usData,1)-1);%./fs;
% plot(t,abs(hilbert(usData(1:end-1,iter))))
% hold all
% %plot(t,abs(hilbert(meanCalib(1:end,1)*0.6)))
% d=usDataMeanFull(1:end-1,1);
% d(3200:3400)=d(3400:3600);
% off=0.0005;%mean(d(3200:3400));
% d=abs(hilbert(d*0.6));
% off=mean(d(3200:3400))
% d=d-off;
% plot(t,d)
% 
% ylabel('Amplitude (V)')
% xlabel('time(\mus)')
% grid on
% 
% 
% %% denoise test
% 
% figure
% iter=1845
% iter=1170 %bom exemplo
% hold all
% t=(10^6)*(1:size(usData,1)-1);%./fs;
% data=usData(1:end-1,iter);
% wdata1=wdenoise(data,1);
% %wdata5=mlptdenoise(data,t,5)
% wdata5=wdenoise(data,5);
% plot(t,data)
% plot(t,wdata1)
% plot(t,wdata5)
% grid on
% legend('or','1','2')
% %%
% hold all
% %plot(t,abs(hilbert(meanCalib(1:end,1)*0.6)))
% d=calib.usData(1:end-1,1);
% d(3200:3400)=d(3400:3600);
% wd=wdenoise(d,4);
% %d=abs(hilbert(d*0.6));
% %off=mean(d(3200:3400))
% %d=d-off;
% plot(t,d)
% plot(t,wd)
% 
% ylabel('Amplitude (V)')
% xlabel('time(\mus)')
% grid on
% %  figure
% %  data=detrend(usData(1:end-1,1243));
% %  plot(abs(data))
% %  grid on
% % %%
% % figure
% % hist(d)
% % 
% % plot(d2)
% % legend('or','m')
% % %%
% % figure
% % hold all
% % intervalo=1:5000;
% % plot(d(intervalo))
% % d2=movmean(d(intervalo),9)
% % plot(d2(intervalo))
% % 
% % grid on
% % legend('original','mov avg 9')
% % %dsadsa
% 
% % %%
% % plot(dist)
% % media=mean(dist);
% % 
% % %%
% % 
% % fprf=1000;
% % t_axis=(0:ndata-1)/fprf;
% % y_axis=1:8193;
% % intervalo=250;
% % vfilt=dist_zeros;%medfilt1(dist_zeros);
% % 
% % 
% % figure
% % intervalo=1:1000
% % yfilt=util.filtro_zeros(vfilt);
% % %plot(t_axis(intervalo),vfilt(intervalo)*1000)
% % plot(vfilt(intervalo))
% % media=mean(vfilt);
% % title(['mean=' num2str(media) ' mm'])
% % xlabel('time(s)')
% % grid on
% %%
% %plotresults_filme1(vfilt, usData,t_axis, y_axis, intervalo,c)
% 
% %plotresults_filme1(yfilt, usData,t_axis, y_axis, intervalo,c)
% % %
% % figure
% % subplot 211
% % %plot(vfilt(1:1000))
% % usDataA=abs(usData(range1(1):range2(end),:));
% % usDataF=flipud(usDataA);
% % 
% % imagesc(usDataF)
% % hold all
% % subplot 212
% % plot(yfilt)
% % %%
% % clc
% % fprf=1000;
% % time=(1:nwaves)/fprf;
% % clear data
% % data=[time',yfilt'];
% % T=table(time',yfilt','VariableNames',{'time','thickness(mm)'});
% % if(~isempty(vel))
% %    % a=strrep(vel{1},'.','');
% %     % b=strrep(vel{2},'.','');
% % else
% %     a='';
% %     b='';
% % end
% 
% % name=['JL' a 'JG' b 'filmCh1.csv'];
% writetable(T,[path name],'Delimiter',';')  
% csvwrite([path name],data);
%%
% save([path filename(1:8) '.mat'],'dn','out','do')

%save([path filename(1:8) '.mat'],'usData')

