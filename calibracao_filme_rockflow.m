clc
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the mean value of a full liquid pipe
%

filepath='E:\dados\rockFlow\us\teste0\'

filename='date_2020_10_28_time_18_35_06_,341_tubo_vazio_2k.xml';
filename='2020-12-17 19_03_53_,372tuboCheio.xml';

file=[filepath filename];


%% read data

channels=1;
d = System.IO.File.GetCreationTime(file);
% to read use the xml path
usObj=Ultrasonic.loadData(file,channels);
usObj.fc=usObj.fc*1e6

%%
% show Ultrasonic Object
usObj

%% crop the region of interest of this experiment
usData=double(usObj.data);

usDataMeanFull=mean(usData,2);
figure
hold all
plot(usDataMeanFull)
%y=usData(:,1);
%plot(y)
%plot(teste(1:2:40000))

%% TUBO VAZIO

filepath='E:\dados\rockFlow\us\teste0\'
filename='2020-12-17 19_05_42_,612tuboQuaseVazio.xml';
file=[filepath filename];


%% read data

channels=1;
d = System.IO.File.GetCreationTime(file);
% to read use the xml path
usObj=Ultrasonic.loadData(file,channels);
usObj.fc=4e6;
usObj.fprf=2000;
usObj.channels=1;
% show Ultrasonic Object
usObj

%% crop the region of interest of this experiment
usData=double(usObj.data);

usDataMeanVoid=mean(usData,2);
figure
hold all
plot(usDataMeanVoid)
plot(usDataMeanFull)
legend('void','full')

figure
hold all
plot(abs(hilbert(usDataMeanVoid)))
plot(abs(hilbert(usDataMeanFull)))
legend('void','full')
%%
save([filepath 'calib'],'usDataMeanFull','usDataMeanVoid')

