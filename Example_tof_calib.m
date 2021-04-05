clear all; close all; clc;

%% Load folder data
folder='E:\dados\rockFlow\oil_calib\us\subida\';

channels=1;% to read use the xml path
filter='*.xml';
usObj=Ultrasonic.loadFolder(folder,channels,filter);

%%
usObj.fprf=0.5;


[tSamples, sSamples]=size(usObj.data)

t=(1:tSamples)/(usObj.fprf*60);
s=1:sSamples;

figure;
imagesc(s,t,abs(usObj.data))
ylabel('minute')


%%
figure
hold all
plot(usObj.data(1,:))

plot(usObj.data(end,:))

%%
range1=120:190;
range2=2000:2800;

figure
plot(usObj.data(1,:))

%%
thr=0.7;
%US=usObj.tofMaxPeak(usObj.data(1,:)',range1,range2,thr,'max')
US=usObj.tofMaxPeak(usObj.data(1,:)',range1,range2,thr,'first_peak')

%% for test

for i=1:tSamples
    usd_vector(i)=usObj.tofMaxPeak(usObj.data(i,:)',range1,range2,thr,'first_peak');
    tof_vector(i)=usd_vector(i).tt;
end

%%
figure

plot(t,tof_vector)
xlabel('minutes')
%%
folder='E:\dados\rockFlow\oil_calib\us\subida\';
filter_type='*.xml'
files = dir([folder filter_type]);

filename = files(1).name(1:end-4)
f = filename(1:19);

aux = strrep(f,'_','-');
dateIni = datetime(aux,'InputFormat','yyyy-MM-dd HH-mm-ss')
dateIni=dateIni-hours(1)
%%

distance=33/100;
c_oil=2*distance./tof_vector;

figure
plot(c_oil)
%%

TTUS=timetable(c_oil','SampleRate',0.5);
TTUS.Time = TTUS.Time + dateIni;


%% LVM calibration

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [23, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["LabVIEWMeasurement", "VarName2", "Temperature", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName9", "EmptyFieldRule", "auto");

% Import the data
TT = readtable("E:\dados\rockFlow\oil_calib\RockingCellNova_21-02-10_1328_001.lvm", opts);
%% Clear temporary variables
clear opts

Y=2021
M=2
D=10
H=13
MI=28
S=45
lvmTimeIni = datetime(Y,M,D,H,MI,S);

TTLvm=timetable(TT.Temperature,'SampleRate',2);
TTLvm.Time=TTLvm.Time+lvmTimeIni;


%% sincronia
LocalTT=synchronize(TTUS,TTLvm,'commonrange','nearest');
LocalTT.Properties.VariableNames = {'C','Temp'};
%LocalTT.Properties.VariableNames{'C'} = 'c oil m/s';
%LocalTT.Properties.VariableNames{'Temp'} = 'temperature ºC';

%%
figure
stackedplot(LocalTT)
grid on
%%
figure

plot(LocalTT.C,LocalTT.Temp)
title('oil c')
xlabel('Sound Speed (m/s)')
ylabel('Temp')
grid on
%%
list=[541,4030,7058,9965,12732,15849,18698,21505]
%%


plot(LocalTT.C(list),LocalTT.Temp(list))
grid on
title('Oil Sound Speed x Temperature')
xlabel('Sound Speed (m/s)')
ylabel('Temperature ºC ')
grid on