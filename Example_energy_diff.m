clear all; close all; clc;

%% Load folder data
folder='E:\dados\celula_pressao_lino\EXP\validacao\Com US\SF6 - 288 K\US\teste3\';

channels=1;% to read use the xml path
usObj=Ultrasonic.loadFolder(folder,channels);

%% verify data
usObj.us_energy_diff(usObj.data,100,100)


