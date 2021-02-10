% us_map_diff=us_energy_diff(usObj,data,sStep,tStep)
% Find amplitude difference in time and space of ultrasonic signals
% 
% data: ultrasonic matrix s x t (space time)
% sStep: differennce spatial step 
% sStep: differennce time step 
% return: diffrece map according to the steps used

%examples:

% 5000 samples (space/distance)
% 1000 waves (100 waves p/sec) (time)
% tStep = 100 , sStep = 1000 
% result map= 5000/1000 x 1000/100 = 5x10 matrix
% each 1000 samples we make a difference between]
% us( 1:1000, 1)  - wave 1 minus
% us (1:1000,101)   wave 101 
% --------------
% mean abs difference of each point. => 1 pixel of 5x10 matrix


function us_map_diff=us_energy_diff(usObj,data,sStep,tStep)

%tStep=50;
%sStep=200;
tTotal=size(data,1);
sTotal=size(data,2);

sWindow=1:sStep;
sCurrent=1
sNext=sCurrent+sWindow;


tCurrent=1;
tNext=tCurrent+tStep;
tIter=1;
sIter=1;
exit=0;

while(sNext(end)<sTotal)
    %get difference in time
    sWindow(end)
    while(tNext<tTotal)
        dif=abs(data(tCurrent,sWindow))-abs(data(tNext,sWindow));
        us_map_diff(sIter,tIter)=mean(dif);
        tCurrent=tCurrent+tStep;
        tNext=tNext+tStep;
        tIter=tIter+1;
    end
    
    if(exit==1)
        break;
    end
    if((sWindow(end)+sStep)>=sTotal)
        sWindow=sWindow+ (sTotal-sWindow);
        exit=1;
    else
        sWindow=sWindow+sStep;
    end    
    tIter=1;
    sIter=sIter+1;
    tCurrent=1;
    tNext=tCurrent+tStep;
   
    
end
%%
figure
subplot 121

sPlot=size(us_map_diff,1)
tPlot=size(us_map_diff,2)

ss=linspace(1,sTotal,sPlot);
tt=linspace(1,tTotal,tPlot);

imagesc(ss,tt,abs(us_map_diff)')
colorbar
xlabel('space samples(depth)')
ylabel('time samples(fprf)')
title('diff plot')

subplot 122
imagesc(abs(data))
colorbar
title('original plot')
xlabel('space samples(depth)')
ylabel('time samples(fprf)')

end
