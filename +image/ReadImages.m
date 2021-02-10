function [TTImage]=ReadImages(path,image_type,frame_rate,varargin)

isCreationTime=1;
if(length(varargin)>0)
    isCreationTime=0;    
end
%clc
 %clear all
 %close all
% %path='O:\posDoc\03092019\set03092019\';
 %path='O:\posDoc\dados\19_50wc_200lh_75bar_4C_10000ml_30glNaCl_00glAA_NG_4\original\pvm\';
 %image_type='bmp';
% load files
if(path(end)~='\')
    path=[path '\'];
end
S=dir([path '*.' image_type]);
%order by creation date
S = S(~[S.isdir]);
[~,idx] = sort([S.datenum]);
S = S(idx);
pathFiles=S;
nfiles=size(pathFiles,1);
% filter images > 3 not black
i=1;
%b=0;

if(nfiles==0)
     warndlg(['No files with extension ' image_type ' found'])
     return;
end

for iterFile=1:nfiles
    files{i}=pathFiles(iterFile).name;
    %filename=[pathFiles(iterFile).folder '\' pathFiles(iterFile).name];
    %image=imread(filename); 
    index(i)=iterFile;
    filedate=pathFiles(iterFile).date;    
    if(iterFile==1)
        if(isCreationTime==1)
            time(i)=datenum(filedate,'dd-mmm-yyyy HH:MM:SS');
        else
            time(i)=varargin{1};
        end
        newtime(i)=time(i);
    else
        if(isCreationTime==1)
            time(i)=datenum(filedate,'dd-mmm-yyyy HH:MM:SS');
            if(time(i-1)==time(i))            
                newtime(i)=newtime(i-1)+(1/(frame_rate*86400));
            else
                newtime(i)=time(i);
            end
        else
            newtime(i)=newtime(i-1)+seconds(varargin{2}/(frame_rate));
        end
    end
    %index{i}=[pathFiles(iterFile).folder '\' pathFiles(iterFile).name];
    i=i+1;
    %else
    %    b=b+1;
    %end
end
str=datestr(newtime, 'YYYY-mm-DD hh:MM:ss.fff');
Time = datetime(str,'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

% dt=correctBR_FR_Timezone(dates(1));
% dates=dates+dt;
TTImage=timetable(Time,files');
%%
%strList=TTPVM.Var1;

%% code to plot the table
% 
% plot(TTPVM.Time,TTPVM.Var2)

