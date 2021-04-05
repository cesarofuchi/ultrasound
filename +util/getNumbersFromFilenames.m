%%%%%%%%%%%%%%%%%%%%%%%%
% Find intengers values in strings
% str= temperature_10_wind_45_angle_23.xml
% in this example, values are in the front
% getNumbersFromFilenames(str,'front','temperature')  
% returns: {cell array 1.0}

function values=getNumbersFromFilenames(str,numberPos,varargin)   
str(str=='.')=[]
nargs=length(varargin);
% find numbers in string
expression = '\d+';
[matchStr,init,fim]= regexpi(str,expression,'match');
values={};
for i=1:nargs
    exp=varargin{i};
    indexExp=findstr(str,exp);
    if(isempty(indexExp))
        warning([exp ' : not found'])        
        return;
    end
    
    switch lower(numberPos)
          case 'front'
              indexExp=indexExp+length(exp);
          case 'back'
              indexBack=init(init<indexExp);
              indexExp=indexBack(end);
          otherwise
              error(['Unexpected option: ' numberPos])
    end    
   
    strExp=cell2mat(matchStr(find(init==indexExp)));
    %% details of each implementation
    if(length(strExp)>1)
        strExp=insertAfter(strExp,strExp(end-1),'.');   
    end
    values=[values strExp];
    

end