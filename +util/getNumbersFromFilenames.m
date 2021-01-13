%%%%%%%%%%%%%%%%%%%%%%%%
%input expressions
%

function values=getNumbersFromFilenames(str,varargin)   

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
    indexExp=indexExp+length(exp);
    strExp=cell2mat(matchStr(find(init==indexExp)));
    %% details of each implementation
    strExp=insertAfter(strExp,1,'.');
    values=[values strExp];
    

end