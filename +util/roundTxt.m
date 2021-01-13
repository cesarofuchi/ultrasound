%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% função para arredondar os valores
%
% function roundTxt(num,casas)
% num:   numero a ser arredondado
% casas: quantas casas serão necessárias    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author(s): Ofuchi, C.Y.
%   Copyright 2017 The MathWorks, Inc.
%   $Revision: 1.0 $  $Date: 2017/03/06 11:47:07 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str=roundTxt(num,casas)
v=round(num*10^casas)/10^casas;
str=num2str(v);
end