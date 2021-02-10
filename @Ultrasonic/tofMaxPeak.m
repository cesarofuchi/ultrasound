function US= tofMaxPeak(usObj,data,range1,range2,c, thr,out,sinalRef)
%tofMaxPeak Cálculo de tempo de trânsito por amplitude

%   data = A-scan signal
%   range1: range where the 1st echo reference is
%   range2: range where the desired echo may appear
%   c : speed of sound at the medium
%   thr: minimum threshould for the desired echo. If lower, return
%        US.valid=0
%   corte: % de corte desejado, caso não encontre nada, retorna zero%
%   US: strutcture with information
%
%   Author(s): Ofuchi, C.Y.            
%   $Revision: 1.0 $  $Date: 02/11/20$

%----------------------------------------------------------
% Info: simple technique to find the time of flight using the max echo peak
% 
%----------------------------------------------------------           
    offset1=range1(1);
    offset2=range2(1);

    % achar os máximos
    %[max1, posicao_eco1] = max(abs(hilbert(dados(range1,1))));
    [max1, posicao_eco1] = max(sinalRef);
    posicao_eco1=posicao_eco1+offset1-1;

    atual=abs(hilbert(data(range2,1)));            
    z=length(atual)-length(sinalRef);
    sinalRef(end+1:end+z)=zeros(1,z);
    sinalSub=atual-sinalRef;
%-debug--------------------------------------------
%      figure
%      hold all
%      plot(atual)
%      plot(sinalRef)
%      plot(sinalSub)
%      legend('atual','ref','sub')
%--------------------------------------------------

    [max2, posicao_eco2] = max(sinalSub);
    US.valid=1;
    if(max2<thr)
        US.valid=0;
    end
    %percent=max2/max1*100;          
    US.pos1=posicao_eco1;
    US.pos2=posicao_eco2+range2(1);
    US.tt=(US.pos2-US.pos1)/usObj.fs;
    US.max2=max2;    
    US.dist=c*US.tt/2;    
    
end