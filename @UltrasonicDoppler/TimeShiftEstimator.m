function CCM = TimeShiftEstimator(usdObj,c)
% Estima a velocidade radial através do atraso de tempo dos pulsos
% d matriz de dados m pontos x n pulsos)
% c velocidade do som no meio
% prf : frequência das emissões de pulsos
% fs : frequencia amostragem
% Ns : número de canais de profundidade
% Nc: número de linhas para realizar a média
% Irs: Indice de confiabilidade da estimativa
% vmax=(c/2)*(prf/fs)*Ns;
% tsmax=vmax*(2/c)*(1/prf);
% lagmax=fs*tsmax;

fprf = usdObj.usObj.fprf;
fs = usdObj.usObj.fs;

Ns = usdObj.ns(k);
Nc = usdObj.nc(j);
ovs = usdObj.ovs;
ovt = usdObj.ovt;

disp('TimeShiftEstimator running');
tic;

% adapta os num de amostras espaciais para ser janelado com ovs
resto = rem(size(usdObj.usObj.samples,1)-Ns,Ns*ovs);
if resto ~= 0
    pad_n = uint8(Ns*ovs-rem(usdObj.usObj.samples-Ns,Ns*ovs));
    pad = zeros(pad_n,size(usdObj.usObj.samples,2));
    d = [d ; pad];
end
% zero pads temporal
pad = rem(size(d,2),Nc);
d = d(:,1:end-pad);

nchannels = floor((size(d,1)-Ns)/(ovs*Ns))+1;
% calcula o número de canais temporais com ovserflow
nchannelst = floor((size(d,2)-Nc)/(ovt*Nc))+1;

% inicializa com zeros as matrizes de resultados de saída do algoritmo
vr = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial
Irs = zeros(nchannels,nchannelst);% mapa espaço-temporal confiabilidade
corr = zeros(2*Ns+1,Nc-1);

for i = 1 : nchannelst % varre todos os canais temporais
    p_i = (i-1)*Nc*ovt+1; % pulso inicial
    p_f = p_i+Nc-1; % pulso final
    
    for j = 1:nchannels % varre todos os canais espaciais
        ps_i = (j-1)*Ns*ovs+1;
        
        if j == nchannels
            dataf = [d(ps_i-Ns:ps_i+Ns-1,p_i:p_f); zeros(Ns,Nc)];
        else
            if j == 1
                dataf = [zeros(Ns,Nc);d(ps_i:ps_i+2*Ns-1,p_i:p_f)];
            else
                dataf = d(ps_i-Ns:ps_i+2*Ns-1,p_i:p_f);
            end
        end
       
        for t = 1 : Nc-1
%             w=1;
%             corr(Ns+1,t)=dataf(1:Ns,t+1)'*data(1:Ns,t);
             for l = 0 : Ns
                if l > 0
                 norm(Ns-l+1) = sqrt(dataf(1+l:end,t+1)'*dataf(1+l:end,t+1)...
                    *(dataf(1:end-l,t)'*dataf(1:end-l,t)));
               end
                 norm(Ns+l+1)=sqrt(dataf(1+l:end,t)'*dataf(1+l:end,t)...
                    *(dataf(1:end-l,t+1)'*dataf(1:end-l,t+1)));
%                 corr(w+Ns+1,t)=dataf(1+lag:lag+Ns,t)'*data(1:Ns,t+1);
%                 w=w+1;
             end
%             [corr(:,t),lag]=xcorr([zeros(Ns,1); data(:,t); zeros(Ns,1)],dataf(:,t+1));
            [corr(:,t),lag] = xcorr(dataf(:,t),dataf(:,t+1),Ns);
            %max(corr(:,t))
            %dataf(1:end+lag(ind),t)
            corr(:,t) = corr(:,t)./norm'; % normaliza
        end
        c_corr = zeros(size(corr,1),1);
        for j2 = 1 : Nc-1
            if max(corr(:,j2)) > 0.9
                c_corr = c_corr+corr(:,j2);
            end
        end
        %c_corr=sum(corr,2)/(Nc-1);

        [Irs(j,i),ind] = max((c_corr));
        if ind > 1 && ind < (2*Ns+1)
            interp = ( c_corr(ind+1)-c_corr(ind-1) ) / ...
             ( (c_corr(ind+1)-2*c_corr(ind)+c_corr(ind-1)) );
            vr(j,i) = -(c/2)*((lag(ind)-interp)/fs)*fprf;
        else
            vr(j,i) = -(c/2)*((lag(ind))/fs)*fprf;
        end  
    end
end

% vr = vr.*media_msk;

CCM.processing_time = toc;
CCM.vel = vr;
CCM.Irs = Irs;

end
