function EAM = ExtAutoCorrelation(usdObj,dbubble,c,dif, NpRange, j, k)
% EAM = Extended Autocorrelation Method
% Estima a velocidade radial através do desvio de fase do sinal e 
% combinado com a correlação cruzada (Time shifting method)
% iq matriz de dados m pontos x n pulsos)
% c velocidade do som no meio
% prf : pulse repetition frequency, frequência das emissões de pulsos
% fc : frequencia central do transdutor
% channel : número de canais de profundidade
% ov: sobreposição de canais 
% debug: plota os resultados para conferência

fprf = usdObj.usObj.fprf;
fs   = usdObj.usObj.fs;
fc   = usdObj.usObj.fc;
iq   = usdObj.iq;
Ns   = usdObj.ns(k);
Nc   = usdObj.nc(j);
ovs  = usdObj.ovs;
ovt  = usdObj.ovt;

RMS = std(iq(:,5));

% zero pads
% adapta os num de amostras espaciais para ser janelado com ov
resto = rem(size(iq,1) - Ns,Ns*ovs);
if resto ~= 0
    pad_n = uint8(Ns*ovs-rem(size(iq,1)-Ns,Ns*ovs));
    pad = zeros(pad_n,size(iq,2));
    iq = [iq ; pad];
    dbubble = [dbubble ; pad]; % update - detecção da bolha
end

% zero pads
% adapta os num de amostras espaciais para ser janelado com ov
% pad_n=Ns*ovs-rem(size(iq,1)-Ns,Ns*ovs);
% pad=zeros(pad_n,size(iq,2));
% 
% iq=[iq ; pad]; % acrescenta zeros na profundidade no fim do vetor

% calcula o numero de canais espaciais já com o pad de zeros no final e
% overflow
if ovs ~= 1
    nchannels = uint16(size(iq,1)/(ovs*Ns)-1);
else
    nchannels = uint16(size(iq,1)/(ovs*Ns));
end
% calcula o número de canais temporais com overflow
nchannelst = floor((size(iq,2)-Nc)/(ovt*Nc))+1;

% inicializa a matriz de velocidade de saída do algoritmo
vr = zeros(nchannels,nchannelst);
desv = zeros(nchannels,nchannelst);

% Np = [-2 -1 0 1 2];
Np = NpRange;

% ccross=zeros(Np,M-2);
% tic % inicia contagem de tempo de execução
for i = 1 : nchannelst % varre todos os canais temporais
    p_i = (i-1)*Nc*ovt+1;
    p_f = p_i+Nc-1;
    
    for j = 1 : nchannels % varre todos os canais espaciais
        ps_i = (j-1)*Ns*ovs+1;
        ps_f = ps_i+Ns-1;
        data = iq(ps_i:ps_f,p_i:p_f);
        [Nsamples,M] = size(data);
        vdata = sum(data,1);
        
        % usado no bifásico (slug-flow) ===================================
        dbbl = dbubble(ps_i:ps_f,p_i:p_f);
        
        media(j,i) = mean(mean(abs(dbbl)));
                      
        if media(j,i) < 80
            media_msk(j,i) = 0;
        else
            media_msk(j,i) = 1;
        end
        % =================================================================
%         index=1;
%         for k=1:Nsamples
            %  Find the proper data
%         if dif == 1
%              vdata = (diff(data(1,:),1));
%         else
%              vdata = (data(1,:));
%         end
            %  Calculate the autocorrelation and the velocity
            if (abs(std(vdata)) > abs(RMS/500))
                desv(j,i) = std(abs(vdata));
                auto  = vdata(2:(M-1)) * vdata(1:(M-2))' ;
%                 phi_est(index) =  atan2(imag(auto),real(auto));
                phi_est = atan2(imag(auto),real(auto));
            else
                phi_est = 0;
            end 
%             index=index+1;
%         end
%         phi_auto=mean(phi_est);

        phi_auto = phi_est;
        % Time shift estimator
        phi_possible = phi_auto+Np*2*pi;
        v_possible = c*fprf/(4*pi*fc) * phi_possible;
        ts_possible = v_possible./((c/2)*fprf);
        lag_possible = ceil(ts_possible.*fs);
        
        if dif == 1
             dataf = diff(data,1,2);
        else
             dataf = data;
        end
        
        for w = 1 : length(Np)
            for t = 1 : M-2
                if lag_possible(w) < 0
                    ccross(w,t) = dataf(1+abs(lag_possible(w)):end,t+1)'*dataf(1:end-abs(lag_possible(w)),t);
                else
                    ccross(w,t) = dataf(1:end-abs(lag_possible(w)),t+1)'*dataf(1+abs(lag_possible(w)):end,t);
                end
            end
            cross_s(w) = sum(ccross(w,:));
        end     
        [~,lag] = max(cross_s');
        phi_true = phi_auto+Np(lag)*2*pi;
        vr(j,i) = -c*fprf/(4*pi*fc) * phi_true;        
    end

end

% vr = vr.*media_msk;

EAM.processing_time = toc;
EAM.vel = vr;
EAM.msk = media_msk;
EAM.desv = desv;

end