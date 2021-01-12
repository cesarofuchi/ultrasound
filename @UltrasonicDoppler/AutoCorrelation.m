function ACM = AutoCorrelation(usdObj,c,subsampling,j,k)
% Autocorrelation Method
%
% Estima a velocidade radial (componente ao longo da direção do feixe do 
% transdutor através do desvio de fase do sinal (phase shift method) 
%
% VARIÁVEIS DE ENTRADA:
% iq  : matriz de dados complexos (dados RF após aplicação de hilbert)
% c   :velocidade do som no meio
% prf : pulse repetition frequency, frequência das emissões de pulsos
% fc  : frequência central do transdutor
% fs  : frequência de amostragem dos dados 
% Ns  : números de amostras espaciais (define a resolução espacial do algoritmo)
% Nc  : números de amostras temporais (define a resolução temporal do algoritmo)
% ovs : define sobreposição de canais espaciais 
% ovt : define sobreposição de canais temporais 
% subsampling :
%
% VARIÁVEIS DE SAÍDA
% vr: mapa espaço-temporal da velocidade radial estimada (componente da 
% velocidade na  direção do feixe do transdutor)
% W: mapa espaço-temporal da potência média estimada do sinal
% var: mapa espaço-temporal variância da velocidade estimada pela autocorrelação
% phi: mapa espaço-temporal dos ângulos estimados
% processing_time: tempo decorrido para efetuar a estimação da velocidade
% vmax: velocidade máxima que pode ser medida com a configuração entrada

% Inicialização de variáveis
fprf = usdObj.usObj.fprf;
fs   = usdObj.usObj.fs;
fc   = usdObj.usObj.fc;
iq   = usdObj.iq;
Ns   = usdObj.ns(k);
Nc   = usdObj.nc(j);
ovs  = usdObj.ovs;
ovt  = usdObj.ovt;
 
% NoisePower=mean(std(iq'))+3*abs(std(std(iq')));
% disp('ACM running')

% adapta os num de amostras espaciais para ser janelado com ov e com o Ns
% por meio de adicção de valores zeros no final do vetor (zeros pads)
pad_n = Ns*ovs-rem(size(iq,1)-Ns,Ns*ovs);
pad = zeros(pad_n,size(iq,2));
iq = [iq ; pad]; % acrescenta zeros na profundidade no fim do vetor

% calcula o numero de canais espaciais já com o pad de zeros no final e
% overflow
nchannels = (size(iq,1)-Ns)/(ovs*Ns);

% calcula o número de canais temporais com overflow
nchannelst = floor((size(iq,2)-Nc)/(ovt*Nc))+1;

%calcula a velocidade máxima 
vmax = (pi)*fprf*c/(4*pi*fc);

%inicializa com zeros as matrizes de resultados de saída do algoritmo
vr = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial
phi = zeros(nchannels,nchannelst);% mapa espaço-temporal dos ângulos estimados
W = zeros(nchannels,nchannelst); % mapa espaço-temporal da potência média
var = zeros(nchannels,nchannelst);% mapa espaço-temporal variância da velocidade estimada pela autocorrelação
f = zeros(nchannels,nchannelst);% mapa espaço-temporal freq central do transdutor
% tic % inicia contagem de tempo de execução

for i = 1 : nchannelst % varre todos os canais temporais
    p_i = (i-1)*Nc*ovt+1;
    p_f = p_i+Nc-1;
    
    for j = 1 : nchannels % varre todos os canais espaciais
        ps_i = (j-1)*Ns*ovs+1;
        ps_f = ps_i+Ns-1;
        
        data = iq(ps_i:ps_f,p_i:p_f); % janela -> (amostras,ondas)
        auto = zeros(1,Ns/subsampling);
        Wt = 0;
        vart = zeros(1,Ns/subsampling);
        
        for nsi = 1 : Ns/subsampling : Ns          
            vdata2 = data(nsi,2:(Nc));
            vdata1 = data(nsi,1:(Nc-1));
            %  Calculate the autocorrelation of R(T)
%             r2=xcorr(data(nsi,1:Nc),1,'unbiased');
            
%             auto(nsi)  = data(nsi,2:(Nc)) * data(nsi,1:(Nc-1))' ;
            auto(nsi) = vdata2 *vdata1';
            
%             auto(nsi)=(r2(3));%+r2(1))/2;

            %  Calculate the average power
            %Wt=Wt+(1/(Nc-1))*(data(nsi,1:Nc-1)*data(nsi,1:Nc-1)');
            vd2 = (vdata1*vdata1');
            Wt = Wt+(1/(Nc-1))*vd2;

            
            %  Calculate the variance of the estimate
            %vart(nsi)=abs(auto(nsi))/(data(nsi,1:Nc-1)*data(nsi,1:Nc-1)');
            vart(nsi) = abs(auto(nsi))/vd2;
            
            %vart(nsi)=abs(auto(nsi))/(vdata(1:Nc)*vdata(1:Nc)');
        end
% apenas se quiser estimar a frequencia de RF        
        ftemp = zeros(1,Nc);
        for nsi = 1 : Nc 
%             ftemp=ftemp+angle(data(2:end,nsi)'*data(1:end-1,nsi))*fs/(2*pi);
            ftemp(nsi) = (data(2:end,nsi)'*data(1:end-1,nsi));
        end
%         f(j,i)=ftemp/Nc;
        f(j,i) = angle(mean(ftemp))*fs/(2*pi);
% %          f(j,i)=angle(max(ftemp))*fs/(2*pi); %estima pela max
% %          f(j,i)=angle(mean(ftemp(abs(ftemp)>0.95)))*fs/(2*pi);
%          %  Calculate the phase angle from the autocorrelation
%         [value,ind]=max(vart); 
        
        %phi_est =  atan2(imag(auto(ind)),real(auto(ind)));
        phi_est = angle(mean(auto));
        phi_auto = phi_est;
        W(j,i) = Wt/subsampling;
        %var(j,i)=(c*prf)/(2*pi*sqrt(2)*fc)*sqrt(abs(1-mean(vart)));
        var(j,i) = ((c/(4*pi*fc))^2)*2*fprf^2*(1-mean(vart));
        %kur(j,i)=mean(kurt);
        
        % calculate the true phase angle by combining the estimate 
        % from autocorrelation with the estimate from cross correlation 
        phi(j,i) = phi_auto;
%        vr(j,i) = -c*prf/(4*pi*(fc+f(j,i))) * phi_auto;
        vr(j,i) = -c*fprf/(4*pi*fc) * phi_auto;
    end

end
processing_time = toc;
% disp(num2str(processing_time));

ACM.vel = vr;
ACM.W = W;
ACM.var = var;
ACM.phi = phi;
ACM.processing_time = processing_time;
ACM.vmax = vmax;
ACM.f = f;

end