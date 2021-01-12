
function STT = StaggeredTriggerTorres(usdObj,c,T1,T2,subsampling,m,n)
% Modificação do Método de autocorrelação ou phase shift estimator
% para estimar a velocidade usando a técnica staggered trigger ou
% dual PRT ou staggered PRT. Ao final mplementa as regras definidas 
% por Torres et al. no artigo "Design, Implementation, and Demonstration 
%of a Staggered PRT Algorithm for the WSR-88D", Journal of Atmospheric and Oceanic
% Technology, 2004 , volume 21, 1389-1399, para realizar o dealiasing da
% velocidade estimada.
% 
% Modificado para ter como entrada os tempos T1 e T2 e ao invés de Ts
% pois dessa forma pode-se saber qual é a sequência correta de Tprf, ou
% seja, se o maior período vem antes do menor e vice versa. Usando apenas o
% Ts não tem como saber...
% Estima a velocidade radial (componente ao longo da direção do feixe do 
% transdutor através do desvio de fase do sinal (phase shift method) 
%
% VARIÁVEIS DE ENTRADA:
% iq : matriz de dados complexos (dados RF após aplicação de hilbert)
% c : velocidade do som no meio
% fc : frequência central do transdutor
% T1 : período de tempo entre a primeira emissão e a segunda - pulse
% repetition time 1 ou PRT1) - T1 precisa ser < que T2
% T2 : período de tempo entre a segunda emissãr e a terceira pulse
% repetition time 2 ou PRT2)
% Ns : números de amostras espaciais (define a resolução espacial do
% algoritmo)
% Nc : números de amostras temporais (define a resolução temporal do
% algoritmo)
% ovs : define sobreposição de canais espaciais 
% ovt : define sobreposição de canais temporais 
% subsampling : pode ser de 1 a Ns, se for 1 apenas a primeira amostra
% espacial de um range gate é autocorrelacionada, se for 2 Ns/2 amostras
% são autocorrelacionadas, padrão é deixar igual a Ns
% m e n: staggered PRT ratio m/n

% VARIÁVEIS DE SAÍDA
% v1d : mapa espaço-temporal da velocidade radial estimada usando o método 
% de Torres et al.(componente da velocidade na  direção do feixe do transdutor) 
% para PRT1 ou T1 dealiased
% v2d : mapa espaço-temporal da velocidade radial estimada usando o método 
% de Torres et al.(componente da velocidade na  direção do feixe do transdutor) 
% para PRT2 ou T2 dealiased 
% v1 : mapa espaço-temporal da velocidade radial estimada usando método convencional
% (componente da  velocidade na  direção do feixe do transdutor) para PRT1 ou T1 
% v2 : mapa espaço-temporal da velocidade radial estimada usando método convencional
% (componente da velocidade na  direção do feixe do transdutor) para PRT2 ou T2 
% va1 : máxima velocidade de Nyquist para PRT1 ou T1
% va2 : máxima velocidade de Nyquista para PRT2 ou T2
% W : mapa espaço-temporal da potência ou energia média estimada do sinal
% var1 : mapa espaço-temporal variância da velocidade relativa ou T1 ou PRT1
% estimada pela autocorrelação
% var2 : mapa espaço-temporal variância da velocidade relativa ou T2 ou PRT2
% estimada pela autocorrelação
% vst : mapa espaço-temporal da velocidade radial estimada usando método
% staggered trigger sem a dealiasing rules de Torres et al. -ou seja, a
% velocidade é estimada fazendo a subtração dos ângulos e dividindo pela
% subtração dos tempos T1 e T2..
% (componente da  velocidade na  direção do feixe do transdutor) para PRT1 ou T1 
% processing_time: tempo decorrido para efetuar a estimação da velocidade

% NoisePower=mean(std(iq'))+3*abs(std(std(iq')));

fprf = usdObj.usObj.fprf;
fs = usdObj.usObj.fs;

fc = usdObj.usObj.fc;
iq = usdObj.iq;
Ns = usdObj.ns;
Nc = usdObj.nc;
ovs = usdObj.ovs;
ovt = usdObj.ovt;

disp('Stt running')
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

% va1=c/(4*T*fc);
% va2=c/(4*(T+Ts)*fc);
va1 = c/(4*T1*fc);
va2 = c/(4*(T2)*fc);
% if (rem(m,2)==0) 
%     %m é par
%     L=(m+n-1)/2;
% else %m é ímpar
%     L=(m+n-2)/2;
% end
% 
% C=zeros(1,L+1);
% 
% for il=1:1:L+1
%     if (rem(il,2)==1)
%         C(il+1)=C(il)+2*va2;
%     else
%         C(il+1)=C(il)-2*va1;
%     end
% end

C0 = 0;
C1 = C0+2*va2;
C2 = C1-2*va1;
C3 = C2+2*va2;
C4 = C3-2*va1;
C5 = C4+2*va2;
C6 = C5-2*va1;

%inicializa com zeros as matrizes de resultados de saída do algoritmo
v1d = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial
v2d = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial

v12 = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial

v1 = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial
v2 = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial

vst = zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade st
% m = zeros(nchannels,nchannelst);

W = zeros(nchannels,nchannelst); % mapa espaço-temporal da potência média

var1 = zeros(nchannels,nchannelst);% mapa espaço-temporal variância da velocidade estimada pela autocorrelação
var2 = zeros(nchannels,nchannelst);% mapa espaço-temporal variância da velocidade estimada pela autocorrelação

% autov  = zeros(1,Ns);
% auto2v = zeros(1,Ns);

tic % inicia contagem de tempo de execução
for i = 1 : nchannelst % varre todos os canais temporais
    p_i = (i-1)*Nc*ovt+1;
    p_f = p_i+Nc-1;
    
    for j = 1 : nchannels % varre todos os canais espaciais
        ps_i = (j-1)*Ns*ovs+1;
        ps_f = ps_i+Ns-1;
        
        data = iq(ps_i:ps_f,p_i:p_f);
        
        auto = 0; auto2 = 0;
        Wt = 0;
        vart = 0; vart2 = 0;
        
        for nsi = 1 : Ns/subsampling : Ns          
            vdata = (data(nsi,:));
            % Calculate the autocorrelation of R(T)
            Wtemp = (vdata(1:Nc)*vdata(1:Nc)');
            autot = vdata(2:2:end)* vdata(1:2:end)';
            auto  = auto+autot;
%             autov(nsi)=autot;
            auto2t = vdata(3:2:end) * vdata(2:2:end-1)';
            auto2  = auto2+auto2t;
%             auto2v(nsi)=auto2t;

            %  Calculate the average power
            Wt = Wt+(1/(Nc-1))*Wtemp;
            
            %  Calculate the variance of the estimate
            %vart = vart+2*abs(auto)/Wtemp;
            vart = vart+(1-abs(autot)/Wtemp);
            %vart2 = vart2+2*abs(auto2)/Wtemp;
            vart2 = vart2+(1-abs(auto2t)/Wtemp);
        end
         %  Calculate the phase angle from the autocorrelation
%          auto=mode(autov);
%          auto2=mode(auto2v);
        phi_est1 = atan2(imag(auto),real(auto));
        phi_est2 = atan2(imag(auto2),real(auto2));
        W(j,i) = Wt/subsampling;
        vart  = vart/subsampling;
        vart2 = vart2/subsampling;
%         var1(j,i)=(c*(1/(T)))/(2*pi*sqrt(2)*fc)*sqrt(vart);
        var1(j,i) = (c*(1/(T1)))/(2*pi*sqrt(2)*fc)*sqrt(vart);
        %var(j,i)=(c^2*(1/(T^2)))/(16*pi^2*(2)*fc^2)*(vart);
%         var2(j,i)=(c*(1/(T+Ts)))/(2*pi*sqrt(2)*fc)*sqrt(vart);
        var2(j,i) = (c*(1/(T2)))/(2*pi*sqrt(2)*fc)*sqrt(vart2);
%         var2(j,i)=(c^2*(1/(T+Ts)^2))/(16*pi^2*(2)*fc^2)*(vart2);

        % calculate the true phase angle by combining the estimate 
        % from autocorrelation with the estimate from cross correlation 
%         v1(j,i) = -c*(1/(T))/(4*pi*fc) *  phi_est1;
%         v2(j,i) = -c*(1/(T+Ts))/(4*pi*fc) *  phi_est2;
        v1(j,i) = -c*(1/(T1))/(4*pi*fc) *  phi_est1;
        v2(j,i) = -c*(1/(T2))/(4*pi*fc) *  phi_est2;
        
%         rpm1=60*v1(j,i)/(2*pi)/0.019;
%         rpm2=60*v2(j,i)/(2*pi)/0.019;
        
        v12(j,i)= v1(j,i)-v2(j,i); 
%         if ((phi_est2<0 && phi_est1<0) )
%                         vst(j,i)=0;
%         else
            vst(j,i)=-c*(1/abs(T1-T2))/(4*pi*fc) * (phi_est2-phi_est1);
%         end
            %         vst(j,i)=-c*(1/Ts)/(4*pi*fc) * atan2(imag(auto2*auto'),real(auto2*auto'));
%         vst(j,i)=-c*(1/Ts)/(4*pi*fc) * angle(auto2*auto');
%         vst(j,i)=-c*(1/Ts)/(4*pi*fc) * angle(auto2/auto);
    end

end

% apply dealising rules
va = m*va1;
dC = 2*va/(m*n);

% m = 0
rule = (v12 > (C0-dC/2)) & (v12<(C0+dC/2));
P = 0; Q = 0;
v1d(rule) = v1(rule)+2*P*va1;
v2d(rule) = v2(rule)+2*Q*va2;

if m >= 1
    rule = v12 > (C1-dC/2);
    P = 0; Q = 1;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%
    rule = v12 < -(C1-dC/2);
    P = 0; Q = -1;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end

if m >= 2
    rule = (v12 > (C2-dC/2)) & (v12 < (C2+dC/2));
    P = 1; Q = 1;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%%
    rule = (v12 < -(C2-dC/2)) & (v12 > -(C2+dC/2));
    P = -1; Q = -1;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end

if m >= 3
    rule = (v12 > (C3-dC/2)) & (v12 < (C3+dC/2));
    P = 1; Q = 2;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%
    rule = (v12 < -(C3-dC/2)) & (v12 > -(C3+dC/2));
    P = -1; Q = -2;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end

if m >= 4
    rule = (v12 > (C4-dC/2)) & (v12 < (C4+dC/2));
    P = 2; Q = 2;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%
    rule = (v12 < -(C4-dC/2)) & (v12 > -(C4+dC/2));
    P = -2; Q = -2;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end

if m >= 5
    rule = (v12 > (C5-dC/2)) & (v12 < (C5+dC/2));
    P = 2; Q = 3;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%
    rule = (v12 < -(C5-dC/2)) & (v12 > -(C5+dC/2));
    P = -2; Q = -3;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end

if m >= 6
    rule = (v12 > (C6-dC/2)) & (v12 < (C6+dC/2));
    P = 3; Q = 3;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
    %%%
    rule = (v12 < -(C6-dC/2)) & (v12 > -(C6+dC/2));
    P = -3; Q = -3;
    v1d(rule) = v1(rule)+2*P*va1;
    v2d(rule) = v2(rule)+2*Q*va2;
end


processing_time = toc;
disp(num2str(processing_time));

% [v1d,v2d,v1,v2,va1,va2,W,var1,var2,vst]

STT.v1d  = v1d;
STT.v2d  = v2d;
STT.v1   = v1;
STT.v2   = v2;
STT.va1  = va1;
STT.va2  = va2;
STT.W    = W;
STT.var1 = var1;
STT.var2 = var2;
STT.vel  = vst;

end