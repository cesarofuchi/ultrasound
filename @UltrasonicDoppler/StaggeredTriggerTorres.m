function STT = StaggeredTriggerTorres(usdObj,diq,c,T1,T2,noisepw,m,n,j,k)
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
% iq: matriz de dados complexos (dados RF após aplicação de hilbert)
% c velocidade do som no meio
% fc : frequência central do transdutor
% T1: período de tempo entre a primeira emissão e a segunda - pulse
% repetition time 1 ou PRT1) - T1 precisa ser < que T2
% T2: período de tempo entre a segunda emissãr e a terceira pulse
% repetition time 2 ou PRT2)
% Ns: números de amostras espaciais (define a resolução espacial do
% algoritmo)
% Nc: números de amostras temporais (define a resolução temporal do
% algoritmo)
% ovs: define sobreposição de canais espaciais 
% ovt: define sobreposição de canais temporais 
% subsampling: pode ser de 1 a Ns, se for 1 apenas a primeira amostra
% espacial de um range gate é autocorrelacionada, se for 2 Ns/2 amostras
% são autocorrelacionadas, padrão é deixar igual a Ns
% m e n: staggered PRT ratio m/n

% VARIÁVEIS DE SAÍDA
% v1d: mapa espaço-temporal da velocidade radial estimada usando o método 
% de Torres et al.(componente da velocidade na  direção do feixe do transdutor) 
% para PRT1 ou T1 dealiased
% v2d: mapa espaço-temporal da velocidade radial estimada usando o método 
% de Torres et al.(componente da velocidade na  direção do feixe do transdutor) 
% para PRT2 ou T2 dealiased 
% v1: mapa espaço-temporal da velocidade radial estimada usando método convencional
% (componente da  velocidade na  direção do feixe do transdutor) para PRT1 ou T1 
% v2: mapa espaço-temporal da velocidade radial estimada usando método convencional
% (componente da velocidade na  direção do feixe do transdutor) para PRT2 ou T2 
% va1: máxima velocidade de Nyquist para PRT1 ou T1
% va2: máxima velocidade de Nyquista para PRT2 ou T2
% W: mapa espaço-temporal da potência ou energia média estimada do sinal
% var1: mapa espaço-temporal variância da velocidade relativa ou T1 ou PRT1
% estimada pela autocorrelação
% var2: mapa espaço-temporal variância da velocidade relativa ou T2 ou PRT2
% estimada pela autocorrelação
% vst: mapa espaço-temporal da velocidade radial estimada usando método
% staggered trigger sem a dealiasing rules de Torres et al. -ou seja, a
% velocidade é estimada fazendo a subtração dos ângulos e dividindo pela
% subtração dos tempos T1 e T2..
% (componente da  velocidade na  direção do feixe do transdutor) para PRT1 ou T1 
% processing_time: tempo decorrido para efetuar a estimação da velocidade

%NoisePower=mean(std(iq'))+3*abs(std(std(iq')));

fprf = usdObj.usObj.fprf;
fs = usdObj.usObj.fs;

fc = usdObj.usObj.fc;
iq = diq;
Ns = usdObj.ns(k);
Nc = usdObj.nc(j);
ovs = usdObj.ovs;
ovt = usdObj.ovt;


%disp('St running')
%adapta os num de amostras espaciais para ser janelado com ov e com o Ns
%por meio de adicção de valores zeros no final do vetor (zeros pads)
pad_n=Ns*ovs-rem(size(iq,1)-Ns,Ns*ovs);
pad=zeros(pad_n,size(iq,2));
iq=[iq ; pad]; % acrescenta zeros na profundidade no fim do vetor

% calcula o numero de canais espaciais já com o pad de zeros no final e
% overflow
nchannels=(size(iq,1)-Ns)/(ovs*Ns);

% calcula o número de canais temporais com overflow
nchannelst=floor((size(iq,2)-Nc)/(ovt*Nc))+1;

va1 = c/(4*T1*fc); % vmax T1
va2 = c/(4*T2*fc); % vmax T2

C0=0;
C1=C0+2*va2;
C2=C1-2*va1;
C3=C2+2*va2;
C4=C3-2*va1;
C5=C4+2*va2;
C6=C5-2*va1;

%inicializa com zeros as matrizes de resultados de saída do algoritmo
v1d=zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade (dealiased) radial
v2d=zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade (dealiased) radial

v12=zeros(nchannels,nchannelst); % mapa espaço-temporal da diferenças das velocidades (aliased) radial

v1=zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial relativa a T1
v2=zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade radial relativa a T2
vst=zeros(nchannels,nchannelst); % mapa espaço-temporal da velocidade st

W=zeros(nchannels,nchannelst); % mapa espaço-temporal da potência média
% W2=zeros(nchannels,nchannelst); % mapa espaço-temporal da potência média

sigmaT1=zeros(nchannels,nchannelst);% mapa espaço-temporal variância da velocidade estimada pela autocorrelação
sigmaT2=zeros(nchannels,nchannelst);% mapa espaço-temporal variância da velocidade estimada pela autocorrelação

%tic %inicia contagem de tempo de execução
for i=1:nchannelst % varre todos os canais temporais
    p_i=(i-1)*Nc*ovt+1;
    p_f=p_i+Nc-1;
    
    for j=1:nchannels % varre todos os canais espaciais
        ps_i=(j-1)*Ns*ovs+1;
        ps_f=ps_i+Ns-1;
        
        data=iq(ps_i:ps_f,p_i:p_f);
        vdata=sum(data,1); % integra todos os valores de um sample range
        
        Nc1 = length(vdata(2:2:end));
        Nc2 = length(vdata(3:2:end));
        
        % de acordo com Loupas1995 página 5 equação 21 em diante
        % Calculate the autocorrelation of R(T1) and R(T2)
        auto1T1= (vdata(2:2:end)* vdata(1:2:end)') / (Nc1);  % R(T1)
        auto1T2= (vdata(3:2:end)* vdata(2:2:end-1)') / (Nc2);% R(T2)
        
        %  Calculate the autocorrelation of RT1(0) and RT2(0)
        auto0T1=(sum(abs(vdata(2:2:end)).^2)+sum(abs(vdata(1:2:end)).^2)) / (2*(length(vdata(2:2:end))));   % R(0) = W
        auto0T2=(sum(abs(vdata(3:2:end)).^2)+sum(abs(vdata(2:2:end-1)).^2)) / (2*(length(vdata(3:2:end)))); % R(0) = W
        
        W(j,i) = (auto0T1+auto0T2)/2;
        
        gamaT1=(abs(auto1T1)^2)/(abs(auto0T1)^2);
        gamaT2=(abs(auto1T2)^2)/(abs(auto0T2)^2);
        
        SNRT1=abs(auto0T1)/abs(noisepw);
        SNRT2=abs(auto0T2)/abs(noisepw);
        
        %  Calculate the standard deviation of the estimate
        sigmaT1(j,i)= (c/(2*fc))*sqrt(((1+1/SNRT1)^2-gamaT1)/(8*Nc1*pi^2*T1^2*gamaT1));
        sigmaT2(j,i)= (c/(2*fc))*sqrt(((1+1/SNRT2)^2-gamaT2)/(8*Nc2*pi^2*T2^2*gamaT2));
       
        %  Calculate the phase angle from the autocorrelation of lag 1
        phi_est1 =  atan2(imag(auto1T1),real(auto1T1));
        phi_est2 =  atan2(imag(auto1T2),real(auto1T2));
        
        snr1(j,i)=10*log10(SNRT1);
        snr2(j,i)=10*log10(SNRT2);
        
        % calculate the axial velocity
        v1(j,i) = -c*(1/(T1))/(4*pi*fc) *  phi_est1;
        v2(j,i) = -c*(1/(T2))/(4*pi*fc) *  phi_est2;       
        v12(j,i)= v1(j,i)-v2(j,i); 
        vst(j,i)=-c*(1/abs(T1-T2))/(4*pi*fc) * (phi_est2-phi_est1);
    end

end
%apply dealising rules acoording Torres et al (2004)
va=m*va1;
dC=2*va/(m*n);

% m=0
rule=(v12>(C0-dC/2))&(v12<(C0+dC/2));
P=0;Q=0;
v1d(rule)=v1(rule)+2*P*va1;
v2d(rule)=v2(rule)+2*Q*va2;


if m>=1
    rule=v12>(C1-dC/2);
    P=0;Q=1;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=v12<-(C1-dC/2);
    P=0;Q=-1;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=2
    rule=(v12>(C2-dC/2))&(v12<(C2+dC/2));
    P=1;Q=1;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%%
    rule=(v12<-(C2-dC/2))&(v12>-(C2+dC/2));
    P=-1;Q=-1;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=3
    rule=(v12>(C3-dC/2))&(v12<(C3+dC/2));
    P=1;Q=2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C3-dC/2))&(v12>-(C3+dC/2));
    P=-1;Q=-2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=4
    rule=(v12>(C4-dC/2))&(v12<(C4+dC/2));
    P=2;Q=2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C4-dC/2))&(v12>-(C4+dC/2));
    P=-2;Q=-2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=5
    rule=(v12>(C5-dC/2))&(v12<(C5+dC/2));
    P=2;Q=3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    % % % % %
    rule=(v12<-(C5-dC/2))&(v12>-(C5+dC/2));
    P=-2;Q=-3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=6
    rule=(v12>(C6-dC/2))&(v12<(C6+dC/2));
    P=3;Q=3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C6-dC/2))&(v12>-(C6+dC/2));
    P=-3;Q=-3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end


% processing_time=toc;
% disp(num2str(processing_time));

STT.v1d  = v1d;
STT.v2d  = v2d;
STT.v1   = v1;
STT.v2   = v2;
STT.va1  = va1;
STT.va2  = va2;
STT.W    = W;
STT.var1 = sigmaT1;
STT.var2 = sigmaT2;
STT.vel  = vst;
STT.snr1 = snr1;
STT.snr2 = snr2;

end