function US= tofMaxPeak(usObj,dados,range1,range2,vel_som, corte,out,sinalRef)
%tofMax Cálculo de tempo de trânsito por amplitude
%   [tt]=US_tt(dados,freq_transd,ciclos_onda,freq_amostragem,periodo_de_intere
%   sse1,periodo_de_interesse2)
%   dados: forma de onda
%   freq_transd: frequência do transdutor
%   ciclos_onda: quantos ciclos possui a forma de onda
%   freq_amostragem: frequência de amostragem dos dados
%   range1: período onde pode estar o primeiro eco
%   range2: período onde pode estar o segundo eco
%   corte: % de corte desejado, caso não encontre nada, retorna zero
%
%   US: retorna diversas variaveis relacionadas ao eco
%
%   See also US_tt, US_lerDados.
%   Author(s): Ofuchi, C.Y.            
%   $Revision: 1.0 $  $Date: 02/11/20$

%----------------------------------------------------------
% Técnica mais simples, apenas para encontrar o máximo
%----------------------------------------------------------           
            
                     
            %dados=dados-repmat(mean(dados)',[1 size(dados,1)])';
            %dados=detrend(dados(1:end-1));
            %size(dados);
            %size(onda_padrao);
            % faço a correlação com a onda padrão, para melhorar a onda que pode estar
            % influenciada por ruídos. A autocorrelação é usada somente para isso!
            offset1=range1(1);
            offset2=range2(1);
            %periodo_de_interesse3=fliplr(periodo_de_interesse2);
            % achar os máximos
            %[max1, posicao_eco1] = max(abs(hilbert(dados(range1,1))));
            [max1, posicao_eco1] = max(sinalRef);
            posicao_eco1=posicao_eco1+offset1-1;
            
            atual=abs(hilbert(dados(range2,1)));            
            z=length(atual)-length(sinalRef);
            sinalRef(end+1:end+z)=zeros(1,z);
            sinalSub=atual-sinalRef;
            % debug
%             figure
%             hold all
%             plot(atual)
%             plot(sinalRef)
%             plot(sinalSub)
%             legend('atual','ref','sub')
            
            
            [max2, posicao_eco2] = max(sinalSub);
            
            %percent=max2/max1*100;          
                
            US.pos1=posicao_eco1;
            US.pos2=posicao_eco2+range2(1);
            US.tt=(US.pos2-US.pos1)/usObj.fs;
            US.max2=max2;

            %tt=(posicao_eco2-posicao_eco1);
            US.dist=vel_som*US.tt/2;
            %end
        end