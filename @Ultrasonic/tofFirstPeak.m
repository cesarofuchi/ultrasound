function US= tofFirstPeak(this,usData,range1,range2,c, corte,out,sinalRef)
            %tofFirstPeak  Cálculo de tempo de trânsito pelo primeiro pico
            %   [tt]=US_tt(dados,freq_transd,ciclos_onda,freq_amostragem,periodo_de_intere
            %   sse1,periodo_de_interesse2)
            %   dados: forma de onda
            %   freq_transd: frequência do transdutor
            %   ciclos_onda: quantos ciclos possui a forma de onda
            %   freq_amostragem: frequência de amostragem dos dados
            %   range1: período onde pode estar o primeiro eco
            %   range2: período onde pode estar o segundo eco
            %   c: velocidade do som no meio
            %   corte: % de corte desejado, caso não encontre nada, retorna zero
            %
            %   tt: tempo de trânsito calculado
            %
            %   See also US_tt, US_lerDados.
            %   Author(s): Ofuchi, C.Y.
            %   Copyright 2012 The MathWorks, Inc.
            %   $Revision: 1.0 $  $Date: 2012/07/20 11:47:07 $
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
            %[max1, posicao_eco1] = max(abs(hilbert(usData(range1,1))));
            [max1, posicao_eco1] = max(sinalRef);
            posicao_eco1=posicao_eco1+offset1-1;
            
            novo=abs(hilbert(usData(range2,1)));            
            z=length(novo)-length(sinalRef);
            sinalRef(end+1:end+z)=zeros(1,z);
            sinalSub=novo-sinalRef;
            
            [max2, posicao_eco2] = max(sinalSub);
            
            percent=max2/max1*100;
            % se for abaixo do corte, ignorar
            %if max2 < percent
            %    tt=0;
            %    dist=0;
            %else %senão 
                %[max1, posicao_eco1] = max(dados(periodo_de_interesse1,1));
                %[max2, posicao_eco2] = max(dados(periodo_de_interesse2,1));
                %[pk,lk]=findpeaks(abs(dados(periodo_de_interesse2,1)),'MinPeakProminence',out);                
             %--------------- procurar pelo pico minimo -----------------------------------------%  
                [pk,lk]=findpeaks(sinalSub,'MinPeakProminence',out);                
                %posicao_eco2=periodo_de_interesse3(posicao_eco2)
                if(out~=0)
                     if(~isempty(lk))
                         if(lk(1)<posicao_eco2-40) % seria o pico
                            posicao_eco2=lk(1)+offset2-1;
                         else
                            posicao_eco2=posicao_eco2+offset2-1;
                         end
                     else
                            posicao_eco2=posicao_eco2+offset2-1;
                     end
                else
                    posicao_eco2=posicao_eco2+offset2-1;
                end
                
                %--------------- procurar pelo pico minimo -----------------------------------------%
                US.pos1=posicao_eco1;
                US.pos2=posicao_eco2;
                %US.pos2=posicao_eco2+periodo_de_interesse2(1);
                US.tt=(US.pos2-US.pos1)/this.fs;
                US.max2=max2;

                %tt=(posicao_eco2-posicao_eco1);
                US.dist=c*US.tt/2;
            %end
        end