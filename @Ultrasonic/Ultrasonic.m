classdef Ultrasonic
    
    properties        
        filename;
        filepath;
        type;
        fs
        position
        initTime
        ts
        fc
        fprf
        cycles
        nwaves
        channels
        samples
        data     % full experimental data
        workData % cropped data
    end
    
    methods (Static)
      
       function usObj = loadData(file,channels,samples)
           
           if file < 0
                msg = ['Error to open file:' file];
                sprintf(msg)
                return
           end
           import sensor.Ultrasonic;
           usObj = Ultrasonic;
           
           %% separe this after
           usObj = usObj.readPXI5752Xml(file);
           
           usObj.cycles = 4;
           usObj.initTime = System.IO.File.GetCreationTime(file); 
           usObj.samples = samples;
           
           fn = usObj.nwaves*channels*samples;
           data = int16(zeros(1,fn));
           
           %%
           usFile = [usObj.filepath '\' usObj.filename];
           f = fopen(usFile);
           if f < 0
                msg = ['Error to open file:' file];
                sprintf(msg)
                return
           end
           % fread(fileID,sizeA,precision,skip,machinefmt) %
           data = fread(f,fn,'int16=>int16',0,'b');           
           usObj.data = reshape(data(1:fn),samples,usObj.nwaves);
           
       end    
    end
    
    methods(Access = private)
        usObj = readPXI5752Xml(usObj,file)
        sound_speed = water_by_temperature(temp) 

    end
    
    methods
         function obj = Ultrasonic()
          
         end
         
         function US = tofMinPeak(this,dados,periodo_de_interesse1,periodo_de_interesse2,vel_som, corte,out)
            %US_tt Cálculo de tempo de trânsito por amplitude
            %   [tt]=US_tt(dados,freq_transd,ciclos_onda,freq_amostragem,periodo_de_intere
            %   sse1,periodo_de_interesse2)
            %   dados: forma de onda
            %   freq_transd: frequência do transdutor
            %   ciclos_onda: quantos ciclos possui a forma de onda
            %   freq_amostragem: frequência de amostragem dos dados
            %   periodo_de_interesse1: período onde pode estar o primeiro eco
            %   periodo_de_interesse2: período onde pode estar o segundo eco
            %   corte: % de corte desejado, caso não encontre nada, retorna zero
            %
            %   tt: tempo de trânsito calculado
            %
            %   See also US_tt, US_lerDados.
            %   Author(s): Ofuchi, C.Y.
            %   Copyright 2012 The MathWorks, Inc.
            %   $Revision: 1.0 $  $Date: 2012/07/20 11:47:07 $
            %dados=dados-repmat(mean(dados)',[1 size(dados,1)])';
            dados=detrend(dados(1:end-1));
            %size(dados);
            %size(onda_padrao);
            % faço a correlação com a onda padrão, para melhorar a onda que pode estar
            % influenciada por ruídos. A autocorrelação é usada somente para isso!
            offset1=periodo_de_interesse1(1);
            offset2=periodo_de_interesse2(1);
            %periodo_de_interesse3=fliplr(periodo_de_interesse2);
            % achar os máximos
            [max1, posicao_eco1] = max(abs(dados(periodo_de_interesse1,1)));
            [max2, posicao_eco2] = max(abs(dados(periodo_de_interesse2,1)));
            percent=max2/max1*100;
            % se for abaixo do corte, ignorar
            %if max2 < percent
            %    tt=0;
            %    dist=0;
            %else %senão 
                [max1, posicao_eco1] = max(dados(periodo_de_interesse1,1));
                [max2, posicao_eco2] = max(dados(periodo_de_interesse2,1));
                
                
                [pk,lk]=findpeaks(abs(dados(periodo_de_interesse2,1)),'MinPeakProminence',out);
                
                
                posicao_eco1=posicao_eco1+offset1-1;
                %posicao_eco2=periodo_de_interesse3(posicao_eco2)
                posicao_eco2=posicao_eco2+offset2-1;
                if(out~=0)
                     if(~isempty(lk))
                         posicao_eco2=lk(1)+offset2-1;
                     end
                end

                US.tt=(posicao_eco2-posicao_eco1)/this.fs;
                US.max2=max2;

                %tt=(posicao_eco2-posicao_eco1);
                US.dist=vel_som*US.tt/2;
            %end
        end
        
        function plotData(obj,v)
             
             plot((v-100)*randn(100,1))
             title('US Sensor')
         end
    end
end