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
        function usObj = loadData(file,channels,varargin)
            
           if file < 0
                msg = ['Error to open file:' file];
                sprintf(msg)
                return
           end
           
           iter=1
           
           
           %import sensor.Ultrasonic;
           usObj = Ultrasonic;
           
           %% separe this after
           usObj = usObj.readPXI5752Xml(file);
           usObj.cycles = 4;
           while ~isempty(varargin)
                switch lower(varargin{1})
                      case 'samples'
                          usObj.samples=varargin{2};
                      case 'fc'
                          usObj.fc=varargin{2};
                      case 'cycles'
                          usObj.cycles=varargin{2};
                      %otherwise
                      %    error(['Unexpected option: ' varargin{1}])
                 end
                 varargin(1) = [];
                  
           end
           
           
           
           usObj.initTime = System.IO.File.GetCreationTime(file);           
           
           fn = usObj.nwaves*channels*usObj.samples;
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
           data=double(data);
           %raw data in ADC 12 bits data. Need to normalize to Volts.
           adc_bits=12;
           data=data./(2^adc_bits);
           
           fclose(f);
           usObj.data = reshape(data(1:fn),usObj.samples,usObj.nwaves);
           usObj.workData=usObj.data;
           
        end 
       
       function usObj = loadFolder(folder,channels,filter_type)
           
           if folder < 0
                msg = ['Error to open file:' folder];
                sprintf(msg)
                return
           end
           %find xml files
%           filter_type='*.xml';
           files = dir([folder filter_type]);
           files.name;
           files([files.isdir]) = [];
           
           if isempty(files)
                msg = ['No files found using this filter:' filter_type ' or this path:' folder];
                sprintf(msg)
                return
           end
           
           
           % first file info          
           usObj = Ultrasonic;
           usObj = usObj.loadData([folder files(1).name],channels)


           %iterate folder 
           
            start = 1; % onde começa efetivamente arquivos com amostras (ver 'files' no Workspace)             
            stop = size(files, 1); % último arquivo com amostras            
            nfiles=length(start : stop);
            data = zeros(nfiles,usObj.samples);
            % loop que salva em uma matriz todos os pacotes de dados adquiridos durante o experimento
            tic
            %import sensor.Ultrasonic;
            usObj = Ultrasonic;
            usObj.cycles = 4;            
            for index_i = start : stop                
                file = [folder files(index_i).name];
                usObj = usObj.loadData(file,channels);                
                data(index_i,:)=usObj.data;
            end
            usObj.data=data;
            toc          
            
       end    
       
      
       
    end
    
    methods(Access = public)
        usObj = readPXI5752Xml(usObj,file)        
        US = tofFirstPeak(usObj,dados,range1,range2,c, corte,out,sinalRef)
        US = tofMaxPeak(usObj,dados,range1,range2,vel_som, corte,out,sinalRef)
        us_map_diff=us_energy_diff(usObj,data,sStep,tStep)
    end
    
    methods
         function obj = Ultrasonic()
          
         end
         
    end
end