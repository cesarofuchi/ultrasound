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
        US = tofFirstPeak(this,dados,periodo_de_interesse1,periodo_de_interesse2,vel_som, corte,out)
        US = tofMaxPeak(this,dados,range1,range2,vel_som, corte,out,sinalRef)
    end
    
    methods
         function obj = Ultrasonic()
          
         end
         
    end
end