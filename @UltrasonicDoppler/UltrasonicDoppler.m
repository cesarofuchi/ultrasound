classdef UltrasonicDoppler 
    
    properties        
        usObj % usObj
        ns    % spatial window -> samples
        nc    % temporal window (channels)-> PRF
        ovs   % overwrite samples
        ovt   % overwrite time
        iq    % IQ data (In phase Quadratura)
        t     % time
        
    end
    methods (Static)
        
       
    end
    
    methods
        % Constructor
        % get Ultrasonic Object and add other parameters
        % create IQ data for velocity 
         function dopplerObj = UltrasonicDoppler(usObj,ns,nc,ovs,ovt)             
             dopplerObj.usObj = usObj; % usObj             
             dopplerObj.ns = ns;
             dopplerObj.nc = nc;
             dopplerObj.ovs = ovs;
             dopplerObj.ovt = ovt;
             
             %% Matched Filter
             impulse_response=sin(2*pi*usObj.fc*(0:1/usObj.fs:usObj.cycles/usObj.fc));
             impulse_response=impulse_response.*hamming(max(size(impulse_response)))';
             b = flipud(impulse_response);
             dataf1=filter(b,1,usObj.workData);
             
             % create demodulated data             
             h = hilbert(dataf1);
             t = (0:1/usObj.fs:size(h,1)*(1/usObj.fs)-(1/usObj.fs))';
             t = repmat(t,1,size(h,2));
             % dopplerObj.iq = h.*exp(-1i*2*pi*(usObj.fc)*t);
             dopplerObj.iq = h.*exp(-1i*2*pi*(usObj.fc)*t);
             dopplerObj.t = t;
         end
        
        
         function plotData(obj,v)
             
            plot((v-100)*randn(100,1))
            title('US Sensor')
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
% %             obj.Property1 = inputArg1 + inputArg2;
         end
        
        function this = load(this,filename,pathname)
            this.filename=filename;
            this.pathname=pathname;
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods%(Access = private)
           ACM = AutoCorrelation(usdObj,c,subsampling,j,k);
           CCM = TimeShiftEstimator(usdObj,d,c,j,k);
           %EAM = ExtAutoCorrelation(usdObj,dbubble,c,dif,NpRange,j,k);
           EAM= ExtAutoCorrelation(usdObj,c,NpRange);
           STT = StaggeredTriggerTorres(usdObj,diq,c,T1,T2,noisepw,m,n,j,k)
    end
end