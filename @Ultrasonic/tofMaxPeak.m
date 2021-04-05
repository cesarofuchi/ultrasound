function US= tofMaxPeak(usObj,data,range1,range2,thr,varargin)
%tofMaxPeak Cálculo de tempo de trânsito por amplitude

%   data = A-scan signal
%   range1: range where the 1st echo reference is
%   range2: range where the desired echo may appear
%   c : speed of sound at the medium
%   thr: minimum threshould for the desired echo. If lower, return
%        US.valid=0

%   US: strutcture with information
%
%   Author(s): Ofuchi, C.Y.            
%   $Revision: 1.0 $  $Date: 02/11/20$

%----------------------------------------------------------
% Info: simple technique to find the time of flight using the max echo peak
% 
%----------------------------------------------------------      
MAX=0;
FIRST_PEAK=1;
US.valid=1; %default is a valid measure
%default type
type=MAX;

while ~isempty(varargin)
    switch lower(varargin{1})
          case 'max'
              type=MAX;
          case 'first_peak'
              type=FIRST_PEAK;
          otherwise
              error(['Unexpected option: ' varargin{1}])
      end
      varargin(1:end) = [];
end



    % find max position of abs hilbert of the first echo
    [max1, pos1] = max(abs(hilbert(data(range1,1))));    
    pos1=pos1+range1(1)-1;
    
    % find position of echo 2
    % if the second echo is too close, use another method
    echo_signal2=abs(hilbert(data(range2,1)));            
    
    %if subtraction of reference ...
        % subtraction method using a reference signal from calibration
%         z=length(echo_signal2)-length(sinalRef);    
%         sinalRef(end+1:end+z)=zeros(1,z);    
%         sinalSub=echo_signal2-sinalRef;
        
%-debug--------------------------------------------
%      figure
%      hold all
%      plot(echo_signal2)
%      plot(sinalRef)
%      plot(sinalSub)
%      legend('echo_signal2','ref','sub')
%--------------------------------------------------
if(type==MAX)
    [max2, pos2] = max(echo_signal2);
    pos2=pos2+range2(1);
    %if signal is too weak, zero
    if(max2<thr)
        US.valid=0;
    end
elseif(type==FIRST_PEAK)
     [peaks,locations]=findpeaks(echo_signal2,'MinPeakProminence',thr); 
     % pks = findpeaks(data) returns a vector with the local maxima (peaks) of the input signal vector, data. 
     % A local peak is a data sample that is either larger than its two neighboring samples or is equal to Inf. Non-Inf signal endpoints are excluded.
     % If a peak is flat, the function returns only the point with the lowest index.
     % [pks,locs] = findpeaks(data) additionally returns the indices at which the peaks occur.
%-debug--------------------------------------------
%       figure
%       hold all
%       plot(echo_signal2)
%       scatter(locations,peaks)

%--------------------------------------------------
     
     if(~isempty(locations))
         %1st peak is locations(1) ... 2nd peak locations(2) ...
         %if(locations(1)<pos2-40) % seria o pico
         pos2=locations(1)+range2(1)-1;
         max2=peaks(1);
         %else
         %   pos2=pos2+offset2-1;
         %end
     else
         pos2=0;
         US.valid=0;
     end
end
       
    %return values
    US.pos1=pos1;
    US.pos2=pos2;
    US.tt=(US.pos2-US.pos1)/usObj.fs;
    US.max2=max2;    
    %US.dist=c*US.tt/2;    
    
end