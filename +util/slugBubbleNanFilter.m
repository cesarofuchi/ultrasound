function filtData=slugBubbleNanFilter(data, positions)
%% ------------------------------------------------------
% after a certain position transform all values to NaN
% ex: data=1,2,3,4,5,6
% position=4;
% filtData=1,2,3,NaN,NaN,NaN
%--------------------------------------------------------

filtData=data;
for i=1:length(positions)
    filtData(positions(i):end,i)=NaN;
end


