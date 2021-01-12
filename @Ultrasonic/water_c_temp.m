%%
%velsom_agua 
% funcao utilizada por Lubber and Graaff
%
function sound_speed = water_by_temperature(temp) 
    
sound_speed=1404.3+4.7*temp-0.04*temp^2;

