%
% Pulse class 
%
classdef Pulse < handle
    
    properties
        time
        f
    end
    
    methods
        function Pu = Pulse(time, f)
            Pu.time = time;
            Pu.f = f;
        end
    end
    
end
