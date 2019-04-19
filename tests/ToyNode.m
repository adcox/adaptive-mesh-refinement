classdef ToyNode < adaptiveMesh.Node
    
    methods(Access = public)
        function metric = getMetric(this)
            % A crescent-moon shape
            r1 = sqrt((this.state(1) - 0.07)^2 + (this.state(2) - 0.07)^2);
            r2 = sqrt(this.state(1)^2 + this.state(2)^2);
            metric = (r1 < pi/4 && r2 > pi/4 - 0.1 ).*1;
        end
    end
    
end