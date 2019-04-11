classdef TestNode < adaptiveMesh.Node
    
    methods(Access = public)
        function metric = getMetric(this)
            r = sqrt(this.state(1)^2 + this.state(2)^2);
            metric = (r < pi/4 && r > pi/5 ).*1;
        end
    end
    
end