% Node - An object that represents one corner of a <a href="matlab:doc('adaptiveMesh.Cell')">Cell</a>.
%
% The Node stores a state and a metric; the state locates the node while
% the metric is the quantity of interest for the mesh refinement.
%
% Author: Andrew Cox
% Version: 9 April 2019
classdef Node < handle & matlab.mixin.Copyable
    
    properties(GetAccess = public, SetAccess = protected)
        % state - the Node state, initially an empty array
        state
    end
    
    properties(Access = protected)
        % metric - the Node metric, initially an empty array
        metric
    end
    
    methods(Access = public)
        
        function this = Node(state)
            % Node constructor
            %
            %   this = Node() constructs a default node with no properties
            %   assigned. Additional functions must be called to set the
            %   state
            %
            %   this = Node(state) constructs a node located at the
            %   specified state (a row vector).
            if(nargin > 0)
                this.state = state;
            else
                this.state = [];
            end
            
            this.metric = [];
        end
        
        function setState(this, state)
            % setState - set the full node state
            %
            %   setState(state)
            this.state = state;
        end
        
        function pos = getPosition(this)
            % getPosition - Get the position of the node
            %
            %   pos = getPosition() returns the first two elements of the
            %   state vector, i.e., the node position in 2D space
            
            if(isempty(this.state))
                error('The state has not been defined');
            end
            
            pos = this.state(1:2);
        end
        
        function metric = getMetric(this)
            % getMetric - get the metric for this node
            %
            %   metric = getMetric() returns the metric for this node.
            metric = this.metric;
        end
    end
end