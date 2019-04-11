% Node - An object that represents one corner of a <a href="matlab:doc('adaptiveMesh.Cell')">Cell</a>.
%
% The Node stores a state and a metric; the state locates the node while
% the metric is the quantity of interest for the mesh refinement.
%
% Author: Andrew Cox
% Version: 9 April 2019
classdef Node < handle & matlab.mixin.Copyable
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties(GetAccess = public, SetAccess = protected)
        % state - the Node state, initially an empty array
        state
    end
    
    properties(Access = protected)
        % metric - the Node metric, initially an empty array
        metric
        
        % key - unique key representing the Node
        key
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            this.key = [];
        end
        
        function setState(this, state)
            % setState - set the full node state
            %
            %   setState(state)
            this.state = state;
            this.metric = []; % reset - must be recomputed
            this.key = []; % reset - must be recomputed
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
        
        function key = getKey(this, mesh)
            if(isempty(this.key))
                this.key = adaptiveMesh.Node.computeKey(mesh, this.getPosition());
            end
            key = this.key;
        end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Public Methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static, Access = public)
        function key = computeKey(mesh, position)
            widthIndex = round( (position(1) - mesh.position(1))/mesh.trueMinCellSize(1), 0);
            heightIndex = round( (position(2) - mesh.position(2))/mesh.trueMinCellSize(2), 0);
            key = sprintf('i%d_%d', widthIndex, heightIndex);
        end
    end
end