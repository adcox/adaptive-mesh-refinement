% Cell - A collection of four <a href="matlab:doc('adaptiveMesh.Node')">Node</a> objects.
% 
% A cell is the basic unit of the mesh. The border of the rectangular cell
% is defined by four <a href="matlab:doc('adapativeGrid.Node')">Node</a>s.
% The cell is considered "uniform" if the <a href="matlab:doc('adapativeGrid.Node/metric')">metric</a>
% values associated with each of the nodes are all identical. A
% "nonuniform" Cell is subdivided via the "subdivide" function.
%
% Author: Andrew Cox
% Version: 9 April 2019
classdef Cell < handle & matlab.mixin.Copyable
    
    properties(GetAccess = public, SetAccess = protected)
        % nodes - a matrix of nodes located at the corners of this cell
        %
        % The nodes are indexed by their location in cartesian space, i.e.,
        %
        %   (1,1) ---- (1,2)
        %     |          |
        %     |          |
        %     |          |
        %   (2,1) ---- (2,2)
        nodes
        
        % cellPosition - Location of the top-left node (1,1)
        cellPosition
        
        % cellSize - the [width, height] of the cell in position
        % coordinates
        cellSize
        
        % Minimum acceptable width
        %
        % The cell wil not subdivide if doing so creates a cell with a
        % width smaller than this size
        minWidth = 1e-4;
        
        % Minimum acceptable height
        %
        % The cell will not subidivide if doing so creates a cell with a
        % height smaller than this size
        minHeight = 1e-4;
        
        % isSubdivided - whether or not this cell has been subdivided
        isSubdivided = false;
    end
    
    methods(Access = public)
        
        function this = Cell(nodes)
            % Cell constructor
            %
            % this = Cell() constructs an empty cell with no nodes
            %
            % this = Cell(nodes) accepts a 2x2 matrix of Node objects,
            % "nodes", with the nodes located in space consistent with
            % their location in the matrix, i.e., nodes(1,1) represents the
            % top-left corner of the Cell, nodes(1,2) represents the
            % top-right corner, nodes(2,1) is the bottom-left corner, and
            % nodes(2,2) is the bottom-right corner.
            
            if(nargin > 0)
                this.setNodes(nodes)
            end
        end
        
        function setNodes(this, nodes)
            % setNodes - set the nodes that define the Cell boundary
            %
            %   setNodes(nodes) accepts a 2x2 matrix of Node objects,
            % "nodes", with the nodes located in space consistent with
            % their location in the matrix, i.e., nodes(1,1) represents the
            % top-left corner of the Cell, nodes(1,2) represents the
            % top-right corner, nodes(2,1) is the bottom-left corner, and
            % nodes(2,2) is the bottom-right corner.
            
            if(~isa(nodes, 'adaptiveMesh.Node'))
                error('Expecting matrix of Node objects');
            end

            if(sum(size(nodes) == [2,2]) ~= 2)
                error('Expecting 2x2 matrix');
            end

            this.nodes = nodes;

            % Compute the size of the node
            r1 = this.nodes(1,1).getPosition();
            r2 = this.nodes(1,2).getPosition();
            r3 = this.nodes(2,1).getPosition();
            this.cellPosition = r3;
            this.cellSize = [r2(1) - r1(1), r1(2) - r3(2)];

            if(this.cellSize(1) < 0 || this.cellSize(2) < 0)
                error('Nodes are not in correct order');
            end
        end
        
        function [newCells, newNodes] = subdivide(this)
            % subdivide - Split the cell into four smaller cells
            %
            %   [newCells, newNodes] = subdivide() subdivides the cell into
            %   four smaller cells as long as the new Cells do not violate
            %   the minimum width or minimum height properties of this
            %   Cell. The new Cells and new Nodes are returned 
            if(isempty(this.nodes))
                error('Nodes matrix is empty');
            end
            
            if(this.cellSize(1)/2 > this.minWidth &&...
                    this.cellSize(2)/2 > this.minHeight)
                
                % Coordinates of the center of the cell
                xMid = this.cellPosition(1) + 0.5*this.cellSize(1);
                yMid = this.cellPosition(2) + 0.5*this.cellSize(2);
                
                % Make copies of one Node and this Cell to preserve their
                % properties in the new Nodes and Cells. The new Nodes and
                % Cells are organized as labeled below, with existing nodes
                % represeted by "o", new nodes represented by numbers, and
                % new cells represented by numbers:
                %
                %   o --- 1 --- o
                %   |  A  |  B  |
                %   3 --- 0 --- 4
                %   |  C  |  D  |
                %   o --- 2 --- o   
                node0 = copy(this.nodes(1,1));
                node1 = copy(node0);
                node2 = copy(node0);
                node3 = copy(node0);
                node4 = copy(node0);
                
                cellA = copy(this);
                cellB = copy(cellA);
                cellC = copy(cellA);
                cellD = copy(cellA);
                
                % Initialize new nodes (#) between the exiting nodes (o):        
                node0.setState([xMid, yMid]);
                node1.setState([xMid, this.cellPosition(2) + this.cellSize(2)]);
                node2.setState([xMid, this.cellPosition(2)]);
                node3.setState([this.cellPosition(1), yMid]);
                node4.setState([this.cellPosition(1) + this.cellSize(1), yMid]);
                
                % Initialize new Cells, labeled in the ascii art above
                cellA.setNodes([this.nodes(1,1), node1; node3, node0]);
                cellB.setNodes([node1, this.nodes(1,2); node0, node4]);
                cellC.setNodes([node3, node0; this.nodes(2,1), node2]);
                cellD.setNodes([node0, node4; node2, this.nodes(2,2)]);

                % Output the new cells and nodes
                this.isSubdivided = true;
                newCells = [cellA, cellB, cellC, cellD];
                newNodes = [node1, node3, node0, node4, node2];
            else
                this.isSubdivided = false;
                newCells = [];
                newNodes = [];
            end
        end
        
        function isUniform = isUniform(this)
            % isUniform - Determine if the cell is uniform
            %
            %   tf = isUniform() returns true if the metrics associated with
            %   all four Nodes are identical and false otherwise
            
            value = this.nodes(1,1).getMetric();
            for r = 1:2
                for c = 1:2
                    if(this.nodes(r,c).getMetric() ~= value)
                        isUniform = false;
                        return;
                    end
                end
            end
            isUniform = true;
        end
        
        function setMinWidth(this, minWidth)
            % setMinWidth - set the minimum acceptable Cell width
            %
            %   setMinWidth(minWidth) defines the minimum cell width, i.e.,
            %   the smallest width a cell can have. This width limits the
            %   "subdivide()" function
            %
            % See also: setMinHeight, subdivide
            this.minWidth = minWidth;
        end
        
        function setMinHeight(this, minHeight)
            % setMinHeight - set the minimum acceptable Cell height
            %
            %   setMinHeight(minHeight) defines the minimum cell height,
            %   i.e., the smallest height a cell can have. This height
            %   limits the "subdivide() function
            %
            % See also: setMinWidth, subdivide
            this.minHeight = minHeight;
        end
    end % End of public methods
    
    methods(Access = protected)
        function cpObj = copyElement(this)
            
            % Make a shallow copy
            cpObj = copyElement@matlab.mixin.Copyable(this);
            
            % Copy the handle REFERENCES; don't want to make copy of the
            % handles themselves to save memory
            cpObj.nodes = this.nodes;
            
        end
    end
end