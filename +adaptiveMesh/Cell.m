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
        
        level
        
        index
        
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
        
        
        % position - Location of the top-left node (1,1)
        position
        
        % size - the [width, height] of the cell in position
        % coordinates
        size
        
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
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public Methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
            this.isSubdivided = false;
            this.index = [0,0];
            this.level = 0;
        end
        
        function delete(this)
            % Cell destructor
            for n = 1:length(this.nodes)
                if(isvalid(this.nodes(n)))
                    delete(this.nodes(n));
                end
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
            this.position = r3;
            this.size = [r2(1) - r1(1), r1(2) - r3(2)];

            if(this.size(1) < 0 || this.size(2) < 0)
                error('Nodes are not in correct order');
            end
        end
        
        function [newCells, newNodes] = subdivide(this, mesh, bCueNeighbors)
            % subdivide - Split the cell into four smaller cells
            %
            %   [newCells, newNodes] = subdivide() subdivides the cell into
            %   four smaller cells as long as the new Cells do not violate
            %   the minimum width or minimum height properties of this
            %   Cell. The new Cells and new Nodes are returned 
            
            newCells = [];
            newNodes = [];
                
            % Skip if this cell has already been subdivided
            if(this.isSubdivided)
                return;
            end
            
            if(isempty(this.nodes))
                error('Nodes matrix is empty');
            end
            
            if(~isa(mesh, 'adaptiveMesh.Mesh'))
                error('Expecting an adaptiveMesh.Mesh');
            end
            
            if(nargin < 3)
                bCueNeighbors = false;
            elseif(~islogical(bCueNeighbors))
                error('bNeighborCatalyst must be logical (a boolean)');
            end
            
            if(this.size(1)/2 > this.minWidth &&...
                    this.size(2)/2 > this.minHeight)
                
                % Coordinates of the center of the cell
                xMid = this.position(1) + 0.5*this.size(1);
                yMid = this.position(2) + 0.5*this.size(2);
                
                % Make copies of one Node and this Cell to preserve their
                % properties in the new Nodes and Cells. The new Nodes and
                % Cells are organized as labeled below, with existing nodes
                % represeted by "o", new nodes represented by numbers, and
                % new cells represented by numbers:
                %
                %   o --- 2 --- o
                %   |  A  |  B  |
                %   4 --- 1 --- 5
                %   |  C  |  D  |
                %   o --- 3 --- o
                newNodes = [copy(this.nodes(1,1)), copy(this.nodes(1,1)), ...
                    copy(this.nodes(1,1)), copy(this.nodes(1,1)), ...
                    copy(this.nodes(1,1))];
                
                newNodePos = [xMid, yMid;
                    xMid, this.position(2) + this.size(2);
                    xMid, this.position(2);
                    this.position(1), yMid
                    this.position(1) + this.size(1), yMid];

                for n = 1:5
                    newNodes(n).setState(newNodePos(n,:));
                end
                
                childIndices = this.getChildIndices();
                
                cellA = copy(this);
                cellA.setLevel(this.level + 1);
                
                cellB = copy(cellA);
                cellC = copy(cellA);
                cellD = copy(cellA);
                
                % Initialize new Cells, labeled in the ascii art above
                cellA.setNodes([this.nodes(1,1), newNodes(2); newNodes(4), newNodes(1)]);
                cellA.setIndex(childIndices{1,1});
                cellB.setNodes([newNodes(2), this.nodes(1,2); newNodes(1), newNodes(5)]);
                cellB.setIndex(childIndices{1,2});
                cellC.setNodes([newNodes(4), newNodes(1); this.nodes(2,1), newNodes(3)]);
                cellC.setIndex(childIndices{2,1});
                cellD.setNodes([newNodes(1), newNodes(5); newNodes(3), this.nodes(2,2)]);
                cellD.setIndex(childIndices{2,2});
                
                % Output the new cells and nodes
                this.isSubdivided = true;
                newCells = [cellA, cellB, cellC, cellD];
                
                % Get the neighbors of this cell that are larger or the
                % same size
                if(bCueNeighbors)
                    neighbors = this.getNeighbors(mesh, true);
                    for n = 1:length(neighbors)
                        [otherCells, otherNodes] = neighbors(n).subdivide(mesh, false);
                        if(~isempty(otherCells))
                            newCells(end+1:end+length(otherCells)) = otherCells;
                        end
                        if(~isempty(otherNodes))
                            newNodes(end+1:end+length(otherNodes)) = otherNodes;
                        end
                    end
                end
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
        
        function setLevel(this, level)            
            this.level = level;
        end
        
        function setIndex(this, index)
            if(size(index,1)*size(index,2) ~= 2)
                error('Expecting 2-element array');
            end
            this.index = index;
        end
        
        function childIndices = getChildIndices(this)
            % getChildIndices - returns the indices of the child cells
            %
            %   childIndices = getChildIndices returns a 2x2 matrix of
            %   2-element arrays; each aray is the index of a child cell,
            %   arranged in spatial order:
            %
            %       o --- o --- o
            %       |  A  |  B  |
            %       o --- o --- o
            %       |  C  |  D  |
            %       o --- o --- o 
            
            if(length(this.index) < 2)
                keyboard;
            end
            
            childIndices = cell(2,2);
            childIndices{1,1} = [2*this.index(1) + 0, 2*this.index(2) + 1];
            childIndices{1,2} = [2*this.index(1) + 1, 2*this.index(2) + 1];
            childIndices{2,1} = [2*this.index(1) + 0, 2*this.index(2) + 0];
            childIndices{2,2} = [2*this.index(1) + 1, 2*this.index(2) + 0];
        end
        
        function parentIndex = getParentIndex(this)
            parentIndex = adaptiveMesh.Cell.computeParentIndex(this.index);
        end
        
        function key = getKey(this)
            key = adaptiveMesh.Cell.computeKey(this.level, this.index);
        end
        
        function neighbors = getNeighbors(this, mesh, bNoSmaller)
            sameLevelNeighbors = [this.index(1) - 1, this.index(2); % left
                this.index(1), this.index(2) + 1; % top
                this.index(1) + 1, this.index(2); % right
                this.index(1), this.index(2) - 1]; % bottom
            
            % The largest index for a cell at the same level
            maxIndex = 2^this.level - 1;
            
            neighbors = this;
            neighbors(1) = [];
            
            for n = 1:4
                if(sum(sameLevelNeighbors(n,:) < 0) == 0 && ...
                        sum(sameLevelNeighbors(n,:) > maxIndex) == 0)
                    neighbor = mesh.getCell(...
                        adaptiveMesh.Cell.computeKey(this.level, sameLevelNeighbors(n,:)));
                else
                    neighbor = [];
                end
                
                if(~isempty(neighbor))
                    if(~bNoSmaller && neighbor.isSubdivided)
                        kidsIx = neighbor.getChildIndices();

                        switch(n)
                            case 1
                                % Neighbor is to the left, get right-side
                                neighbors(end+1:end+3) =...
                                    [mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(1,2))),...
                                    mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(2,2)))];
                            case 2
                                % Neighbor is above, get bottom-side
                                neighbors(end+1:end+3) =...
                                    [mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(2,1))),...
                                    mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(2,2)))];
                            case 3
                                % Neighbor is to the right, get left-side
                                neighbors(end+1:end+3) =...
                                	[mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(1,1))),...
                                    mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(2,1)))];
                            case 4
                                % Neighbor is below, get top-side
                                neighbors(end+1:end+3) =...
                                    [mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(1,1))),...
                                    mesh.getCell(adaptiveMesh.Cell.computeKey(this.level+1, kidsIx(1,2)))];
                        end
                        
                    else
                        neighbors(end+1) = neighbor;
                    end
                else
                    L = this.level - 1;
                    myIndex = this.index;
                    neighborIndex = sameLevelNeighbors(n,:);
                    while(L > 0)
                        maxParentIndex = 2^L - 1;
                        
                        % Try the parant of the non-existent neighboring
                        % cell
                        myParentIndex = this.computeParentIndex(myIndex);
                        neighborParentIndex = this.computeParentIndex(neighborIndex);
                        
                        % If neighbor has the same parent as me
                        if(sum(myParentIndex == neighborParentIndex) == 2)
                            break;
                        end
                        
                        % Make sure the neighbor parent index is in bounds
                        if(sum(neighborParentIndex < 0) == 0 &&...
                                sum(neighborParentIndex > maxParentIndex) == 0)
                            
                            neighborParent = mesh.getCell(...
                                this.computeKey(L, neighborParentIndex));
                            
                            if(~isempty(neighborParent))
                                % If the lower-level cell exists, output it
                                neighbors(end+1) = neighborParent;
                                % exit the loop since we found a neighbor;
                                % if no neighbor is found, continue looking
                                % up at higher levels
                                break;
                            end
                        end
                        myIndex = myParentIndex;
                        neighborIndex = neighborParentIndex;
                        L = L - 1;
                    end
                end
            end
            
            % Debugging: plot neighbors:
%             figure(); hold on;
%             rectangle('position', [this.position, this.size], 'edgecolor', 'b',...
%                 'linewidth', 2);
%             for n = 1:length(neighbors)
%                 rectangle('position', [neighbors(n).position, neighbors(n).size], ...
%                     'edgecolor', 'r', 'linewidth', 1);
%             end
%             keyboard;
        end % End of getNeighbors()
    end % End of public methods
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static Public Methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static, Access = public)
        function key = computeKey(level, index)
            key = 0;
            for l = 0:level-1
                key = key + (2^l * 2^l);
            end
            key = int32(key + index(1)*2^level + index(2));
        end
        
        function parentIndex = computeParentIndex(childIndex)
            parentIndex = [floor(childIndex(1)/2), floor(childIndex(2)/2)];
        end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected Methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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