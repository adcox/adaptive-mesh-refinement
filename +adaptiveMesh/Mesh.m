% Mesh - a rectangular mesh that may be refined
classdef Mesh < handle
    
    properties(GetAccess = public, SetAccess = protected)
        
        % minLevel - the minimum level for refinement
        %
        % The level increases as the mesh is refined; thus, the minimum
        % level specifies the minimum amount of refinement. The mesh always
        % begins at level 0 with a single cell, so the minimum number of
        % intermediate points is 2^minLevel - 1
        minLevel
        
        % minCellSize - the minimum allowable cell size
        %
        % A 2-entry array specifying the minimum [width, height] of a cell.
        % If a subdivision would reduce a cell to a size smaller than one
        % or both of these dimensions, the subdivision is not conducted.
        minCellSize = [1e-4, 1e-4];
        
        % cellMap - maps an integer key to a Cell object
        %
        % Each cell is represented by a unique key, available from the
        % Cell.computeKey() or cell.getKey() functions
        cellMap
        
        % nodeMap - maps a character array to a Node object
        %
        % Since the mesh always begins with a single cell and every
        % subdivision divides an existing cell into four equally-sized
        % "sub-cells," the location of a node is directly available from a
        % set of indices (ix, iy). The maximum values of ix and iy are
        % computed from the quotient of the mesh size and the true minimum
        % cell size. Accordingly, the key is a string "%d,%d" where the
        % first integer is ix and the second integer is iy.
        nodeMap
        
        % vis - a structure of options for visualization
        vis = struct('bVisualize', true,...
            'hFig', [],...
            'hAx', [],...
            'colors', [],...
            'cLim', []);
        
        % size - the size = [width, height] of the mesh
        size
        
        % position - the position [x0, y0] of the mesh
        position
        
        % trueMinCellSize - the true minimum cell size
        %
        % The "true" minimum cell size is generally larger than
        % "minCellSize" and is available from size./[2^NX, 2^NY] where NX
        % and NY are the maximum number of subdivisions that can be
        % completed before a resulting "sub-cell" is smaller than the
        % "minCellSize."
        trueMinCellSize
    end
    
    methods(Access = public)
        
        function this = Mesh(minLevel)
            % Mesh - Construct a mesh
            %
            %   this = Mesh() constructs a mesh with a default minimum
            %   level of 3
            %
            %   this = Mesh(minLevel) constructs a mesh with the specified
            %   minimum level
            %
            % See also: minLevel
            
            if(nargin == 0)
                this.minLevel = 3;
            else
                this.minLevel =  minLevel;
            end
            
            this.cellMap = containers.Map('KeyType', 'int32', 'ValueType', 'any');
            this.nodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            if(this.vis.bVisualize)
                this.vis.hFig = figure();
                this.vis.hAx = gca;
                this.vis.colors = lines(7);
                colormap(this.vis.colors);
                title(this.vis.hAx, 'Mesh Visualization');
                hold(this.vis.hAx, 'on');
                colorbar(this.vis.hAx);
            end
        end
        
        function delete(this)
            % delete - destructor: delete all the nodes and cells
            keys = this.cellMap.keys;
            for k = 1:length(keys)
                cell = this.cellMap(keys{k});
                if(isvalid(cell))
                    delete(cell);
                end
            end
            delete(this.cellMap);
            
            keys = this.nodeMap.keys;
            for k = 1:length(keys)
                node = this.nodeMap(keys{k});
                if(isvalid(node))
                    delete(node);
                end
            end
            delete(this.nodeMap);
        end
        
        function initMesh(this, bounds, templateNode, templateCell)
            % initMesh - initialize the mesh
            %
            %   initMesh(bounds, templateNode) initializes the mesh to the
            %   specified bounds = [minX, maxX, minY, maxY] using the
            %   specified "templateNode" to construct all nodes added to 
            %   the mesh.
            %
            %   initMesh(bounds, templateNode, templateCell) initializes 
            %   the mesh and uses the specified "templateCell" as the 
            %   template for each cell added to the mesh. This templateCell
            %   must be derived from adaptiveMesh.Cell
            
            if(length(bounds) ~= 4)
                error('Bounds must have four elements');
            end
            
            if(~isa(templateNode, 'adaptiveMesh.Node'))
                error('Template Node must be derived from adaptiveMesh.Node');
            end
            
            if(nargin < 4)
                % If no template cell is supplied, use the base class
                templateCell = adaptiveMesh.Cell();
                templateCell.setMinSize(this.minCellSize);
            else
                % If a template cell IS supplied, check its class
                if(~isa(templateCell, 'adaptiveMesh.Cell'))
                    error('Template Cell must be derived from adaptiveMesh.Cell');
                end
            end
                   
            % Compute size of mesh: [width, height]
            this.size = [bounds(2) - bounds(1), bounds(4) - bounds(3)];
            % Compute mesh position (bottom-left): [x0, y0]
            this.position = [bounds(1), bounds(3)];
            % Get maximum number of subidivisions of width and height
            nMaxWDivisions = max(this.minLevel,...
                floor( log(this.size(1)/this.minCellSize(1))/log(2) ));
            nMaxHDivisions = max(this.minLevel,...
                floor(log(this.size(2)/this.minCellSize(2))/log(2) ));
            nMaxPtsX = 2^nMaxWDivisions;
            nMaxPtsY = 2^nMaxHDivisions;
            
            % Compute true minimum width and height
            this.trueMinCellSize = this.minCellSize./[nMaxPtsX, nMaxPtsY];
            
            % Initialize four boundary nodes
            nodeGrid(2,2) = copy(templateNode);
            nodeGrid(1,1) = copy(templateNode);
            nodeGrid(1,2) = copy(templateNode);
            nodeGrid(2,1) = copy(templateNode);
            
            nodeGrid(1,1).setState(bounds([1,3]));   % bottom-left
            nodeGrid(2,1).setState(bounds([2,3]));   % bottom-right
            nodeGrid(1,2).setState(bounds([1,4]));   % top-left
            nodeGrid(2,2).setState(bounds([2,4]));   % top-right
            
            % Initialize level-0 cell
            cell = copy(templateCell);
            cell.setLevel(0);
            cell.setIndex([0,0]);
            cell.setNodes(nodeGrid);
            this.cellMap(cell.getKey()) = cell;
            
            if(this.vis.bVisualize)
                rectangle(this.vis.hAx,...
                    'position', [cell.position, cell.size],...
                    'edgecolor', 'k', 'linewidth', 1);
            end
            
            
            % Store nodes
            for r = 1:2
                for c = 1:2
                    this.nodeMap(nodeGrid(r,c).getKey(this)) = nodeGrid(r,c);
                end
            end
            
            if(this.vis.bVisualize)
                q = zeros(4, 3);
                n = 1;
                for r = 1:2
                    for c = 1:2
                        q(n,1:2) = nodeGrid(r,c).getPosition();
                        q(n,3) = nodeGrid(r,c).getMetric();
                        n = n + 1;
                    end
                end
                this.vis.cLim = [min(q(:,3)), max(q(:,3))];
                scatter(q(:,1), q(:,2), 25, q(:,3), 'filled');
                caxis(this.vis.hAx, this.vis.cLim);
            end
        end
        
        function refine(this)
            % refine - Refine the mesh
            %
            %   refine - refines the mesh until either the minimum level
            %   has been reached or until no further subdivision is
            %   possible and/or necessary.
            
            level = 0;
            bContinue = true;
            while(bContinue)
                bContinue = level < this.minLevel;
                keys = this.cellMap.keys;
                
                for i = 1:length(keys)
                    
                    if(this.cellMap(keys{i}).isSubdivided)
                        % No need to check an already subdivided cell
                        continue;
                    end
                    
                    if(level < this.minLevel || ~this.cellMap(keys{i}).isUniform())
                        
                        % Subdivide the non-uniform cell
                        [newCells, newNodes] = ...
                            this.cellMap(keys{i}).subdivide(this,...
                            level >= this.minLevel);
                        
                        % If the subdivision succeeded, save the results;
                        % The most common reason that newCells and newNodes
                        % are empty is that the minimum cell size is
                        % reached and the subdivision does not occur.
                        if(~isempty(newCells) && ~isempty(newNodes))
                       
                            % Keep refining - there is work to do!
                            bContinue = true;
                            
                            % Plot debugging stuff
                            if(this.vis.bVisualize)
                                for c = 1:length(newCells)
                                    rectangle(this.vis.hAx, 'position',...
                                        [newCells(c).position, newCells(c).size],...
                                        'edgecolor', 'k', 'linewidth', 1);
                                end
                                q = zeros(length(newNodes), 3);
                                for n = 1:length(newNodes)
                                    q(n,:) = [newNodes(n).state(1:2), newNodes(n).getMetric()];
                                end
                                this.vis.cLim(1) = min(min(q(:,3)), this.vis.cLim(1));
                                this.vis.cLim(2) = max(max(q(:,3)), this.vis.cLim(2));
                                scatter(this.vis.hAx, q(:,1), q(:,2), 25, q(:,3), 'filled');
                                caxis(this.vis.hAx, this.vis.cLim);
                            end

                            % Store the new Nodes
                            for n = 1:length(newNodes)
                                key = newNodes(n).getKey(this);
                                if(this.nodeMap.isKey(key))
                                    % Reassign "new" node to the existing node
                                    newNodes(n) = this.nodeMap(key);
                                else
                                    % Store the new node
                                    this.nodeMap(key) = newNodes(n);
                                end
                            end

                            % Store the new Cells
                            for c = 1:length(newCells)
                                key = newCells(c).getKey();
                                if(this.cellMap.isKey(key))
                                    error('Duplicate key: %d', key);
                                end
                                this.cellMap(key) = newCells(c);
                            end
                        end
                    end
                end % End of loop through cells
                
                level = level + 1;
            end
        end
        
        function setMinCellSize(this, minSize)
            % setMinCellSize - set the minimum cell size
            %
            %   setMinCellSize(minSize) - sets the minimum cell size to the
            %   specified "minSize" vector. This must be called before
            %   "initMesh"
            
            if(length(minSize) < 2)
                error('minSize must have at least two elements: [width, heigh]');
            end
            
            this.minCellSize = minSize;
        end
        
        function cell = getCell(this, key)
            % getCell - retrieve a cell via its key
            %
            %   cell = getCell(key) returns a cell identified by the "key".
            %   If no cell matching that key exists, the "cell" return
            %   value is an empty array.
            
            if(~isa(key, 'int32'))
                error('Key must be an int32');
            end
            
            if(this.cellMap.isKey(key))
                cell = this.cellMap(key);
            else
                cell = [];
            end
        end
        
        function node = getNode(this, key)
            % getNode - get a node via its key
            %
            %   node = getNode(key) returns a key identified by the "key".
            %   If no node matching that key exists, the "node" return
            %   value is an empty array.
            
            if(~isa(key, 'char'))
                error('Key must be a char array');
            end
            
            if(this.nodeMap.isKey(key))
                node = this.nodeMap(key);
            else
                node = [];
            end
        end
        
        function node = getNodeFromPos(this, pos)
            % getNodeFromPos - get a node given its (x,y) position
            %
            %   node = getNodeFromPos(pos) returns a node located at the
            %   (x,y) position, "pos." If no node is located at the
            %   specified position, the return value, "node", is an empty
            %   array.
            key = adaptiveMesh.Node.computeKey(this, pos);
            node = this.getNode(key);
        end
    end
    
end