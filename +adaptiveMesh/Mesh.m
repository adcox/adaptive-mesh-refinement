% Mesh - a rectangular mesh that may be refined
classdef Mesh < handle
    
    properties(GetAccess = public, SetAccess = protected)
        
        minDepth
        
        minCellSize = [1e-4, 1e-4];
        
        cellMap
        
        nodeMap
        
        vis = struct('bVisualize', false,...
            'hFig', [],...
            'hAx', [],...
            'colors', [],...
            'cLim', []);
        
        size
        position
        trueMinCellSize
    end
    
    methods(Access = public)
        
        function this = Mesh(minDepth)
            if(nargin == 0)
                this.minDepth = 3;
            else
                this.minDepth =  minDepth;
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
            % Delete all the cells
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
            %   initMesh(bounds, templateNode)
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
                templateCell.setMinWidth(this.minCellSize(1));
                templateCell.setMinHeight(this.minCellSize(2));
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
            nMaxWDivisions = max(this.minDepth,...
                floor( log(this.size(1)/this.minCellSize(1))/log(2) ));
            nMaxHDivisions = max(this.minDepth,...
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
            
            nodeGrid(1,1).setState(bounds([1,4]));   % top-left
            nodeGrid(2,1).setState(bounds([1,3]));   % bottom-left
            nodeGrid(1,2).setState(bounds([2,4]));   % top-right
            nodeGrid(2,2).setState(bounds([2,3]));   % bottom-right
            
            % Initialize level-0 cell
            cell = copy(templateCell);
            cell.setLevel(0);
            cell.setIndex([0,0]);
            cell.setNodes(nodeGrid);
            this.cellMap(cell.getKey()) = cell;
            
            if(this.vis.bVisualize)
                rectangle(this.vis.hAx,...
                    'position', [cell.position, cell.size],...
                    'edgecolor', 'k');
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
            % refine - Recursively refine the mesh
            %
            %   refine - refines the mesh using recursion, ending once the
            %   minimum cell size is reached or when all cells are uniform.
            level = 0;
            bContinue = true;
            while(bContinue)
                bContinue = level < this.minDepth;
                keys = this.cellMap.keys;
                
                for i = 1:length(keys)
                    if(level < this.minDepth || ~this.cellMap(keys{i}).isUniform())
                        [newCells, newNodes] = ...
                            this.cellMap(keys{i}).subdivide(this,...
                            level >= this.minDepth);
                        
                        if(~isempty(newCells) && ~isempty(newNodes))
                       
                            bContinue = true;
                            
                            if(this.vis.bVisualize)
                                for c = 1:length(newCells)
                                    rectangle(this.vis.hAx, 'position',...
                                        [newCells(c).position, newCells(c).size],...
                                        'edgecolor', 'k');
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
%             
%             for i = 1:length(cells)
%                 if(depth < this.minDepth || ~cells(i).isUniform())
%                     [newCells, newNodes] = cells(i).subdivide(this, depth >= this.minDepth);
% 
%                     if(~isempty(newCells) && ~isempty(newNodes))
%                        
%                         if(this.vis.bVisualize)
%                             for c = 1:length(newCells)
%                                 rectangle(this.vis.hAx, 'position',...
%                                     [newCells(c).position, newCells(c).size],...
%                                     'edgecolor', 'k');
%                             end
%                             q = zeros(length(newNodes), 3);
%                             for n = 1:length(newNodes)
%                                 q(n,:) = [newNodes(n).state(1:2), newNodes(n).getMetric()];
%                             end
%                             this.vis.cLim(1) = min(min(q(:,3)), this.vis.cLim(1));
%                             this.vis.cLim(2) = max(max(q(:,3)), this.vis.cLim(2));
%                             scatter(this.vis.hAx, q(:,1), q(:,2), 25, q(:,3), 'filled');
%                             caxis(this.vis.hAx, this.vis.cLim);
%                         end
%                         
%                         % Store the new Nodes
%                         for n = 1:length(newNodes)
%                             key = newNodes(n).getKey(this);
%                             if(this.nodeMap.isKey(key))
%                                 % Reassign "new" node to the existing node
%                                 newNodes(n) = this.nodeMap(key);
%                             else
%                                 % Store the new node
%                                 this.nodeMap(key) = newNodes(n);
%                             end
%                         end
%                         
%                         % Store the new Cells
%                         for c = 1:length(newCells)
%                             key = newCells(c).getKey();
%                             if(this.cellMap.isKey(key))
%                                 error('Duplicate key: %d', key);
%                             end
%                             this.cellMap(key) = newCells(c);
%                         end
%                         
%                         % Use recursion to continue; newCells and newNodes
%                         % will be empty once the minimum cell size is
%                         % reached, ending the recursion
%                         this.refine(newCells, depth + 1);
%                     end
%                 end
%             end
        end
        
        function setMinCellSize(this, minSize)
            % setMinCellSize
            %
            %   setMinCellSize(minSize)
            
            if(length(minSize) < 2)
                error('minSize must have at least two elements: [width, heigh]');
            end
            
            this.minCellSize = minSize;
        end
        
        function cell = getCell(this, key)
            
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
            key = adaptiveMesh.Node.computeKey(this, pos);
            node = this.getNode(key);
        end
    end
    
end