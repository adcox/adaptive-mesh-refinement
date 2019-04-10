% Mesh - a rectangular mesh that may be refined
classdef Mesh < handle
    
    properties(GetAccess = public, SetAccess = protected)
        minCellSize = [1e-4, 1e-4];
        
        allNodes
        
        cellMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    
    methods(Access = public)
        
        function this = Mesh()
            
        end
        
        function initMesh(this, nodeMesh, templateCell)
            % initMesh - initialize the mesh
            %
            %   initMesh(nodeMesh) initializes the mesh from a rectangular
            %   grid of nodes, "nodeMesh." This grid, represented as a 2D
            %   array (i.e., a matrix), must be oriented consistent with
            %   the positions of the nodes: nodes(1,1) is the top-left
            %   corner of the mesh, nodes(1,2) has a larger x-coordinate
            %   but equal y-coordinate as nodes(1,1), and nodes(2,1) has a
            %   smaller y-coordinate but equal x-coordinate as nodes(1,1).
            %   Rectangular cells are created for every set of four nodes.
            %
            %   initMesh(nodeMesh, templateCell) initializes the mesh and
            %   uses the specified "templateCell" as the template for each
            %   cell added to the mesh. This templateCell must be derived
            %   from adaptiveMesh.Cell
            
            if(~isa(nodeMesh, 'adaptiveMesh.Node'))
                error(['Node mesh must contain adaptiveMesh.Node objects,',...
                    ' or derivatives thereof']);
            end
            
            S = size(nodeMesh);
            
            if(S(1)*S(2) < 4)
                error('Must have at least four nodes in the mesh');
            end
            
            % Create a template cell with the specified properties
            if(nargin < 3)
                templateCell = adaptiveMesh.Cell();
                templateCell.setMinWidth(this.minCellSize(1));
                templateCell.setMinHeight(this.minCellSize(2));
            else
                if(~isa(templateCell, 'adaptiveMesh.Cell'))
                    error('Template Cell must be derived from adaptiveMesh.Cell');
                end
            end
                        
            % Populate cell map
            for r = 1:(S(1) - 1)
                for c = 1:(S(2) - 1)
                    % Copy template to preserve all properties
                    cell = copy(templateCell);
                    
                    % Assume the node mesh is oriented properly in position
                    % space
                    cell.setNodes(nodeMesh(r:r+1, c:c+1));
                    cell.setIndex([S(1) - 1 - r, c - 1 ]);
                    
                    this.cellMap(cell.getKey()) = cell;
                end
            end
            
            % Store nodes
            this.allNodes = reshape(nodeMesh, 1, S(1)*S(2));
        end
        
        function refine(this, cells)
            % refine - Recursively refine the mesh
            %
            %   refine - refines the mesh using recursion, ending once the
            %   minimum cell size is reached or when all cells are uniform.
            %
            %   refine(cells) - refines only the specified cell array of
            %   Cell objects
            
            % When called without arguments, run on all the cells in the
            % mesh
            if(nargin == 1)
                cells_cellArray = this.cellMap.values;
                cells(length(cells_cellArray)) = cells_cellArray{end};
                for c = 1:length(cells_cellArray)
                    cells(c) = cells_cellArray{c};
                end
            else
                if(~isa(cells, 'adaptiveMesh.Cell'))
                    error('Cells must be derived from adpativeMesh.Cell');
                end
            end
            
            for i = 1:length(cells)
                if(~cells(i).isUniform())
                    [newCells, newNodes] = cells(i).subdivide(this);

                    if(~isempty(newCells) && ~isempty(newNodes))
                                               
                        % Save the new Nodes and Cells
                        this.allNodes(end+1:end+length(newNodes)) = newNodes;
                        
                        for c = 1:length(newCells)
                            this.cellMap(newCells(c).getKey()) = newCells(c);
                        end
                        % Use recursion to continue; newCells and newNodes
                        % will be empty once the minimum cell size is
                        % reached, ending the recursion
                        this.refine(newCells);
                    end
                end
            end
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
            
            if(~isa(key, 'char'))
                error('Key must be a character array');
            end
            
            if(this.cellMap.isKey(key))
                cell = this.cellMap(key);
            else
                cell = [];
            end
        end
    end
    
end