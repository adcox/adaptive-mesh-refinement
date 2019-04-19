clear; clc; close all;
addpath tests;

bounds = [-1, 1, -1, 1];
templateNode = ToyNode();

%% Create the mesh and refine it
mapMesh = adaptiveMesh.Mesh();
mapMesh.setMinCellSize([1e-2, 1e-2]);
mapMesh.initMesh(bounds, templateNode);
mapMesh.refine();

%% Extract data
keys = fieldnames(mapMesh.nodeMap);
N = length(keys);
mapData = struct();
mapData.ics = zeros(N,2);
mapData.metric = zeros(N,1);
for r = 1:N
    mapData.ics(r,:) = mapMesh.nodeMap.(keys{r}).state;
    mapData.metric(r) = mapMesh.nodeMap.(keys{r}).getMetric();
end

%% Check Data
cellKeys = fieldnames(mapMesh.cellMap);
failingKeys = {};
for c = 1:length(cellKeys)   
    cell = mapMesh.cellMap.(cellKeys{c});
        
    if(~cell.isSubdivided)
        neighbors = cell.getNeighbors(mapMesh);

        for n = 1:length(neighbors)
            if(abs(neighbors(n).level - cell.level) > 1)
                failingKeys{end+1} = cellKeys{c};
                break;
            end
    %         assert(abs(neighbors(n).level - cell.level) <= 1)
        end
    end
end

%% Plot Results
colors = lines(length(unique(mapData.metric)));

hFig = figure(); hold on;
colormap(colors);
hcb = colorbar;

keys = fieldnames(mapMesh.cellMap);
for c = 1:length(keys)
    % Plot cells
    cell = mapMesh.cellMap.(keys{c});
    if(~cell.isSubdivided)
        rectangle('position', [cell.position, cell.size],...
            'edgeColor', 'k', 'linewidth', 0.5);
        if(sum(cellfun(@(a)strcmp(keys{c}, a), failingKeys)) > 0)
            plot([cell.position(1), cell.position(1) + cell.size(1)], ...
                [cell.position(2), cell.position(2) + cell.size(2)], 'r');
        end
    end
end
scatter(mapData.ics(:,1), mapData.ics(:,2), 16, mapData.metric, 'filled');
hold off; grid on; axis equal;

%% Get Image
figure();
colormap(colors);
hcb = colorbar;
[xData, yData, CData] = mapMesh.toImage();
image(xData, yData, CData, 'CDataMapping', 'scaled');
axis equal;
set(gca, 'XDir', 'norma', 'YDir', 'normal');