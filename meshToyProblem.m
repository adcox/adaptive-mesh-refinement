clear; clc; close all;

bounds = [-1, 1, -1, 1];
templateNode = TestNode();

%% Create the mesh and refine it
mapMesh = adaptiveMesh.Mesh();
mapMesh.setMinCellSize([5e-2, 5e-2]);
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

%% Plot Results
colors = lines(length(unique(mapData.metric)));

figure(); hold on;
colormap(colors);
keys = fieldnames(mapMesh.cellMap);
for c = 1:length(keys)
    % Plot cells
    cell = mapMesh.cellMap.(keys{c});
    if(~cell.isSubdivided)
        rectangle('position', [cell.position, cell.size],...
            'edgeColor', 'k', 'linewidth', 0.5);
    end
end
scatter(mapData.ics(:,1), mapData.ics(:,2), 16, mapData.metric, 'filled');
hold off; grid on; axis equal;