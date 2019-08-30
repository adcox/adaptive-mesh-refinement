clear; clc; close all;
addpath tests;  % The ToyNode class is located here

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
meshData = struct();
meshData.pos = zeros(N,2);
meshData.metric = zeros(N,1);
for r = 1:N
    meshData.pos(r,:) = mapMesh.nodeMap.(keys{r}).state;
    meshData.metric(r) = mapMesh.nodeMap.(keys{r}).getMetric();
end

%% Plot Results
colors = lines(length(unique(meshData.metric)));

hFig = figure(); hold on;
colormap(colors);
colorbar;

keys = fieldnames(mapMesh.cellMap);
for c = 1:length(keys)
    % Plot cells
    cell = mapMesh.cellMap.(keys{c});
    if(~cell.isSubdivided)
        rectangle('position', [cell.position, cell.size],...
            'edgeColor', 'k', 'linewidth', 0.5);
    end
end
% Plot the nodes, colored by their metric
scatter(meshData.pos(:,1), meshData.pos(:,2), 16, meshData.metric, 'filled');
hold off; grid on; axis equal;
xlabel('x');
ylabel('y');

%% Plot the mesh as an image
figure();
colormap(colors);
hcb = colorbar;
[xData, yData, CData] = mapMesh.toImage();
image(xData, yData, CData, 'CDataMapping', 'scaled');
axis equal;
set(gca, 'XDir', 'normal', 'YDir', 'normal');