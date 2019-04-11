%% Main function to generate tests
function results = testCell
    results = functiontests(localfunctions);
end

%% Fixtures
function setupOnce(testCase)
    xGrid = linspace(-1, 1, 2);
    yGrid = linspace(1, -1, 2);

    [XX, YY] = meshgrid(xGrid,  yGrid);
    
    templateNode = TestNode();
    nodes(2, 2) = templateNode;
    for r = 1:2
        for c = 1:2
            % Make a copy of templateNode and then set the state

            nodes(r,c) = copy(templateNode);
            nodes(r,c).setState([XX(r,c), YY(r,c)]);
        end
    end
    testCase.TestData.nodeMesh = nodes;
    
    mapMesh = adaptiveMesh.Mesh();
    mapMesh.setMinCellSize([5e-2, 5e-2]);
    mapMesh.initMesh(nodes);
    
    testCase.TestData.mapMesh = mapMesh;
end

function setup(testCase)
    % Nothing to do
end

function teardown(testCase)
    % Nothing to do
end

function teardownOnce(testCase)
    % Nothing to do
end

%% Test Functions

function testGetNeighbor(testCase)
    keyboard
end