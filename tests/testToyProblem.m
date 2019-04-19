%% Main function to generate tests
function results = testToyProblem
    results = functiontests(localfunctions);
end

%% Fixtures
function setupOnce(testCase)
    addpath ../
    
    testCase.TestData.mapMesh = adaptiveMesh.Mesh();
    testCase.TestData.mapMesh.setMinCellSize([2e-2, 2e-2]);
    testCase.TestData.mapMesh.initMesh([-1, 1, -1, 1], ToyNode());
end

function setup(testCase)
    % Nothing to do
end

function teardown(testCase)
    % Nothing to do
end

function teardownOnce(testCase)
    % Nothing to do
    delete(testCase.TestData.mapMesh);
    rmpath ../
end

%% Test Functions
function testMeshRefinement(testCase)
    mesh = testCase.TestData.mapMesh;
    mesh.refine();
    
    cellKeys = fieldnames(mesh.cellMap);
    for k = 1:length(cellKeys)
        cell = mesh.cellMap.(cellKeys{k});
        
        if(~cell.isSubdivided)
            neighbors = cell.getNeighbors(mesh);
            
            for n = 1:length(neighbors)
                assertTrue(testCase, abs(neighbors(n).level - cell.level) <= 1);
            end
        end
    end
end