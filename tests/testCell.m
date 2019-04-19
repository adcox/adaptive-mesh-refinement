%% Main function to generate tests
function results = testCell
    results = functiontests(localfunctions);
end

%% Fixtures
function setupOnce(testCase)
    addpath ../
    
    testCase.TestData.nodes = [...
        adaptiveMesh.Node([-1, -1]), adaptiveMesh.Node([1, -1]),...
        adaptiveMesh.Node([-1, 1]), adaptiveMesh.Node([1, 1])];
end

function setup(testCase)
    % Nothing to do
end

function teardown(testCase)
    % Nothing to do
end

function teardownOnce(testCase)
    % Nothing to do
    rmpath ../
end

%% Test Functions

function testSetNodes(testCase)
    
    nodes = reshape(testCase.TestData.nodes, 2, 2);
    nodesWrongOrder = nodes';
    cell = adaptiveMesh.Cell();
    
    % No input arguments throws an error
    assertError(testCase, @()cell.setNodes(), 'MATLAB:minrhs');
    
    % Passing in doubles is the wrong type
    assertError(testCase, @()cell.setNodes([1, 2; 3, 4]), 'adaptiveMesh:Cell:type');
    
    % Passing in too few nodes
    assertError(testCase, @()cell.setNodes(nodes(1:2)), 'adaptiveMesh:Cell:size');
    
    % Passing in the wrong order
    assertError(testCase, @()cell.setNodes(nodesWrongOrder), 'adaptiveMesh:Cell');
    
    try
        cell.setNodes(nodes);
    catch err
        fprintf('Failure: %s', err.message);
        assertFail(testCase);
    end
end