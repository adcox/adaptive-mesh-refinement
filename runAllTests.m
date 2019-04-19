clear; clc; close all;

testScripts = dir('tests/test*.m');
files = cell(length(testScripts), 1);
for f = 1:length(testScripts)
    files{f} = [testScripts(f).folder, filesep, testScripts(f).name];
end
results = runtests(files);
disp(results);