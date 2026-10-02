classdef LoadOriginalDiceResultsTest < matlab.unittest.TestCase
    %LoadOriginalDiceResultsTest - Test published benchmark result loading
    %   LoadOriginalDiceResultsTest verifies the published GAMS scenarios
    %   are available with their expected dimensions.
    %
    %   Copyright 2024-2026 The MathWorks, Inc.

    methods (TestClassSetup)
        function addSourcePath(testCase)
            projectFolder = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(fullfile(projectFolder, "dice")));
        end
    end

    methods (Test)
        function testReturnsExpectedBenchmarkTables(testCase)
            results = loadOriginalDiceResults();

            testCase.verifyTrue(isfield(results, "tb"));
            testCase.verifyTrue(isfield(results, "tb_20"));
            testCase.verifyTrue(isfield(results, "tb_15"));
            testCase.verifyTrue(isfield(results, "altdam"));
            testCase.verifyTrue(isfield(results, "paris"));
            testCase.verifyTrue(isfield(results, "disc1"));
            testCase.verifySize(results.tb, [60 81]);
        end
    end
end
