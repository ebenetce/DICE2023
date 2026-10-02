classdef RecoverAllVarsTest < matlab.unittest.TestCase
    %RecoverAllVarsTest - Test DICE variable recovery
    %   RecoverAllVarsTest verifies recovery with fixed control variables.
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
        function testFixedControlContextRestoresMissingControl(testCase)
            params = LoadParams(3);
            controls = struct( ...
                "MIU", params.miuup(1:3), ...
                "S", 0.25*ones(3,1), ...
                "alpha", [params.a0; 0.4; 0.4]);
            expected = recoverAllVars(controls, params);
            partialSolution = rmfield(controls, "MIU");
            fixedControls = struct("MIU", controls.MIU);

            actual = recoverAllVars(partialSolution, params, ...
                FixedControls=fixedControls);

            testCase.verifyEqual(actual.Consumption, expected.Consumption, ...
                AbsTol=1e-12);
        end
    end
end
