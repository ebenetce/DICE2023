classdef DiceTrajectoryTest < matlab.unittest.TestCase
    %DiceTrajectoryTest - Test DICE trajectory calculation
    %   DiceTrajectoryTest verifies finite trajectories and expected output
    %   sizes for DICE control variables.
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
        function testReturnsFiniteTrajectory(testCase)
            params = LoadParams(4);
            miu = params.miuup(1:4);
            savings = 0.25*ones(4,1);
            alpha = [params.a0; 0.4*ones(3,1)];

            [utility, consumption, capital] = ...
                diceTrajectory(params, 4, miu, savings, alpha);

            testCase.verifyTrue(isfinite(utility));
            testCase.verifySize(consumption, [4 1]);
            testCase.verifySize(capital, [4 1]);
            testCase.verifyTrue(all(isfinite(consumption)));
            testCase.verifyTrue(all(isfinite(capital)));
        end
    end
end
