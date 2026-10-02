classdef LoadParamsTest < matlab.unittest.TestCase

    methods (TestClassSetup)
        function addSourcePath(testCase)
            projectFolder = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(fullfile(projectFolder, "dice")));
        end
    end

    methods (Test)
        function testDefaultFiveYearCalibration(testCase)
            params = LoadParams();

            testCase.verifyEqual(params.tstep, 5);
            testCase.verifySize(params.L, [81 1]);
            testCase.verifyEqual(params.miuup(1:2), [0.05; 0.10], AbsTol=eps);
        end

        function testOneYearStepProducesAnnualParameters(testCase)
            params = LoadParams(21, tstep=1);

            testCase.verifyEqual(params.tstep, 1);
            testCase.verifySize(params.L, [21 1]);
            testCase.verifyTrue(all(isfinite(params.miuup)));
        end
    end
end
