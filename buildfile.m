function plan = buildfile
%buildfile - Configure the DICE-2023 build tasks
%   PLAN = buildfile creates code-check and test tasks for the project.
%
%   Copyright 2024-2026 The MathWorks, Inc.

import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask

plan = buildplan(localfunctions);

plan("check") = CodeIssuesTask("dice", WarningThreshold=0);

plan("test") = TestTask("tests", SourceFiles="dice");

plan.DefaultTasks = ["check", "test"];
end
