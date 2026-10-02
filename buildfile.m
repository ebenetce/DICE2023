function plan = buildfile
import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask

plan = buildplan(localfunctions);

plan("check") = CodeIssuesTask("dice", WarningThreshold=0);

plan("test") = TestTask("tests", SourceFiles="dice");

plan.DefaultTasks = ["check", "test"];
end
