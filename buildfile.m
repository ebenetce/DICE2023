function plan = buildfile
import matlab.buildtool.tasks.CodeIssuesTask

% Create a plan from task functions
plan = buildplan(localfunctions);

% Add the "check" task to identify code issues
plan("check") = CodeIssuesTask(WarningThreshold=0);

% Make the "archive" task the default task in the plan
plan.DefaultTasks = "archive";

% Make the "archive" task dependent on the "check" and "test" tasks
plan("archive").Dependencies = "check";
end

function archiveTask(~, version)

opts = matlab.addons.toolbox.ToolboxOptions('dice', "92916d00-6626-4c09-b4d7-1cc087b13a6f"); %pick any identifier you want, but keep it constant.
opts.AuthorCompany = 'MathWorks';
opts.AuthorEmail = 'ebenetce@mathworks.com';
opts.AuthorName = 'Eduard Benet Cerda';
opts.Description = 'MATLAB Implementation of the DICE 2023 model from William Nordhaus.';
opts.OutputFile = 'DICE2023.mltbx';
opts.Summary = 'DICE Climate Model 2023';
opts.ToolboxName = 'DICE2023';
opts.ToolboxVersion = version;

matlab.addons.toolbox.packageToolbox(opts);
                
end