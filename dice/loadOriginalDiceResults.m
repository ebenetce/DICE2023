function results = loadOriginalDiceResults()
    %loadOriginalDiceResults - Load published DICE-2023 GAMS results
    %   RESULTS = loadOriginalDiceResults() loads the local benchmark or
    %   recreates it from the published DICE-2023 workbook.
    %
    %   Copyright 2024-2026 The MathWorks, Inc.

    projectFolder = fileparts(fileparts(mfilename("fullpath")));
    resultsPath = fullfile(projectFolder, "GAMSresults.mat");

    if ~isfile(resultsPath)
        workbookURL = "https://www.dicemodel.org/_files/ugd/66d8d1_bff32c2165564c7c94a2f54bbe03cc6f.xlsx?dn=DICE2023-Excel-b-4-3-10-v18.3.xlsx";
        workbookPath = string(tempname) + ".xlsx";
        cleanup = onCleanup(@() delete(workbookPath));

        try
            websave(workbookPath, workbookURL);
        catch exception
            error("DICE2023:BenchmarkDownloadFailed", ...
                "Could not download the published DICE-2023 workbook. %s", ...
                exception.message);
        end

        scenarioRows = string(readcell(workbookPath, Sheet="GAMS", Range="A1:A800"));
        sourceScenarios = ["SCENARIO: OPTIMAL" "SCENARIO: T < 2%" ...
            "SCENARIO: T < 1.5" "SCENARIO: ALTDAM" ...
            "SCENARIO: PARIS-UPDATE" "SCENARIO: R = 1%"];
        resultNames = ["tb" "tb_20" "tb_15" "altdam" "paris" "disc1"];
        results = struct;

        for scenarioIndex = 1:numel(sourceScenarios)
            scenarioRow = find(scenarioRows == sourceScenarios(scenarioIndex), 1);
            if isempty(scenarioRow)
                error("DICE2023:MissingBenchmarkScenario", ...
                    "The published workbook does not contain the %s scenario.", ...
                    sourceScenarios(scenarioIndex));
            end

            headerRow = scenarioRow + 6;
            raw = readcell(workbookPath, Sheet="GAMS", ...
                Range=sprintf("A%d:CD%d", headerRow, headerRow + 59));
            rowNames = string(raw(:,1));
            rowNames = replace(rowNames, ...
                ["Atmospheric temperaturer (deg c above preind)" "Abatement/0utput"], ...
                ["Atmospheric temperature (deg c above preind)" "Abatement/Output"]);
            years = cell2mat(raw(2,2:end));
            values = nan(size(raw,1), size(raw,2)-1);

            for row = 1:size(values,1)
                for column = 1:size(values,2)
                    value = raw{row,column+1};
                    if isnumeric(value)
                        values(row,column) = value;
                    elseif ischar(value) || isstring(value)
                        values(row,column) = str2double(value);
                    end
                end
            end

            if ~isequal(rowNames(1:2), ["Period"; "Year"]) || ...
                    numel(years) ~= 81 || any(diff(years) ~= 5)
                error("DICE2023:UnexpectedBenchmarkFormat", ...
                    "The published workbook has an unexpected GAMS scenario layout.");
            end

            results.(resultNames(scenarioIndex)) = array2table(values, ...
                VariableNames="Period" + string(years), RowNames=rowNames);
        end

        save(resultsPath, "-struct", "results");
    end

    results = load(resultsPath);
end
