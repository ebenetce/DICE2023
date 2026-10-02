# DICE 2023 — MATLAB

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ebenetce/DICE2023&project=DICE2023.prj&file=Main.m)

MATLAB&reg; implementation of the DICE-2023 integrated assessment model described by [Barrage and Nordhaus (2024)](https://doi.org/10.1073/pnas.2312030121).

The published DICE model, GAMS implementation, documentation, and supporting files are available from [DICE Model](https://www.dicemodel.org/).

## Run the model

Open the project, then create and solve the default 81-period DICE problem:

Requirements: MATLAB&reg; and Optimization Toolbox&trade;. Parallel Computing Toolbox&trade; is optional.

```matlab
openProject("DICE2023.prj")
```

Load the parameters:

```matlab
np = 81;
params = LoadParams(np);
```

Create the optimization problem:

```matlab
prob = diceFunc(np, params);
```

Define an initial point:

```matlab
x0.MIU = params.miuup(1:np) - eps;
x0.S = 0.25*ones(np, 1);
x0.alpha = [params.a0; 0.4*ones(np-1, 1)];
```

Solve and post-process the results:

```matlab
opts = optimoptions("fmincon", Display="iter", UseParallel=false);
[sol, fval] = solve(prob, x0, Options=opts);
allVars = recoverAllVars(sol, params);
```

## Published benchmark results

`loadOriginalResults` loads the local `GAMSresults.mat` benchmark file. If the file is absent, it downloads the [published DICE-2023 Excel workbook](https://www.dicemodel.org/_files/ugd/66d8d1_bff32c2165564c7c94a2f54bbe03cc6f.xlsx?dn=DICE2023-Excel-b-4-3-10-v18.3.xlsx), recreates `GAMSresults.mat`, and returns the six GAMS scenarios used by the examples.

For questions, contact ebenetce@mathworks.com.
