This repository contains the code used to conduct the analysis in Setterholm et
al. (2024). Both the parametric models and image reconstruction codes are
written in the [Julia programming language](https://julialang.org)[^1], however
each of these two analyses uses a different version of the `OIFITS.jl` library,
so they are seperated into their own directories.

To run the code, simply `cd` into either the `ImageReconstruction` or
`ParametricModelFitting` directory, start a Julia REPL session, press the
<kbd>]</kbd> key to enter the package manager mode, and enter in
```activate .```
to install/load the relevant library dependencies. Once this is complete, press
the <kbd>Backspace</kbd> key to exit package manager mode and then use the
`include` function to load the .jl files.

For the parametric modeling code, modify the paths/analysis as needed in the
`example.jl` file for your system before including it. Once you have issued the
`include("example.jl")` command in the Julia REPL, the model fitting will begin
automatically.

For image reconstruction, you will need to use the `mkimage` function from the
Julia REPL directly passing in the path to the raw (merged) OIFITS file that you
would like to image.

[^1]: Specifically, Julia v1.9.4 was used in the article but there should be no
issue running this code on newer 1.x versions of Julia