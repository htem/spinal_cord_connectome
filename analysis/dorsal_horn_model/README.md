Connectome-constrained model of the dorsal horn microcircuit.
Copyright (MIT) 2026 Ramin Khajeh


SUMMARY:
This code implements the model in: W. Xiang, R. Khajeh, K. Delgado, J. Liu, A. Holtrup, Q. Zhang, S. Cameron, T. Nguyen, J.A. Bae, A. Halageri, N. Kemnitz, T. Macrina, B.M. Llamas, A. Schildkamp, M.F. Mazri, J. Ni, R. Roberts, L.S. Capdevila, R. Wang, A. Dsouza, A. Norton, S. Wang, D. Zhang, P.L. Tang, W-C.A. Lee, D.D. Ginty. A synaptic wiring diagram underlying modality-specific ascending pathways of the superficial spinal dorsal horn.


REQUIREMENTS:
This program requires Julia (julialang.org) and tested on Julia 1.10.2. Install all dependencies with `julia install.jl`


USAGE:
To produce the main figure in the manuscript, inside the Julia REPL, do
```julia
include("dh_model.jl")
fig = v.plot_main_figure()
```

To produce the optimal-input panel, do
```julia
include("dh_model_optimInput.jl")
results = v.compute_optimal_h_multiple()
fig = v.plot_optimal(results)
```

To produce the panels in the corresponding supplementary figure, do
```julia
include("dh_model.jl")
results = v.compute_spearman_heatmap()
fig = v.plot_spearman_heatmap(results)
```

And
```julia
include("dh_model.jl")
fig = v.plot_selectivity_combined()
```

To plot the circuit schematic, do
```julia
include("dh_circuit_schematic.jl")
fig = v.plot_circuit()
```

Variable names in the code are as consistent as possible with those that appear in the manuscript;  see Methods for more information.

CONTACT:
Ramin Khajeh
ramin_khajeh at hms dot harvard dot edu
