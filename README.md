# autodiff-drc

This repo contains the code to calculate degree of rate control for surface chemical reactions.
The usage of the code is pretty simple:
```console
julia run_drc.jl input_path [--output_path] [--delsa] [--delsa_sample]
```
where 
- `input_path` points to a input yaml file that contains reactions, reaction conditions, time intervals for calculation, initial conditions and the kinetic parameters. Examples (Propylene partial oxidation) for DRC calculation w/wo DELSA are provided. (see PO-input.yaml and PO-input-delsa.yaml).
- `output_path` is the folder to store the results. If it does not exit, will be created.
- `delsa` is the flag of bool type. If set, delsa will be used to estimate the distribution of the degree of rate control. The input file should contain the lower and upper bound of the kinetic parameters.
- `delsa_sample`: int flag specifies the number of samples generated over the uncertain range of the kinetic parameters.
