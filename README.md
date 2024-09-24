This repository contains the code for the article titled "Formation and Regulation of Calcium Sparks on a Spatial Network of Ryanodine Receptors."

### Usage

The main file is `main.jl`. To run the code, use the following command:

```julia
julia main.jl
```

Parameters can be found in `p.csv`, and you can choose the parameter set using the onff option.

### File Descriptions

- `cluster.jl`: Constructs the spatial network.
  
- `sys.jl`: Handles time evolution.
  
- `plot.jl`: Generates figures for the results located in the `results` directory.
  
- `plotpaper.jl`: Creates figures for the paper.
  
The files in the `Data` folder are solely for the figure regarding the transition rate.
