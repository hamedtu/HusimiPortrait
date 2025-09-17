## Pre-Floquet States – Data and Reproduction Code

![Irises Illustration](irises.jpeg)

This folder contains the data and scripts used to generate figures for the preprint:

- "Pre-Floquet states facilitating coherent subharmonic response of periodically driven many-body systems" by Steffen Seligmann, Hamed Koochakikelardeh, and Martin Holthaus (DOI: https://doi.org/10.1103/vgpm-s6tv).
  
The codebase is primarily MATLAB for figure generation and C++ (MPI + Eigen) for generating Floquet eigenvectors used in some figures.


### Citation
If you use this dataset or code, please cite the manuscript above. For questions, contact **Dr. Hamed Koochaki Kelardeh** (`hamed.koochakikelardeh@uni-oldenburg.de`).

### Contents
- `Main.m`: Entry point to reproduce figures. Select a figure via `Fig_choice` and run.
- `Coherence.m`: Computes coherence measure `eta` for Floquet states.
- `Husimiplot.m`: Computes and plots Husimi distributions over phase space and overlays a phase-space portrait.
- `pSpacePortrait.m`: Generates the Poincaré surface of section used by `Husimiplot`.
- `main.cpp`: Time-slicing solver (MPI-parallel) to compute Floquet eigenvectors; writes the `TS_eigenvectors_*.dat` files used in Fig. 3 and Fig. 6.

### Dataset
- Text files `Fig*.txt` hold precomputed data arrays used directly by `Main.m` to render figures (return probability, Poincaré maps, occupation probabilities, etc.).
- `TS_eigenvectors_*.dat` contain real and imaginary parts of Floquet eigenvectors computed by `main.cpp` for parameters `(alpha, mu, omega, N, L)` encoded in the filename.

### Task and Flow
1. Reproduce figures from the paper by running MATLAB:
   - Open `Main.m`.
   - Set `Fig_choice` to one of `{1, 2, 3, 5, 6, 7, 8, 9}` as indicated in the script.
   - Run the script to generate the selected figure.
2. For figures 3 and 6: if the `TS_eigenvectors_*.dat` files matching the parameters are missing, generate them with the C++ program `main.cpp` (see Build & Run below). Place the resulting `.dat` files in the same folder and re-run `Main.m`.

### MATLAB Details
- `Coherence(u_n, N)`: Computes the coherence measure `eta` for each Floquet state (columns of `u_n`) for particle number `N`, returning sorted values and indices.
- `Husimiplot(V, V_ind, param, N, gridSize)`: Builds Husimi distributions on a `(p, φ)` grid, normalizes per-state maxima for visibility when summing multiple states, and overlays the phase-space portrait via `pSpacePortrait(param...)`.
- `pSpacePortrait(param1, param2, param3)`: Integrates the mean-field equations over many periods and scatters section points to show the Poincaré map.

Required MATLAB toolboxes and notes:
- MATLAB R2022b or newer is recommended.
- Uses `ode45` (base MATLAB ODE suite) and `normpdf` (Statistics and Machine Learning Toolbox) in `Husimiplot`.
- LaTeX axis labels are used; ensure LaTeX interpreter support is available in your MATLAB.

### C++ Generator (optional, for eigenvectors)
`main.cpp` computes Floquet eigenvectors using a time-slicing method and saves them to `TS_eigenvectors_<alpha>_<mu>_<omega>_<N>_<L>_{real,imag}.dat`.

Dependencies:
- C++17 compiler
- MPI implementation (e.g., MPICH or OpenMPI)
- Eigen 3.x (header-only)

Example build (adjust include paths as needed):
```bash
mpic++ -O3 -std=c++17 -I /path/to/eigen main.cpp -o timeslicing_solver
```

Example run (parameters: `alpha mu omega N L`):
```bash
mpirun -np 8 ./timeslicing_solver 0.92 0.40 1.90 10000 64
```
This will produce files matching those already included in this folder for Fig. 3 and Fig. 6.

Performance tip:
- Large `N` and `L` significantly increase cost. Use MPI with enough ranks and ensure Eigen is compiled with optimizations (`-O3`).

### Figure Mapping
- `Fig_choice == 1`: Return probability (`Fig1.txt`).
- `Fig_choice == 2`: Poincaré map from precomputed data (`Fig2.txt`).
- `Fig_choice == 3`: Color-coded Husimi projections from `TS_eigenvectors_*.dat` at `(alpha=0.92, mu=0.40, ω=1.90, N=10000, L=64)`.
- `Fig_choice == 5`: Two color plots from `Fig5a.txt` and `Fig5b.txt`.
- `Fig_choice == 6`: Husimi projection for a selected state, using `TS_eigenvectors_*.dat`.
- `Fig_choice == 7`: Occupation probabilities heatmap from `Fig7.txt`.
- `Fig_choice == 8`: Zoomed Poincaré map from `Fig8.txt`.
- `Fig_choice == 9`: Return probability (two system sizes) from `Fig9a.txt` and `Fig9b.txt`.

### How to Reproduce a Figure (quick start)
1. Start MATLAB.
2. Open `Main.m`, set `Fig_choice` to the desired figure number.
3. Ensure required input files for that figure are present in this folder.
4. Run.

### MATLAB quick-start snippets per figure
Copy one of the snippets below into the MATLAB console or adapt in `Main.m`.

```matlab
% Fig. 1 – Return probability
clc; clear; close all
Fig_choice = 1; % uses Fig1.txt
run('Main.m');
```

```matlab
% Fig. 2 – Poincaré map (precomputed)
clc; clear; close all
Fig_choice = 2; % uses Fig2.txt
run('Main.m');
```

```matlab
% Fig. 3 – Husimi projections from Floquet eigenvectors
clc; clear; close all
Fig_choice = 3; % uses TS_eigenvectors_0.92_0.40_1.90_10000_64_{real,imag}.dat
run('Main.m');
```

```matlab
% Fig. 5 – Two color plots
clc; clear; close all
Fig_choice = 5; % uses Fig5a.txt and Fig5b.txt
run('Main.m');
```

```matlab
% Fig. 6 – Husimi projection (single state)
clc; clear; close all
Fig_choice = 6; % uses TS_eigenvectors_0.92_0.40_1.90_10000_64_{real,imag}.dat
run('Main.m');
```

```matlab
% Fig. 7 – Occupation probabilities heatmap
clc; clear; close all
Fig_choice = 7; % uses Fig7.txt
run('Main.m');
```

```matlab
% Fig. 8 – Zoomed Poincaré map
clc; clear; close all
Fig_choice = 8; % uses Fig8.txt
run('Main.m');
```

```matlab
% Fig. 9 – Return probability (two system sizes)
clc; clear; close all
Fig_choice = 9; % uses Fig9a.txt and Fig9b.txt
run('Main.m');
```




