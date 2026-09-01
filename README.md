# HR-GGM

HR-GGM is a Gaussian graphical modeling method for estimating conditional
associations among observed variables while incorporating a biologically
specified but unmeasured hidden master regulatory variable and a hidden global
variable representing broad variation. The current implementation is written
in C++17 and is built with CMake and OpenMP.

The source retains legacy identifiers such as `Nlocal`, `Nhidden`, `Nglobal`,
`LHG`, `LL`, `LH`, and `LG`. In the manuscript terminology used here, these
refer to the observed ($O$), hidden master regulatory ($H$), and hidden global
($G$) components.

## Requirements

- 64-bit Linux
- CMake 3.16 or later
- a C++17-compatible compiler, such as GCC
- OpenMP
- MATLAB, used with
  [`Write_Data_whole_bin.m`](https://github.com/songtaow-git/HRGGM_Model/blob/main/HRGGM_Model/Source_code/Write_Data_whole_bin.m)
  to convert a MATLAB `double` matrix into `Data_whole.bin`

## Input files

Before running HR-GGM, prepare a directory named exactly `Data` with the
following structure:

```text
Data/
├── Data_whole.bin
├── Idx_h_list.txt
├── Idx_g_list.txt
├── Lambda_list.txt
├── Alpha_list.txt
└── Theta_full.txt        # simulation GT only; omit for real data
```

The first five files are required by the C++ program. `Theta_full.txt` is an
optional ground-truth precision matrix used only for evaluating simulated
data; it is not read by `Main.cpp`.

Use
[`Write_Data_whole_bin.m`](https://github.com/songtaow-git/HRGGM_Model/blob/main/HRGGM_Model/Source_code/Write_Data_whole_bin.m)
to write the observed-variable matrix as `Data_whole.bin`. The input must be a
MATLAB `double` matrix with:

```text
rows    = observed variables
columns = samples
```

Before calling the function, make sure its `Source_code` directory is the
current MATLAB folder or is included on the MATLAB search path.

For a matrix named `Data`, write the file directly to the prepared input
directory with:

```matlab
Write_Data_whole_bin(Data, fullfile('build', 'Data', 'Data_whole.bin'));
```

Alternatively, `Write_Data_whole_bin(Data)` writes `Data_whole.bin` to the
current MATLAB folder. The function only converts an existing matrix to the
required binary format; data preprocessing and feature selection must be
completed beforehand.

`Idx_h_list.txt` contains C++ zero-based indices of the prior-guided observed
variables used to initialize the hidden master regulatory variable; these
variables remain in the observed data, and the prior set may be incomplete or
noisy.

`Idx_g_list.txt` contains the observed-variable indices used to initialize the
hidden global variable; under the current model formulation, it should include
all observed rows (`0` through `Nlocal-1`) and may overlap with
`Idx_h_list.txt`. Each file stores one index per line. For example, C++ index
`15` selects MATLAB row 16.

`Lambda_list.txt` and `Alpha_list.txt` contain the candidate values of
$\lambda$ and $\alpha$, respectively.

See the [HR-GGM User Manual](docs/HRGGM_User_Manual.md) for the complete binary
format and detailed input-file specifications.

## Configure the model

Edit the adjustable settings near the beginning of `Main.cpp`, including the
number of hidden master regulatory and hidden global variables, the number of
cross-validation folds, the OpenMP thread count, $\gamma$, and the convergence
settings. `Nlocal` is the legacy source identifier for the number of observed
variables; it is read automatically from the number of rows in
`Data_whole.bin` and should not be hard-coded.

## Build and run

From the project root, configure and compile a release build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8
```

After building, place the prepared `Data` directory inside `build`:

```text
build/
├── C_GRN
└── Data/
    ├── Data_whole.bin
    ├── Idx_h_list.txt
    ├── Idx_g_list.txt
    ├── Lambda_list.txt
    ├── Alpha_list.txt
    └── Theta_full.txt    # simulation GT only
```

The build process does not copy the input files automatically. Run HR-GGM from
the project root with:

```bash
./build/C_GRN 2>&1 | tee build/HRGGM_run.log
```

If only `Main.cpp` or another source file is subsequently modified, rebuild
with:

```bash
cmake --build build --parallel 8
```

## Output

The program creates `build/Result/`. The final fitted precision matrix is:

```text
build/Result/Theta_final.txt
```

Its rows and columns are ordered as observed variables, hidden master
regulatory variables, and hidden global variables.

Parameter-selection diagnostics are written to:

```text
build/Result/lambda_ratio/
build/Result/Alpha_metric/
```

For iterative selection outputs, the filename suffix indicates the
selection-loop number (e.g., `1`, `2`, ...).

## Documentation

The complete [HR-GGM User Manual](docs/HRGGM_User_Manual.md) provides:

- detailed descriptions of the source files;
- exact formats of all model inputs;
- the MATLAB-to-C++ binary data layout;
- the C++ zero-based index convention;
- all user-adjustable settings in `Main.cpp`;
- the parameter-selection and final-fitting workflow;
- output-matrix block interpretation;
- computational guidance and troubleshooting.

## Citation

HR-GGM paper: the complete citation will be added upon publication.

Eigen: HR-GGM uses Eigen for dense and sparse matrix computations in its C++
implementation. Please cite:
Gaël Guennebaud, Benoît Jacob, and others. (2010). Eigen.
https://libeigen.gitlab.io.
