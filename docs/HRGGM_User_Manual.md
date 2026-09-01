# HR-GGM User Manual

## 1. Overview

HR-GGM is a Gaussian graphical modeling method for estimating conditional
associations among observed variables while explicitly incorporating a
biologically specified but unmeasured hidden master regulatory variable and a
hidden global variable that represents broad variation.

The variables are ordered as:

1. observed variables, denoted by $O$;
2. hidden master regulatory variables, denoted by $H$;
3. hidden global variables, denoted by $G$.

Let $N_O$, $N_H$, and $N_G$ denote their respective numbers. The analyses in
the current manuscript use $N_H=N_G=1$, while the implementation retains
`Nhidden` and `Nglobal` as configurable counts.

The estimated precision matrix therefore has the block structure

$$
\Theta =
\begin{bmatrix}
\Theta_O & \Theta_{OH} & \Theta_{OG} \\
\Theta_{HO} & \Theta_H & \Theta_{HG} \\
\Theta_{GO} & \Theta_{GH} & \Theta_G
\end{bmatrix}.
$$

The source code retains legacy implementation identifiers such as `Nlocal`,
`Nhidden`, `Nglobal`, `LHG`, `LL`, `LH`, and `LG`. In this manual, `local` or
`L` in a literal code identifier maps to an observed variable or the observed
block $O$; `hidden` or `H` maps to the hidden master regulatory variable; and
`global` or `G` maps to the hidden global variable. These identifiers must not
be renamed unless the source code is updated consistently.

The current implementation is a C++17 project. Users place the required input
files in the **Data** directory, edit the adjustable settings at the beginning
of **Main.cpp**, build the project with CMake, and run the resulting executable.

## 2. Software requirements

The following software is required:

* a 64-bit Linux system;
* CMake 3.16 or later;
* a C++17-compatible compiler such as GCC;
* OpenMP support;
* MATLAB, used to convert an arbitrary MATLAB `double` matrix into
`Data_whole.bin` with the provided binary-writer code.

On Ubuntu or Debian, the basic compilation tools can be installed with:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake
```

On a managed computing server, use the compiler and CMake modules provided by
the system administrator when system-wide installation is unavailable.

## 3. Project files

The source directory contains the HR-GGM implementation, the bundled Eigen
headers, and the CMake build definition. In general, each `.h` file provides
shared type definitions and function declarations, while the corresponding
`.cpp` file provides the implementations. The principal files are:

|File|Purpose|
|-|-|
|`Main.cpp`|Program entry point and user-adjustable settings. Locates `Data/`, reads all inputs, constructs cross-validation folds, selects lambda and alpha, fits the final HR-GGM model, and writes results.|
|`Class_define.h`|Defines the shared matrix aliases and declares `Input_format`, update-index helpers, KKT-residual containers, fold/z-score structures, and lambda/alpha metric structures.|
|`Class_define.cpp`|Implements constructors, accessors, diagnostic printing, marked-update handling, KKT reporting, and the lambda- and alpha-selection score calculations declared in `Class_define.h`.|
|`LHG_update.h`|Declares data standardization, covariance calculation, prior-guided PCA construction, updates of the hidden master regulatory and hidden global variables, fold-specific data readers, extended covariance construction, alternating model fitting, and validation metrics. `LHG` is a legacy source identifier.|
|`LHG_update.cpp`|Implements initialization and updating of the observed, hidden-regulatory, and global components, block-wise covariance and PCA operations, K-fold preprocessing, train/test construction, alternating precision-matrix/data estimation, and likelihood, entropy, and stability metrics.|
|`Matrix_compute.h`|Declares numerical routines for log determinants, traces, likelihoods, block penalties, sparse-matrix inversion, active update sets, coordinate updates, gradients, KKT convergence checks, and precision-matrix expansion.|
|`Matrix_compute.cpp`|Implements the numerical matrix operations, observed--observed, observed--hidden-regulatory, and observed--global penalty calculations, coordinate-wise update formulas, sparse Cholesky-based inversion, gradient calculations, and KKT residual checks declared in `Matrix_compute.h`.|
|`Matrix_io.h`|Declares the binary and text input/output interfaces, including data loading, parameter/index reading, matrix/vector writing, and timestamp generation.|
|`Matrix_io.cpp`|Implements loading of `Data_whole.bin`, reading parameter and C++ index files, and writing dense matrices, sparse matrices, metrics, vectors, and update pairs.|
|`Newton_method_offdiag.h`|Declares the two precision-matrix optimization interfaces: the observed-only solver and the full extended HR-GGM solver.|
|`Newton_method_offdiag_l.cpp`|Implements the observed-only sparse precision-matrix optimizer, including observed-block coordinate updates, positive-definiteness-preserving line search, and KKT-based convergence.|
|`Newton_method_offdiag_lhg.cpp`|Implements the full HR-GGM precision-matrix optimizer with separate observed--observed, observed--hidden-regulatory, and observed--global updates, the corresponding block penalties, Cholesky-based line search, and KKT-based convergence.|
|`Parameter_select_parallel.h`|Declares input-copy helpers and the fold-wise covariance, lambda-selection, and alpha-selection routines used during parallel cross-validation.|
|`Parameter_select_parallel.cpp`|Implements OpenMP-parallel fold-wise computation and evaluation of lambda and alpha candidates, including model warm starts and collection of validation metrics.|
|`CMakeLists.txt`|Defines the C++17 `C_GRN` executable, lists the implementation files, adds the bundled Eigen include path, links OpenMP, and enables optimized compiler settings.|
|`Eigen/`|Bundled header-only Eigen linear-algebra dependency used by the HR-GGM source code.|
|`README.md`|Repository overview, quick-start commands, and link to the complete user manual.|
|`docs/HRGGM_User_Manual.md`|Complete model build, configuration, input, execution, output, and troubleshooting documentation.|

Users normally edit only the adjustable parameter block at the beginning of
`Main.cpp`. The remaining source files implement the method.

## 4. Model input files

### 4.1 Input location and directory structure

HR-GGM reads all model inputs from a directory named exactly `Data`. The
executable starts from its own directory and searches upward for the nearest
`Data/` directory. For a standard Linux build, the following layout is
recommended:

```text
HRGGM/
├── build/
│   ├── C_GRN
│   └── Data/
│       ├── Data_whole.bin
│       ├── Idx_h_list.txt
│       ├── Idx_g_list.txt
│       ├── Lambda_list.txt
│       ├── Alpha_list.txt
│       └── Theta_full.txt    # simulation data only
├── Main.cpp
├── CMakeLists.txt
└── other source files
```

The first five files are required and are read directly by the C++ program.
`Theta_full.txt` is optional and appears only in simulation examples that have
a known ground-truth precision matrix. File and directory names are
case-sensitive on Linux.

### 4.2 `Data_whole.bin`

The provided MATLAB binary-writer code converts an arbitrary MATLAB `double`
matrix into `Data_whole.bin`. This conversion code only writes the matrix in
the format expected by C++; dataset-specific simulation or preprocessing code
can be supplied separately for each analysis.

Let `Data` denote the matrix passed to the MATLAB writer. Its orientation must
be:

```text
rows    = observed variables
columns = samples
```

Thus, `Data(i,j)` is the value of observed variable `i` in sample `j`. For
example, 200 observed variables measured in 10,000 samples must be stored as a
`200 x 10000` MATLAB `double` matrix. The C++ program obtains `Nlocal` from the
number of rows and the sample count from the number of columns; users do not
set `Nlocal` manually in `Main.cpp`. Here, `Nlocal` is the legacy C++ name for
the number of observed variables.

The binary file stores, in order:

1. the number of rows as an unsigned 64-bit integer;
2. the number of columns as an unsigned 64-bit integer;
3. all values as 64-bit floating-point numbers in row-major order.

The row order must remain unchanged because `Idx_h_list.txt` and
`Idx_g_list.txt` refer directly to these rows. Complete normalization,
transformation, quality control, and feature selection before creating the
binary file. HR-GGM then performs row-wise z-score standardization internally.
The final matrix must contain only finite values and no constant rows.

### 4.3 C++ index convention

Both index files contain one **C++ zero-based row index** per line. The first
row of `Data` is index `0` in C++ but row `1` in MATLAB:

$$
\text{MATLAB row number} = \text{C++ index} + 1.
$$

|Stored C++ index|Selected feature|MATLAB expression|
|-:|-|-|
|`0`|first observed variable|`Data(1,:)`|
|`15`|16th observed variable|`Data(16,:)`|
|`28`|29th observed variable|`Data(29,:)`|
|`183`|184th observed variable|`Data(184,:)`|
|`199`|200th observed variable|`Data(200,:)`|

The stored value and its line number are different. For example, line 16 of
`Idx_g_list.txt` stores `15`, while line 168 stores `183`. Here, `16` and `168`
are positions within the index file; the stored C++ indices `15` and `183`
select MATLAB rows 16 and 184, respectively. Do not add one to values written
to the index files.

### 4.4 `Idx_h_list.txt`

`Idx_h_list.txt` contains the C++ indices of the prior-guided observed
variables used to initialize the hidden master regulatory variable. These
variables remain part of the observed data; the index file does not remove or
reclassify them. The example contains 16 indices, shown compactly below; the
actual file stores one index per line:

```text
28, 29, 40, 42, 62, 64, 83, 84, 86, 89, 111, 112, 170, 172, 180, 181
```

These values correspond to MATLAB rows 29, 30, 41, 43, 63, 65, 84, 85, 87,
90, 112, 113, 171, 173, 181, and 182. When `Nhidden = 1`, each sample's score
on the first principal component of these selected rows initializes the hidden
master regulatory variable. The number `16` is the number of entries, not a
stored C++ index.

### 4.5 `Idx_g_list.txt`

`Idx_g_list.txt` contains the C++ indices of the observed variables used to
initialize the hidden global variable, with one index per line. For the 200-variable example, this file contains can indice `0` through `199`. When`Nglobal = 1`, each sample's score on the first principal component of these rows initializes the hidden global variable.

### 4.6 `Lambda_list.txt`

`Lambda_list.txt` stores the candidate lambda values in one comma-separated
row. The example contains:

```text
0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06
```

Lambda controls sparsity in the observed-variable block $\Theta_O$ and
contributes to the observed--hidden-regulatory penalty.

### 4.7 `Alpha_list.txt`

`Alpha_list.txt` stores the candidate alpha values in one comma-separated row.
The example contains:

```text
0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2
```

Alpha controls the squared $\ell_2$ penalty on the observed--global block and
the squared $\ell_2$ component of the observed--hidden-regulatory penalty.

### 4.8 `Theta_full.txt` for simulation data only

`Theta_full.txt` is the **ground-truth (GT) precision matrix** and is included
only with simulated data for evaluation. It is not a required HR-GGM input,
and `Main.cpp` does not read it during parameter selection or final fitting.
Real-data analyses normally do not contain this file.

For a simulation with 200 observed variables, one hidden master regulatory
variable, and one hidden global variable, the GT matrix has dimension:

```text
(200 + 1 + 1) x (200 + 1 + 1) = 202 x 202
```

## 5. Configure `Main.cpp`

Open `Main.cpp` and edit the adjustable parameter block near the beginning of
`main()`:

```cpp
const int Max_iter = 10000;
const double Converge_thre = 1e-3;
const int N_thread = 8;
const int Blockcols = 4096;
const int Fold_k = 3;
const int Nhidden = 1;
const int Nglobal = 1;
const double Gamma = 0.9;
const double Alpha_retention_tau = 0.20;
```

|Setting|Meaning|Example value|
|-|-|-:|
|`Max_iter`|Maximum iterations in each precision-matrix optimization call|`10000`|
|`Converge_thre`|KKT convergence tolerance|`1e-3`|
|`N_thread`|OpenMP workers used during parameter selection|`8`|
|`Blockcols`|Number of sample columns processed per data block|`4096`|
|`Fold_k`|Number of cross-validation folds|`3`|
|`Nhidden`|Number of modeled hidden master regulatory variables|`1`|
|`Nglobal`|Number of modeled hidden global variables|`1`|
|`Gamma`|Mixing weight for observed--hidden-regulatory regularization|`0.9`|
|`Alpha_retention_tau`|Retention threshold used in alpha selection|`0.20`|

For the supplied 200-observed-variable example, use:

```cpp
const int Nhidden = 1;
const int Nglobal = 1;
```

`Nlocal` is read automatically from `Data_whole.bin`. It is the legacy source
identifier for the number of observed variables and must not be hard-coded.

After modifying `Main.cpp`, save the file before building the project.

## 6. Build HR-GGM on Linux

From the top-level source directory, configure a release build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Compile the project:

```bash
cmake --build build --parallel 8
```

After the build is complete, place the prepared `Data` folder inside the
`build` directory. Its contents and file formats must follow Section 4:

```text
build/
├── C_GRN
└── Data/
    ├── Data_whole.bin
    ├── Idx_h_list.txt
    ├── Idx_g_list.txt
    ├── Lambda_list.txt
    ├── Alpha_list.txt
    └── Theta_full.txt        # simulation GT only; omit for real data
```

**The program must be able to locate `build/Data/` when it starts. Building the
executable does not copy the input files automatically.**

The Linux executable is normally created at:

```text
build/C_GRN
```

If only `Main.cpp` or another source file is modified after the initial CMake
configuration, rebuild with:

```bash
cmake --build build --parallel 8
```

It is not necessary to rerun the CMake configuration command after every
parameter edit unless the build configuration or source-file list changes.

## 7. Run HR-GGM

From the project root, run:

```bash
./build/C_GRN
```

To retain the complete terminal output in a log file, use:

```bash
./build/C_GRN 2>&1 | tee build/HRGGM_run.log
```

During execution, HR-GGM reports:

* fold-specific covariance construction;
* the initial $\lambda$ search;
* alternating $\alpha$ and $\lambda$ searches;
* elapsed time for parameter/fold combinations;
* the selected $\lambda$, $\alpha$, and $\gamma$;
* completion of the final full-data fit.

The program creates `Result` beside the located `Data` directory.

## 8. Computational workflow

### 8.1 Cross-validation folds

The samples are randomly permuted using the default fold-plan seed in the C++
implementation. Equal-sized folds are constructed. If the number of samples is
not divisible by `Fold_k`, the remainder is not used in cross-validation.

Training and test data are standardized within their respective fold-specific
data subsets.

### 8.2 Initialization of the unobserved variables

For every fold:

1. the prior-guided observed-variable rows listed in `Idx_h_list.txt` are used
   for PCA;
2. the leading component scores initialize the hidden master regulatory
   variable;
3. the observed-variable rows listed in `Idx_g_list.txt` are used in a
   separate PCA;
4. the leading component scores initialize the hidden global variable;
5. the initialized variables are standardized and combined with the observed
   data to form the extended data matrix.

With `Nhidden = 1` and `Nglobal = 1`, the top component from each index subset
is retained. In the source orientation (variables by samples), the centered
extended covariance matrix is calculated with the sample-covariance
normalization factor $1/(n-1)$.

### 8.3 Initial $\lambda$ selection

The initial $\lambda$ search fits observed-only graphical models to the
training and test subsets. For each candidate $\lambda$, HR-GGM measures the
shared same-sign edge strength and the difference between the fitted
observed-variable precision matrices. The scores are averaged across folds and
normalized.

The initial $\lambda$ is the candidate with the largest value in:

```text
Result/lambda_ratio/Score_norm_0.txt
```

### 8.4 Alternating $\alpha$ and $\lambda$ selection

After the initial observed-only $\lambda$ is selected, the precision matrix is
expanded to include the hidden master regulatory and hidden global variables.

$\alpha$ selection uses:

* held-out Gaussian loss;
* normalized entropy of the observed--global block;
* the retained relative magnitude of the global contribution.

Entropy measures how broadly the observed--global associations are distributed,
whereas the magnitude term prevents a diffuse but negligible global effect
from receiving a high adjusted-entropy score. At outer loop `<loop>`, the
$\alpha$-selection quantities are stored under
`Alpha_metric/` using the same loop index. The combined score is
`Score_vec_<loop>.txt`, and the candidate $\alpha$ with the smallest score is
selected.

The subsequent $\lambda$ search evaluates the stability of the
observed-variable blocks from full HR-GGM models fitted to the training and
test data. At each outer loop, its outputs are saved
under `lambda_ratio/` as `Score_norm<loop>.txt` and
`Score_ratio<loop>.txt`. The suffix indicates the selection-loop number (e.g., 1, 2, ...).

The candidate $\lambda$ with the largest normalized score is selected.

$\alpha$ and $\lambda$ selection alternate until the selected value is unchanged or
the maximum number of outer selection loops is reached.

### 8.5 Final fit

The selected $\lambda$ and $\alpha$, together with the user-specified
$\gamma$, are used to fit HR-GGM to the complete input matrix. Precision-matrix
optimization and estimation of the unobserved columns of the extended data
matrix alternate until the stopping rule is satisfied or the maximum number of
data-update loops is reached. The covariance matrix is recalculated after each
extended-data update.

## 9. Output files

The output directory has the structure:

```text
Result/
├── Theta_final.txt
├── lambda_ratio/
│   ├── Score_norm_0.txt
│   ├── Score_ratio_0.txt
│   ├── Score_norm<loop>.txt
│   └── Score_ratio<loop>.txt
└── Alpha_metric/
    ├── Likelihood_mat_<loop>.txt
    ├── Entropy_mat_<loop>.txt
    ├── Global_magnitude_mat_<loop>.txt
    ├── Retention_vec_<loop>.txt
    ├── Adjusted_entropy_vec_<loop>.txt
    └── Score_vec_<loop>.txt
```

### 9.1 `Theta_final.txt`

`Theta_final.txt` is the final precision matrix fitted to all samples.

Its dimension is:

$$
(N_O+N_H+N_G)
\times
(N_O+N_H+N_G).
$$

For the example:

```text
(200 + 1 + 1) x (200 + 1 + 1) = 202 x 202
```

The C++ matrix-index ranges are:

|Variable type|C++ index range|Text-file row/column numbers|
|-|-|-|
|Observed variable|`0:199`|`1:200`|
|Hidden master regulatory variable|`200`|`201`|
|Hidden global variable|`201`|`202`|

Therefore:

* `Theta_final.txt(1:200,1:200)` is the observed-variable block when loaded in MATLAB;
* `Theta_final.txt(1:200,201)` is the observed--hidden-regulatory block;
* `Theta_final.txt(1:200,202)` is the observed--global block;
* row/column 201 represents the hidden master regulatory variable;
* row/column 202 represents the hidden global variable.

Load the matrix in MATLAB with:

```matlab
Theta = readmatrix('Result/Theta_final.txt');

Nlocal = 200;  % legacy code name: number of observed variables
Nhidden = 1;   % number of hidden master regulatory variables
Nglobal = 1;   % number of hidden global variables

Theta_O = Theta(1:Nlocal, 1:Nlocal);
Theta_OH = Theta(1:Nlocal, Nlocal + (1:Nhidden));
Theta_OG = Theta(1:Nlocal, Nlocal + Nhidden + (1:Nglobal));
```

### 9.2 Conditional-association interpretation

An off-diagonal precision-matrix value represents a conditional association
under the fitted Gaussian graphical model. A corresponding partial correlation
can be calculated as:

$$
\rho_{ij\mid -ij}
=
-\frac{\Theta_{ij}}
{\sqrt{\Theta_{ii}\Theta_{jj}}}.
$$

In MATLAB:

```matlab
d = sqrt(diag(Theta));
PartialCorr = -Theta ./ (d * d.');
PartialCorr(1:size(Theta,1)+1:end) = 1;
```

The estimated graph is an undirected conditional-association network. The
precision matrix alone does not determine causal direction.

### 9.3 Parameter-selection files

The files under `lambda_ratio/` and `Alpha_metric/` retain the fold- and
candidate-specific quantities used during parameter selection. The selected
$\lambda$, $\alpha$, and $\gamma$ are also printed to the terminal and are retained when
the run is executed with `tee` as shown above.

## 10. Running a user-provided dataset

For a new dataset:

1. complete normalization, transformation, quality control, and feature
selection outside HR-GGM;
2. arrange the final matrix as observed variables by samples;
3. preserve the observed-variable order;
4. identify the prior-guided observed-variable subset for initializing the
   hidden master regulatory variable;
5. convert any one-based row numbers to C++ indices by subtracting one;
6. list all observed-variable indices in `Idx_g_list.txt` under the current
   hidden-global initialization strategy;
7. prepare the required files in the `Data` directory;
8. edit `Nhidden`, `Nglobal`, `Gamma`, thread count, and other settings in
   `Main.cpp`;
9. rebuild the executable;
10. place the new `Data` directory beside the executable;
11. run HR-GGM and retain the terminal log;
12. interpret the final precision matrix using the
    observed--hidden-regulatory--global ordering.

For example, if the prior set for the hidden master regulator contains MATLAB
rows 3, 10, and 25:

```matlab
idx_h_matlab = [3; 10; 25];
idx_h_cpp = idx_h_matlab - 1;
```

The index file must contain:

```text
2
9
24
```

## 11. Computational considerations

The implementation forms dense covariance and inverse-related matrices during
optimization. Memory and runtime increase substantially with the number of
observed variables. The number of samples also affects covariance construction
and unobserved-data updating.

Recommendations:

* use a release build;
* select an appropriate number of OpenMP threads for the available hardware;
* avoid requesting more threads than the server allocation;
* keep `Blockcols` at a value supported by available memory;
* begin with the supplied example before running a large dataset;
* retain the terminal log for long-running jobs;
* use a persistent terminal session such as `tmux` or the cluster scheduler for
server runs.

## 12. Troubleshooting

### `Could not locate 'Data' directory`

The executable could not find a directory named exactly `Data` in its own
directory or any parent directory. Place `Data` at `build/Data` when running
`build/C_GRN`.

### `open failed: .../Data_whole.bin`

Confirm that `Data_whole.bin` exists in the located `Data` directory and that
its name and capitalization are unchanged.

### OpenMP is not found during CMake configuration

Load or install a compiler with OpenMP support. On Linux, GCC is the recommended
compiler for the supplied build configuration.

### Cholesky or positive-definiteness failure

Check that:

* the input contains only finite values;
* constant observed variables were removed;
* the parameter grids contain suitable positive values;
* the matrix dimensions and index files are consistent;
* the initialization matrix corresponds to the selected model dimensions.

### Invalid or unexpected unobserved-variable estimates

Verify the C++ indices against the original MATLAB matrix:

```matlab
cpp_idx = readmatrix('Data/Idx_h_list.txt');
matlab_rows = cpp_idx + 1;
Prior_data = Data(matlab_rows, :);
```

Do not use the stored C++ indices directly as MATLAB subscripts without adding
one.

### The executable was rebuilt but an old version is running

Confirm the executable path:

```bash
realpath ./build/C_GRN
ls -l ./build/C_GRN
```

Run the executable from the same build directory that was rebuilt.

## 13. Reproducibility checklist

Record the following information for every analysis:

* HR-GGM code version or Git commit;
* dimensions and preprocessing of `Data_whole.bin`;
* ordered observed-variable list;
* complete `Idx_h_list.txt` and `Idx_g_list.txt`;
* lambda and alpha candidate grids;
* all adjustable values in `Main.cpp`;
* compiler, CMake version, and number of OpenMP threads;
* selected lambda, alpha, and gamma from the terminal log;
* final `Theta_final.txt`;
* parameter-selection metric directories.

## 14. Citation and contact

HR-GGM paper: the complete citation will be added upon publication.

Eigen: HR-GGM uses Eigen for dense and sparse matrix computations in its C++
implementation. Please cite:
Gaël Guennebaud, Benoît Jacob, and others. (2010). Eigen.
https://libeigen.gitlab.io.

Contact us: songtaow@vt.edu, yzwang@vt.edu