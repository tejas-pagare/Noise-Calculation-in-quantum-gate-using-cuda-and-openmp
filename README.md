# HPC Quantum Noise Benchmarking with QuEST

## Overview

This project benchmarks quantum algorithms using the [QuEST](https://quest.qtechtheory.org/) (Quantum Exact Simulation Toolkit) framework. A key contribution is the integration of a **tracked quantum noise simulation function** using **OpenMP-based parallel programming** to measure decoherence error efficiently.

## Project Structure

```
├── HPC/
│   ├── Codes/                  # Main source code files
│   │   ├── HEA_*.cpp           # Hardware Efficient Ansatz implementations
│   │   └── ...
│   ├── QuEST-main/             # QuEST framework with modifications
│   └── QuEST-main_Code/        # Additional QuEST code variants
│
├── OUTPUT/
│   ├── Density/                # Density matrix simulation results
│   │   └── final graphs/       # Generated visualizations
│   └── Pure state/             # Pure state simulation results
│       └── final graphs/       # Generated visualizations
│
└── README.md
```

## Key Features

- **Pure State Simulations**: Efficient state vector simulations
- **Density Matrix Simulations**: Full density matrix approach for noise modeling
- **Noise Benchmarking**: Analysis of noise vs error relationships
- **Parallel Performance**: OpenMP parallelization for multi-threaded execution
- **Comprehensive Metrics**: Threads vs time, qubits vs error analysis

## Prerequisites

### Linux
```bash
sudo apt update
sudo apt install g++ cmake libomp-dev
```

### macOS
```bash
brew install cmake llvm libomp
```

## Build Instructions

1. Navigate to the QuEST directory:
```bash
cd HPC/QuEST-main
mkdir build && cd build
```

2. Configure with CMake:
```bash
cmake .. \
  -DCOMPILE_OPENMP=1 \
  -DCOMPILE_MPI=0 \
  -DCOMPILE_CUDA=0 \
  -DCOMPILE_CUQUANTUM=0 \
  -DCOMPILE_ZE=0 \
  -DCOMPILE_KOKKOS=0
```

3. Build:
```bash
make -j
```

4. Compile your benchmark:
```bash
g++ ../simulate_bell.cpp \
    -I../quest/include \
    -L. -lQuEST -lm \
    -DCOMPILE_OPENMP=1 -DCOMPILE_MPI=0 -DCOMPILE_CUDA=0 -DCOMPILE_CUQUANTUM=0 \
    -o output
```

5. Run:
```bash
LD_LIBRARY_PATH=. ./output
```

## Results

Output CSV files are generated in the `OUTPUT/` directory with analysis data for:
- Noise probability vs measured error
- Number of qubits vs error
- Thread count vs execution time
- Sequential vs parallel performance comparisons

## License

See [LICENCE.txt](HPC/QuEST-main/LICENCE.txt) for QuEST licensing information.

## Authors

See [AUTHORS.txt](HPC/QuEST-main/AUTHORS.txt) for contributor information.
