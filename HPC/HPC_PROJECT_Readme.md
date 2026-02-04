# HPC Quantum Noise Benchmarking with QuEST

## Overview

This project involves benchmarking a quantum algorithm (Bell state creation) using the [QuEST](https://quest.qtechtheory.org/) framework. A key contribution of the project is the integration of a **tracked quantum noise simulation function** using **OpenMP-based parallel programming** to measure decoherence error more efficiently. The results of the simulation are stored in a CSV file for analysis.

The file structure follows the standard QuEST repository, with additional files and modifications relevant to this project. This repository is submitted as a zip file. Upon unzipping, it creates a top-level folder named `QuEST-main`.

---

## Folder Structure

```bash
QuEST-main/
├── build/ # (Generated after building) Contains compiled libQuEST.so and output
├── quest/
│ ├── include/ # Header files of QuEST
│ └── src/ # Source files (API, CPU, GPU, MPI communication logic)
├── simulate_bell.cpp # Main benchmarking script (calls new noise-tracked API)
├── CMakeLists.txt # CMake build script
├── README.md # This file
```

---

## Prerequisites (Linux)

Before building and running the simulation, ensure the following dependencies are installed:

```bash
# Update packages
sudo apt update

# Install C++ compiler
sudo apt install g++

# Install OpenMP (usually bundled with g++)
sudo apt install libomp-dev

# Install CMake (for building QuEST)
sudo apt install cmake

# (Optional) Check CMake and g++ versions
cmake --version
g++ --version
```

---

## Build and Execution Steps

#### 💡 IMPORTANT: First, unzip the submitted file to create a folder QuEST-main. Then follow the steps below from inside that folder.

### 1. Create a build directory and navigate into it
```bash
mkdir build
cd build
```
---
### 2. Build QuEST using CMake
```bash
cmake .. \
  -DCOMPILE_OPENMP=1 \
  -DCOMPILE_MPI=0 \
  -DCOMPILE_CUDA=0 \
  -DCOMPILE_CUQUANTUM=0 \
  -DCOMPILE_ZE=0 \
  -DCOMPILE_KOKKOS=0
```
```bash
make -j
```
This will generate the static/shared QuEST library (libQuEST.so) in the build/ directory.

---

### 3. Compile the benchmark script
```bash
g++ ../simulate_bell.cpp \
    -I../quest/include \
    -L. -lQuEST -lm \
    -DCOMPILE_OPENMP=1 -DCOMPILE_MPI=0 -DCOMPILE_CUDA=0 -DCOMPILE_CUQUANTUM=0 \
    -o output1
```
This compiles simulate_bell.cpp using the built libQuEST.so.

---

### 4. Run the simulation
```bash
LD_LIBRARY_PATH=. ./output1
```
This will execute the Bell state benchmarking with and without noise tracking, and generate an output file:
```bash
    results.csv: Contains a row of the form
    Bell,<noise_probability>,<useTracked_flag>,<measured_error>
```
Example Output (in results.csv)
```bash
Bell,0.3,0,0
Bell,0.3,1,0.84
```
This output shows how much decoherence error was measured using the tracked function (with OpenMP), compared to the untracked version.

---

## Notes

You can rerun the simulation by deleting the existing results.csv file in the build/ folder.

You can adjust the noise probability and toggle tracking within the simulate_bell.cpp file.

---

## License

This project builds upon the QuEST framework which is open-source under the MIT license. This code is submitted solely for academic evaluation purposes.