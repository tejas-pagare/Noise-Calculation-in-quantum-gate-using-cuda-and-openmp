# CUDA HEA Simulation - Build and Run Instructions

## Prerequisites

Before running, ensure you have the following installed on your target system:

1. **NVIDIA CUDA Toolkit** (version 11.0 or later recommended)
   ```bash
   # Check if CUDA is installed
   nvcc --version
   
   # Check if GPU is available
   nvidia-smi
   ```

2. **CMake** (version 3.21 or later)
   ```bash
   cmake --version
   ```

3. **C++ Compiler** (g++ 9+ or clang++)
   ```bash
   g++ --version
   ```

---

## Step-by-Step Build Instructions

### Step 1: Navigate to the QuEST directory
```bash
cd /path/to/HPC/QuEST-main
```

### Step 2: Create a CUDA build directory
```bash
mkdir -p build_cuda
cd build_cuda
```

### Step 3: Configure CMake with CUDA enabled
```bash
# First, find your GPU compute capability
nvidia-smi --query-gpu=compute_cap --format=csv

# Then configure with the appropriate architecture (e.g., 75 for RTX 20xx/GTX 16xx)
cmake .. \
    -DENABLE_CUDA=ON \
    -DENABLE_MULTITHREADING=OFF \
    -DENABLE_DISTRIBUTION=OFF \
    -DCMAKE_CUDA_ARCHITECTURES=75 \
    -DCMAKE_BUILD_TYPE=Release
```

**Important CMake Options:**
| Option | Value | Description |
|--------|-------|-------------|
| `ENABLE_CUDA` | ON | Enable NVIDIA GPU acceleration |
| `ENABLE_MULTITHREADING` | OFF | Disable OpenMP (single thread) |
| `ENABLE_DISTRIBUTION` | OFF | Disable MPI |
| `CMAKE_CUDA_ARCHITECTURES` | 75 | **REQUIRED** - Your GPU compute capability |
| `CMAKE_BUILD_TYPE` | Release | Optimized build |

**CUDA Architecture Reference:**
| GPU Series | Architecture Code |
|------------|------------------|
| GTX 10xx (Pascal) | 61 |
| GTX 16xx, RTX 20xx (Turing) | 75 |
| RTX 30xx (Ampere) | 86 |
| RTX 40xx (Ada Lovelace) | 89 |
| Tesla V100 | 70 |
| A100 | 80 |

**Or use "native" for auto-detection (CMake 3.24+):**
```bash
cmake .. -DENABLE_CUDA=ON -DENABLE_MULTITHREADING=OFF -DENABLE_DISTRIBUTION=OFF -DCMAKE_CUDA_ARCHITECTURES=native
```

### Step 4: Build the QuEST library
```bash
make -j$(nproc)
```

### Step 5: Copy the CUDA source file to build directory
```bash
cp ../Codes/HEA_cuda_single_thread.cu .
# OR if file is elsewhere:
cp /path/to/HEA_cuda_single_thread.cu .
```

### Step 6: Compile the CUDA application
```bash
nvcc -o hea_cuda_single HEA_cuda_single_thread.cu \
    -I../quest/include \
    -I.. \
    -L./quest -lQuEST \
    -std=c++17 \
    --expt-relaxed-constexpr
```

**Alternative compilation using g++ (if QuEST was built with CUDA):**
```bash
g++ -o hea_cuda_single HEA_cuda_single_thread.cu \
    -I../quest/include \
    -I.. \
    -L./quest -lQuEST \
    -std=c++17 \
    -lcudart
```

### Step 7: Set library path and run
```bash
# Set the library path
export LD_LIBRARY_PATH=$PWD/quest:$LD_LIBRARY_PATH

# Run the simulation
./hea_cuda_single
```

---

## Quick One-Liner Script

Save this as `build_and_run_cuda.sh`:

```bash
#!/bin/bash
set -e

QUEST_DIR="$(pwd)"

# Create and enter build directory
mkdir -p build_cuda && cd build_cuda

# Configure with CUDA
cmake .. -DENABLE_CUDA=ON -DENABLE_MULTITHREADING=OFF -DENABLE_DISTRIBUTION=OFF -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Copy source file
cp ../Codes/HEA_cuda_single_thread.cu . 2>/dev/null || cp ../../Codes/HEA_cuda_single_thread.cu .

# Compile application
nvcc -o hea_cuda_single HEA_cuda_single_thread.cu \
    -I../quest/include -I.. -L./quest -lQuEST \
    -std=c++17 --expt-relaxed-constexpr

# Run
export LD_LIBRARY_PATH=$PWD/quest:$LD_LIBRARY_PATH
./hea_cuda_single
```

Run with:
```bash
chmod +x build_and_run_cuda.sh
./build_and_run_cuda.sh
```

---

## Troubleshooting

### Issue: "CUDA not found"
```bash
# Add CUDA to PATH
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

### Issue: "No GPU found"
```bash
# Check GPU availability
nvidia-smi
# If no output, ensure NVIDIA drivers are installed
```

### Issue: "Library not found"
```bash
# Ensure LD_LIBRARY_PATH includes the QuEST library
export LD_LIBRARY_PATH=/path/to/build_cuda/quest:$LD_LIBRARY_PATH
```

### Issue: CMake version too old
```bash
# Install newer CMake
pip install cmake --upgrade
# Or download from https://cmake.org/download/
```

---

## Configuration Parameters

You can modify these constants in `HEA_cuda_single_thread.cu`:

```cpp
const int NUM_QUBITS = 12;      // Increase for larger simulations
const int NUM_LAYERS = 30;      // Number of HEA circuit layers
const qreal NOISE_LEVEL = 0.2;  // Noise probability (0.0 to 1.0)
```

**Note:** Density matrix simulations require memory of size 2^(2*NUM_QUBITS). 
- 12 qubits → ~128 MB GPU memory
- 14 qubits → ~2 GB GPU memory
- 16 qubits → ~32 GB GPU memory

---

## Expected Output

```
========================================
  CUDA Single-Thread HEA Simulation
========================================
Configuration:
  Qubits: 12
  Layers: 30
  Noise Level: 0.2
========================================

Starting CUDA-accelerated simulation...

========================================
  Results
========================================
Simulation Time: X.XXX seconds
========================================
```
