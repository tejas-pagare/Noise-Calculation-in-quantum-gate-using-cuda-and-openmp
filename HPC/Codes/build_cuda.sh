#!/bin/bash
# Build and Run CUDA HEA Simulation Script
# Usage: ./build_cuda.sh [path_to_quest_main]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  CUDA HEA Build Script${NC}"
echo -e "${GREEN}========================================${NC}"

# Check prerequisites
echo -e "\n${YELLOW}Checking prerequisites...${NC}"

# Check CUDA
if ! command -v nvcc &> /dev/null; then
    echo -e "${RED}ERROR: nvcc (CUDA compiler) not found!${NC}"
    echo "Please install CUDA Toolkit or add it to PATH:"
    echo "  export PATH=/usr/local/cuda/bin:\$PATH"
    exit 1
fi
echo -e "${GREEN}✓ CUDA found:${NC} $(nvcc --version | grep release)"

# Check CMake
if ! command -v cmake &> /dev/null; then
    echo -e "${RED}ERROR: cmake not found!${NC}"
    exit 1
fi
CMAKE_VERSION=$(cmake --version | head -n1)
echo -e "${GREEN}✓ CMake found:${NC} $CMAKE_VERSION"

# Check GPU
if command -v nvidia-smi &> /dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n1)
    echo -e "${GREEN}✓ GPU found:${NC} $GPU_NAME"
else
    echo -e "${YELLOW}WARNING: nvidia-smi not found. GPU may not be available.${NC}"
fi

# Determine QuEST directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUEST_DIR="${1:-$(dirname "$SCRIPT_DIR")/QuEST-main}"

if [ ! -f "$QUEST_DIR/CMakeLists.txt" ]; then
    echo -e "${RED}ERROR: QuEST not found at $QUEST_DIR${NC}"
    echo "Usage: $0 [path_to_quest_main]"
    exit 1
fi
echo -e "${GREEN}✓ QuEST found at:${NC} $QUEST_DIR"

# Create build directory
BUILD_DIR="$QUEST_DIR/build_cuda"
echo -e "\n${YELLOW}Creating build directory: $BUILD_DIR${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure CMake
echo -e "\n${YELLOW}Configuring CMake with CUDA...${NC}"
cmake .. \
    -DENABLE_CUDA=ON \
    -DENABLE_MULTITHREADING=OFF \
    -DENABLE_DISTRIBUTION=OFF \
    -DCMAKE_BUILD_TYPE=Release

# Build QuEST library
echo -e "\n${YELLOW}Building QuEST library...${NC}"
make -j$(nproc)

# Copy and compile the CUDA source
CUDA_SOURCE="$SCRIPT_DIR/HEA_cuda_single_thread.cu"
if [ ! -f "$CUDA_SOURCE" ]; then
    CUDA_SOURCE="$SCRIPT_DIR/../Codes/HEA_cuda_single_thread.cu"
fi

if [ -f "$CUDA_SOURCE" ]; then
    echo -e "\n${YELLOW}Compiling CUDA application...${NC}"
    cp "$CUDA_SOURCE" .
    
    nvcc -o hea_cuda_single HEA_cuda_single_thread.cu \
        -I../quest/include \
        -I.. \
        -L./quest -lQuEST \
        -std=c++17 \
        --expt-relaxed-constexpr \
        -Wno-deprecated-gpu-targets
    
    echo -e "\n${GREEN}========================================${NC}"
    echo -e "${GREEN}  Build Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "\nTo run the simulation:"
    echo -e "  cd $BUILD_DIR"
    echo -e "  export LD_LIBRARY_PATH=\$PWD/quest:\$LD_LIBRARY_PATH"
    echo -e "  ./hea_cuda_single"
    
    # Ask to run
    echo -e "\n${YELLOW}Run now? (y/n)${NC}"
    read -r response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        export LD_LIBRARY_PATH=$PWD/quest:$LD_LIBRARY_PATH
        ./hea_cuda_single
    fi
else
    echo -e "${RED}WARNING: CUDA source file not found at $CUDA_SOURCE${NC}"
    echo "Please copy HEA_cuda_single_thread.cu to $BUILD_DIR and compile manually"
fi
