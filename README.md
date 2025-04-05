# SmartPIR

SmartPIR is an efficient Index PIR (Private Information Retrieval) system developed via a protocol and hardware co-design approach. The project now includes three main implementations:  
1. **CPU-based SmartPIR** (relying on the Microsoft SEAL library)  
2. **GPU-based SmartPIR** (template for GPU acceleration, to be updated)  
3. **CSD-based SmartPIR** (leveraging Samsung SmartSSD and Xilinx FPGA for acceleration)

---

## Features

- **Flexible PIR Protocol**  
  - Efficiently handles variable-length entries, maintaining high performance for large-scale datasets.  
  - Implements a sophisticated BFV (Brakerski-Fan-Vaikuntanathan) homomorphic encryption scheme for secure data retrieval.

- **CPU-based Implementation**  
  - Utilizes [Microsoft SEAL](https://github.com/microsoft/SEAL) for homomorphic operations.  
  - Built with standard C++ tools (CMake, etc.) for easy development and deployment on typical CPU-based environments.

- **GPU-based Implementation** (Template)  
  - (Planned) Accelerates PIR using GPU resources.  
  - Will use specialized libraries and frameworks to offload homomorphic operations onto GPUs.  
  - To be updated in future releases.

- **CSD-based Implementation**  
  - Employs a Samsung SmartSSD (with an onboard Xilinx FPGA) to execute homomorphic encryption kernels near the storage.  
  - Reduces data movement by offloading significant computation onto the FPGA, enabling higher throughput and lower latency.

---

## Repository Structure

This repository is organized into the following high-level folders (plus supporting documentation):

1. **CPU-based SmartPIR/**  
   - Contains source code for the CPU-only PIR solution, relying on Microsoft SEAL.

2. **GPU-based SmartPIR/** (Template)  
   - Intended for a future GPU-accelerated PIR implementation.  
   - Currently a placeholder for code structure and build scripts (to be added).

3. **CSD-based SmartPIR/**  
   - Host-side code (CPU) and FPGA kernel code for SmartSSD-based BFV acceleration.  
   - Includes High-Level Synthesis (HLS) code (`kernel/`) and host orchestration code (`host/`).

4. **README** (this file)  
   - Provides an overview of SmartPIR, features, directory organization, and build instructions.

---

## Build & Run Instructions

The instructions below outline how to build and run each of the three implementations.

### 1. CPU-based SmartPIR

This version uses [Microsoft SEAL](https://github.com/microsoft/SEAL) for homomorphic operations.

1. **Prerequisites**  
   - A C++17 (or later) compiler.  
   - CMake (3.10+ recommended).  
   - The Microsoft SEAL library (you can build it from source or install via package managers).  

2. **Build and Install**  
   1. Navigate to the `CPU-based SmartPIR/` directory:
      ```bash
      cd CPU-based\ SmartPIR
      ```
   2. Build and install the project:
      ```bash
      cmake -S . -B build
      cmake --build build
      sudo cmake --install build
      ```
      > Note: You may specify a custom install prefix by adding `-DCMAKE_INSTALL_PREFIX=<path>` during the initial CMake configuration.

3. **Run Examples**  
   After installation, you can run sample programs (such as `sealexamples`):
   ```bash
   cd ./native/examples
   cmake -S . -B build
   cmake --build build
   cd build/bin
   ./sealexamples 1
