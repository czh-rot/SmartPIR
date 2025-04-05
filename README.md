<<<<<<< HEAD
# SmartPIR

SmartPIR is an efficient Index PIR (Private Information Retrieval) system designed through protocol and architecture co-design. By optimizing both the PIR protocol and hardware architecture, SmartPIR delivers high-performance and privacy-preserving solutions for large-scale data retrieval.

------

### Features

- **Protocol**: A sophisticated PIR protocol capable of handling variable-length entries efficiently, ensuring minimal overhead and maximum performance for data retrieval.
- **Architecture**: Utilizes SmartSSDs (Smart Solid-State Drives) to enhance both computational and I/O performance. The architecture optimizes data flow between the CPU and SmartSSD, significantly reducing latency and improving throughput.

------

### Code Organization

The repository is organized into two primary directories:

```
/host/          CPU+SmartSSD Coordination
/kernel/        BFV kernel implemented in FPGA HLS code
```

- **/host/**: Contains the code responsible for coordinating tasks between the CPU and SmartSSD, including system-level orchestration.
- **/kernel/**: Implements the BFV (Brakerski-Fan-Vaikuntanathan) cryptographic kernel on FPGA using High-Level Synthesis (HLS) to accelerate encryption and decryption operations.

------

### Build Instructions

Before running any tests, you need to construct a 3-D matrix and store it on the SmartSSD. This process ensures that the data is properly prepared for the PIR protocol.

After that, you can build and run the system with the following command:

```
make run TARGET=hw PLATFORM=xilinx_u2_gen3x4_xdma_gc_2_202110_1
```

This command compiles and executes the program targeting the specified hardware platform (`xilinx_u2_gen3x4_xdma_gc_2_202110_1`), which is an FPGA-based setup for high-performance computation.

------

### License

This project is licensed under the MIT License.
=======
# SmartPIR

SmartPIR is an efficient Index PIR (Private Information Retrieval) system developed via a protocol and hardware co-design approach. The project now includes three main implementations:  

1. **CPU-based SmartPIR** (relying on the Microsoft SEAL library)  
2. **GPU-based SmartPIR** (inspired by HEonGPU)  
3. **CSD-based SmartPIR** (leveraging Samsung SmartSSD and Xilinx FPGA for acceleration)

---

## Features

- **Flexible PIR Protocol**  
  - Efficiently handles variable-length entries, maintaining high performance for large-scale datasets.  
  - Implements a sophisticated BFV (Brakerski-Fan-Vaikuntanathan) homomorphic encryption scheme for secure data retrieval.

- **CPU-based Implementation**  
  - Utilizes [Microsoft SEAL](https://github.com/microsoft/SEAL) for homomorphic operations.  
  - Built with standard C++ tools (CMake, etc.) for straightforward development and deployment on typical CPU-based environments.

- **GPU-based Implementation**  
  - Inspired by the [HEonGPU project](https://github.com/Alisah-Ozcan/HEonGPU), aiming to accelerate PIR queries using GPU resources.  
  - Offloads homomorphic operations onto GPUs by leveraging CUDA/OpenCL and specialized libraries.  
  - Provides significantly faster PIR operations on systems equipped with high-performance GPUs.

- **CSD-based Implementation**  
  - Employs a Samsung SmartSSD (with an onboard Xilinx FPGA) to execute homomorphic encryption kernels near the storage.  
  - Reduces data movement by offloading a significant portion of computation onto the FPGA, enabling higher throughput and lower latency.

---

## Repository Structure

This repository is organized into the following high-level folders (plus supporting documentation):

1. **CPU-based SmartPIR/**  
   - Contains source code for the CPU-only PIR solution, relying on Microsoft SEAL.

2. **GPU-based SmartPIR/**  
   - Based on the core ideas of [HEonGPU](https://github.com/Alisah-Ozcan/HEonGPU) for GPU acceleration.  
   - Includes build scripts, CUDA/OpenCL kernels, and sample usage for running BFV-based PIR on GPUs.

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
   ```

### 2. GPU-based SmartPIR

This version aims to accelerate BFV homomorphic operations by offloading them to GPUs, inspired by the [HEonGPU](https://github.com/Alisah-Ozcan/HEonGPU) approach.

1. **Prerequisites**  

   - A GPU with CUDA or OpenCL support.  
   - C++17 (or later), CMake (3.10+), and any required GPU libraries (e.g., CUDA Toolkit or OpenCL headers).  

2. **Build Instructions**  

   - Navigate to the `GPU-based SmartPIR/` directory:

     ```bash
     cd GPU-based\ SmartPIR
     ```

   - Configure and build (example with CUDA):

     ```bash
     $ cmake -S . -D HEonGPU_BUILD_BENCHMARKS=ON -D CMAKE_CUDA_ARCHITECTURES=89 -B build
     $ cmake --build ./build/
     
     $ ./build/bin/benchmark/<...>
     $ Example: ./build/bin/benchmark/smartPIR_benchmark
     ```

     > Adjust flags or environment variables based on your GPU and toolkit.

3. **Usage**  

   - (Placeholder) After building, you can run the GPU-accelerated PIR demos:

     ```bash
     cd build/bin
     ./gpu_pir_demo
     ```

   - Future updates will provide example configurations and performance benchmarks.

### 3. CSD-based SmartPIR

This version targets a Samsung SmartSSD (with a Xilinx FPGA) to run the BFV kernel.

1. **Prerequisites**  

   - Xilinx toolchain (e.g., Vitis, Vivado) compatible with your FPGA platform.  
   - The specific hardware platform (e.g., `xilinx_u2_gen3x4_xdma_gc_2_202110_1`) or equivalent.  

2. **Directory Structure**  

   - **`host/`**: Contains the CPU-side code that communicates with the SmartSSD.  
   - **`kernel/`**: Includes Xilinx High-Level Synthesis (HLS) code to implement BFV in hardware.

3. **Build and Run**  

   - Navigate to `CSD-based SmartPIR/` and run:

     ```bash
     make run TARGET=hw PLATFORM=xilinx_u2_gen3x4_xdma_gc_2_202110_1
     ```

     > Replace `PLATFORM` with your target hardware platform as necessary.

4. **Data Preparation**  

   - Before running, construct the data (e.g., a 3D matrix) and store it on the SmartSSD.  
   - Project-specific scripts or documentation will guide you on how to load the data into the SmartSSD.

---

## License

SmartPIR is licensed under the [MIT License](./LICENSE). You are free to use, modify, and distribute this project. For academic use, please cite the original SmartPIR publication or mention this repository.
>>>>>>> cd4fcd7ec03807baf29b38fa1bad0309733ec26c
