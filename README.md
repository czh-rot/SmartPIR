# SmartPIR

SmartPIR is an efficient Index PIR (Private Information Retrieval) system designed through protocol and architecture co-design. By optimizing both the PIR protocol and hardware architecture, SmartPIR delivers high-performance and privacy-preserving solutions for large-scale data retrieval.

------

### Features

- **Protocol**: A sophisticated PIR protocol capable of handling variable-length entries efficiently, ensuring minimal overhead and maximum performance for data retrieval.
- **Architecture**: Utilizes SAMSUNG SmartSSDs to enhance both computational and I/O performance. The architecture optimizes data flow between the CPU and SmartSSD, significantly reducing latency and improving throughput.

------

### Code Organization

The repository is organized into two primary directories:

```
/host/          CPU+SmartSSD Coordination
/kernel/        BFV kernel implemented in FPGA HLS code
```

- **/host/**: Contains the code responsible for coordinating tasks between the CPU and SmartSSD, including system-level orchestration.
- **/kernel/**: Implements the BFV (Brakerski-Fan-Vaikuntanathan) cryptographic kernel on FPGA using High-Level Synthesis (HLS) to accelerate Homomorphic operations.

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
