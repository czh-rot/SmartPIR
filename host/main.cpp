#include "utils.hpp"
#include "answer.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>
#include <math.h>

#include <iomanip>
#include <iosfwd>
#include <unistd.h>
#include <fcntl.h>


void load_data(std::string filepath,unsigned long long* data_in, int size){
    char line[size] = {0};
    int index = 0;

    std::fstream filestream(filepath.c_str(), std::ios::in);
    if(!filestream){
        std::cout << "Error : " << filepath << " file doesn't exist !" << std::endl;
        exit(1);
    }
    while(filestream.getline(line, sizeof(line), '\n')){
        std::stringstream ss(line);
        ss >> data_in[index];
        index++;
    }
}

void estimated_cryptonets(std::string xclbin_path) {
    xf::common::utils_sw::Logger logger(std::cout, std::cerr);
    cl_int err;

    // timing checker
    struct timeval start_time;

    // Getting Device
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[2];
    std::string devName = device.getInfo<CL_DEVICE_NAME>();
    printf("Found Device=%s\n", devName.c_str());
    devices.resize(1);
    
    // Creating Context and Command Queue for selected Device
    cl::Context context(device, NULL, NULL, NULL, &err);  // create opencl context
    logger.logCreateContext(err);
    cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err); // create cq, non-sync performance analysis
    logger.logCreateCommandQueue(err);

    // Import Binary File
    std::cout<<"binaries xclbins read to import"<<std::endl;
    cl::Program::Binaries xclBins = xcl::import_binary_file(xclbin_path);
    if(xclBins.size() == 0) {
        printf("Error: Failed to import Xclbin");
        exit(1);
    }

    // Create Program and Kernel
    std::vector<cl::Device> contextDevices;
    contextDevices.push_back(device);
    cl::Program program(context, contextDevices, xclBins, NULL, &err);
    logger.logCreateProgram(err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to create program");
    }
    cl::Kernel cryptonets;
    cryptonets = cl::Kernel(program, "cryptonets", &err); // apply the kernel name
    logger.logCreateKernel(err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to create program");
    } else {
        std::cout << "kernel has been created" << std::endl;
    }
 

    unsigned long long averagetime = 0;
    unsigned long long averagetime2 = 0;

    // worksize: the fix space in SmartSSD DRAM
    int worksize = sizeof(unsigned long long) * 16384
    size_t vector_size_bytes = worksize;
    
    cl_mem_ext_ptr_t data_in_ssd;    // define a OpenCL extension pointer
    data_in_ssd = {XCL_MEM_EXT_P2P_BUFFER, nullptr, 0};  // flag, ass obj, ass param
    //data_in_ssd.banks = XCL_MEM_DDR_BANK0;
    cl::Buffer data_in_buf(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX, vector_size_bytes, &data_in_ssd);

    // enqueuteMapBuffer map the data_in_buf to host
    unsigned long long* data_in_ptr = (unsigned long long*)q.enqueueMapBuffer(data_in_buf,CL_TRUE,CL_MAP_READ | CL_MAP_WRITE,0,vector_size_bytes,nullptr, nullptr);
    q.finish();

    // 3-D Matrix store
    std::string data_in_filepath="/mnt/smartssd/Matrix.txt"; // file path

    std::chrono::high_resolution_clock::time_point p2pStart = std::chrono::high_resolution_clock::now();
    // P2P transfer from NVMe to DDR using pread
    int fd = open(data_in_filepath.c_str(), O_RDONLY|O_DIRECT);  // Open file in read-only mode
    if (fd == -1) {
        // Handle error if the file cannot be opened
        perror("Error opening file");
        return -1;
    }

    // Read 16384 bytes from the file starting at offset 0
    ssize_t bytesRead = pread(fd, data_in_ptr, worksize, 0);
    if (bytesRead == -1) {
        // Handle error if read fails
        perror("Error reading file");
        close(fd);
        return -1;
    }

    close(fd);
    q.finish();
    std::chrono::high_resolution_clock::time_point p2pEnd = std::chrono::high_resolution_clock::now();
    cl_ulong p2pTime = std::chrono::duration_cast<std::chrono::microseconds>(p2pEnd - p2pStart).count();
    averagetime+=p2pTime;

    int j = 0;
    // 4 param. allocate these param
    cryptonets.setArg(j++, data_in_buf);

    std::chrono::high_resolution_clock::time_point p2pStart2 = std::chrono::high_resolution_clock::now();
    
    q.enqueueTask(cryptonets);
    q.finish();

    std::chrono::high_resolution_clock::time_point p2pEnd2 = std::chrono::high_resolution_clock::now();
    cl_ulong p2pTime2 = std::chrono::duration_cast<std::chrono::microseconds>(p2pEnd2 - p2pStart2).count();
    averagetime2+=p2pTime2;
    
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "INFO: Kernel execution average time: " << averagetime2 << " us\n";
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "verify the data transfer method" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "INFO: P2P transfer average time: " << averagetime << " us\n";
    std::cout << "-------------------------------------------------------" << std::endl;
}

int main(int argc, const char* argv[]) {
    std::cout << "\n-------------------estimated diameter----------------\n";
    //Command Line Parser
    ArgParser parser(argc, argv);
    
    // print the argument
    std::cout << "\n-------------------print param----------------\n";
    std::cout<< argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "argv[" << i << "]:" << argv[i] << std::endl;
    }

    std::cout<<"\n-------------------parser finished----------------\n"<<std::endl;

    std::string tmpStr;
    if (!parser.getCmdOption("--xclbin", tmpStr)) {
        std::cout << "ERROR:xclbin path is not set!\n";
        return 1;
    }
    std::string xclbin_path = tmpStr;
    std::cout << xclbin_path << std::endl;    

    struct timeval start_time, end_time;
    gettimeofday(&start_time, 0);
    std::cout<< "estimated_cryptonets ready"<<std::endl;

    estimated_cryptonets(xclbin_path);
  
    gettimeofday(&end_time, 0);
    std::cout << "User function execution time is: " << tvdiff(&start_time, &end_time) / 1000UL << " ms" << std::endl;

    return 0;
}

