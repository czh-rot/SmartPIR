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


void estimated_cryptonets(std::string xclbin_path,unsigned long long *ntt_out) {
    xf::common::utils_sw::Logger logger(std::cout, std::cerr);
    cl_int err;

    // timing checker
    struct timeval start_time;

    // Getting Device
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];
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
    
    // Create Program and Kernel
    cl::Program program(context, devices, xclBins, NULL, &err);
    logger.logCreateProgram(err);
    cl::Kernel cryptonets;
    cryptonets = cl::Kernel(program, "cryptonets", &err); // apply the kernel name
    logger.logCreateKernel(err);
    std::cout << "kernel has been created" << std::endl;


    unsigned long long averagetime = 0;
    unsigned long long averagetime2 = 0;
    std::vector<unsigned long long> hostBuffer(sizeof(unsigned long long) * 16384);
    const int loop = 1;
    for (int i = 0; i < loop; i++) {
        size_t vector_size_bytes=sizeof(unsigned long long) * 16384;
		
        cl_mem_ext_ptr_t ntt_in_ext;    // define a OpenCL extension pointer
        ntt_in_ext = {XCL_MEM_EXT_P2P_BUFFER, nullptr, 0};  // flag, ass obj, ass param
        //ntt_in_ext.banks = XCL_MEM_DDR_BANK0;
        cl::Buffer ntt_in_buf(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX, vector_size_bytes, &ntt_in_ext);

        cl_mem_ext_ptr_t rp_ext;
        rp_ext = {XCL_MEM_EXT_P2P_BUFFER, nullptr, 0};
        //rp_ext.banks = XCL_MEM_DDR_BANK1;
        cl::Buffer rp_buf(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX, vector_size_bytes, &rp_ext);

        cl_mem_ext_ptr_t srp_ext;
        srp_ext = {XCL_MEM_EXT_P2P_BUFFER, nullptr, 0};
        //srp_ext.banks = XCL_MEM_DDR_BANK2;
        cl::Buffer srp_buf(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX, vector_size_bytes, &srp_ext);

        //std::cout << "\nMap P2P device buffers to host access pointers\n" << std::endl;
        // enqueuteMapBuffer map the ntt_in_buf to host
        unsigned long long* ntt_in_ptr = (unsigned long long*)q.enqueueMapBuffer(ntt_in_buf,CL_TRUE,CL_MAP_READ | CL_MAP_WRITE,0,vector_size_bytes,nullptr, nullptr);
        unsigned long long* rp_ptr = (unsigned long long*)q.enqueueMapBuffer(rp_buf,CL_TRUE,CL_MAP_READ | CL_MAP_WRITE,0,vector_size_bytes,nullptr, nullptr);
        unsigned long long* srp_ptr = (unsigned long long*)q.enqueueMapBuffer(srp_buf,CL_TRUE,CL_MAP_READ | CL_MAP_WRITE,0,vector_size_bytes,nullptr, nullptr);
        q.finish();

        std::string ntt_in_filepath="/media/embedded-415/f29fd339-0d3e-4eb2-ae93-016cd9fcee27/ntt/data_16384/ntt_in.dat"; // file path
        std::string rp_filepath = "/media/embedded-415/f29fd339-0d3e-4eb2-ae93-016cd9fcee27/ntt/data_16384/rp.dat";
        std::string srp_filepath = "/media/embedded-415/f29fd339-0d3e-4eb2-ae93-016cd9fcee27/ntt/data_16384/srp.dat";
        q.finish();


        std::chrono::high_resolution_clock::time_point p2pStart = std::chrono::high_resolution_clock::now();

        //P2P transfer from NVMe to DDR
        load_data(ntt_in_filepath, ntt_in_ptr, 16384);
        load_data(rp_filepath, rp_ptr, 16384);
        load_data(srp_filepath, srp_ptr, 16384);
        q.finish();
        std::chrono::high_resolution_clock::time_point p2pEnd = std::chrono::high_resolution_clock::now();
        cl_ulong p2pTime = std::chrono::duration_cast<std::chrono::microseconds>(p2pEnd - p2pStart).count();
        averagetime+=p2pTime;
        // std::cout << "-------------------------------------------------------" << std::endl;
        // std::cout << "INFO: P2P transfer average time: " << averagetime / 100.0 << " us\n";
        // std::cout << "-------------------------------------------------------" << std::endl;
        // my test ntt_in_buf
        q.enqueueReadBuffer(ntt_in_buf, CL_TRUE, 0, vector_size_bytes, hostBuffer.data());
        // std::cout << "verify the data transfer method" << std::endl;
        // std::cout << hostBuffer[0] << " " << hostBuffer[1] << std::endl;
        


        // std::chrono::high_resolution_clock::time_point p2pEnd = std::chrono::high_resolution_clock::now();
        // cl_ulong p2pTime = std::chrono::duration_cast<std::chrono::microseconds>(p2pEnd - p2pStart).count();
        // averagetime+=p2pTime;

        cl_mem_ext_ptr_t ntt_out_ext={0};
        ntt_out_ext.obj=ntt_out;
        //ntt_out_ext.banks=XCL_MEM_DDR_BANK3;
        cl::Buffer ntt_out_buf(context, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,vector_size_bytes, &ntt_out_ext);

        int j = 0;
        // 4 param. allocate these param
        cryptonets.setArg(j++, ntt_in_buf);
        cryptonets.setArg(j++, ntt_out_buf);
        cryptonets.setArg(j++, rp_buf);
        cryptonets.setArg(j++, srp_buf);
        // cryptonets.setArg(j++, ntt_out_buf);

        std::chrono::high_resolution_clock::time_point p2pStart2 = std::chrono::high_resolution_clock::now();
        
        q.enqueueTask(cryptonets);
        q.finish();
        q.enqueueMigrateMemObjects({ntt_out_buf}, 1); // migrate the obj from (1: device to host)
        q.finish();
    
        std::chrono::high_resolution_clock::time_point p2pEnd2 = std::chrono::high_resolution_clock::now();
        cl_ulong p2pTime2 = std::chrono::duration_cast<std::chrono::microseconds>(p2pEnd2 - p2pStart2).count();
        averagetime2+=p2pTime2;
    }
    
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "INFO: Kernel execution average time: " << averagetime2/double(loop) << " us\n";
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "verify the data transfer method" << std::endl;
    std::cout << static_cast<unsigned long long>(hostBuffer[0]) << " " << static_cast<unsigned long long>(hostBuffer[1])  << " " << static_cast<unsigned long long>(hostBuffer[2]) << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "INFO: P2P transfer average time: " << averagetime/double(loop) << " us\n";
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
 
    std::string ntt_out_filepath = "/media/embedded-415/f29fd339-0d3e-4eb2-ae93-016cd9fcee27/ntt/data_16384/ntt_out.dat";
    unsigned long long* ntt_out = aligned_alloc<unsigned long long>(16384);     // res
    unsigned long long* ntt_out_testbench = aligned_alloc<unsigned long long>(16384); // test_bench res
    load_data(ntt_out_filepath,ntt_out_testbench,16384);
    

    struct timeval start_time, end_time;
    gettimeofday(&start_time, 0);
    std::cout<< "estimated_cryptonets ready"<<std::endl;

    estimated_cryptonets(xclbin_path,ntt_out);
  
    gettimeofday(&end_time, 0);
    std::cout << "User function execution time is: " << tvdiff(&start_time, &end_time) / 1000UL << " ms" << std::endl;

    
    std::cout<<"---test example---"<<std::endl;
    std::cout<<"ntt_out_testbench[5] is "<<ntt_out_testbench[5]<<std::endl;
    std::cout<<"ntt_out[5] is "<<ntt_out_testbench[5]<<std::endl;
    std::cout<<"---test example---"<<std::endl;
    

    unsigned errnum = 0;
    for (int i = 0; i < 16384; i++) {
       if (ntt_out_testbench[i] != ntt_out[i]) {
            errnum++;
        }
    }
	std::cout << "error: " << errnum << std::endl;
    return errnum;
}

