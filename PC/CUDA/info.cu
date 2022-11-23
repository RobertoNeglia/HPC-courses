#include <cuda.h>
#include <iostream>

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    cout << "Device number: " << i << endl;
    cout << "\tDevice name: " << prop.name << endl;
    cout << "\tMax block per multiprocessor: "
         << prop.maxBlocksPerMultiProcessor << endl;
    cout << "\tMax threads per multiprocessor: "
         << prop.maxThreadsPerMultiProcessor << endl;
    cout << "\tMax threads per blocks: " << prop.maxThreadsPerBlock << endl;
    cout << "\tNumber of streaming multiprocessors: "
         << prop.multiProcessorCount << endl;
    cout << "\tNumber of bytes of shared memory per block: "
         << prop.sharedMemPerBlock << endl;
    cout << "\tNumber of bytes of shared memory per multiprocessor: "
         << prop.sharedMemPerMultiprocessor << endl;
    cout << "\tMemory clock rate: " << prop.memoryClockRate << endl;
    cout << "\tMemory bus width: " << prop.memoryBusWidth << endl;
  }
}