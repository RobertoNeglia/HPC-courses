#include <cuda.h>
#include <iostream>

#define MATRIX_SIZE 8192
#define MAX_BLOCK_SIZE 32
#define MASK_SIZE 5

__global__ void convolution_2D_basic_kernel(unsigned char *in,
                                            unsigned char *mask,
                                            unsigned char *out, int maskwidth,
                                            int w, int h) {
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;

  if (col < w && row < h) {
    int pix_val = 0;

    int n_start_col = col - maskwidth / 2;
    int n_start_row = row - maskwidth / 2;

    // Get surrounding block
    for (int j = 0; j < maskwidth; j++) {
      for (int k = 0; k < maskwidth; k++) {
        int cur_row = n_start_row + j;
        int cur_col = n_start_col + k;

        // Verify the image pixel is valid
        if (cur_row > -1 && cur_row < h && cur_col > -1 && cur_col < w)
          pix_val += in[cur_row * w + cur_col] * mask[j * maskwidth + k];
      }
    }
    // Write new pixel out
    out[row * w + col] = (unsigned char)(pix_val);
  }
}

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;

  cout << "Size of matrix: " << MATRIX_SIZE << "x" << MATRIX_SIZE << endl;
  cout << "Size of mask: " << MASK_SIZE << "x" << MASK_SIZE << endl;

  for (int block_size = 8; block_size <= MAX_BLOCK_SIZE; block_size *= 2) {
    unsigned char *n, *m, *p;

    cudaMallocManaged((void **)&n,
                      sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);
    cudaMallocManaged((void **)&m,
                      sizeof(unsigned char) * MASK_SIZE * MASK_SIZE);
    cudaMallocManaged((void **)&p,
                      sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);

    int add;
    int val = 1;
    // Instantiation of the matrix
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        add = i * MATRIX_SIZE + j;
        n[add] = 2 * val;
      }
      val++;
    }

    val = 1;
    // Instantiation of the mask
    for (int i = 0; i < MASK_SIZE; i++) {
      for (int j = 0; j < MASK_SIZE; j++) {
        add = i * MASK_SIZE + j;
        m[add] = val;
        if (j < MASK_SIZE / 2)
          val++;
        else
          val--;
      }
      if (i < MASK_SIZE / 2)
        val += 2;
    }

    // NO TILING KERNEL EXECUTION

    // static parameters to calculate execution time
    float no_tiling_gpu_elapsed_time_ms;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    unsigned int grid_rows = (MATRIX_SIZE + block_size - 1) / block_size;
    unsigned int grid_cols = (MATRIX_SIZE + block_size - 1) / block_size;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(block_size, block_size);

    cout << "Launching a kernel on the GPU with " << grid_rows << " rows and "
         << grid_cols << " columns, and " << block_size << "x" << block_size
         << " threads per block" << endl;

    cudaEventRecord(start, 0);
    convolution_2D_basic_kernel<<<dimGrid, dimBlock>>>(
        n, m, p, MASK_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&no_tiling_gpu_elapsed_time_ms, start, stop);

    cout << "Time elapsed without tiling: " << no_tiling_gpu_elapsed_time_ms
         << "ms" << endl;

    cudaFree(n);
    cudaFree(m);
    cudaFree(p);
  }
  return 0;
}