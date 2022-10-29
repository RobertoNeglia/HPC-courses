#include <cuda.h>
#include <iostream>

#define VERBOSE 1
#define MATRIX_SIZE 16
#define BLOCK_SIZE (TILE_WIDTH)
#define MASK_SIZE 5
#define O_TILE_WIDTH 4
#define TILE_WIDTH (O_TILE_WIDTH + MASK_SIZE - 1)

__global__ void convolution_2D_basic_kernel(unsigned char *in,
                                            unsigned char *mask,
                                            unsigned char *out, int maskwidth,
                                            int w, int h) {

  // each thread computes a single cell in the final matrix
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

__global__ void convolution_2D_tiling_kernel(unsigned char *in,
                                             unsigned char *mask,
                                             unsigned char *out, int maskwidth,
                                             int w, int h) {
  // Create shared memory block tile
  __shared__ unsigned char *d_in[TILE_WIDTH][TILE_WIDTH];

  int bx = blockIdx.x, by = blockIdx.y;
  int tx = threadIdx.x, ty = threadIdx.y;

  int col = bx * blockDim.x + tx;
  int row = by * blockDim.y + ty;

  int pix_val = 0;

  // Loop over input tiles required to compute the output element
  for (int o = 0; o < (MATRIX_SIZE - 1) / TILE_WIDTH + 1; o++) {
    // Collaborative loading of input tiles into shared memory
    if (o * TILE_WIDTH + ty < MATRIX_SIZE && o * TILE_WIDTH + tx < MATRIX_SIZE)
      d_in[ty][tx] = in[];

    __syncthreads();

    if (row < MATRIX_SIZE && col < MATRIX_SIZE) {
      for (int i = 0; i < MASK_SIZE; i++)
    }
  }

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

  unsigned char *n, *m, *p;

  cudaMallocManaged((void **)&n,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);
  cudaMallocManaged((void **)&m, sizeof(unsigned char) * MASK_SIZE * MASK_SIZE);
  cudaMallocManaged((void **)&p,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);

  int val = 1;
  int add;
  for (int i = 0; i < MATRIX_SIZE; i++) {
    for (int j = 0; j < MATRIX_SIZE; j++) {
      add = i * MATRIX_SIZE + j;
      n[add] = val;
      if (j < MATRIX_SIZE / 2)
        val++;
      else
        val--;
    }
    if (i < MATRIX_SIZE / 2)
      val += 2;
  }

  if (VERBOSE) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)n[i * MATRIX_SIZE + j] << "\t";
      }
      cout << endl;
    }
  }

  val = 1;
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

  if (VERBOSE) {
    for (int i = 0; i < MASK_SIZE; i++) {
      for (int j = 0; j < MASK_SIZE; j++) {
        cout << (int)m[i * MASK_SIZE + j] << "\t";
      }
      cout << endl;
    }
  }

  // static parameters to calculate execution time
  float no_tiling_gpu_elapsed_time_ms;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  unsigned int grid_rows = (MATRIX_SIZE + BLOCK_SIZE - 1) / BLOCK_SIZE;
  unsigned int grid_cols = (MATRIX_SIZE + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 dimGrid(grid_cols, grid_rows);
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

  cudaEventRecord(start, 0);
  convolution_2D_basic_kernel<<<dimGrid, dimBlock>>>(n, m, p, MASK_SIZE,
                                                     MATRIX_SIZE, MATRIX_SIZE);

  cudaDeviceSynchronize();
  cudaEventRecord(stop, 0);

  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&no_tiling_gpu_elapsed_time_ms, start, stop);

  if (VERBOSE) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)p[i * MATRIX_SIZE + j] << "\t";
      }
      cout << endl;
    }
  }

  cout << "Time elapsed without tiling: " << no_tiling_gpu_elapsed_time_ms
       << endl;

  cudaFree(n);
  cudaFree(m);
  cudaFree(p);
  return 0;
}