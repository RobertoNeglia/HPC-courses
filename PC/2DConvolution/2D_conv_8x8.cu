#include <cuda.h>
#include <iostream>

#define DEBUG 0
#define MATRIX_SIZE 8192
#define BLOCK_SIZE 8
#define MASK_SIZE 5
#define O_TILE_SIZE (BLOCK_SIZE - MASK_SIZE + 1)

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

__global__ void convolution_2D_tiling_kernel(unsigned char *in,
                                             unsigned char *mask,
                                             unsigned char *out, int maskwidth,
                                             int w, int h) {
  // Shared memory block tile
  __shared__ unsigned char d_in[BLOCK_SIZE][BLOCK_SIZE];

  int bx = blockIdx.x, by = blockIdx.y;
  int tx = threadIdx.x, ty = threadIdx.y;

  // Calculate the offset from the element position in the output tile to the
  // element position in the input tile
  int offset = maskwidth / 2;

  // Position of the thread in the GPU grid, output tile POV
  int o_row = by * O_TILE_SIZE + ty;
  int o_col = bx * O_TILE_SIZE + tx;

  // Position of the thread in the GPU grid, input tile POV
  int i_row = o_row - offset;
  int i_col = o_col - offset;

  // Accumulation variable
  int output;

  if ((i_row > -1 && i_row < h) && (i_col > -1 && i_col < w))
    d_in[ty][tx] = in[i_row * w + i_col];
  else
    d_in[ty][tx] = 0;
  // Wait for all threads in the block to load the tile
  __syncthreads();

  // Compute the ouput element for each thread
  if (ty < O_TILE_SIZE && tx < O_TILE_SIZE) {
    output = 0;
    for (int i = 0; i < maskwidth; i++) {
      for (int j = 0; j < maskwidth; j++) {
        output += mask[i * maskwidth + j] * d_in[ty + i][tx + j];
      }
    }

    if (o_row < h && o_col < w)
      out[o_row * w + o_col] = output;
  }
}

__global__ void matrix_diff(unsigned char *m, unsigned char *n,
                            unsigned char *diff, int w, int h) {
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;

  if (row < h && col < w)
    diff[row * w + col] = m[row * w + col] - n[row * w + col];
}

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;
  unsigned char *n, *m, *p_no_tiling, *p_tiling, *diff;

  cudaMallocManaged((void **)&n,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);
  cudaMallocManaged((void **)&m, sizeof(unsigned char) * MASK_SIZE * MASK_SIZE);
  cudaMallocManaged((void **)&p_no_tiling,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);
  cudaMallocManaged((void **)&p_tiling,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);
  cudaMallocManaged((void **)&diff,
                    sizeof(unsigned char) * MATRIX_SIZE * MATRIX_SIZE);

  int add;
  int val = 1;
  for (int i = 0; i < MATRIX_SIZE; i++) {
    for (int j = 0; j < MATRIX_SIZE; j++) {
      add = i * MATRIX_SIZE + j;
      n[add] = 2 * val;
    }
    val++;
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

  // NO TILING KERNEL EXECUTION

  // static parameters to calculate execution time
  float no_tiling_gpu_elapsed_time_ms;
  cudaEvent_t start_no_tiling, stop_no_tiling;
  cudaEventCreate(&start_no_tiling);
  cudaEventCreate(&stop_no_tiling);

  unsigned int grid_rows = (MATRIX_SIZE + O_TILE_SIZE - 1) / O_TILE_SIZE;
  unsigned int grid_cols = (MATRIX_SIZE + O_TILE_SIZE - 1) / O_TILE_SIZE;
  dim3 dimGrid(grid_cols, grid_rows);
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

  cout << "Size of matrix: " << MATRIX_SIZE << "x" << MATRIX_SIZE << endl;
  cout << "Size of mask: " << MASK_SIZE << "x" << MASK_SIZE << endl;

  cout << "Launching a kernel on the GPU with " << grid_rows << " rows and "
       << grid_cols << " columns, and " << BLOCK_SIZE << "x" << BLOCK_SIZE
       << " threads per block" << endl;

  cudaEventRecord(start_no_tiling, 0);
  convolution_2D_basic_kernel<<<dimGrid, dimBlock>>>(
      n, m, p_no_tiling, MASK_SIZE, MATRIX_SIZE, MATRIX_SIZE);
  cudaDeviceSynchronize();
  cudaEventRecord(stop_no_tiling, 0);
  cudaEventSynchronize(stop_no_tiling);

  cudaEventElapsedTime(&no_tiling_gpu_elapsed_time_ms, start_no_tiling,
                       stop_no_tiling);

  cout << "Time elapsed without tiling: " << no_tiling_gpu_elapsed_time_ms
       << "ms" << endl;

  // TILING KERNEL EXECUTION

  // static parameters to calculate execution time
  float tiling_gpu_elapsed_time_ms;
  cudaEvent_t start_tiling, stop_tiling;
  cudaEventCreate(&start_tiling);
  cudaEventCreate(&stop_tiling);

  cudaEventRecord(start_tiling, 0);
  convolution_2D_tiling_kernel<<<dimGrid, dimBlock>>>(n, m, p_tiling, MASK_SIZE,
                                                      MATRIX_SIZE, MATRIX_SIZE);
  cudaDeviceSynchronize();
  cudaEventRecord(stop_tiling, 0);
  cudaEventSynchronize(stop_tiling);

  cudaEventElapsedTime(&tiling_gpu_elapsed_time_ms, start_tiling, stop_tiling);
  cout << "Time elapsed with tiling: " << tiling_gpu_elapsed_time_ms << "ms"
       << endl;

  matrix_diff<<<dimGrid, dimBlock>>>(p_no_tiling, p_tiling, diff, MATRIX_SIZE,
                                     MATRIX_SIZE);

  bool equal = true;
  for (int i = 0; i < MATRIX_SIZE; i++) {
    for (int j = 0; j < MATRIX_SIZE; j++) {
      if ((int)diff[i * MATRIX_SIZE + j] != 0) {
        equal = false;
        break;
      }
    }
  }

  if (equal)
    cout << "Convolution matrices with tiling and without tiling are equal"
         << endl;
  else
    cout << "Convolution matrices with tiling and without tiling are NOT equal"
         << endl;

  // DEBUG
  if (DEBUG) {
    cout << "---- N -----------" << endl;
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)n[i * MATRIX_SIZE + j] << "--";
      }
      cout << endl;
    }

    cout << "---- TILING -----------" << endl;
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)p_tiling[i * MATRIX_SIZE + j] << "--";
      }
      cout << "__" << endl;
    }

    cout << "---- NO TILING --------" << endl;
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)p_no_tiling[i * MATRIX_SIZE + j] << "--";
      }
      cout << endl;
    }

    cout << "---- DIFF --------" << endl;
    for (int i = 0; i < MATRIX_SIZE; i++) {
      for (int j = 0; j < MATRIX_SIZE; j++) {
        cout << (int)diff[i * MATRIX_SIZE + j] << "    ";
      }
      cout << endl;
    }
  }

  cudaFree(n);
  cudaFree(m);
  cudaFree(p_no_tiling);
  cudaFree(p_tiling);
  cudaFree(diff);
  return 0;
}