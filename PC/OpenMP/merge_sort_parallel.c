#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG 0

/*** The size of the array to be processed */
unsigned int array_size = 0;

/**
 * Merge two arrays already sorted
 * @param input_data is the array containing the two arrays to be merged
 * @param starting_cell is the index of the first cell of the first array
 * @param size is the sum of the sizes of the two arrays
 * @param output_data is where result has to be stored
 */
void
bottom_up_merge(float *input_data, int starting_cell, int size, float *output_data) {
  if (starting_cell > array_size)
    return;

  /*The last position to be written */
  const int last_cell =
    (starting_cell + size <= array_size) ? starting_cell + size : array_size;

  /*The position in the output data to be written */
  int index = starting_cell;

  /*The position in the left part of input data */
  int left_index = starting_cell;

  /*The position in the right part of input data */
  int right_index = starting_cell + size / 2;

  /*The last position in the left part to be read*/
  const int last_left = (right_index < last_cell) ? right_index : last_cell;

    for (index = starting_cell; index < last_cell; index++) {
        if (left_index < last_left &&
            (right_index >= last_cell ||
             input_data[left_index] <= input_data[right_index])) {
          output_data[index] = input_data[left_index];
          left_index++;
        } else {
          output_data[index] = input_data[right_index];
          right_index++;
        }
    }
}

bool
check_sort(const float *data, int size) {
    if (size == 0 || size == 1) {
      return true;
  }

    for (int i = 1; i < size; i++) {
      if (data[i - 1] > data[i])
        return false;
    }
  return true;
}

int
main(int argc, char **argv) {
    if (argc != 3) {
      printf("Wrong number of parameters\n");
      return 0;
  }
  array_size               = (unsigned int)atoi(argv[1]);
  unsigned int num_threads = (unsigned int)atoi(argv[2]);

  size_t size = array_size * sizeof(float);

  float *even_data = (float *)malloc(size);
  float *odd_data  = (float *)malloc(size);

  int          iteration = 0;
  unsigned int index     = 0;
  unsigned int width     = 2;
  unsigned int sequence_number;

  /* Initialize data in a random way */
    for (index = 0; index < array_size; index++) {
      unsigned int seed                = (index + 1) * 2;
      odd_data[array_size - index - 1] = rand_r(&seed) / (double)RAND_MAX;
    }

    if (DEBUG) {
      /* Print the initial array */
      printf("unordered array: \n");
        for (index = 0; index < array_size; index++) {
          printf("%f ", odd_data[index]);
        }
      printf("\n");
  }

#ifdef _OPENMP
  const double start = omp_get_wtime();
#else
  const clock_t start = clock();
#endif

    while (width / 2 < array_size) {
      sequence_number = array_size / width + (array_size % width != 0 ? 1 : 0);
#pragma omp parallel for num_threads(num_threads)
        for (unsigned int sequence = 0; sequence < sequence_number; sequence++) {
            /* Even iteration: the result is stored in even_data */
            if (iteration % 2 == 0) {
              bottom_up_merge(odd_data, sequence * width, width, even_data);
            } else {
              bottom_up_merge(even_data, sequence * width, width, odd_data);
            }
        }
      iteration++;
      width = width * 2;
    }

#ifdef _OPENMP
  const double end = omp_get_wtime();
  const double dt  = (end - start) * 1000;
#else
  const clock_t end   = clock();
  const double  dt    = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
#endif

  printf("Time elapsed: %f[ms]\n", dt);

  const float *final_data = iteration % 2 == 0 ? odd_data : even_data;

    if (check_sort(final_data, array_size)) {
      printf("Array has been correctly sorted\n");
    } else {
      printf("Array has not been correctly sorted\n");
    }

    if (DEBUG) {
      /* Print the final result */
      printf("ordered array: \n");
        for (index = 0; index < array_size; index++) {
          printf("%f ", final_data[index]);
        }
      printf("\n");
  }

  free(odd_data);
  free(even_data);

  return 0;
}
