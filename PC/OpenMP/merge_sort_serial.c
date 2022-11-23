#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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

int
main(int argc, char **argv) {
    if (argc != 3) {
      printf("Wrong number of parameters\n");
      return 0;
  }
  array_size               = (unsigned int)atoi(argv[1]);
  unsigned int num_threads = (unsigned int)atoi(argv[2]);

  float even_data[array_size];
  float odd_data[array_size];

  int          iteration = 0;
  int          index     = 0;
  unsigned int width     = 2;
  unsigned int sequence_number;

    /* Initialize data in a random way */
    for (index = 0; index < array_size; index++) {
      int seed        = index;
      odd_data[index] = rand_r(&seed) / (double)RAND_MAX;
    }

    while (width / 2 < array_size) {
      sequence_number = array_size / width + (array_size % width != 0 ? 1 : 0);
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

  const float *final_data = iteration % 2 == 0 ? odd_data : even_data;

    /* Print the final result */
    for (index = 0; index < array_size; index++) {
      printf("%f ", final_data[index]);
    }
  printf("\n");
  return 0;
}
