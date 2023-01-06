#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <functional>

#define DEBUG 0

auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

/*** The size of the array to be processed */
unsigned long int array_size = 0;

/**
 * Merge two arrays already sorted
 * @param input_data is the array containing the two arrays to be merged
 * @param starting_cell is the index of the first cell of the first array
 * @param size is the sum of the sizes of the two arrays
 * @param output_data is where result has to be stored
 */
void
bottom_up_merge(float            *input_data,
                unsigned long int starting_cell,
                unsigned long int size,
                float            *output_data) {
  if (starting_cell > array_size)
    return;

  /*The last position to be written */
  const unsigned long int last_cell =
    (starting_cell + size <= array_size) ? starting_cell + size : array_size;

  /*The position in the output data to be written */
  unsigned long int index = starting_cell;

  /*The position in the left part of input data */
  unsigned long int left_index = starting_cell;

  /*The position in the right part of input data */
  unsigned long int right_index = starting_cell + size / 2;

  /*The last position in the left part to be read*/
  const unsigned long int last_left = (right_index < last_cell) ? right_index : last_cell;

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
  array_size                    = (unsigned long int)atoi(argv[1]);
  unsigned long int num_threads = (unsigned long int)atoi(argv[2]);

  float        even_data[array_size];
  float        odd_data[array_size];
  const float *final_data;

  int               iteration = 0;
  unsigned long int index     = 0;
  unsigned long int width     = 2;
  unsigned long int sequence_number;

  /* Initialize data in a random way */
    for (index = 0; index < array_size; index++) {
      unsigned int seed = index;
      odd_data[index]   = rand_r(&seed) / (double)RAND_MAX;
    }

    if (DEBUG) {
      printf("Unordered array: ");
        for (index = 0; index < array_size; index++) {
          printf("%f ", odd_data[index]);
        }
      printf("\n");
  }

  const auto dt = timeit([&]() {
    while (width / 2 < array_size) {
      sequence_number = array_size / width + (array_size % width != 0 ? 1 : 0);
        for (unsigned long int sequence = 0; sequence < sequence_number; sequence++) {
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

  final_data = iteration % 2 == 0 ? odd_data : even_data;
  });

  printf("Elapsed time: %ld[ms]\n", dt);

    if (DEBUG) {
      printf("Ordered array: ");
        /* Print the final result */
        for (index = 0; index < array_size; index++) {
          printf("%f ", final_data[index]);
        }
      printf("\n");
  }
  return 0;
}
