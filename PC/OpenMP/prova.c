#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

void
vector_sum(const double *input1, const double *input2, double *output, int size) {
  int thread_id = -1;
#ifdef _OPENMP
  thread_id = omp_get_thread_num();
#endif

#pragma omp critical
  printf("Hello from thread: %d from inside the function\n", thread_id);
#pragma omp barrier

#pragma omp for
  for (int i = 0; i < size; i++)
    output[i] = input1[i] + input2[i];
}

bool
vector_equal_to(const double val, const double *v, int size) {
    for (int i = 0; i < size; i++) {
      if (v[i] != val)
        return false;
    }
  return true;
}

int
main() {
  const int size = 10;
  double   *input1, *input2, *output;

  const size_t array_size = size * sizeof(double);

  input1 = (double *)malloc(array_size);
  input2 = (double *)malloc(array_size);
  output = (double *)malloc(array_size);

#ifdef _OPENMP
  const double start = omp_get_wtime();
#endif

#pragma omp parallel num_threads(4)
  {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif

      for (int i = 0; i < size; i++) {
        input1[i] = 3;
        input2[i] = 4;
        output[i] = 0;
      }

#pragma omp critical
    printf("Hello from thread: %d from outside the function\n", thread_id);

    vector_sum(input1, input2, output, size);
  }

  if (vector_equal_to(7, output, size))
    printf("Vector sum correct\n");
    else {
      printf("Vector sum not correct\n");
        for (int i = 0; i < size; i++) {
          printf("%f - ", output[i]);
        }
      printf("\n");
    }

#ifdef _OPENMP
  const double end = omp_get_wtime();
  const double dt  = (end - start) * 1000;
  printf("Time elapsed: %f\n", dt);
#endif

  free(input1);
  free(input2);
  free(output);

  return 0;
}