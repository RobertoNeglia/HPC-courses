#include <mpi.h>

#include <iostream>
#include <string>

int
main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int size, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string who = (argc > 1) ? argv[1] : "world";
  std::string message =
    "Hello, " + who + ", from rank " + std::to_string(rank) + " of " + std::to_string(size);

  // messaging management

    if (rank) { // I'm not the master
      // I have to send my message to him
      const unsigned length = message.size();
      MPI_Send(&length, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&message[0], length, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    } else

    { // I'm the master
      // first print my message
      std::cout << message << std::endl;
        // then print the other message
        for (int i = 1; i < size; i++) {
          unsigned    length;
          std::string message;
          MPI_Recv(&length, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          message.resize(length);
          MPI_Recv(&message[0], length, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          std::cout << message << std::endl;
        }
    }

  MPI_Finalize();
}