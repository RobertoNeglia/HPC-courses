#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <regex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include "horner.hpp"

int main(int argc, char** argv) {
    bool verbose = 0;

    // int degree;
    // std::cout << "Polynomial degree: " << std::endl;
    // std::cout << "=> ";
    // std::cin >> degree;

    const auto param = parse_parameters(get_file_contents("params.dat"));

    const auto degree = static_cast<unsigned int>(param.at("degree"));
    const auto n = static_cast<unsigned int>(param.at("n_points"));
    const auto x_0 = param.at("x_0");
    const auto x_f = param.at("x_f");

    std::cout << "Coefficients are computed automatically (with 2sin(2k) function)" << std::endl;
    if (verbose)
        std::cout << "Equation: " << std::endl;
    std::vector<double> coeff(degree + 1);
    for (int i = 0; i <= degree; i++) {
        coeff[i] = 2 * std::sin(2 * i);
        if (verbose) {
            if (i == 0)
                std::cout << coeff[i];
            else
                std::cout << " + (" << coeff[i] << ") * x^" << i;
        }
    }

    std::cout << std::endl;


    // // static parameters
    // const double x_0 = 0.0;
    // const double x_f = 2.0;
    // const unsigned int n = 1000000;

    // compute evaluation points and put them in a vector
    std::vector<double> points(n);
    const double h = (x_f - x_0) / n;

    for (int i = 0; i < n; i++)
        points[i] = x_0 + i * h;

    std::vector<double> naive_results(n);
    std::vector<double> recursive_horner_results(n);
    std::vector<double> iterative_horner_results(n);

    if (verbose) {

        std::cout << "Computing " << n << " evaluations of polynomial with standard formula" << std::endl;
        {
            // we use the namespace in a limitedd scope thanks to the curly brackets
            using namespace std::chrono;
            // get the current time
            const auto start = high_resolution_clock::now();
            naive_results = evaluate_poly(points, coeff, eval_x);
            const auto end = high_resolution_clock::now();

            const auto diff = duration_cast<milliseconds>(end - start).count();

            std::cout << "Time elapsed: " << diff << " [ms]" << std::endl;
        }

        std::cout << "Computing " << n << " evaluations of polynomial with horner recursive formula" << std::endl;
        {
            // we use the namespace in a limitedd scope thanks to the curly brackets
            using namespace std::chrono;
            // get the current time
            const auto start = high_resolution_clock::now();
            recursive_horner_results = evaluate_poly(points, coeff, eval_horner_x_recursive);
            const auto end = high_resolution_clock::now();

            const auto diff = duration_cast<milliseconds>(end - start).count();

            std::cout << "Time elapsed: " << diff << " [ms]" << std::endl;
        }

        std::cout << "Computing " << n << " evaluations of polynomial with horner iterative formula" << std::endl;
        {
            // we use the namespace in a limitedd scope thanks to the curly brackets
            using namespace std::chrono;
            // get the current time
            const auto start = high_resolution_clock::now();
            iterative_horner_results = evaluate_poly(points, coeff, eval_horner_x_iterative);
            const auto end = high_resolution_clock::now();

            const auto diff = duration_cast<milliseconds>(end - start).count();

            std::cout << "Time elapsed: " << diff << " [ms]" << std::endl;
        }
    }
    else {
        std::cout << "Computing " << n << " evaluations of polynomial with standard formula" << std::endl;
        std::cout << "Time elapsed: " << timeit([&]() {
            evaluate_poly(points, coeff, eval_x);
            }) << "[ms]" << std::endl;

        std::cout << "Computing " << n << " evaluations of polynomial with horner recursive formula" << std::endl;
        std::cout << "Time elapsed: " << timeit([&]() {
            evaluate_poly(points, coeff, eval_horner_x_recursive);
            }) << "[ms]" << std::endl;

        std::cout << "Computing " << n << " evaluations of polynomial with horner iterative formula" << std::endl;
        std::cout << "Time elapsed: " << timeit([&]() {
            evaluate_poly(points, coeff, eval_horner_x_iterative);
            }) << "[ms]" << std::endl;
    }

    if (verbose) {
        std::cout << "Coefficients: " << std::endl;
        print_vector(coeff);

        std::cout << "Points: " << std::endl;
        print_vector(points);

        std::cout << "Naive results: " << std::endl;
        print_vector(naive_results);
        std::cout << "Horner (recursive) results: " << std::endl;
        print_vector(recursive_horner_results);
        std::cout << "Horner (iterative) results: " << std::endl;
        print_vector(iterative_horner_results);

        std::cout << "Diff vector norm: " << vector_norm(abs_diff_vector(naive_results, recursive_horner_results)) << std::endl;


    }
    std::cout << "Sanity check: " << std::endl << "naive results == Horner (recursive)? " << compare_vectors(naive_results, recursive_horner_results) << std::endl;
    std::cout << "Sanity check: " << std::endl << "naive results == Horner (iterative)? " << compare_vectors(naive_results, iterative_horner_results) << std::endl;
}