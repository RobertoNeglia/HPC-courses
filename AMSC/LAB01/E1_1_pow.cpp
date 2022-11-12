#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

template<typename T>
using pow_method_t = T (*)(T, unsigned int);

template<typename T>
T pow_iterative(T base, unsigned int exp) {
    T res = 1;
    
    while(exp > 0) {
        if(exp & 1)
            res *= base;
        base *= base;
        exp >>= 1;
    }

    return res;
}

template<typename T>
T pow_recursive(T base, unsigned int exp) {
    if(exp == 0)
        return 1;
    
    if(exp == 1)
        return base;

    if(exp % 2)
        return base * pow_recursive(base * base, (exp - 1) / 2);
    else
        return pow_recursive(base * base, exp / 2);
}

double rel_error(double a, double b) {
    if(a == 0 && b == 0)
        return 0;
    return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

void test_error(double err, double toll) {
    std::cout << err << " " << (err < toll ? "PASSED" : "FAILED") << std::endl;
}

int main(int argc, char** argv) {
    const std::vector<unsigned int> exponents = {0, 1, 2, 3, 5, 8, 13};
    const std::vector<unsigned long long int> base_uint = {0, 1, 2, 3, 5, 8, 13};
    const std::vector<double> base_double = { 0, 1, 2.3, 3.7, 5.4, 7.9, 8.1 };

    std::cout << "---- unsigned int tests ----------------------------------------" << std::endl;
    for(const auto& e : exponents) {
        for(const auto& b : base_uint) {
            std::cout << pow_recursive(b, e) - std::pow(b, e) << std::endl;
            std::cout << pow_iterative(b, e) - std::pow(b, e) << std::endl;
        }
    }

    std::cout << "---- double tests ----------------------------------------------" << std::endl;
    constexpr auto TOLL = std::numeric_limits<double>::epsilon() * 10;
    for(const auto& e : exponents) {
        for(const auto& b : base_double) {
            test_error(rel_error(pow_recursive(b, e), std::pow(b, e)), TOLL);
            test_error(rel_error(pow_iterative(b, e), std::pow(b, e)), TOLL);
        }
    }

    return 0;
}