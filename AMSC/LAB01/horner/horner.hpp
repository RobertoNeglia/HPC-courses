#ifndef HORNER_H
#define HORNER_H

using eval_method_t = std::function<double(const std::vector<double>&, const double&)>;
// using eval_method_t = double (*)(const std::vector<double>&, const double&);

std::string get_file_contents(const char* filename) {
    std::ifstream in(filename, std::ios::in);
    if (in) {
        std::string content;
        in.seekg(0, std::ios::end);
        content.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&content[0], content.size());
        in.close();
        return content;
    }
    throw(errno);
}

std::unordered_map<std::string, double> parse_parameters(const std::string& s) {
    const std::regex parameter_pattern(R"((\w+)=(-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?)(\n|\r\n)?)");
    std::cmatch base_match;
    std::unordered_map<std::string, double> parameters;
    auto cs = s.c_str();
    while (std::regex_search(cs, base_match, parameter_pattern)) {
        cs += base_match.length();
        // convert int to double is safe until 2^53-1
        parameters.emplace(base_match[1], std::stod(base_match[2]));
    }

    return parameters;
}

auto timeit(const std::function<void()>& f) {
    using namespace std::chrono;
    const auto start = high_resolution_clock::now();
    f();
    const auto end = high_resolution_clock::now();
    return duration_cast<milliseconds>(end - start).count();
}

double pow_iterative(double base, unsigned int exp) {
    double r = 1;

    while (exp > 0) {
        if (exp & 1)
            r *= base;
        base *= base;
        exp >>= 1;
    }

    return r;
}

std::vector<double> evaluate_poly(const std::vector<double>& points, const std::vector<double>& coeff, eval_method_t method) {
    std::vector<double> results;
    std::size_t n = points.size();
    results.reserve(n);
    for (std::size_t i = 0; i < n; i++)
        results.emplace_back(method(coeff, points[i]));
    return results;
}

double eval_x(const std::vector<double>& coeff, const double& x) {
    double r = 0;

    for (int i = 0; i < coeff.size(); i++) {
        r += coeff[i] * pow_iterative(x, i);
    }

    return r;
}

double horner_recursive(const std::vector<double>& coeff, const double& x, int i) {
    if (i == coeff.size())
        return coeff[i - 1];

    return coeff[i - 1] + x * horner_recursive(coeff, x, i + 1);
}

double eval_horner_x_recursive(const std::vector<double>& coeff, const double& x) {
    double r = 0;
    r = horner_recursive(coeff, x, 1);
    return r;
}

double eval_horner_x_iterative(const std::vector<double>& coeff, const double& x) {
    double result = coeff.back();
    for (auto i = coeff.crbegin() + 1; i != coeff.crend(); i++) {
        result = x * result + (*i);
    }

    return result;
}

// double eval_horner_x_iterative(const std::vector<double>& coeff, const double& x) {
//     int deg = coeff.size() - 1;
//     double r = coeff[deg];
//     for(int i = deg; i > 0; i--) {
//         r = coeff[i - 1] + x * r;
//     }

//     return r;
// }

std::vector<double> abs_diff_vector(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); i++)
        c[i] = std::abs(a[i] - b[i]);
    return c;
}

void print_vector(const std::vector<double>& a) {
    for (const auto& i : a) std::cout << i << " - ";
    std::cout << std::endl;
}

double vector_norm(const std::vector<double>& a) {
    double n = 0;
    for (int i = 0; i < a.size(); i++)
        n += pow_iterative(a[i], 2);
    return std::sqrt(n);
}

bool compare_vectors(const std::vector<double>& a, const std::vector<double>& b) {
    return vector_norm(abs_diff_vector(a, b)) < 1.e-6;
}

#endif // HORNER_H