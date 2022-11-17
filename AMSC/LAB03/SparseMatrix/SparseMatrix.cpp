#include <vector>
#include <iostream>
#include <cassert>
#include <map>
#include <unordered_map>
#include <chrono>

using elem_t = double;
using Vector = std::vector<elem_t>;

// Step 1: abstract matrix class
class SparseMatrix { // abstract class
public:
    SparseMatrix() : rows(0), cols(0), nnz(0) {}
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    size_t getNnz() const { return nnz; }

    void print(std::ostream& os = std::cout) const {
        os << "Numbero of rows: " << rows << std::endl;
        os << "Number of columns: " << cols << std::endl;
        os << "Number of non-zero elements: " << nnz << std::endl;
        _print(os);
    }
    virtual Vector vmult(const Vector& v) const = 0;
    // element access in read only (const version)
    virtual const elem_t& operator()(size_t i, size_t j) const = 0;
    // element access in write (returns non-const reference)
    virtual elem_t& operator()(size_t i, size_t j) = 0;
    // virtual destructor
    virtual ~SparseMatrix() = default;

protected: // protected because need to be accessible to children!
    size_t rows;
    size_t cols;
    size_t nnz;
    virtual void _print(std::ostream& os) const = 0;
};

class MapMatrix : public SparseMatrix {
public:
    virtual Vector vmult(const Vector& v) const override {
        // checks that matrix-vector multiplication col-row dimension are compatible
        assert(v.size() == cols);
        Vector result(v.size());

        // for every row i...

        for (size_t i = 0; i < data.size(); i++) {
            // ...and for every column element key...
            for (const auto [key, value] : data[i]) {
                // ... multiply the value in the vector with the value in the matrix
                // !!RMK¡¡: only elements that are present in the matrix are considered, others are zero (not how the key value is used also as input vector index)
                result[i] += v[key] * value;
            }
        }
        return result;
    }

    // element access in read only (const version)
    virtual const elem_t& operator()(size_t i, size_t j) const override {
        return data.at(i).at(j);
    }

    // element access in write (returns non-const reference)
    virtual elem_t& operator()(size_t i, size_t j) override {
        // Row boundary check
        if (i + 1 > data.size()) { // must use this notation and not i > data.size() - 1 because it returns a size_t, that cant be negative
            // if accessed element is outside row boundaries, increment the number of rows
            data.resize(i + 1);
            rows = i + 1;
        }
        // Now i have to access the column
        // find the element in the map at position j
        const auto it = data[i].find(j);
        if (it == data[i].end()) // if the element in (i, j) is not present
        {
            // update number of columns, if it's outside the column boundaries
            cols = std::max(cols, j + 1);
            // increment number of non-zero elements
            nnz++;

            // find the map in position i and insert the value mapped through key j.
            // emplace returns pair of an iterator to the newly created pair in the map, and a bool that denotes if the operation was successfull:
            // so i take the pair pointed by the iterator (so (*data|...|.first) and the return the reference to the mapped value through the key)
            return (*data.at(i).emplace(j, 0).first).second;
        }
        // if the element is present, just return the reference to the mapped value through the key
        return (*it).second;
    }

    void printDataSize(std::ostream& os = std::cout) {
        os << data.size() << std::endl;
    }
    virtual ~MapMatrix() override = default;

protected:
    virtual void _print(std::ostream& os) const override {
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto [key, value] : data[i]) {
                os << "(" << i << ", " << key << ") -> " << value << std::endl;
            }
        }
    }

private:
    std::vector<std::map<size_t, elem_t>> data;
};

class UnMapMatrix : public SparseMatrix {
public:
    virtual Vector vmult(const Vector& v) const override {
        // checks that matrix-vector multiplication col-row dimension are compatible
        assert(v.size() == cols);
        Vector result(v.size());

        // for every row i...

        for (size_t i = 0; i < data.size(); i++) {
            // ...and for every column element key...
            for (const auto [key, value] : data[i]) {
                // ... multiply the value in the vector with the value in the matrix
                // !!RMK¡¡: only elements that are present in the matrix are considered, others are zero (not how the key value is used also as input vector index)
                result[i] += v[key] * value;
            }
        }
        return result;
    }

    // element access in read only (const version)
    virtual const elem_t& operator()(size_t i, size_t j) const override {
        return data.at(i).at(j);
    }

    // element access in write (returns non-const reference)
    virtual elem_t& operator()(size_t i, size_t j) override {
        // Row boundary check
        if (i + 1 > data.size()) { // must use this notation and not i > data.size() - 1 because it returns a size_t, that cant be negative
            // if accessed element is outside row boundaries, increment the number of rows
            data.resize(i + 1);
            rows = i + 1;
        }
        // Now i have to access the column
        // find the element in the map at position j
        const auto it = data[i].find(j);
        if (it == data[i].end()) // if the element in (i, j) is not present
        {
            // update number of columns, if it's outside the column boundaries
            cols = std::max(cols, j + 1);
            // increment number of non-zero elements
            nnz++;

            // find the map in position i and insert the value mapped through key j.
            // emplace returns pair of an iterator to the newly created pair in the map, and a bool that denotes if the operation was successfull:
            // so i take the pair pointed by the iterator (so (*data|...|.first) and the return the reference to the mapped value through the key)
            return (*data.at(i).emplace(j, 0).first).second;
        }
        // if the element is present, just return the reference to the mapped value through the key
        return (*it).second;
    }

    void printDataSize(std::ostream& os = std::cout) {
        os << data.size() << std::endl;
    }
    virtual ~UnMapMatrix() override = default;

protected:
    virtual void _print(std::ostream& os) const override {
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto [key, value] : data[i]) {
                os << "(" << i << ", " << key << ") -> " << value << std::endl;
            }
        }
    }
private:
    std::vector<std::unordered_map<size_t, elem_t>> data;
};

void fill_matrix(SparseMatrix& m, size_t n) {
    for (size_t i = 0; i < n; i++) {
        m(i, i) = -2;
        if (i > 0) m(i, i - 1) = 1;
        if (i < n - 1) m(i, i + 1) = 1;
    }
}

template <typename T>
bool eq(const T& lhs, const T& rhs) {
    if (lhs.size() != rhs.size())
        return false;
    for (auto i = lhs.begin(), j = rhs.begin(); i != lhs.end() && j != rhs.end(); i++, j++) {
        if (*i != *j)
            return false;
    }
    return true;
}

void printVector(const Vector& v) {
    std::cout << "| ";
    for (size_t i = 0; i < v.size(); i++) {
        std::cout << v[i] << " | ";
    }
    std::cout << std::endl;
}

int main() {
    size_t n = 1000;
    MapMatrix m;
    UnMapMatrix um;
    Vector v(n), e_res(n), res(n), ures(n);


    {
        // elapsed time for ordered map
        using namespace std::chrono;

        const auto start = high_resolution_clock::now();
        fill_matrix(m, n);
        const auto end = high_resolution_clock::now();

        const auto diff = duration_cast<milliseconds>(end - start).count();

        std::cout << "Time elapsed to populate ordered matrix: " << diff << std::endl;
    }
    {
        // elapsed time for ordered map
        using namespace std::chrono;

        const auto start = high_resolution_clock::now();
        fill_matrix(um, n);
        const auto end = high_resolution_clock::now();

        const auto diff = duration_cast<milliseconds>(end - start).count();

        std::cout << "Time elapsed to populate unordered matrix: " << diff << std::endl;
    }
    for (size_t i = 0; i < n; i++) {
        v[i] = i;
    }
    e_res[0] = 1;
    e_res[n - 1] = -static_cast<int>(n);
    for (size_t i = 1; i < n - 1; i++) {
        e_res[i] = 0;
    }

    // m.printDataSize();
    // m.print();
    // um.print();

    {
        // elapsed time for ordered map
        using namespace std::chrono;

        const auto start = high_resolution_clock::now();
        res = m.vmult(v);
        const auto end = high_resolution_clock::now();

        const auto diff = duration_cast<milliseconds>(end - start).count();

        std::cout << "Time elapsed to compute matrix-vector multiplication with ordered matrix: " << diff << std::endl;
    }

    {
        // elapsed time for ordered map
        using namespace std::chrono;

        const auto start = high_resolution_clock::now();
        ures = um.vmult(v);
        const auto end = high_resolution_clock::now();

        const auto diff = duration_cast<milliseconds>(end - start).count();

        std::cout << "Time elapsed to compute matrix-vector multiplication with unordered matrix: " << diff << std::endl;
    }

    std::cout << "correct mult (map)? " << eq(e_res, res) << std::endl;
    std::cout << "correct mult (umap)? " << eq(e_res, ures) << std::endl;

    return 0;
}
