#include <vector>
#include <iostream>
#include <cassert>
#include <map>

using elem_t = double;
using Vector = std::vector<elem_t>;

// Step 1: abstract matrix class
class SparseMatrix // abstract class
{
public:
    SparseMatrix() : rows(0), cols(0), nnz(0) {}
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    size_t getNnz() const { return nnz; }

    void print(std::ostream &os = std::cout) const
    {
        os << "Numbero of rows: " << rows << std::endl;
        os << "Number of columns: " << cols << std::endl;
        os << "Number of non-zero elements: " << nnz << std::endl;
        _print(os);
    }
    virtual Vector vmult(const Vector &v) const = 0;
    // element access in read only (const version)
    virtual const elem_t &operator()(size_t i, size_t j) const = 0;
    // element access in write (returns non-const reference)
    virtual elem_t &operator()(size_t i, size_t j) = 0;
    // virtual destructor
    virtual ~SparseMatrix() = default;

protected: // protected because need to be accessible to children!
    size_t rows;
    size_t cols;
    size_t nnz;
    virtual void _print(std::ostream &os) const = 0;
};

class MapMatrix : SparseMatrix
{
public:
    virtual Vector vmult(const Vector &v) const override
    {
        assert(v.size() == cols);
        Vector result(v.size());
    }

    // element access in read only (const version)
    virtual const elem_t &operator()(size_t i, size_t j) const override
    {
        return data.at(i).at(j);
    }

    // element access in write (returns non-const reference)
    virtual elem_t &operator()(size_t i, size_t j) override
    {
        // Row boundary check
        if (i > data.size() - 1) // if accessed element is outside row boundaries, increment the number of rows
        {
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
            // emplace returns an iterator to the newly inserted pair:
            // so i get the memory position pointed by that operator (*data{...}.first)
            return (*data.at(i).emplace(j, 0).first).second;
        }
        return (*it).second;
    }
    virtual ~MapMatrix() override = default;

protected:
    virtual void _print(std::ostream &os) const
    {
        for (size_t i = 0; i < data.size(); i++)
        {
            for (const auto &[j, v] : data[i])
            {
                os << i << ", " << j << ", " << v << std::endl;
            }
        }
    }

private:
    std::vector<std::map<size_t, elem_t>> data;
};

int main()
{
    MapMatrix m;

    return 0;
}
