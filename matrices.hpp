#pragma once

#include <cstddef>

class matrix {
private:
    size_t N_, M_;
    double* data_;

public:
    matrix() noexcept;
    matrix(size_t N, size_t M);
    matrix(const matrix& m);
    matrix(matrix&& m) noexcept;
    ~matrix();

    size_t get_rows() const noexcept;
    size_t get_columns() const noexcept;

    matrix& resize(size_t N, size_t M);

    matrix& operator=(matrix&& m) noexcept;
    matrix& operator=(const matrix& m);

    double& operator()(const size_t i, const size_t j);
    double operator()(const size_t i, const size_t j) const;

    matrix operator*(const matrix& m) const;
    matrix& operator*=(const matrix& m);

    matrix operator+(const matrix& m) const;
    matrix& operator+=(const matrix& m);
};

std::istream& operator>>(std::istream& stream, matrix& m);
std::ostream& operator<<(std::ostream& stream, const matrix& m);
