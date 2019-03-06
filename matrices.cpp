#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include "matrices.hpp"

matrix::matrix() noexcept : N_(0), M_(0), data_(nullptr) {
}

matrix::matrix(size_t N, size_t M) : N_(N), M_(M), data_(nullptr) {
    if ((N == 0) ^ (M == 0))
        throw std::invalid_argument("matrix sizes should be positive");

    if (N * M / M != N)
        throw std::invalid_argument("matrix size is too large (size_t overflow)");

    data_ = new double[N * M];
}

matrix& matrix::resize(size_t N, size_t M) {
    if (N == 0 || M == 0)
        throw std::invalid_argument("matrix sizes should be positive!");

    if (N * M / M != N)
        throw std::invalid_argument("matrix size is too large (size_t overflow)");

    delete[] data_;

    N_ = N;
    M_ = M;
    data_ = new double[M * N];

    return *this;
}

matrix::matrix(const matrix& m) : N_(m.N_), M_(m.M_) {
    if (m.data_ == nullptr)
        return;

    data_ = new double[N_ * M_];

    std::copy(data_, data_ + N_ * M_, m.data_);
}

matrix::matrix(matrix&& m) noexcept : N_(m.N_), M_(m.N_), data_(m.data_) {
    m.data_ = nullptr;
}

matrix& matrix::operator=(matrix&& m) noexcept {
    if (this == &m)
        return *this;

    delete[] data_;

    N_ = m.N_;
    M_ = m.M_;
    data_ = m.data_;

    m.data_ = nullptr;
    m.N_ = 0;
    m.M_ = 0;

    return *this;
}

matrix& matrix::operator=(const matrix& m) {
    if (this == &m)
        return *this;

    delete[] data_;
    N_ = m.N_;
    M_ = m.M_;

    if (m.data_ == nullptr) {
        data_ = nullptr;
        return *this;
    }
    data_ = new double[N_ * M_];

    std::copy(data_, data_ + N_ * M_, m.data_);

    return *this;
}

double& matrix::operator()(const size_t i, const size_t j) {
    if (data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix");

    if (i > N_ || j > M_)
        throw std::out_of_range("incorrect matrix indexing");

    return data_[i * M_ + j];
}

double matrix::operator()(const size_t i, const size_t j) const {
    if (data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix");

    if (i > N_ || j > M_)
        throw std::out_of_range("incorrect matrix indexing");
    return data_[i * M_ + j];
}

matrix matrix::operator*(const matrix& m) const {
    if (data_ == nullptr || m.data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix in multiplication");

    if (M_ != m.N_)
        throw std::invalid_argument("matrix dimensions in multiplication are incompatible");

    matrix res(N_, m.M_);

    for (size_t i = 0; i < res.N_; i++)
        for (size_t j = 0; j < res.M_; j++) {
            double sum = 0;
            for (size_t k = 0; k < M_; k++) {
                sum += (*this)(i, k) * m(k, j);
            }
            res(i, j) = sum;
        }

    return res;
}


matrix& matrix::operator*=(const matrix& m) {
    if (data_ == nullptr || m.data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix in multiplication");

    if (M_ != m.N_)
        throw std::invalid_argument("matrix dimensions in multiplication are incompatible");

    matrix res(N_, m.M_);

    for (size_t i = 0; i < res.N_; i++)
        for (size_t j = 0; j < res.M_; j++) {
            double val = 0;
            for (size_t k = 0; k < M_; k++) {
                val += (*this)(i, k) * m(k, j);
            }
            res(i, j) = val;
        }

    *this = std::move(res);
    return *this;
}

matrix matrix::operator+(const matrix& m) const {
    if (data_ == nullptr || m.data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix in add");

    if (N_ != m.N_ || M_ != m.M_)
        throw std::invalid_argument("matrix dimensions in add are incompatible");

    auto res = *this;

    res += m;

    return res;
}

matrix& matrix::operator+=(const matrix& m) {
    if (data_ == nullptr || m.data_ == nullptr)
        throw std::runtime_error("use of uninitialized matrix in addition");

    if (N_ != m.N_ || M_ != m.M_)
        throw std::invalid_argument("matrix dimensions in addition are incompatible");

    for (size_t i = 0; i < N_; i++)
        for (size_t j = 0; j < M_; j++) {
            (*this)(i, j) += m(i, j);
        }

    return *this;
}

matrix::~matrix() {
    delete[] data_;
}

size_t matrix::get_rows() const noexcept {
    return N_;
}

size_t matrix::get_columns() const noexcept {
    return M_;
}

std::istream& operator>>(std::istream& stream, matrix& m) {
    size_t N, M;

    stream >> N >> M;

    m.resize(N, M);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < M; j++)
            stream >> m(i, j);

    return stream;
}


std::ostream& operator<<(std::ostream& stream, const matrix& m) {
    if (m.get_rows() == 0)
        return stream;

    stream << m.get_rows() << " " << m.get_columns() << std::endl;

    for (size_t i = 0; i < m.get_rows(); i++) {
        stream << m(i, 0);
        for (size_t j = 1; j < m.get_columns(); j++)
            stream << " " << m(i, j);
        stream << std::endl;
    }

    return stream;
}
