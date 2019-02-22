/* Copyright (c) 2019 Abdullah Güntepe, <abdullah@guentepe.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/* runtime-size array view */

#include <hbrs/mpl/dt/Range.hpp>
#include <hbrs/mpl/dt/CVector.hpp>
#include <hbrs/mpl/dt/RVector.hpp>
#include <algorithm>
#include <iostream>
#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

class Matrix {
    double* matrix_;
    std::size_t m_, n_;

public:
    Matrix(std::size_t m, std::size_t n)
        : matrix_ {new double[m * n]}, m_ {m}, n_ {n} {
        BOOST_ASSERT(m * n >= 0);
    }
    Matrix(Matrix const& M)
        : matrix_ {new double[M.m_ * M.n_]}, m_ {M.m_}, n_ {M.n_} {
        std::copy(M.matrix_, M.matrix_ + M.m_ * M.n_, this->matrix_);
    }
    Matrix(std::initializer_list<double> const& l, std::size_t const rows)
        : matrix_ {new double[l.size()]}, m_ {rows}, n_ {l.size() / rows} {
        double const* begin{l.begin()};
        for (std::size_t i{0}; i < m_; ++i) {
            for (std::size_t j{0}; j < n_; ++j) {
                at(i,j) = *(begin + i * n_ + j);
            }
        }
    }
    explicit Matrix(CVector<double> const& v)
        : matrix_ {new double[v.m()]}, m_ {v.m()}, n_ {1} {
        std::copy(v.vector(), v.vector() + v.m(), matrix_);
    }
    explicit Matrix(RVector const& v)
        : matrix_ {new double[v.n()]}, m_ {1}, n_ {v.n()} {
        std::copy(v.vector(), v.vector() + v.n(), matrix_);
    }
    ~Matrix() noexcept {
        delete[] matrix_;
    }
    void swap(Matrix& M) noexcept {
        std::swap(this->matrix_, M.matrix_);
        std::swap(this->m_, M.m_);
        std::swap(this->n_, M.n_);
    }
    Matrix& operator= (Matrix M) {
        swap(M);
        return *this;
    }

    double& operator[](std::size_t i) {
        return this->matrix_[i];
    }
    auto m() const {
        return this->m_;
    }
    auto n() const {
        return this->n_;
    }
    
    // Several operator() functions to return submatrices.
    CVector<double> operator()(Range const& rows, std::size_t const column) const {
        std::size_t const length {rows.end - rows.begin + 1};
        CVector<double> v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = at(i + rows.begin, column);
        }
        return v;
    }
    RVector operator()(std::size_t const& row, Range const& columns) const {
        std::size_t const length {columns.end - columns.begin + 1};
        RVector v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = at(row, i + columns.begin);
        }
        return v;
    }
    Matrix operator()(Range const& rows, Range const& columns) const {
        BOOST_ASSERT(rows.end    < m_);
        BOOST_ASSERT(columns.end < n_);
        std::size_t const rLength {rows.end - rows.begin + 1};
        std::size_t const cLength {columns.end - columns.begin + 1};
        Matrix M {rLength, cLength};
        for (std::size_t i {0}; i < rLength; ++i) {
            for (std::size_t k {0}; k < cLength; ++k) {
                M.at(i, k) = at(i + rows.begin, k + columns.begin);
            }
        }
        return M;
    }
    Matrix operator()(CVector<std::size_t> const& rows, Range const& columns) const {
        std::size_t const cLength {columns.end - columns.begin + 1};
        Matrix M {rows.m(), cLength};
        for (auto const& i : rows) {
            for (std::size_t k {0}; k < cLength; ++k) {
                M.at(i, k) = at(i, k + columns.begin);
            }
        }
        return M;
    }
    Matrix operator()(Range const& rows, CVector<std::size_t> const& columns) const {
        std::size_t const rLength {rows.end - rows.begin + 1};
        Matrix M {rLength, columns.m()};
        for (std::size_t i {0}; i < rLength; ++i) {
            for (auto const k : columns) {
                M.at(i, k) = at(i + rows.begin, k);
            }
        }
        return M;
    }
    double& at(std::size_t const x, std::size_t const y) const {
        BOOST_ASSERT(x < m_);
        BOOST_ASSERT(y < n_);
        return this->matrix_[x * n_ + y];
    }
    friend Matrix& fill(Matrix& M, double const value) {
        std::fill(M.matrix_, M.matrix_ + M.m_ * M.n_, value);
        return M;
    }
};

inline void swap(Matrix& m1, Matrix& m2) noexcept {
    m1.swap(m2);
}

inline std::ostream& operator<< (std::ostream& os, Matrix const& M) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < M.m(); ++i) {
        for (std::size_t j {0}; j < M.n(); ++j)
            os << M.at(i, j) << "\t\t";
        os << std::endl;
    }
    return os << '-' << std::endl;
}

// overwrites rows and columns of M1 with M2
inline void overwrite(Matrix& M1, Range const& rows, Range const& columns, Matrix const& M2) {
    for (std::size_t i {rows.begin}; i <= rows.end; ++i) {
        for (std::size_t j {columns.begin}; j <= columns.end; ++j) {
            M1.at(i,j) = M2.at(i - rows.begin, j - columns.begin);
        }
    }
}
inline void overwrite(Matrix& M, Range const& rows, double const column, CVector<double> const& v) {
    for (std::size_t i {rows.begin}; i <= rows.end; ++i) {
        M.at(i,column) = v.at(i - rows.begin);
    }
}
inline void overwrite(Matrix& M, double const row, Range const& columns, RVector const& v) {
    for (std::size_t j {columns.begin}; j <= columns.end; ++j) {
        M.at(row,j) = v.at(j - columns.begin);
    }
}
inline void overwrite(Matrix& M1, CVector<std::size_t> const rows, Range const& columns, Matrix const& M2){
    for (auto const i : rows) {
        for (std::size_t j {columns.begin}; j <= columns.end; ++j) {
            M1.at(i,j) = M2.at(i, j - columns.begin);
        }
    }
}
inline void overwrite(Matrix& M1, Range const& rows, CVector<std::size_t> const columns, Matrix const& M2){
    for (std::size_t i {rows.begin}; i <= rows.end; ++i) {
        for (auto const j : columns) {
            M1.at(i,j) = M2.at(i - rows.begin, j);
        }
    }
}

// returns a square Identity Matrix with size amount of rows and columns
inline Matrix identity(std::size_t const size) {
    Matrix result {size, size};
    for (std::size_t i {0}; i < size; ++i) {
        for (std::size_t j {0}; j < size; ++j) {
            result.at(i, j) = 0;
        }
        result.at(i, i) = 1;
    }
    return result;
}

inline Matrix operator* (Matrix const& M1, Matrix const& M2) {
    BOOST_ASSERT(M1.n() == M2.m());

    Matrix result {M1.m(), M2.n()};
    for (std::size_t i {0}; i < result.m(); ++i) {
        for (std::size_t j {0}; j < result.n(); ++j) {
            result.at(i, j) = M1(i, Range(0, M1.n() - 1)) * M2(Range(0, M2.m() - 1), j);
        }
    }
    return result;
}
inline RVector operator* (RVector const& v, Matrix const& M) {
    BOOST_ASSERT(v.n() == M.m());

    RVector result(M.n());
    for (std::size_t i {0}; i < result.n(); ++i) {
        result.at(i) = v * M(Range(0, M.m() - 1), i);
    }
    return result;
}
inline CVector<double> operator* (Matrix const& M, CVector<double> const& v) {
    BOOST_ASSERT(M.n() == v.m());

    CVector<double> result(M.m());
    for (std::size_t i {0}; i < result.m(); ++i) {
        result.at(i) = M(i, Range(0, M.n() - 1)) * v;
    }
    return result;
}
inline Matrix operator* (CVector<double> const& v1, RVector const& v2) {
    Matrix result {v1.m(), v2.n()};
    for (std::size_t i {0}; i < result.m(); ++i) {
        for (std::size_t j {0}; j < result.n(); ++j) {
            result.at(i, j) = v1.at(i) * v2.at(j);
        }
    }
    return result;
}
inline Matrix operator* (Matrix M, double const d) {
    for (std::size_t i {0}; i < M.m(); ++i) {
        for (std::size_t j {0}; j < M.n(); ++j) {
            M.at(i,j) = M.at(i,j) * d;
        }
    }
    return M;
}

inline Matrix operator* (double const d, Matrix M) {
    return M * d;
}

inline Matrix operator- (Matrix const& M1, Matrix const& M2) {
    BOOST_ASSERT(M1.m() == M2.m());
    BOOST_ASSERT(M1.n() == M2.n());

    Matrix result {M1.m(), M1.n()};
    for (std::size_t i {0}; i < M1.m(); ++i) {
        for (std::size_t j {0}; j < M1.n(); ++j) {
            result.at(i, j) = M1.at(i, j) - M2.at(i, j);
        }
    }
    return result;
}

inline bool operator== (Matrix const& M1, Matrix const& M2) {
    if (M1.m() != M2.m() || M1.n() != M2.n()) {
        return false;
    }
    for (std::size_t i {0}; i < M1.m(); ++i) {
        for (std::size_t j {0}; j < M1.n(); ++j) {
            if (!(std::abs(M1.at(i, j) - M2.at(i, j)) <= std::numeric_limits<double>::epsilon() * 10000000000)) { // FIXME use proper epsilon
                return false;
            }
        }
    }
    return true;
}

inline bool operator== (Matrix const& M, const std::initializer_list< double >& l) {
    if (M.m() * M.n() != l.size()) {
        return false;
    }
    for (std::size_t i {0}; i < M.m(); ++i) {
        for (std::size_t j {0}; j < M.n(); ++j) {
            if (M.at(i, j) != * (l.begin() + i * M.n() + j)) {
                return false;
            }
        }
    }
    return true;
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

/* template <typename T> */
/* struct tag_of< hbrs::mpl::rtsav<T> > { */
/* 	using type = hbrs::mpl::rtsav_tag; */
/* }; */

/* template <> */
/* struct make_impl<hbrs::mpl::rtsav_tag> { */
/* 	template <typename T> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T* & a, std::size_t l) { */
/* 		return {a, l}; */
/* 	} */
	
/* 	template <typename T> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T* const& a, std::size_t l) { */
/* 		return {a, l}; */
/* 	} */
	
/* 	template <typename T, std::size_t Length> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T (&a) [Length]) { */
/* 		return {a, Length}; */
/* 	} */
/* }; */

/* namespace hana */ } /* namespace boost */ }
