/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_DT_RTSAM_HPP
#define HBRS_MPL_DT_RTSAM_HPP

#include <hbrs/mpl/fwd/dt/rtsam.hpp>

#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/smr.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/fn/minus.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/detail/translate_index.hpp>
#include <hbrs/mpl/dt/exception.hpp>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>
#include <vector>

HBRS_MPL_NAMESPACE_BEGIN

template<
	typename /* type of matrix entries */ Ring,
	storage_order Order
>
struct rtsam {
	rtsam(matrix_size<std::size_t, std::size_t> sz) : rtsam(sz.m(), sz.n()) {}
	
	rtsam(std::vector<Ring> data, matrix_size<std::size_t, std::size_t> sz) : data_{data}, size_{sz} {
		if (sz.m() * sz.n() != data.size()) {
			/* BOOST_THROW_EXCEPTION( */
			/* 	incompatible_matrix_sequence_exception{} */
			/* 	<< errinfo_matrix_size{sz} */
			/* 	<< errinfo_sequence_size{data.size()} */
			/* ); */
		}
	}
	
	rtsam(std::size_t m, std::size_t n) : data_(m * n, Ring{0}), size_{m,n} {
		BOOST_ASSERT(m * n >= 0);
	}

	rtsam(std::initializer_list<Ring> const& l, std::size_t const rows) : data_(l.size()), size_{rows, l.size() / rows} {
		decltype(auto) begin{l.begin()};
		for (std::size_t i{0}; i < size_.m(); ++i) {
			for (std::size_t j{0}; j < size_.n(); ++j) {
				at(make_matrix_index(i,j)) = *(begin + i * size_.n() + j);
			}
		}
	}

    explicit rtsam(rtsacv<Ring> const& v) : data_ {v.m()}, size_{v.m(), 1} {
        std::copy(v.vector().begin(), v.vector().end(), data_.begin());
    }
    explicit rtsam(rtsarv<Ring> const& v) : data_ {v.n()}, size_{1, v.n()} {
        std::copy(v.vector().begin(), v.vector().end(), data_.begin());
    }
	
	rtsam(rtsam const&) = default;
	rtsam(rtsam &&) = default;
	
	rtsam&
	operator=(rtsam const&) = default;
	rtsam&
	operator=(rtsam &&) = default;
	
	auto const&
	size() const { return size_; };
	
	decltype(auto)
	order() const { return storage_order_c<Order>; }

	auto m() const {
		return size_.m();
	}

	auto n() const {
		return size_.n();
	}
	
	template<typename Index>
	decltype(auto)
	at(Index && i) {
		return data_[
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		];
	}
	
	template<typename Index>
	decltype(auto)
	at(Index && i) const {
		return data_[
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		];
	}
	
	template<typename Index>
	auto
	operator[](Index && i) & { return smr<rtsam &, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	auto
	operator[](Index && i) const& { return smr<rtsam const&, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	auto
	operator[](Index && i) && { return make_smr(std::move(*this), HBRS_MPL_FWD(i)); }
	
    // Several operator() functions to return submatrices.
	auto
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const {
        std::size_t const length {rows.last() - rows.first() + 1};
        rtsacv<Ring> v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = at(make_matrix_index(i + rows.first(), column));
        }
        return v;
    }
	auto
	operator()(std::size_t const& row, range<std::size_t,std::size_t> const& columns) const {
        std::size_t const length {columns.last() - columns.first() + 1};
        rtsarv<Ring> v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = at(make_matrix_index(row, i + columns.first()));
        }
        return v;
    }
	auto
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) const {
        BOOST_ASSERT(rows.last()    < size_.m());
        BOOST_ASSERT(columns.last() < size_.n());
        std::size_t const rLength {rows.last() - rows.first() + 1};
        std::size_t const cLength {columns.last() - columns.first() + 1};
        rtsam<Ring,Order> M {rLength, cLength};
        for (std::size_t i {0}; i < rLength; ++i) {
            for (std::size_t k {0}; k < cLength; ++k) {
                M.at(make_matrix_index(i, k)) = at(make_matrix_index(i + rows.first(), k + columns.first()));
            }
        }
        return M;
    }
	auto
	operator()(rtsacv<std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) const {
        std::size_t const cLength {columns.last() - columns.first() + 1};
        rtsam<Ring,Order> M {rows.m(), cLength};
        for (auto const& i : rows) {
            for (std::size_t k {0}; k < cLength; ++k) {
                M.at(make_matrix_index(i, k)) = at(make_matrix_index(i, k + columns.first()));
            }
        }
        return M;
    }
	auto
	operator()(range<std::size_t,std::size_t> const& rows, rtsacv<std::size_t> const& columns) const {
        std::size_t const rLength {rows.last() - rows.first() + 1};
        rtsam<Ring,Order> M {rLength, columns.m()};
        for (std::size_t i {0}; i < rLength; ++i) {
            for (auto const k : columns) {
                M.at(make_matrix_index(i, k)) = at(make_matrix_index(i + rows.first(), k));
            }
        }
        return M;
    }

private:
	std::vector<Ring> data_;
	matrix_size<std::size_t, std::size_t> size_;
};

template<
	typename Ring,
	storage_order Order
>
std::ostream&
operator<< (std::ostream& os, rtsam<Ring,Order> const& M) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < M.m(); ++i) {
        for (std::size_t j {0}; j < M.n(); ++j)
            os << M.at(make_matrix_index(i, j)) << "\t\t";
        os << std::endl;
    }
    return os << '-' << std::endl;
}

// overwrites rows and columns of M1 with M2
template<
	typename Ring,
	storage_order Order
>
void
overwrite(rtsam<Ring,Order>& M1, range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns, rtsam<Ring,Order> const& M2) {
    for (std::size_t i {rows.first()}; i <= rows.last(); ++i) {
        for (std::size_t j {columns.first()}; j <= columns.last(); ++j) {
            M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i - rows.first(), j - columns.first()));
        }
    }
}

template<
	typename Ring,
	storage_order Order
>
void
overwrite(rtsam<Ring,Order>& M, range<std::size_t,std::size_t> const& rows, double const column, rtsacv<Ring> const& v) {
    for (std::size_t i {rows.first()}; i <= rows.last(); ++i) {
        M.at(make_matrix_index(i,column)) = v.at(i - rows.first());
    }
}

template<
	typename Ring,
	storage_order Order
>
void
overwrite(rtsam<Ring,Order>& M, double const row, range<std::size_t,std::size_t> const& columns, rtsarv<Ring> const& v) {
    for (std::size_t j {columns.first()}; j <= columns.last(); ++j) {
        M.at(make_matrix_index(row,j)) = v.at(j - columns.first());
    }
}

template<
	typename Ring,
	storage_order Order
>
void
overwrite(rtsam<Ring,Order>& M1, rtsacv<std::size_t> const rows, range<std::size_t,std::size_t> const& columns, rtsam<Ring,Order> const& M2){
    for (auto const i : rows) {
        for (std::size_t j {columns.first()}; j <= columns.last(); ++j) {
            M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i, j - columns.first()));
        }
    }
}

template<
	typename Ring,
	storage_order Order
>
void
overwrite(rtsam<Ring,Order>& M1, range<std::size_t,std::size_t> const& rows, rtsacv<std::size_t> const columns, rtsam<Ring,Order> const& M2){
    for (std::size_t i {rows.first()}; i <= rows.last(); ++i) {
        for (auto const j : columns) {
            M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i - rows.first(), j));
        }
    }
}

// returns a square Identity Matrix with size amount of rows and columns
template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
identity(std::size_t const size) {
    rtsam<Ring,Order> result {size, size};
    for (std::size_t i {0}; i < size; ++i) {
        for (std::size_t j {0}; j < size; ++j) {
            result.at(make_matrix_index(i, j)) = 0;
        }
        result.at(make_matrix_index(i, i)) = 1;
    }
    return result;
}

template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
operator* (rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) {
	return multiply(M1,M2);
}

template<
	typename Ring,
	storage_order Order
>
rtsarv<Ring>
operator* (rtsarv<Ring> const& v, rtsam<Ring,Order> const& M) {
	return multiply(v,M);
}

template<
	typename Ring,
	storage_order Order
>
rtsacv<Ring>
operator* (rtsam<Ring,Order> const& M, rtsacv<Ring> const& v) {
	return multiply(M,v);
}

template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
operator* (rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) {
	return multiply(v1,v2);
}

template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
operator* (rtsam<Ring,Order> M, Ring const& d) {
	return multiply(M,d);
}

template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
operator* (Ring const& d, rtsam<Ring,Order> M) {
	return multiply(d,M);
}

template<
	typename Ring,
	storage_order Order
>
rtsam<Ring,Order>
operator- (rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) {
	return minus(M1,M2);
}

template<
	typename Ring,
	storage_order Order
>
bool
operator== (rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) {
	return equal(M1,M2);
}

template<
	typename Ring,
	storage_order Order
>
bool
operator== (rtsam<Ring,Order> const& M, const std::initializer_list< double >& l) {
	return equal(M,l);
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <
	typename Ring,
	hbrs::mpl::storage_order Order
>
struct tag_of< hbrs::mpl::rtsam<Ring, Order> > {
	using type = hbrs::mpl::rtsam_tag;
};

template <>
struct make_impl<hbrs::mpl::rtsam_tag> {
	template <
		typename Ring,
		hbrs::mpl::storage_order Order
	>
	static hbrs::mpl::rtsam<Ring, Order>
	apply(hana::basic_type<Ring>, hbrs::mpl::matrix_size<std::size_t, std::size_t> sz, hbrs::mpl::storage_order_<Order>) {
		return {sz};
	}
	
	template <
		typename Ring,
		hbrs::mpl::storage_order Order
	>
	static hbrs::mpl::rtsam<std::remove_const_t<Ring>, Order>
	apply(std::vector<Ring> data, hbrs::mpl::matrix_size<std::size_t, std::size_t> sz, hbrs::mpl::storage_order_<Order>) {
		return {data, sz};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_RTSAM_HPP
