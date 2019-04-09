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
#include <hbrs/mpl/detail/translate_index.hpp>
#include <vector>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>

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

    explicit rtsam(rtsacv<Ring> const& v) : data_ {v.data()}, size_{v.m(), 1} {}

    explicit rtsam(rtsarv<Ring> const& v) : data_ {v.data()}, size_{1, v.n()} {}
	
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
		return data_.at(
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		);
	}
	
	template<typename Index>
	decltype(auto)
	at(Index && i) const {
		return data_.at(
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		);
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
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const;
	decltype(auto)
	operator()(std::size_t const row, range<std::size_t,std::size_t> const& columns) const;
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) &;
	decltype(auto)
	operator()(std::vector<std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) const {
        std::size_t const cLength {columns.last() - columns.first() + 1};
        rtsam<Ring,Order> M {rows.size(), cLength};
        for (auto const& i : rows) {
            for (std::size_t k {0}; k < cLength; ++k) {
                M.at(make_matrix_index(i, k)) = at(make_matrix_index(i, k + columns.first()));
            }
        }
        return M;
    }
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::vector<std::size_t> const& columns) const {
        std::size_t const rLength {rows.last() - rows.first() + 1};
        rtsam<Ring,Order> M {rLength, columns.size()};
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
/* template< */
/* 	typename Ring, */
/* 	storage_order Order */
/* > */
/* void */
/* overwrite(rtsam<Ring,Order>& M1, range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns, rtsam<Ring,Order> const& M2) { */
/*     for (std::size_t i {rows.first()}; i <= rows.last(); ++i) { */
/*         for (std::size_t j {columns.first()}; j <= columns.last(); ++j) { */
/*             M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i - rows.first(), j - columns.first())); */
/*         } */
/*     } */
/* } */

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

/* template< */
/* 	typename Ring, */
/* 	storage_order Order */
/* > */
/* void */
/* overwrite(rtsam<Ring,Order>& M1, std::vector<std::size_t> const rows, range<std::size_t,std::size_t> const& columns, rtsam<Ring,Order> const& M2){ */
/*     for (auto const i : rows) { */
/*         for (std::size_t j {columns.first()}; j <= columns.last(); ++j) { */
/*             M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i, j - columns.first())); */
/*         } */
/*     } */
/* } */

/* template< */
/* 	typename Ring, */
/* 	storage_order Order */
/* > */
/* void */
/* overwrite(rtsam<Ring,Order>& M1, range<std::size_t,std::size_t> const& rows, std::vector<std::size_t> const columns, rtsam<Ring,Order> const& M2){ */
/*     for (std::size_t i {rows.first()}; i <= rows.last(); ++i) { */
/*         for (auto const j : columns) { */
/*             M1.at(make_matrix_index(i,j)) = M2.at(make_matrix_index(i - rows.first(), j)); */
/*         } */
/*     } */
/* } */

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
