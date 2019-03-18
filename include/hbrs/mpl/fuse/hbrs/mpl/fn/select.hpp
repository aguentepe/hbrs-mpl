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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_SELECT_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_SELECT_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/ssm.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>
#include <utility>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<typename Ring, storage_order Order>
decltype(auto)
rtsam<Ring,Order>::operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const {
	return select(*this, std::make_pair(rows, column));
}

template<typename Ring, storage_order Order>
decltype(auto)
rtsam<Ring,Order>::operator()(std::size_t const& row, range<std::size_t,std::size_t> const& columns) const {
	return select(*this, std::make_pair(row, columns));
}

template<typename Ring, storage_order Order>
decltype(auto)
rtsam<Ring,Order>::operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) {
	return select(*this, std::make_pair(rows, columns));
}

template<typename Matrix, typename Index1, typename Index2>
decltype(auto)
ssm<Matrix, Index1, Index2>::operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const {
	return select(*this, std::make_pair(rows, column));
}

template<typename Matrix, typename Index1, typename Index2>
decltype(auto)
ssm<Matrix, Index1, Index2>::operator()(std::size_t const& row, range<std::size_t,std::size_t> const& columns) const {
	return select(*this, std::make_pair(row, columns));
}

template<typename Matrix, typename Index1, typename Index2>
decltype(auto)
ssm<Matrix, Index1, Index2>::operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) {
	return select(*this, std::make_pair(rows, columns));
}

namespace detail {

struct select_impl_rtsam_range_size {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const {
		decltype(auto) rows {range.first};
		decltype(auto) column {range.second};
        std::size_t const length {rows.last() - rows.first() + 1};
        rtsacv<Ring> v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = M.at(make_matrix_index(i + rows.first(), column));
        }
        return v;
    }
	template<
		typename Matrix,
		typename Index1,
		typename Index2
	>
	decltype(auto)
	operator()(ssm<Matrix,Index1,Index2> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const {
	}
};

struct select_impl_rtsam_size_range {
	template<
		typename Matrix,
		typename std::enable_if_t< std::is_same_v< hana::tag_of_t<Matrix>, rtsam_tag > || std::is_same_v< hana::tag_of_t<Matrix>, ssm_tag > >* = nullptr
	>
	/* constexpr */ 
	decltype(auto)
	operator()(Matrix const& M, std::pair<std::size_t, range<std::size_t,std::size_t>> const& range) const {
		decltype(auto) row {range.first};
		decltype(auto) columns {range.second};
        std::size_t const length {columns.last() - columns.first() + 1};
        rtsarv<typename Matrix::type> v(length);
        for (std::size_t i {0}; i < length; ++i) {
            v.at(i) = M.at(make_matrix_index(row, i + columns.first()));
        }
        return v;
    }
};

struct select_impl_rtsam_range_range {
	/* template< */
	/* 	typename Matrix, */
	/* 	typename std::enable_if_t< std::is_same_v< hana::tag_of_t<Matrix>, rtsam_tag > || std::is_same_v< hana::tag_of_t<Matrix>, ssm_tag > >* = nullptr */
	/* > */
	template<typename Ring, storage_order Order>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order>& M, std::pair<range<std::size_t, std::size_t>, range<std::size_t, std::size_t>> const& ranges) const {
		decltype(auto) rows    {ranges.first};
		decltype(auto) columns {ranges.second};
        BOOST_ASSERT(   rows.last() < M.m());
        BOOST_ASSERT(columns.last() < M.n());
		return ssm<rtsam<Ring,Order>&, matrix_index<std::size_t, std::size_t>, matrix_index<std::size_t, std::size_t>>{M, make_matrix_index(rows.first(), columns.first()), make_matrix_index(rows.last(), columns.last())};
    }
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_SELECT_IMPLS boost::hana::make_tuple(                                              \
		hbrs::mpl::detail::select_impl_rtsam_range_size{},                                                           \
		hbrs::mpl::detail::select_impl_rtsam_size_range{},                                                           \
		hbrs::mpl::detail::select_impl_rtsam_range_range{}                                                           \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_SELECT_HPP
