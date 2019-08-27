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

#ifndef HBRS_MPL_FN_SELECT_IMPL_HBRS_MPL_HPP
#define HBRS_MPL_FN_SELECT_IMPL_HBRS_MPL_HPP

#include "../fwd/hbrs_mpl.hpp"

#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

namespace detail {

template<
	typename Ring,
	storage_order Order
>
decltype(auto)
select_impl_matrix_range_size::operator()(rtsam<Ring,Order> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const {
	return impl(M, range, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
select_impl_matrix_range_size::operator()(submatrix<rtsam<Ring,Order>&,Offset,Size> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const {
	return impl(M, range, hana::type_c<Ring>);
}

template<typename Matrix, typename Range, typename Ring>
decltype(auto)
select_impl_matrix_range_size::impl(Matrix const& M, Range const& range, hana::basic_type<Ring>) const {
	auto rows   {range.first };
	auto column {range.second};

	std::size_t const length {rows.last() - rows.first() + 1};
	rtsacv<Ring> v(length);
	for (std::size_t i {0}; i < length; ++i) {
		v.at(i) = HBRS_MPL_FWD(M).at(make_matrix_index(i + rows.first(), column));
	}
	return v;
}

template<
	typename Ring,
	storage_order Order
>
decltype(auto)
select_impl_matrix_size_range::operator()(rtsam<Ring,Order> const& M, std::pair<std::size_t, range<std::size_t,std::size_t>> const& range) const {
	return impl(M,range, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
select_impl_matrix_size_range::operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M, std::pair<std::size_t, range<std::size_t,std::size_t>> const& range) const {
	return impl(M,range, hana::type_c<Ring>);
}

template<
	typename Ring,
	typename Matrix,
	typename Range
>
decltype(auto)
select_impl_matrix_size_range::impl(Matrix const& M, Range const& range, hana::basic_type<Ring>) const {
	auto row     {range.first };
	auto columns {range.second};

	std::size_t const length {columns.last() - columns.first() + 1};
	rtsarv<Ring> v(length);
	for (std::size_t i {0}; i < length; ++i) {
		v.at(i) = HBRS_MPL_FWD(M).at(make_matrix_index(row, i + columns.first()));
	}
	return v;
}

template<typename Ring, storage_order Order>
decltype(auto)
select_impl_rtsam_range_range::operator()(rtsam<Ring,Order>& M, std::pair<range<std::size_t, std::size_t>, range<std::size_t, std::size_t>> const& ranges) const {
	auto rows    {ranges.first };
	auto columns {ranges.second};
	BOOST_ASSERT(   rows.last() < M.size().m());
	BOOST_ASSERT(columns.last() < M.size().n());
	return submatrix<rtsam<Ring,Order>&, matrix_index<std::size_t, std::size_t>, matrix_size<std::size_t, std::size_t>>
		{M, make_matrix_index(rows.first(), columns.first()), make_matrix_size(rows.last() - rows.first() + 1, columns.last() - columns.first() + 1)};
}
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END


#endif // !HBRS_MPL_FN_SELECT_IMPL_HBRS_MPL_HPP
