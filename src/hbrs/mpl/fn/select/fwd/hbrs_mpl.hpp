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

#ifndef HBRS_MPL_FN_SELECT_FWD_HBRS_MPL_HPP
#define HBRS_MPL_FN_SELECT_FWD_HBRS_MPL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/submatrix.hpp>
#include <hbrs/mpl/dt/range.hpp>

#include <boost/hana/ext/std/pair.hpp>
#include <boost/hana/core/tag_of.hpp>

#include <type_traits>
#include <utility>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct select_impl_matrix_range_size {
	template<
		typename Ring,
		storage_order Order
	>
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const;
	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&,Offset,Size> const& M, std::pair<range<std::size_t,std::size_t>, std::size_t> const& range) const;

private:
	template<typename Matrix, typename Range, typename Ring>
	decltype(auto)
	impl(Matrix const& M, Range const& range, hana::basic_type<Ring>) const;
};

struct select_impl_matrix_size_range {
	template<
		typename Ring,
		storage_order Order
	>
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M, std::pair<std::size_t, range<std::size_t,std::size_t>> const& range) const;

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M, std::pair<std::size_t, range<std::size_t,std::size_t>> const& range) const;

private:
	template<
		typename Ring,
		typename Matrix,
		typename Range
	>
	decltype(auto)
	impl(Matrix const& M, Range const& range, hana::basic_type<Ring>) const;
};

struct select_impl_rtsam_range_range {
	template<typename Ring, storage_order Order>
	decltype(auto)
	operator()(rtsam<Ring,Order>& M, std::pair<range<std::size_t, std::size_t>, range<std::size_t, std::size_t>> const& ranges) const;
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_SELECT_IMPLS_HBRS_MPL boost::hana::make_tuple(                                                     \
		hbrs::mpl::detail::select_impl_matrix_range_size{},                                                            \
		hbrs::mpl::detail::select_impl_matrix_size_range{},                                                            \
		hbrs::mpl::detail::select_impl_rtsam_range_range{}                                                             \
	)

#endif // !HBRS_MPL_FN_SELECT_FWD_HBRS_MPL_HPP
