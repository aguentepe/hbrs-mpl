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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_EQUAL_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_EQUAL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/fwd/fn/equal.hpp>
#include <hbrs/mpl/detail/function_object.hpp>
#include <hbrs/mpl/fwd/dt/matrix_size.hpp>
#include <hbrs/mpl/fwd/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>

#include <boost/hana/core/tag_of.hpp>
#include <boost/hana/tuple.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag >
	>* = nullptr
>
bool
operator==(T1 && t1, T2 && t2) {
	return equal(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

namespace detail {

HBRS_MPL_DEF_FO_TRY_OPERATOR(equal_impl_matrix_size , matrix_size_tag , HBRS_MPL_OPERATOR_EQUAL, equal)
HBRS_MPL_DEF_FO_TRY_OPERATOR(equal_impl_matrix_index, matrix_index_tag, HBRS_MPL_OPERATOR_EQUAL, equal)

struct equal_impl_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	constexpr bool 
	operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
		if (M1.m() != M2.m() || M1.n() != M2.n()) {
			return false;
		}
		for (std::size_t i {0}; i < M1.m(); ++i) {
			for (std::size_t j {0}; j < M1.n(); ++j) {
				if (M1.at(make_matrix_index(i, j)) != M2.at(make_matrix_index(i, j))) {
					return false;
				}
			}
		}
		return true;
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_EQUAL_IMPLS boost::hana::make_tuple(                                                 \
		hbrs::mpl::detail::equal_impl_matrix_size{},                                                                   \
		hbrs::mpl::detail::equal_impl_matrix_index{},                                                                  \
		hbrs::mpl::detail::equal_impl_rtsam{}                                                                          \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_EQUAL_HPP
