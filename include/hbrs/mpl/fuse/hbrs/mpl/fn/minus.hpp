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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_MINUS_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_MINUS_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag >
	>* = nullptr
>
decltype(auto)
operator-(T1 && t1, T2 && t2) {
	return minus(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

struct minus_impl_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
		BOOST_ASSERT(M1.m() == M2.m());
		BOOST_ASSERT(M1.n() == M2.n());

		rtsam<Ring,Order> result {M1.m(), M1.n()};
		for (std::size_t i {0}; i < M1.m(); ++i) {
			for (std::size_t j {0}; j < M1.n(); ++j) {
				result.at(make_matrix_index(i, j)) = M1.at(make_matrix_index(i, j)) - M2.at(make_matrix_index(i, j));
			}
		}
		return result;
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_MINUS_IMPLS boost::hana::make_tuple(\
		hbrs::mpl::detail::minus_impl_rtsam{}\
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_MINUS_HPP
