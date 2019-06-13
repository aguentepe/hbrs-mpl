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

#ifndef HBRS_MPL_FN_ROWS_IMPL_ELEMENTAL_HPP
#define HBRS_MPL_FN_ROWS_IMPL_ELEMENTAL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/dt/el_matrix.hpp>
#include <El.hpp>
#include <hbrs/mpl/dt/smrs.hpp>
#include <boost/hana/tuple.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

struct rows_impl_matrix {
	template <
		typename Matrix,
		typename std::enable_if_t< 
			std::is_same< hana::tag_of_t<Matrix>, matrix_tag >::value
		>* = nullptr
	>
	constexpr auto
	operator()(Matrix && a) const {
		return mpl::smrs<Matrix>{HBRS_MPL_FWD(a)};
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_ROWS_IMPLS_ELEMENTAL boost::hana::make_tuple(                                                          \
		elemental::detail::rows_impl_matrix{}                                                                          \
	)

#endif // !HBRS_MPL_FN_ROWS_IMPL_ELEMENTAL_HPP
