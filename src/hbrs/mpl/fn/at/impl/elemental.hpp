/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_AT_IMPL_ELEMENTAL_HPP
#define HBRS_MPL_FN_AT_IMPL_ELEMENTAL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/detail/function_object.hpp>

#include <hbrs/mpl/dt/matrix_index/fwd.hpp>
#include <hbrs/mpl/dt/el_matrix/fwd.hpp>
#include <hbrs/mpl/dt/el_dist_matrix/fwd.hpp>
#include <hbrs/mpl/dt/el_vector/fwd.hpp>
#include <hbrs/mpl/dt/el_dist_vector/fwd.hpp>

#include <hbrs/mpl/dt/smr.hpp>
#include <El.hpp>

#include <boost/hana/tuple.hpp>
#include <boost/hana/core/tag_of.hpp>
#include <boost/assert.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace mpl = hbrs::mpl;

namespace detail {

HBRS_MPL_DEF_FO_TRY_METHOD(at_impl_matrix, matrix_tag, at)

struct at_impl_matrix_smr {
	template <
		typename Matrix,
		typename std::enable_if_t< 
			std::is_same< hana::tag_of_t<Matrix>, matrix_tag >::value
		>* = nullptr
	>
	decltype(auto)
	operator()(Matrix && m, El::Int i) const {
		return mpl::smr<Matrix, El::Int>{HBRS_MPL_FWD(m), i};
	}
};

HBRS_MPL_DEF_FO_TRY_METHOD(at_impl_column_vector, column_vector_tag, at)
HBRS_MPL_DEF_FO_TRY_METHOD(at_impl_row_vector, row_vector_tag, at)

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_AT_IMPLS_ELEMENTAL boost::hana::make_tuple(                                                            \
		elemental::detail::at_impl_column_vector{},                                                                    \
		elemental::detail::at_impl_row_vector{},                                                                       \
		elemental::detail::at_impl_matrix{},                                                                           \
		elemental::detail::at_impl_matrix_smr{}                                                                        \
	)

#endif // !HBRS_MPL_FN_AT_IMPL_ELEMENTAL_HPP
