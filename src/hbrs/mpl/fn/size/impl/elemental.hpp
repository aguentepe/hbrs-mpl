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

#ifndef HBRS_MPL_FN_SIZE_IMPL_ELEMENTAL_HPP
#define HBRS_MPL_FN_SIZE_IMPL_ELEMENTAL_HPP

#include "../fwd/elemental.hpp"
#ifdef HBRS_MPL_ENABLE_ELEMENTAL

#include <hbrs/mpl/detail/function_object.hpp>
#include <hbrs/mpl/dt/el_matrix.hpp>
#include <hbrs/mpl/dt/el_dist_matrix.hpp>
#include <hbrs/mpl/dt/el_vector.hpp>
#include <hbrs/mpl/dt/el_dist_vector.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_el_matrix, el_matrix_tag, size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_el_dist_matrix, el_dist_matrix_tag, size)

HBRS_MPL_DEF_FO_TRY_METHOD(length_impl_el_column_vector, el_column_vector_tag, length)
HBRS_MPL_DEF_FO_TRY_METHOD(length_impl_el_row_vector, el_row_vector_tag, length)

HBRS_MPL_DEF_FO_TRY_METHOD(length_impl_el_dist_column_vector, el_dist_column_vector_tag, length)
HBRS_MPL_DEF_FO_TRY_METHOD(length_impl_el_dist_row_vector, el_dist_row_vector_tag, length)

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_ENABLE_ELEMENTAL
#endif // !HBRS_MPL_FN_SIZE_IMPL_ELEMENTAL_HPP
