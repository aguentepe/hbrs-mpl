/* Copyright (c) 2016-2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_DETAIL_GATHER_IMPL_HPP
#define HBRS_MPL_DETAIL_GATHER_IMPL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>

#ifdef HBRS_MPL_ENABLE_MATLAB
	#include <hbrs/mpl/dt/ml_matrix.hpp>
	#include <hbrs/mpl/dt/ml_vector.hpp>
#endif

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	#include <hbrs/mpl/dt/el_matrix.hpp>
	#include <hbrs/mpl/dt/el_dist_matrix.hpp>
	#include <hbrs/mpl/dt/el_vector.hpp>
	#include <hbrs/mpl/dt/el_dist_vector.hpp>
	#include <El.hpp>
#endif

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/smr.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/submatrix.hpp>
#include <hbrs/mpl/dt/subsequence.hpp>
#include <hbrs/mpl/dt/scv.hpp>
#include <hbrs/mpl/dt/srv.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/ext/std/vector.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

template <
	typename T,
	typename std::enable_if_t< 
		#ifdef HBRS_MPL_ENABLE_MATLAB
			std::is_same< hana::tag_of_t<T>, ml_matrix_tag >::value ||
			std::is_same< hana::tag_of_t<T>, ml_column_vector_tag >::value ||
			std::is_same< hana::tag_of_t<T>, ml_row_vector_tag >::value ||
		#endif
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			std::is_same< hana::tag_of_t<T>, el_matrix_tag >::value ||
			std::is_same< hana::tag_of_t<T>, el_column_vector_tag >::value ||
			std::is_same< hana::tag_of_t<T>, el_row_vector_tag >::value ||
		#endif
		std::is_same< hana::tag_of_t<T>, sm_tag >::value ||
		std::is_same< hana::tag_of_t<T>, smr_tag >::value ||
		std::is_same< hana::tag_of_t<T>, ctsam_tag >::value ||
		std::is_same< hana::tag_of_t<T>, rtsam_tag >::value ||
		std::is_same< hana::tag_of_t<T>, submatrix_tag >::value ||
		std::is_same< hana::tag_of_t<T>, hana::ext::std::array_tag>::value ||
		std::is_same< hana::tag_of_t<T>, hana::ext::std::vector_tag>::value ||
		std::is_same< hana::tag_of_t<T>, subsequence_tag >::value ||
		std::is_same< hana::tag_of_t<T>, scv_tag >::value ||
		std::is_same< hana::tag_of_t<T>, srv_tag >::value
	>* = nullptr
>
constexpr decltype(auto)
gather(T && t) {
	return HBRS_MPL_FWD(t);
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	#define _DEF_GATHER_DIST(kind)                                                                                     \
		template <typename Ring, El::Dist Columnwise, El::Dist Rowwise, El::DistWrap Wrapping>                         \
		constexpr auto                                                                                                 \
		gather(el_dist_ ## kind<Ring, Columnwise, Rowwise, Wrapping> const& t) {                                       \
			typedef std::decay_t<Ring> _Ring_;                                                                         \
			El::DistMatrix<_Ring_, El::CIRC, El::CIRC, Wrapping> dmat{t.data()};                                       \
			                                                                                                           \
			if (dmat.Grid().Rank() == 0) {                                                                             \
				return make_el_ ## kind(dmat.Matrix());                                                                \
			} else {                                                                                                   \
				return make_el_ ## kind(El::Matrix<_Ring_>{dmat.Height(), dmat.Width()});                              \
			}                                                                                                          \
		}
	
	_DEF_GATHER_DIST(column_vector)
	_DEF_GATHER_DIST(row_vector)
	_DEF_GATHER_DIST(matrix)
	
	#undef _DEF_GATHER_DIST
#endif

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_DETAIL_GATHER_IMPL_HPP
