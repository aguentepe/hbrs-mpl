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

#ifndef HBRS_MPL_FN_MINUS_IMPL_HBRS_MPL_HPP
#define HBRS_MPL_FN_MINUS_IMPL_HBRS_MPL_HPP

#include "../fwd/hbrs_mpl.hpp"

#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/submatrix.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<
	typename Ring,
	storage_order Order
>
decltype(auto)
minus_impl_rtsam::operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
	return impl(M1,M2,hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
minus_impl_rtsam::operator()(submatrix<rtsam<Ring,Order>&, Offset, Size> const& M1, rtsam<Ring,Order> const& M2) const {
	return impl(M1,M2,hana::type_c<Ring>);
}

template<
	typename Ring,
	typename Matrix1,
	typename Matrix2
>
decltype(auto)
minus_impl_rtsam::impl(Matrix1 const& M1, Matrix2 const& M2, hana::basic_type<Ring>) const {
	BOOST_ASSERT((*m)((*size)(M1)) == (*m)((*size)(M2)));
	BOOST_ASSERT((*n)((*size)(M1)) == (*n)((*size)(M2)));

	typedef std::decay_t<Ring> _Ring_;
	rtsam<_Ring_, storage_order::row_major> result {M1.size().m(), M1.size().n()};
	for (std::size_t i {0}; i < M1.size().m(); ++i) {
		for (std::size_t j {0}; j < M1.size().n(); ++j) {
			result.at(make_matrix_index(i, j)) = M1.at(make_matrix_index(i, j)) - M2.at(make_matrix_index(i, j));
		}
	}
	return result;
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FN_MINUS_IMPL_HBRS_MPL_HPP
