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

#ifndef HBRS_MPL_FN_ALMOST_EQUAL_IMPL_HBRS_MPL_HPP
#define HBRS_MPL_FN_ALMOST_EQUAL_IMPL_HBRS_MPL_HPP

#include "../fwd/hbrs_mpl.hpp"

#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/fn/almost_equal.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/size.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

bool
almost_equal_impl_rtsam::operator()(rtsam<double,storage_order::row_major> const& M1, rtsam<double,storage_order::row_major> const& M2) const {
	if (m(size(M1)) != m(size(M2)) || n(size(M1)) != n(size(M2))) {
		return false;
	}
	for (std::size_t i {0}; i < m(size(M1)); ++i) {
		for (std::size_t j {0}; j < n(size(M1)); ++j) {
			if (!almost_equal(M1.at(make_matrix_index(i, j)), M2.at(make_matrix_index(i, j)))) {
				return false;
			}
		}
	}
	return true;
}

bool
almost_equal_impl_rtsam::operator()(rtsam<double,storage_order::column_major> const& M1, rtsam<double,storage_order::column_major> const& M2) const {
	if (m(size(M1)) != m(size(M2)) || n(size(M1)) != n(size(M2))) {
		return false;
	}
	for (std::size_t j {0}; j < m(size(M1)); ++j) {
		for (std::size_t i {0}; i < n(size(M1)); ++i) {
			if (!almost_equal(M1.at(make_matrix_index(i, j)), M2.at(make_matrix_index(i, j)))) {
				return false;
			}
		}
	}
	return true;
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FN_ALMOST_EQUAL_IMPL_HBRS_MPL_HPP
