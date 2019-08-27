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

#ifndef HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP
#define HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP

#include "../fwd/hbrs_mpl.hpp"

#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<typename Ring>
auto
multiply_impl_rtsacv_ring::operator()(rtsacv<Ring> v, Ring const& s) const {
	for (std::size_t i {0}; i < v.length(); ++i) {
		v.at(i) *= s;
	}
	return v;
}


template<typename Ring>
decltype(auto)
multiply_impl_ring_rtsacv::operator()(Ring const& s, rtsacv<Ring> v) const {
	return multiply(v,s);
}

template<typename Ring>
decltype(auto)
multiply_impl_rtsarv_rtsacv::operator()(rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) const {
	typedef std::decay_t<Ring> _Ring_;
	_Ring_ sum {0};
	for (std::size_t i {0}; i < v1.length(); ++i)
		sum += v1.at(i) * v2.at(i);
	return sum;
}


template<typename Ring>
decltype(auto)
multiply_impl_rtsacv_rtsarv::operator()(rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) const {
	typedef std::decay_t<Ring> _Ring_;
	//TODO: Implement code to choose best storage_order
	rtsam<_Ring_,storage_order::row_major> result {v1.length(), v2.length()};
	for (std::size_t i {0}; i < (*m)((*size)(result)); ++i) {
		for (std::size_t j {0}; j < (*n)((*size)(result)); ++j) {
			result.at(make_matrix_index(i, j)) = v1.at(i) * v2.at(j);
		}
	}
	return result;
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END
#endif // !HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP
