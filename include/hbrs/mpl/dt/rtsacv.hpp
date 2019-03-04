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

#ifndef HBRS_MPL_DT_RTSACV_HPP
#define HBRS_MPL_DT_RTSACV_HPP

#include <hbrs/mpl/fwd/dt/rtsacv.hpp>

#include <hbrs/mpl/preprocessor/core.hpp>
/* #include <hbrs/mpl/dt/matrix_index.hpp> */
/* #include <hbrs/mpl/dt/matrix_size.hpp> */
/* #include <hbrs/mpl/dt/smr.hpp> */
/* #include <hbrs/mpl/dt/storage_order.hpp> */
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>

#include <vector>

HBRS_MPL_NAMESPACE_BEGIN


template<typename Ring>
rtsacv<Ring>::rtsacv(std::initializer_list<Ring> const& l) : data_(l.size()) {
	decltype(auto) begin{l.first()};
	for (std::size_t i{0}; i < l.size(); ++i) {
		at(i) = *(begin + i);
	}
}
	
template<typename Ring>
auto
rtsacv<Ring>::operator() (range<std::size_t,std::size_t> const& r) const {
	rtsacv v (r.last() - r.first() + 1);
	for (std::size_t i = 0; i < v.m(); ++i) {
		v.at(i) = at(i + r.first());
	}
	return v;
}

template<typename Ring>
std::ostream& operator<< (std::ostream& os, rtsacv<Ring> const& v) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < v.m(); ++i)
        os << v.at(i) << std::endl;
    return os << '-' << std::endl;
}

template<typename Ring>
Ring operator* (rtsacv<Ring> const& v1, rtsacv<Ring> const& v2){
    Ring sum {0};
    for (std::size_t i {1}; i < v1.m(); ++i)
        sum += v1.at(i) * v2.at(i);
    return sum;
}

template<typename Ring>
rtsacv<Ring> operator* (Ring const s, rtsacv<Ring> v){
    for (std::size_t i {0}; i < v.m(); ++i)
        v.at(i) *= s;
    return v;
}

template<typename Ring>
rtsacv<Ring> operator* (rtsacv<Ring> const& v, Ring const s) {
    return s * v;
}

template<typename Ring>
rtsacv<Ring> operator+ (Ring const d, rtsacv<Ring> const& v){
    rtsacv<Ring> result(v.m() + 1);
    result.at (0) = d;
    for (std::size_t i = 1; i < v.m() + 1; ++i)
        result.at(i) = v.at(i - 1);
    return result;
}

template<typename Ring>
rtsacv<Ring> operator/ (rtsacv<Ring> const& v, Ring const d) {
    return 1. / d * v;
}
  
template<typename Ring>
bool operator== (rtsacv<Ring> const& v, std::initializer_list<Ring> const& l){
    if (v.m() != l.size())
        return false;
    for (std::size_t i {0}; i < v.m(); ++i)
        if (v.at (i) != * (l.first() + i))
            return false;
    return true;
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename Ring>
struct tag_of< hbrs::mpl::rtsacv<Ring> > {
	using type = hbrs::mpl::rtsacv_tag;
};

template <>
struct make_impl<hbrs::mpl::rtsacv_tag> {
	template <typename Ring>
	static hbrs::mpl::rtsacv<Ring>
	apply(hana::basic_type<Ring>, std::size_t size) {
		return {size};
	}
	
	template <typename Ring>
	static hbrs::mpl::rtsacv<std::remove_const_t<Ring>>
	apply(std::vector<Ring> data, std::size_t size) {
		return {data, size};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_RTSACV_HPP
