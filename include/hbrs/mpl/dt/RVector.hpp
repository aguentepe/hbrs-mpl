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

/* runtime-size array view */

#include <hbrs/mpl/dt/Range.hpp>
#include <hbrs/mpl/dt/CVector.hpp>
#include <algorithm>
#include <iostream>
//#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

class RVector
{
    CVector<double> vector_;
public:
    RVector(std::size_t const m)
        : vector_ (m) {}
    explicit RVector(CVector<double> const& v)
        : vector_ {v} {}
    RVector(std::initializer_list<double> const& l)
        : vector_ {l} {}

    std::size_t n() const {
        return this->vector_.m();
    }
    double* vector() const {
        return vector_.vector();
    }
    
    double& at(std::size_t const x) const {
        return this->vector_.at(x);
    }

    CVector<double> transpose() const {
        return vector_;
    }
    
    friend RVector operator* (RVector const& v, double const s) {
        return RVector(v.vector_ * s);
    }
    friend RVector operator* (double const s, RVector const& v) {
        return RVector(s * v.vector_);
    }

    friend RVector& fill(RVector& v, double const value) {
        fill(v.vector_, value);
        return v;
    }
    friend RVector& fill(RVector&& v, double const value) {
        fill(v.vector_, value);
        return v;
    }
};

double operator* (RVector const&, CVector<double> const&);

std::ostream& operator<< (std::ostream&, RVector const&);

inline std::ostream& operator<< (std::ostream& os, RVector const& v) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < v.n(); ++i)
        os << v.at(i) << "\t";
    os << std::endl;
    return os << '-' << std::endl;
}

inline double operator* (RVector const& v1, CVector<double> const& v2) {
    double sum {0};
    for (std::size_t i {0}; i < v1.n(); ++i)
        sum += v1.at(i) * v2.at(i);
    return sum;
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

/* template <typename T> */
/* struct tag_of< hbrs::mpl::rtsav<T> > { */
/* 	using type = hbrs::mpl::rtsav_tag; */
/* }; */

/* template <> */
/* struct make_impl<hbrs::mpl::rtsav_tag> { */
/* 	template <typename T> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T* & a, std::size_t l) { */
/* 		return {a, l}; */
/* 	} */
	
/* 	template <typename T> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T* const& a, std::size_t l) { */
/* 		return {a, l}; */
/* 	} */
	
/* 	template <typename T, std::size_t Length> */
/* 	static constexpr hbrs::mpl::rtsav<T> */
/* 	apply(T (&a) [Length]) { */
/* 		return {a, Length}; */
/* 	} */
/* }; */

/* namespace hana */ } /* namespace boost */ }
