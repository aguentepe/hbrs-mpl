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
#include <algorithm>
#include <iostream>
//#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<typename Element>
class CVector {
    Element* vector_;
    std::size_t m_;

public:
    CVector (std::size_t const m)
        : vector_ { new Element[m] }, m_ {m}
    {}
    CVector (CVector const& v)
        :   vector_ { new Element[v.m_] }, m_ {v.m_} {
        std::copy (v.vector_, v.vector_ + v.m_, this->vector_);
    }
    CVector (std::initializer_list<Element> const& l)
        : vector_ { new Element[l.size()] }, m_ {l.size() } {
        Element const* const begin{l.begin() };
        for (std::size_t i{0}; i < m_; ++i)
            at (i) = * (begin + i);
    }
    ~CVector() {
        delete[] vector_;
    }
    void swap (CVector& v) noexcept {
        std::swap (this->vector_, v.vector_);
        std::swap (this->m_, v.m_);
    }
    CVector& operator= (CVector v) {
        swap (v);
        return *this;
    }

    std::size_t m() const {
        return this->m_;
    }
    Element* vector() const {
        return vector_;
    }
    
    /* friend Element& at (CVector const& v, std::size_t const index) { */
    /*     return v.vector_[index]; */
    /* } */
    Element const* begin() const {
        return &at(0);
    }
    Element const* end() const {
        return &at(m_);
    }
    
    CVector operator() (Range const& r) const {
        CVector v (r.end - r.begin + 1);
        for (std::size_t i = 0; i < v.m(); ++i) {
            v.at (i) = this->at (i + r.begin);
        }
        return v;
    }
    Element& at (std::size_t const x) const {
        return this->vector_[x];
    }
    friend CVector& fill(CVector& v, Element const value) {
        std::fill(v.vector_, v.vector_ + v.m(), value);
        return v;
    }
    friend CVector& fill(CVector&& v, Element const value) {
        std::fill(v.vector_, v.vector_ + v.m(), value);
        return v;
    }
};

template<typename Element>
void swap (CVector<Element>& v1, CVector<Element>& v2) noexcept {
    v1.swap (v2);
}

template<typename Element>
std::ostream& operator<< (std::ostream& os, CVector<Element> const& v) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < v.m(); ++i)
        os << v.at (i) << std::endl;
    return os << '-' << std::endl;
}

template<typename Element>
Element operator* (CVector<Element> const& v1, CVector<Element> const& v2){
    Element sum {0};
    for (std::size_t i {1}; i < v1.m(); ++i)
        sum += v1.at (i) * v2.at (i);
    return sum;
}

template<typename Element>
CVector<Element> operator* (Element const s, CVector<Element> const v){
    for (std::size_t i {0}; i < v.m(); ++i)
        v.at (i) *= s;
    return v;
}

template<typename Element>
CVector<Element> operator* (CVector<Element> const& v, Element const s) {
    return s * v;
}

template<typename Element>
CVector<Element> operator+ (Element const d, CVector<Element> const& v){
    CVector<Element> result (v.m() + 1);
    result.at (0) = d;
    for (std::size_t i = 1; i < v.m() + 1; ++i)
        result.at (i) = v.at (i - 1);
    return result;
}

template<typename Element>
CVector<Element> operator/ (CVector<Element> const& v, Element const d) {
    return 1. / d * v;
}
  
template<typename Element>
bool operator== (CVector<Element> const& v, std::initializer_list<Element> const& l){
    if (v.m() != l.size())
        return false;
    for (std::size_t i {0}; i < v.m(); ++i)
        if (v.at (i) != * (l.begin() + i))
            return false;
    return true;
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
