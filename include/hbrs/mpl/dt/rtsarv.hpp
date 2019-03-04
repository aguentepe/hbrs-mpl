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

#ifndef HBRS_MPL_DT_RTSARV_HPP
#define HBRS_MPL_DT_RTSARV_HPP

#include <hbrs/mpl/fwd/dt/rtsarv.hpp>

#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
/* #include <hbrs/mpl/fn/multiply.hpp> */

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>
#include <vector>

HBRS_MPL_NAMESPACE_BEGIN

template<typename /* type of vector entries */ Ring>
struct rtsarv {
	explicit rtsarv(rtsacv<Ring> data) : vector_{data} {}
	
	explicit rtsarv(std::size_t size) : vector_(size) {}

	explicit rtsarv(std::initializer_list<Ring> const& l) : vector_(l) {}
	
	rtsarv(rtsarv const&) = default;
	rtsarv(rtsarv &&) = default;
	
	rtsarv&
	operator=(rtsarv const&) = default;
	rtsarv&
	operator=(rtsarv &&) = default;
	
	auto const
	size() const {
		return vector_.size();
	}

	decltype(auto)
	vector() const {
		return vector_;
	}

	decltype(auto)
	transpose() const {
		return vector_;
	}

	auto
	n() const {
		return vector_.size();
	}

	decltype(auto)
	at(std::size_t const& i) {
		return vector_.at(i);
	}
	
	decltype(auto)
	at(std::size_t const& i) const {
		return vector_.at(i);
	}
	
	auto
	operator() (range<std::size_t,std::size_t> const& r) const {
		return vector_(r);
	}

private:
	rtsacv<Ring> vector_;
};

template<typename Ring>
std::ostream&
operator<< (std::ostream& os, rtsarv<Ring> const& v) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < v.n(); ++i)
        os << v.at(i) << "\t";
    os << std::endl;
    return os << '-' << std::endl;
}

template<typename Ring>
Ring
operator* (rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) {
	/* return multiply(v1,v2); */
	return 0.;
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename Ring>
struct tag_of< hbrs::mpl::rtsarv<Ring> > {
	using type = hbrs::mpl::rtsarv_tag;
};

template <>
struct make_impl<hbrs::mpl::rtsarv_tag> {
	template <typename Ring>
	static hbrs::mpl::rtsarv<Ring>
	apply(hana::basic_type<Ring>, std::size_t size) {
		return {size};
	}
	
	template <typename Ring>
	static hbrs::mpl::rtsarv<std::remove_const_t<Ring>>
	apply(std::vector<Ring> data, std::size_t size) {
		return {data, size};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_RTSARV_HPP
