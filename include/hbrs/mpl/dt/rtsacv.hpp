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
#include <hbrs/mpl/dt/range.hpp>
#include <vector>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename /* type of vector entries */ Ring>
struct rtsacv {
	explicit rtsacv(std::vector<Ring> data) : data_{data} {}
	
	explicit rtsacv(std::size_t size) : data_(size, Ring{0}) {
		BOOST_ASSERT(size >= 0);
	}

	explicit rtsacv(std::initializer_list<Ring> const& l);
	
	rtsacv(rtsacv const&) = default;
	rtsacv(rtsacv &&) = default;
	
	rtsacv&
	operator=(rtsacv const&) = default;
	rtsacv&
	operator=(rtsacv &&) = default;
	
	auto
	size() const {
		return data_.size();
	}

	decltype(auto)
	data() const {
		return data_;
	}

	decltype(auto)
	at(std::size_t const& i) {
		return data_[i];
	}
	
	decltype(auto)
	at(std::size_t const& i) const {
		return data_[i];
	}

	auto
	operator() (range<std::size_t,std::size_t> const& r) const {
		rtsacv v (r.last() - r.first() + 1);
		for (std::size_t i = 0; i < v.size(); ++i) {
			v.at(i) = at(i + r.first());
		}
		return v;
	}

private:
	std::vector<Ring> data_;
};

template<typename Ring>
std::ostream& operator<< (std::ostream& os, rtsacv<Ring> const& v) {
    os << '-' << std::endl;
    for (std::size_t i {0}; i < v.size(); ++i)
        os << v.at(i) << std::endl;
    return os << '-' << std::endl;
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
