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

#ifndef HBRS_MPL_DT_RTSACV_IMPL_HPP
#define HBRS_MPL_DT_RTSACV_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/fn/divide.hpp>
#include <vector>
#include <type_traits>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename /* type of vector entries */ Ring>
//TODO: Merge with rtsarv
struct rtsacv {
	explicit rtsacv(std::vector<Ring> data) : data_{data} {}
	
	rtsacv(std::size_t size) : data_(size, Ring{0}) {}

	rtsacv(rtsacv const&) = default;
	rtsacv(rtsacv &&) = default;
	
	rtsacv&
	operator=(rtsacv const&) = default;
	rtsacv&
	operator=(rtsacv &&) = default;
	
	auto
	length() const {
		return data_.size();
	}

	decltype(auto)
	data() const {
		return data_;
	}

	decltype(auto)
	at(std::size_t const i) {
		return data_.at(i);
	}
	
	decltype(auto)
	at(std::size_t const i) const {
		return data_.at(i);
	}

	auto
	operator() (range<std::size_t,std::size_t> const& r) const {
		rtsacv<Ring> v (r.last() - r.first() + 1);
		for (std::size_t i{0}; i < v.length(); ++i) {
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
    for (std::size_t i {0}; i < v.length(); ++i)
        os << v.at(i) << std::endl;
    return os << '-' << std::endl;
}

template<typename Ring>
decltype(auto)
operator*(Ring const& s, rtsacv<Ring> const& v) {
	return multiply(v,s);
}

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		std::is_same_v< hana::tag_of_t<T1>, submatrix_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag >
	>* = nullptr
>
decltype(auto)
operator*(T1 && t1, T2 && t2) {
	return multiply(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

template<typename Ring>
decltype(auto)
operator/(rtsacv<Ring> const& v, Ring const& d) {
	return divide(v,d);
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

#endif // !HBRS_MPL_DT_RTSACV_IMPL_HPP
