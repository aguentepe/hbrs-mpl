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

#ifndef HBRS_MPL_DT_RTSARV_IMPL_HPP
#define HBRS_MPL_DT_RTSARV_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <vector>
#include <type_traits>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename /* type of vector entries */ Ring>
//TODO: Merge with rtsacv
struct rtsarv {
	explicit rtsarv(rtsacv<Ring> data) : vector_{data} {}
	
	explicit rtsarv(std::size_t size) : vector_(size) {}

	rtsarv(rtsarv const&) = default;
	rtsarv(rtsarv &&) = default;
	
	rtsarv&
	operator=(rtsarv const&) = default;
	rtsarv&
	operator=(rtsarv &&) = default;
	
	auto const
	length() const {
		return vector_.length();
	}

	decltype(auto)
	data() const {
		return vector_;
	}

	decltype(auto)
	transpose() const {
		return vector_;
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
		return rtsarv{vector_(r)};
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

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		(std::is_same_v< hana::tag_of_t<T1>, rtsacv_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsarv_tag    >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag    >)
	>* = nullptr
>
decltype(auto)
operator*(T1 && t1, T2 && t2) {
	return multiply(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
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

#endif // !HBRS_MPL_DT_RTSARV_IMPL_HPP
