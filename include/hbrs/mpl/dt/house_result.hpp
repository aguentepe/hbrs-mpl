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

#ifndef HBRS_MPL_DT_HOUSE_RESULT_HPP
#define HBRS_MPL_DT_HOUSE_RESULT_HPP

#include <hbrs/mpl/fwd/dt/house_result.hpp>
#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN

template<typename Ni, typename Beta>
struct house_result {
	template<
		typename Ni_ = Ni,
		typename Beta_ = Beta,
		typename std::enable_if_t<
			std::is_default_constructible<Ni_>::value &&
			std::is_default_constructible<Beta_>::value
		>* = nullptr
	>
	constexpr
	house_result() {}
	
	template<typename Ni_, typename Beta_>
	constexpr 
	house_result(Ni_ && ni, Beta_ && beta) 
	: ni_{HBRS_MPL_FWD(ni)}, beta_{HBRS_MPL_FWD(beta)}
	{}
	
	constexpr 
	house_result(house_result const&) = default;
	constexpr 
	house_result(house_result &&) = default;
	
	constexpr house_result&
	operator=(house_result const&) = default;
	constexpr house_result&
	operator=(house_result &&) = default;
	
	constexpr decltype(auto)
	ni() & { return (ni_); };
	
	constexpr decltype(auto)
	ni() const& { return (ni_); };
	
	constexpr decltype(auto)
	ni() && { return HBRS_MPL_FWD(ni_); };
	
	constexpr decltype(auto)
	beta() & { return (beta_); };
	
	constexpr decltype(auto)
	beta() const& { return (beta_); };
	
	constexpr decltype(auto)
	beta() && { return HBRS_MPL_FWD(beta_); };

private:
	Ni ni_;
	Beta beta_;
};

struct house_ni{};
struct house_beta{};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename Ni, typename Beta>
struct tag_of< hbrs::mpl::house_result<Ni, Beta> > {
	using type = hbrs::mpl::house_result_tag;
};

template <>
struct make_impl<hbrs::mpl::house_result_tag> {
	template <typename Ni, typename Beta>
	static constexpr hbrs::mpl::house_result<
		std::decay_t<Ni>, 
		std::decay_t<Beta>
	>
	apply(Ni && ni, Beta && beta) {
		return {
			HBRS_MPL_FWD(ni),
			HBRS_MPL_FWD(beta)
		};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_HOUSE_RESULT_HPP
