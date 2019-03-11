
/* Copyright (c) 2016 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_DT_SVD_RESULT_HPP
#define HBRS_MPL_DT_SVD_RESULT_HPP

#include <hbrs/mpl/fwd/dt/svd_result.hpp>
#include <hbrs/mpl/dt/bidiag_result.hpp>
#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN

template<typename U, typename S, typename V>
struct svd_result {
	template<
		typename U_ = U,
		typename S_ = S,
		typename V_ = V,
		typename std::enable_if_t<
			std::is_default_constructible<U_>::value &&
			std::is_default_constructible<S_>::value &&
			std::is_default_constructible<V_>::value
		>* = nullptr
	>
	constexpr
	svd_result() {}
	
	template<typename U_, typename S_, typename V_>
	constexpr 
	svd_result(U_ && u, S_ && s, V_ && v) 
	: u_{HBRS_MPL_FWD(u)}, s_{HBRS_MPL_FWD(s)}, v_{HBRS_MPL_FWD(v)}
	{}
	
	constexpr 
	svd_result(svd_result const&) = default;
	constexpr 
	svd_result(svd_result &&) = default;
	
	constexpr svd_result&
	operator=(svd_result const&) = default;
	constexpr svd_result&
	operator=(svd_result &&) = default;
	
	constexpr decltype(auto)
	u() & { return (u_); };
	
	constexpr decltype(auto)
	u() const& { return (u_); };
	
	constexpr decltype(auto)
	u() && { return HBRS_MPL_FWD(u_); };
	
	constexpr decltype(auto)
	s() & { return (s_); };
	
	constexpr decltype(auto)
	s() const& { return (s_); };
	
	constexpr decltype(auto)
	s() && { return HBRS_MPL_FWD(s_); };
	
	constexpr decltype(auto)
	v() & { return (v_); };
	
	constexpr decltype(auto)
	v() const& { return (v_); };
	
	constexpr decltype(auto)
	v() && { return HBRS_MPL_FWD(v_); };
	
private:
	U u_;
	S s_;
	V v_;
};

struct svd_u{};
struct svd_s{};
struct svd_v{};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename U, typename S, typename V>
struct tag_of< hbrs::mpl::svd_result<U, S, V> > {
	using type = hbrs::mpl::svd_result_tag;
};

template <>
struct make_impl<hbrs::mpl::svd_result_tag> {
	template <typename U, typename S, typename V>
	static constexpr hbrs::mpl::svd_result<
		std::decay_t<U>, 
		std::decay_t<S>, 
		std::decay_t<V>
	>
	apply(U && u, S && s, V && v) {
		return {
			HBRS_MPL_FWD(u),
			HBRS_MPL_FWD(s),
			HBRS_MPL_FWD(v)
		};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_SVD_RESULT_HPP
