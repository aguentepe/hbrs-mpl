/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_TRANSFORM_FWD_STD_HPP
#define HBRS_MPL_FN_TRANSFORM_FWD_STD_HPP

#include <hbrs/mpl/config.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/ext/std/vector.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <array>
#include <vector>
#include <tuple>
#include <utility>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct transform_impl_std_array {
	template <
		typename T, std::size_t N,
		typename F
		//TODO: Add invokable check for F?
	>
	constexpr decltype(auto)
	operator()(std::array<T,N> const& s, F && f) const;
	
	template <
		typename T, std::size_t N,
		typename F
		//TODO: Add invokable check for F?
	>
	constexpr decltype(auto)
	operator()(std::array<T,N> & s, F && f) const;
};

struct transform_impl_std_vector {	
	template<
		typename T, typename Allocator,
		typename F
		//TODO: Add invokable check for F?
	>
	constexpr auto
	operator()(std::vector<T, Allocator> const& s, F && f) const;
};

struct transform_impl_std_tuple {
	template<
		typename S,
		typename F,
		typename std::enable_if_t<
			std::is_same< hana::tag_of_t<S>, hana::ext::std::tuple_tag>::value
			//TODO: Add invokable check for F?
		>* = nullptr
	>
	constexpr decltype(auto)
	operator()(S && s, F && f) const;
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_TRANSFORM_IMPLS_STD boost::hana::make_tuple(                                                       \
		hbrs::mpl::detail::transform_impl_std_array{},                                                                 \
		hbrs::mpl::detail::transform_impl_std_vector{},                                                                \
		hbrs::mpl::detail::transform_impl_std_tuple{}                                                                  \
	)

#endif // !HBRS_MPL_FN_TRANSFORM_FWD_STD_HPP
