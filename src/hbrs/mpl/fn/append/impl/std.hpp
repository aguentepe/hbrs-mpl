/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_APPEND_IMPL_STD_HPP
#define HBRS_MPL_FN_APPEND_IMPL_STD_HPP

#include "../fwd/std.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/append.hpp>
#include <boost/hana/ext/std/vector.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/mpl/if.hpp>
#include <vector>
#include <numeric>
#include <tuple>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<
	typename S,
	typename E,
	typename std::enable_if_t<
		boost::mpl::if_c<
			std::is_same< hana::tag_of_t<S>, hana::ext::std::vector_tag>::value,
			std::is_convertible<E&&, typename std::remove_reference_t<S>::value_type>,
			std::false_type
		>::type::value
	>*
>
constexpr decltype(auto)
append_impl_std_vector::operator()(S && s, E && e) const {
	s.push_back(HBRS_MPL_FWD(e));
	return HBRS_MPL_FWD(s);
}

template<
	typename S,
	typename E,
	typename std::enable_if_t< 
		boost::mpl::if_c<
			std::is_same< hana::tag_of_t<S>, hana::ext::std::vector_tag>::value,
			std::is_convertible<E&&, typename std::remove_reference_t<S const>::value_type>,
			std::false_type
		>::type::value
	>*
>
constexpr decltype(auto)
append_impl_std_vector::operator()(S const& s, E && e) const {
	S  v(s.size() + 1ul);
	std::copy (s.begin(), s.end(), v.begin());
	v[s.size()] = e;
	
	return v;
}


template<
	typename S, 
	typename E, 
	typename std::enable_if_t< std::is_same< hana::tag_of_t<S>, hana::ext::std::tuple_tag >::value >*
>
constexpr decltype(auto)
append_impl_std_tuple::operator()(S && s, E && e) const {
	return hana::append(HBRS_MPL_FWD(s), HBRS_MPL_FWD(e));
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FN_APPEND_IMPL_STD_HPP
