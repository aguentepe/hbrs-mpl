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

#ifndef HBRS_MPL_FN_ZIP_IMPL_STD_HPP
#define HBRS_MPL_FN_ZIP_IMPL_STD_HPP

#include "../fwd/std.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/core/evaluate.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/mpl/dt/zas.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <array>
#include <vector>
#include <tuple>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<
	typename S1,
	typename S2,
	typename std::enable_if_t<
		(
			std::is_same< hana::tag_of_t<S1>, hana::ext::std::array_tag>::value ||
			std::is_same< hana::tag_of_t<S1>, hana::ext::std::vector_tag>::value ||
			std::is_same< hana::tag_of_t<S1>, hana::ext::boost::integer_range_tag>::value
		) &&
		(
			std::is_same< hana::tag_of_t<S2>, hana::ext::std::array_tag>::value ||
			std::is_same< hana::tag_of_t<S2>, hana::ext::std::vector_tag>::value ||
			std::is_same< hana::tag_of_t<S2>, hana::ext::boost::integer_range_tag>::value
		)
	>*
>
constexpr decltype(auto)
zip_impl_std_array_vector_irange::operator()(S1&& s1, S2&& s2) const {
	auto s1sz = (*size)(s1);
	auto s2sz = (*size)(s2);
	
	if (s1sz != s2sz) {
		BOOST_THROW_EXCEPTION((
			incompatible_sequences_exception{} 
			<< errinfo_sequences_sizes{ {s1sz, s2sz} }
		));
	}
	
	return make_zas(HBRS_MPL_FWD(s1), HBRS_MPL_FWD(s2));
}

template<
	typename S1,
	typename S2,
	typename std::enable_if_t<
		std::is_same< hana::tag_of_t<S1>, hana::ext::std::vector_tag>::value &&
		std::is_same< hana::tag_of_t<S2>, hana::ext::std::vector_tag>::value
	>*
>
constexpr auto
zip_impl_std_tuple_vector::operator()(S1&& s1, S2&& s2) const {
	auto s1sz = (*size)(s1);
	auto s2sz = (*size)(s2);
	
	if (s1sz != s2sz) {
		BOOST_THROW_EXCEPTION((
			incompatible_sequences_exception{} 
			<< errinfo_sequences_sizes{ {s1sz, s2sz} }
		));
	}
	
	typedef typename std::remove_reference_t<S1>::value_type T1;
	typedef typename std::remove_reference_t<S2>::value_type T2;
	std::vector<std::tuple<T1, T2>> zipped;
	zipped.reserve(s1.size());
	
	for(std::size_t i = 0; i < s1.size(); ++i) {
		zipped.push_back({
			HBRS_MPL_FWD(s1)[i],
			HBRS_MPL_FWD(s2)[i]
		});
	}
	
	return zipped;
}

template<
	typename S1,
	typename S2,
	typename std::enable_if_t<
		std::is_same< hana::tag_of_t<S1>, hana::ext::std::tuple_tag>::value &&
		std::is_same< hana::tag_of_t<S2>, hana::ext::std::tuple_tag>::value &&
		std::tuple_size<std::remove_reference_t<S1>>::value == std::tuple_size<std::remove_reference_t<S2>>::value
	>*
>
constexpr decltype(auto)
zip_impl_std_tuple::operator()(S1&& s1, S2&& s2) const {
	return make_zas(HBRS_MPL_FWD(s1), HBRS_MPL_FWD(s2));
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FN_ZIP_IMPL_STD_HPP
