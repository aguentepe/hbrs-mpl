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

#ifndef HBRS_MPL_FN_TRANSFORM_IMPL_BOOST_HANA_HPP
#define HBRS_MPL_FN_TRANSFORM_IMPL_BOOST_HANA_HPP

#include "../fwd/boost_hana.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <boost/hana/functional/demux.hpp>
#include <hbrs/mpl/core/evaluate.hpp>
#include <boost/hana/transform.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<
	typename S,
	typename F,
	typename std::enable_if_t<
		std::is_same< hana::tag_of_t<S>, hana::tuple_tag>::value
		//TODO: Add invokable check for F?
	>*
>
constexpr decltype(auto)
transform_impl_hana_tuple::operator()(S && s, F && f) const {
	return hana::transform(HBRS_MPL_FWD(s), hana::demux(evaluate)(HBRS_MPL_FWD(f)));
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FN_TRANSFORM_IMPL_BOOST_HANA_HPP
