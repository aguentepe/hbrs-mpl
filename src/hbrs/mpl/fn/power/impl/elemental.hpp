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

#ifndef HBRS_MPL_FN_POWER_IMPL_ELEMENTAL_HPP
#define HBRS_MPL_FN_POWER_IMPL_ELEMENTAL_HPP

#include "../fwd/elemental.hpp"
#ifdef HBRS_MPL_ENABLE_ELEMENTAL

#include <El.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

template <
	typename F,
	typename T,
	typename std::enable_if_t< 
		pow_is_invokable_trait<F const& , T const&>::value
	>*
>
auto
power_impl_el::operator()(F const& a, T const& b) const {
	return El::Pow(a, b);
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_ENABLE_ELEMENTAL
#endif // !HBRS_MPL_FN_POWER_IMPL_ELEMENTAL_HPP
