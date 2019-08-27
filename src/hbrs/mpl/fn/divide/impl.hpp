/* Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_DIVIDE_IMPL_HPP
#define HBRS_MPL_FN_DIVIDE_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/mpl/dt/function.hpp>

HBRS_MPL_NAMESPACE_BEGIN
HBRS_MPL_DEF_F2(divide, divide_t)
HBRS_MPL_NAMESPACE_END

#include "impl/hbrs_mpl.hpp"
#include "impl/boost_hana.hpp"
#include "impl/elemental.hpp"
#include "impl/std.hpp"

#endif // !HBRS_MPL_FN_DIVIDE_IMPL_HPP
