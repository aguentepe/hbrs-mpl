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

#ifndef HBRS_MPL_FN_MULTIPLY_FWD_HPP
#define HBRS_MPL_FN_MULTIPLY_FWD_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/dt/function/fwd.hpp>
#include <hbrs/mpl/core/implementations_of.hpp>

HBRS_MPL_NAMESPACE_BEGIN
HBRS_MPL_DEC_F2(multiply, multiply_t)
HBRS_MPL_NAMESPACE_END

#include "fwd/hbrs_mpl.hpp"
#include "fwd/boost_hana.hpp"
#include "fwd/elemental.hpp"
#include "fwd/matlab.hpp"
#include "fwd/std.hpp"

HBRS_MPL_MAP_IMPLS(multiply_t, HBRS_MPL_FN_MULTIPLY_IMPLS_HBRS_MPL, HBRS_MPL_FN_MULTIPLY_IMPLS_BOOST_HANA, HBRS_MPL_FN_MULTIPLY_IMPLS_ELEMENTAL, HBRS_MPL_FN_MULTIPLY_IMPLS_MATLAB, HBRS_MPL_FN_MULTIPLY_IMPLS_STD)

#endif // !HBRS_MPL_FN_MULTIPLY_FWD_HPP
