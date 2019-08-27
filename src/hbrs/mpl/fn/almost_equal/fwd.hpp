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

#ifndef HBRS_MPL_FN_ALMOST_EQUAL_FWD_HPP
#define HBRS_MPL_FN_ALMOST_EQUAL_FWD_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/dt/function/fwd.hpp>

HBRS_MPL_NAMESPACE_BEGIN
HBRS_MPL_DEC_F2(almost_equal, almost_equal_t)
HBRS_MPL_NAMESPACE_END

#include "fwd/std.hpp"

HBRS_MPL_MAP_IMPLS(almost_equal_t, HBRS_MPL_FN_ALMOST_EQUAL_IMPLS_STD)

#endif // !HBRS_MPL_FN_ALMOST_EQUAL_FWD_HPP
