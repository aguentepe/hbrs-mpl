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

#include <hbrs/mpl/fwd/fn/house.hpp>
#include <hbrs/mpl/dt/function.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

HBRS_MPL_NAMESPACE_BEGIN
HBRS_MPL_DEF_F1(house, house_t)
HBRS_MPL_NAMESPACE_END

#include <hbrs/mpl/fuse/fn/house.hpp>
#include <hbrs/mpl/core/implementations_of.hpp>
HBRS_MPL_MAP_IMPLS(house_t, HBRS_MPL_FUSE_FN_HOUSE_IMPLS)
