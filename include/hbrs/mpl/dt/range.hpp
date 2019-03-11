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

#pragma once

#ifndef HBRS_MPL_DT_RANGE_HPP
#define HBRS_MPL_DT_RANGE_HPP

#include <hbrs/mpl/fwd/dt/range.hpp>
#include <hbrs/mpl/detail/axis.hpp>

HBRS_MPL_DEFINE_AXIS_2D(range, Range, first, First, last, Last)

/* #include <hbrs/mpl/fn/equal.hpp> */
#include <hbrs/mpl/fn/and.hpp>
#include <hbrs/mpl/fn/not.hpp>
HBRS_MPL_DEFINE_AXIS_2D_EQUALITY_OPERATORS(range, Range, first, First, last, Last)

#endif // !HBRS_MPL_DT_RANGE_HPP
