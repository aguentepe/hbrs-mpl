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

#ifndef HBRS_MPL_FWD_DT_GIVENS_ROTATION_HPP
#define HBRS_MPL_FWD_DT_GIVENS_ROTATION_HPP

#include <hbrs/mpl/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>
#include <hbrs/mpl/fwd/dt/storage_order.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<typename /* type of matrix entries */ Ring>
struct givens_rotation;
struct givens_rotation_tag{};
constexpr auto make_givens_rotation = hana::make<givens_rotation_tag>;
constexpr auto to_givens_rotation = hana::to<givens_rotation_tag>;

HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FWD_DT_GIVENS_ROTATION_HPP
