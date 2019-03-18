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

#ifndef HBRS_MPL_FWD_DT_SSM_HPP
#define HBRS_MPL_FWD_DT_SSM_HPP

#include <hbrs/mpl/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

/* sequential sub matrix */
template<typename Matrix, typename Index1, typename Index2>
struct ssm;
struct ssm_tag{};
constexpr auto make_ssm = hana::make<ssm_tag>;
constexpr auto to_ssm = hana::to<ssm_tag>;

HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FWD_DT_SSM_HPP
