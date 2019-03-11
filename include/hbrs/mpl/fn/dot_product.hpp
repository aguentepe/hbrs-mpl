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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_DOT_PRODUCT_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_DOT_PRODUCT_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

template<typename Ring>
Ring
dot_product (rtsacv<Ring> const& v1, rtsacv<Ring> const& v2){
    Ring sum {0};
    for (std::size_t i {1}; i < v1.m(); ++i)
        sum += v1.at(i) * v2.at(i);
    return sum;
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_DOT_PRODUCT_HPP
