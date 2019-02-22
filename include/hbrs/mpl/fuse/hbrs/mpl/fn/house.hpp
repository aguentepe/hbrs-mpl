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

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/CVector.hpp>
#include <hbrs/mpl/dt/RVector.hpp>
#include <hbrs/mpl/fn/house.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/integral_constant.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
using namespace hana::literals;
namespace detail {

struct house_impl {
        /* constexpr */ 
        decltype(auto)
        operator()(CVector<double> const& x) {
                auto const m {x.m()}; // copy for readability

                /* x2m is a temporary which is written x(2:m) in the book and is equivalent to x(Range(1,m-1)) in this code */
                auto const x2m {x(Range(1, m-1))};
                auto const sigma {transpose(x2m) * x2m};

                /* the tuple 'result' stores the vector 'ni' and the double 'beta' */
                /* The vector ni is the vector x with the value 1 in its first row */
                hana::tuple<CVector<double>, double> result {{x}, {0.}};
                result[0_c].at(0) = 1;

                if (sigma == 0 && x.at(0) >= 0)
                        result[1_c] = 0;
                else if (sigma == 0 && x.at(0) < 0)
                        /* In the book beta is set to -2 but in the Errata it says that it should be +2 */
                        result[1_c] = 2;
                else {
                        auto const mi {std::sqrt(x.at(0) * x.at(0) + sigma)};
                        if (x.at(0) <= 0)
                            result[0_c].at(0) = x.at(0) - mi;
                        else
                            result[0_c].at(0) = -sigma / (x.at(0) + mi);
                        auto const nisq = result[0_c].at(0) * result[0_c].at(0); // sqare of first element of ni
                        result[1_c] = 2 * nisq / (sigma + nisq);
                        result[0_c] = result[0_c] / result[0_c].at(0);
                }
                return result;
        }
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_HOUSE_IMPLS boost::hana::make_tuple(                                             \
		hbrs::mpl::detail::house_impl{}                                                                       \
	)
