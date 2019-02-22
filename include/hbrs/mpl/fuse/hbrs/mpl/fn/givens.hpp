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
#include <hbrs/mpl/dt/srv.hpp>
#include <hbrs/mpl/dt/CVector.hpp>
#include <hbrs/mpl/dt/scv.hpp>
#include <hbrs/mpl/fn/givens.hpp>
#include <boost/hana/tuple.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct givens_impl {
        /* template <typename Vector> */
        /* constexpr */
        decltype(auto)
        operator()(double const a, double const b) const {
                CVector<double> cs(2);
                if (b == 0) {
                        cs.at(0) = 1;
                        cs.at(1) = 0;
                } else {
                        if (std::abs(b) > std::abs(a)) {
                                auto const tau = -a / b;
                                cs.at(1) = 1 / std::sqrt(1 + tau * tau);
                                cs.at(0) = cs.at(1) * tau;
                        } else {
                                auto const tau = -b / a;
                                cs.at(0) = 1 / std::sqrt(1 + tau * tau);
                                cs.at(1) = cs.at(0) * tau;
                        }
                }
                return cs;
        }
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_GIVENS_IMPLS boost::hana::make_tuple(                                             \
		hbrs::mpl::detail::givens_impl{}                                                                       \
	)
