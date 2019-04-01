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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_HOUSE_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_HOUSE_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/house_result.hpp>
#include <hbrs/mpl/fn/house.hpp>
#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/fn/divide.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

/*
 * Algorithm 5.1.1 (Householder Vector) on page 236
 * Given the real column vector x, this function computes the real
 * column vector v with v(0)=1 and the real number beta.
 */
struct house_impl {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsacv<Ring> const& x) {
		auto const m {x.size()}; // copy for readability

		/* x2m is a temporary which is written x(2:m) in the book and is equivalent to x(range(1,m-1)) in this code */
		auto const x2m {x(range<std::size_t,std::size_t>(std::size_t{1}, m-1))};
		auto const sigma {transpose(x2m) * x2m};

		/* The vector ni is the vector x with the value 1 in its first row */
		auto result{make_house_result(x, 0.)};
		auto& ni  {result.ni  ()};
		auto& beta{result.beta()};

		ni.at(0) = 1;

		if (sigma == 0 && x.at(0) >= 0) {
			beta = 0;
		} else if (sigma == 0 && x.at(0) < 0) {
			/* In the book beta is set to -2 but in the Errata it says that it should be +2 */
			beta = 2;
		} else {
			auto const mi {std::sqrt(x.at(0) * x.at(0) + sigma)};
			if (x.at(0) <= 0) {
				ni.at(0) = x.at(0) - mi;
			} else {
				ni.at(0) = -sigma / (x.at(0) + mi);
			}
			auto const nisq = ni.at(0) * ni.at(0); // sqare of first element of ni
			beta = 2 * nisq / (sigma + nisq);
			ni = ni / ni.at(0);
		}
		return result;

	}
};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_HOUSE_IMPLS boost::hana::make_tuple(                                                 \
		hbrs::mpl::detail::house_impl{}                                                                                \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_HOUSE_HPP
