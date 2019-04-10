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

#ifndef HBRS_MPL_FUSE_STD_FN_ALMOST_EQUAL_HPP
#define HBRS_MPL_FUSE_STD_FN_ALMOST_EQUAL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/fn/almost_equal.hpp>
#include <limits>
#include <cmath>

#include <boost/hana/tuple.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

/*
 * Read about ULP here:
 * https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
struct almost_equal_impl_double {
	bool
	operator()(double const a, double const b) const {
		static int const max_ulps = 10;

		if (almost_zero(a)) {
			if (almost_zero(b)) {
				return true;
			} else {
				return false;
			}
		} else if (almost_zero(b)) {
			return false;
		}

		union {
			double ad;
			long   al;
		};
		union {
			double bd;
			long   bl;
		};
		ad = a;
		bd = b;

		if ((al < 0) != (bl < 0)) {
			if (a == b) {
				return true;
			} else {
				return false;
			}
		}

		auto ulps_diff{std::abs(al - bl)};
		
		if (ulps_diff <= max_ulps) {
			return true;
		} else {
			return false;
		}
	}

private:
	bool
	almost_zero(double const x) const {
		static int const precision = 2;

		return std::abs(x) <= std::numeric_limits<double>::epsilon() * precision;
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_STD_FN_ALMOST_EQUAL_IMPLS boost::hana::make_tuple(                                               \
		hbrs::mpl::detail::almost_equal_impl_double{}                                                                  \
	)

#endif // !HBRS_MPL_FUSE_STD_FN_ALMOST_EQUAL_HPP
