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

#ifndef HBRS_MPL_FN_ALMOST_EQUAL_FWD_STD_HPP
#define HBRS_MPL_FN_ALMOST_EQUAL_FWD_STD_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>

#include <boost/hana/tuple.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

struct almost_equal_impl_double {
	bool
	operator()(double const a, double const b) const;

private:
	bool
	almost_zero(double const x) const;
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_ALMOST_EQUAL_IMPLS_STD boost::hana::make_tuple(                                                    \
		hbrs::mpl::detail::almost_equal_impl_double{}                                                                  \
	)

#endif // !HBRS_MPL_FN_ALMOST_EQUAL_FWD_STD_HPP
