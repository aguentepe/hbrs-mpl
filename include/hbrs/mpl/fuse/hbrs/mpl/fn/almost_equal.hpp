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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_ALMOST_EQUAL_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_ALMOST_EQUAL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/fn/almost_equal.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>

#include <boost/hana/tuple.hpp>
#include <cstdio>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct almost_equal_impl_rtsam {

	/* constexpr bool */ 
	bool
	operator()(rtsam<double,storage_order::row_major> const& M1, rtsam<double,storage_order::row_major> const& M2) const {
		/* return impl(M1,M2); */
		if (M1.m() != M2.m() || M1.n() != M2.n()) {
			return false;
		}
		for (std::size_t i {0}; i < M1.m(); ++i) {
			for (std::size_t j {0}; j < M1.n(); ++j) {
				//TODO: almost_equal_double durch almost_equal ersetzen
				if (!almost_equal_double(M1.at(make_matrix_index(i, j)), M2.at(make_matrix_index(i, j)))) {
					return false;
				}
			}
		}
		return true;
	}

	/* constexpr bool */ 
	/* operator()(rtsam<double,storage_order::column_major> const& M1, rtsam<double,storage_order::column_major> const& M2) const { */
	/* 	/1* return impl(M1,M2); *1/ */
	/* 	if (M1.m() != M2.m() || M1.n() != M2.n()) { */
	/* 		return false; */
	/* 	} */
	/* 	for (std::size_t i {0}; i < M1.m(); ++i) { */
	/* 		for (std::size_t j {0}; j < M1.n(); ++j) { */
	/* 			if (!almost_equal(M1.at(make_matrix_index(i, j)), M2.at(make_matrix_index(i, j)))) { */
	/* 				return false; */
	/* 			} */
	/* 		} */
	/* 	} */
	/* 	return true; */
	/* } */

/* private: */
/* 	template<typename Matrix> */
/* 	constexpr bool */
/* 	impl(Matrix const& M1, Matrix const& M2) const { */
/* 		if (M1.m() != M2.m() || M1.n() != M2.n()) { */
/* 			return false; */
/* 		} */
/* 		for (std::size_t i {0}; i < M1.m(); ++i) { */
/* 			for (std::size_t j {0}; j < M1.n(); ++j) { */
/* 				if (!almost_equal(M1.at(make_matrix_index(i, j)), M2.at(make_matrix_index(i, j)))) { */
/* 					return false; */
/* 				} */
/* 			} */
/* 		} */
/* 		return true; */
/* 	} */

/*
 * These two function below belong in the std subdirectory.
 * For the function above to be able to use it they're here as workaround.
 */
private:
	bool
	almost_equal_double(double const a, double const b) const {
		static int const max_ulps = 147;

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
		std::printf("ulp: %ld ; a: %.96lf b: %.96lf\n", ulps_diff, a, b);
		
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

#define HBRS_MPL_FUSE_HBRS_MPL_FN_ALMOST_EQUAL_IMPLS boost::hana::make_tuple(                                                 \
		hbrs::mpl::detail::almost_equal_impl_rtsam{}                                                                          \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_ALMOST_EQUAL_HPP
