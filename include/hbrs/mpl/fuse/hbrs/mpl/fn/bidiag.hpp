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

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_BIDIAG_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_BIDIAG_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/bidiag_result.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/house.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

/*
 * Algorithm 5.4.2 (Householder Bidiagonalization) on page 284
 * Given the real m-by-n matrix A with m>=n, the following algorithm
 * creates B,U' and V with U'AV=B where B is upper bidiagonal and U'
 * and V are products of householder matrices.
 *
 * In comparison with the book:
 *      In the book the Algorithm directly overwrites A. Instead we
 *      make a copy of A, and call it B, and we overwrite B. We don't
 *      overwrite A. And instead of storing U' and V together with B
 *      inside A, we return a type that contains U', B and V called
 *      BidiagResult.
 */
struct bidiag_impl_householder {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& A, int econ) {
			BOOST_ASSERT(A.m() >= A.n());

			// Copy A to result.b()
			auto result {make_bidiag_result(identity<Ring,Order>(A.m()), A, identity<Ring,Order>(A.n()))};

			{ /* This block is here so the reference A in the next line can be
				created for readability */
				auto& U {result.u()};
				auto& A {result.b()};
				auto& V {result.v()};

				auto const m {A.m()};
				auto const n {A.n()};

				range<std::size_t,std::size_t> const rm {std::size_t{0}, m-1};

				for (std::size_t j {0}; j < n - 1; ++j) {
					range<std::size_t,std::size_t> const jm     {j             , m-1};
					range<std::size_t,std::size_t> const jn     {j	           , n-1};
					range<std::size_t,std::size_t> const j1n    {j + 1         , n-1};
					range<std::size_t,std::size_t> const column {std::size_t{0}, m-1}; // a full column

					//Use Algorithm 5.1.1 to compute the householder vector.
					auto const h {house(A(jm, j))};
					auto& ni {h.ni()};
					auto& beta {h.beta()};

					/*
					 * Here the first commented out line is the one on page 285
					 * in Algorithm 5.4.2. But as explained on page 236 in
					 * Chapter 5.1.4 it should never be implemented that way.
					 * So we implemented the way it is suggested on that page.
					 */

					/* Ajmjn is equivalent to A(j:m,j:n) in the book */
					auto const Ajmjn {A(jm, jn)};
					/* overwrite(A, jm, jn, (identity<Ring,Order>(m - j) - beta * ni * transpose(ni)) * Ajmjn); // mathematical notation */
					overwrite(A, jm, jn, Ajmjn -  operator*<Ring,Order>((beta * ni), (transpose(ni) * Ajmjn)));
					/*
					 * In the book here the householder vector would be saved
					 * inside of A. But since we compute and return U we don't
					 * do that.
					 */

					auto Ui {identity<Ring,Order>(m)};
					overwrite(Ui, jm, jm, identity<Ring,Order>(m - j) - beta * operator*<Ring,Order>(ni, transpose(ni)));
					U = Ui * U;

					if (j + 1 <= n - 2) {
						auto const h {house(transpose(A(j, j1n)))};
						auto& ni {h.ni()};
						auto& beta {h.beta()};

						auto const Ajmj1n {A(jm, j1n)}; // equivalent to A(j:m,j+1:n)
						/* overwrite(A, jm, j1n, Ajmj1n * (identity<Ring,Order>(n-1 - j) - beta * ni * transpose(ni))); // mathematical notation */
						overwrite(A, jm, j1n, Ajmj1n - operator*<Ring,Order>((Ajmj1n * ni), transpose(beta * ni)));
						/*
						 * In the book here the householder vector would be saved
						 * inside of A. But since we compute and return V we don't
						 * do that.
						 */

						auto const Vjmj1n {V(jn, j1n)};
						overwrite(V, jn, j1n, Vjmj1n - operator*<Ring,Order>((Vjmj1n * ni), transpose(beta * ni)));
					}
				}
				U = transpose(U);
			}
			return result;
	}

};
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_BIDIAG_IMPLS boost::hana::make_tuple(                                             \
		hbrs::mpl::detail::bidiag_impl_householder{}                                                                       \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_BIDIAG_HPP
