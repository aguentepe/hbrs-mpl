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
#include <hbrs/mpl/dt/Matrix.hpp>
#include <hbrs/mpl/fn/house.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/integral_constant.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
using namespace hana::literals;
namespace detail {

struct bidiag_impl_householder {
        /* constexpr */ 
        decltype(auto)
        operator()(Matrix const& A) {
                BOOST_ASSERT(A.m() >= A.n());

                // Copy A to result.B
                /* BidiagResult result {A}; */
                hana::tuple<Matrix, Matrix, Matrix> result{{identity(A.m())}, {A}, {identity(A.n())}};

                { /* This block is here so the reference A in the next line can be
                    created for readability */
                    auto& U {result[0_c]};
                    auto& A {result[1_c]};
                    auto& V {result[2_c]};

                    auto const m {A.m()};
                    auto const n {A.n()};

                    Range const rm {0, m-1};

                    for (std::size_t j {0}; j < n - 1; ++j) {
                        Range const jm     {j    , m-1};
                        Range const jn     {j    , n-1};
                        Range const j1n    {j + 1, n-1};
                        Range const column {0    , m-1}; // a full column

                        //Use Algorithm 5.1.1 to compute the householder vector.
                        auto const h {house(A(jm, j))};
                        auto& ni {h[0_c]};
                        auto& beta {h[1_c]};

                        /*
                         * Here the first commented out line is the one on page 285
                         * in Algorithm 5.4.2. But as explained on page 236 in
                         * Chapter 5.1.4 it should never be implemented that way.
                         * So we implemented the way it is suggested on that page.
                         */

                        /* Ajmjn is equivalent to A(j:m,j:n) in the book */
                        auto const Ajmjn {A(jm, jn)};
                        /* overwrite(A, jm, jn, (identity(m - j) - beta * ni * transpose(ni)) * Ajmjn); // mathematical notation */
                        overwrite(A, jm, jn, Ajmjn - (beta * ni) * (transpose(ni) * Ajmjn));
                        /*
                         * In the book here the householder vector would be saved
                         * inside of A. But since we compute and return U we don't
                         * do that.
                         */

                        auto Ui {identity(m)};
                        overwrite(Ui, jm, jm, identity(m - j) - beta * (ni * transpose(ni)));
                        U = Ui * U;

                        if (j + 1 <= n - 2) {
                            auto const h {house(transpose(A(j, j1n)))};
                            auto& ni {h[0_c]};
                            auto& beta {h[1_c]};

                            auto const Ajmj1n {A(jm, j1n)}; // equivalent to A(j:m,j+1:n)
                            /* overwrite(A, jm, j1n, Ajmj1n * (identity(n-1 - j) - beta * ni * transpose(ni))); // mathematical notation */
                            overwrite(A, jm, j1n, Ajmj1n - (Ajmj1n * ni) * transpose(beta * ni));
                            /*
                             * In the book here the householder vector would be saved
                             * inside of A. But since we compute and return V we don't
                             * do that.
                             */

                            auto const Vjmj1n {V(jn, j1n)};
                            overwrite(V, jn, j1n, Vjmj1n - (Vjmj1n * ni) * transpose(beta * ni));
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
