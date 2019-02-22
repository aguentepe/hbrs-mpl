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


#define BOOST_TEST_MODULE bidiag_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/dt/zas.hpp>
#include <hbrs/mpl/dt/rtsav.hpp>
#include <hbrs/mpl/fn/bidiag.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <boost/hana/integral_constant.hpp>

#include <boost/hana/tuple.hpp>
#include <boost/hana/first.hpp>
#include <boost/hana/second.hpp>
#include <array>
#include <tuple>
#include <vector>

BOOST_AUTO_TEST_SUITE(bidiag_test)

BOOST_AUTO_TEST_CASE(bidiag_test1) {
        using namespace hbrs::mpl;

        Matrix A (4, 4);
        A.at(0, 0) = 16;
        A.at(0, 1) = 2;
        A.at(0, 2) = 3;
        A.at(0, 3) = 13;

        A.at(1, 0) = 5;
        A.at(1, 1) = 11;
        A.at(1, 2) = 10;
        A.at(1, 3) = 8;

        A.at(2, 0) = 9;
        A.at(2, 1) = 7;
        A.at(2, 2) = 6;
        A.at(2, 3) = 12;

        A.at(3, 0) = 4;
        A.at(3, 1) = 14;
        A.at(3, 2) = 15;
        A.at(3, 3) = 1;
        
        Matrix A2 {{ 2,  2,3,
                     9,  8,1,
                    15,100,7,
                    99,  1,2,
                     5,  7,3},5};

        auto B {bidiag(A)};
        auto B2 {bidiag(A)};
        auto B3 {bidiag(A2)};

        auto C { B[0_c] * B[1_c] * transpose(B[2_c]) };
        auto C2 { B2[0_c] * B2[1_c] * transpose(B2[2_c]) };
        auto C3 { B3[0_c] * B3[1_c] * transpose(B3[2_c]) };
        
        BOOST_TEST(C == A);
        BOOST_TEST(C2 == A);
        BOOST_TEST(C3 == A2);
}

BOOST_AUTO_TEST_SUITE_END()
