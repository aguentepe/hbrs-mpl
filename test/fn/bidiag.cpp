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

#include <hbrs/mpl/fn/bidiag.hpp>
#include <hbrs/mpl/fn/equal.hpp>

BOOST_AUTO_TEST_SUITE(bidiag_test)

BOOST_AUTO_TEST_CASE(bidiag_test1) {
        using namespace hbrs::mpl;

        rtsam<double, storage_order::row_major> A (4, 4);
        A.at(make_matrix_index(0, 0)) = 16;
        A.at(make_matrix_index(0, 1)) = 2;
        A.at(make_matrix_index(0, 2)) = 3;
        A.at(make_matrix_index(0, 3)) = 13;

        A.at(make_matrix_index(1, 0)) = 5;
        A.at(make_matrix_index(1, 1)) = 11;
        A.at(make_matrix_index(1, 2)) = 10;
        A.at(make_matrix_index(1, 3)) = 8;

        A.at(make_matrix_index(2, 0)) = 9;
        A.at(make_matrix_index(2, 1)) = 7;
        A.at(make_matrix_index(2, 2)) = 6;
        A.at(make_matrix_index(2, 3)) = 12;

        A.at(make_matrix_index(3, 0)) = 4;
        A.at(make_matrix_index(3, 1)) = 14;
        A.at(make_matrix_index(3, 2)) = 15;
        A.at(make_matrix_index(3, 3)) = 1;
        
        rtsam<double, storage_order::row_major> A2 {
			{ 2,  2,3,
              9,  8,1,
             15,100,7,
			 99,  1,2,
			  5,  7,3},5};

        auto B {bidiag(A,0)};
        auto B2 {bidiag(A,0)};
        auto B3 {bidiag(A2,0)};

        auto C { B.u() * B.b() * transpose(B.v()) };
        auto C2 { B2.u() * B2.b() * transpose(B2.v()) };
        auto C3 { B3.u() * B3.b() * transpose(B3.v()) };
        
        BOOST_TEST(C == A);
        BOOST_TEST(C2 == A);
        BOOST_TEST(C3 == A2);
}

BOOST_AUTO_TEST_SUITE_END()
