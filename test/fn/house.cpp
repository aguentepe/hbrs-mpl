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


#define BOOST_TEST_MODULE givens_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/fn/house.hpp>

BOOST_AUTO_TEST_SUITE(house_test)

BOOST_AUTO_TEST_CASE(house_test1) {
        using namespace hbrs::mpl;

        rtsacv<double> x (3);
        x.at(0) = 5;
        x.at(1) = 2;
        x.at(2) = 3;

        auto nibeta { house(x) };
}

BOOST_AUTO_TEST_SUITE_END()
