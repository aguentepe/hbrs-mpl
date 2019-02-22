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

#include <hbrs/mpl/dt/zas.hpp>
#include <hbrs/mpl/dt/rtsav.hpp>
#include <hbrs/mpl/fn/givens.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <boost/hana/integral_constant.hpp>

#include <boost/hana/tuple.hpp>
#include <boost/hana/first.hpp>
#include <boost/hana/second.hpp>
#include <array>
#include <tuple>
#include <vector>

BOOST_AUTO_TEST_SUITE(givens_test)

BOOST_AUTO_TEST_CASE(given_test1) {
        using namespace hbrs::mpl;

        auto cs { givens(1,1) };
	BOOST_TEST( cs.at(0) == 0.707107 );
	BOOST_TEST( cs.at(1) == -0.707107 );
}

BOOST_AUTO_TEST_SUITE_END()
