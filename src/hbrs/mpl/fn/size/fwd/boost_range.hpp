/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_SIZE_FWD_BOOST_RANGE_HPP
#define HBRS_MPL_FN_SIZE_FWD_BOOST_RANGE_HPP

#include <hbrs/mpl/config.hpp>
#include <boost/range/irange.hpp>
#include <boost/hana/tuple.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct size_impl_range_integer_range {
	template<typename Integer>
	constexpr decltype(auto)
	operator()(boost::integer_range<Integer> const& s) const;
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_SIZE_IMPLS_BOOST_RANGE boost::hana::make_tuple(                                                    \
		hbrs::mpl::detail::size_impl_range_integer_range{}                                                             \
	)

#endif // !HBRS_MPL_FN_SIZE_FWD_BOOST_RANGE_HPP
