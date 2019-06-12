/* Copyright (c) 2016 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_DT_NO_CONTEXT_HPP
#define HBRS_MPL_DT_NO_CONTEXT_HPP

#include <hbrs/mpl/fwd/dt/no_context.hpp>
#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
struct no_context{};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <>
struct tag_of< hbrs::mpl::no_context > {
	using type = hbrs::mpl::no_context_tag;
};

template <>
struct make_impl<hbrs::mpl::no_context_tag> {
	
	static constexpr hbrs::mpl::no_context
	apply() {
		return {};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_NO_CONTEXT_HPP