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

#ifndef HBRS_MPL_FN_MULTIPLY_FWD_HBRS_MPL_HPP
#define HBRS_MPL_FN_MULTIPLY_FWD_HBRS_MPL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/rtsacv/fwd.hpp>
#include <hbrs/mpl/dt/rtsarv/fwd.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct multiply_impl_rtsacv_ring {
	template<typename Ring>
	auto
	operator()(rtsacv<Ring> v, Ring const& s) const;
};

struct multiply_impl_ring_rtsacv {
	template<typename Ring>
	decltype(auto)
	operator()(Ring const& s, rtsacv<Ring> v) const;
};

struct multiply_impl_rtsarv_rtsacv {
	template<typename Ring>
	decltype(auto)
	operator()(rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) const;
};

struct multiply_impl_rtsacv_rtsarv {
	template<typename Ring>
	decltype(auto)
	operator()(rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) const;
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_MULTIPLY_IMPLS_HBRS_MPL boost::hana::make_tuple(                                                   \
		hbrs::mpl::detail::multiply_impl_rtsacv_ring{},                                                                \
		hbrs::mpl::detail::multiply_impl_ring_rtsacv{},                                                                \
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsacv{},                                                              \
		hbrs::mpl::detail::multiply_impl_rtsacv_rtsarv{}                                                               \
	)

#endif // !HBRS_MPL_FN_MULTIPLY_FWD_HBRS_MPL_HPP
