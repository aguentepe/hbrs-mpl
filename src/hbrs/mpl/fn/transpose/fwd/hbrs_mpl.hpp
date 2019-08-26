/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_MPL_FN_TRANSPOSE_FWD_HBRS_MPL_HPP
#define HBRS_MPL_FN_TRANSPOSE_FWD_HBRS_MPL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/dt/srv/fwd.hpp>
#include <hbrs/mpl/dt/scv/fwd.hpp>
#include <hbrs/mpl/dt/rtsam/fwd.hpp>
#include <hbrs/mpl/dt/submatrix/fwd.hpp>
#include <hbrs/mpl/dt/rtsacv/fwd.hpp>
#include <hbrs/mpl/dt/rtsarv/fwd.hpp>

#include <boost/hana/tuple.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct transpose_impl_srv {
	template <
		typename Vector,
		typename std::enable_if_t<
			std::is_same< hana::tag_of_t<Vector>, srv_tag>::value
		>* = nullptr
	>
	constexpr auto
	operator()(Vector && v) const;
};

struct transpose_impl_scv {
	template <
		typename Vector,
		typename std::enable_if_t<
			std::is_same< hana::tag_of_t<Vector>, scv_tag>::value
		>* = nullptr
	>
	constexpr auto 
	operator()(Vector && v) const;
};

struct transpose_impl_matrix {
//TODO: Rename to e.g. transpose_impl_rtsam
	template<
		typename Ring,
		storage_order Order
	>
    decltype(auto)
    operator()(rtsam<Ring,Order> const& M) const;

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M) const;

private:
	template<typename Ring, typename Matrix>
	decltype(auto)
	impl(Matrix const& Mp, hana::basic_type<Ring>) const;
};

struct transpose_impl_rtsacv {
	template<typename Ring>
	decltype(auto)
    operator()(rtsacv<Ring> const& v) const;
};

struct transpose_impl_rtsarv {
	template<typename Ring>
	decltype(auto)
    operator()(rtsarv<Ring> const& v) const;
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_TRANSPOSE_IMPLS_HBRS_MPL boost::hana::make_tuple(                                                  \
		hbrs::mpl::detail::transpose_impl_srv{},                                                                       \
		hbrs::mpl::detail::transpose_impl_scv{},                                                                       \
		hbrs::mpl::detail::transpose_impl_matrix{},                                                                    \
		hbrs::mpl::detail::transpose_impl_rtsacv{},                                                                    \
		hbrs::mpl::detail::transpose_impl_rtsarv{}                                                                     \
	)

#endif // !HBRS_MPL_FN_TRANSPOSE_FWD_HBRS_MPL_HPP
