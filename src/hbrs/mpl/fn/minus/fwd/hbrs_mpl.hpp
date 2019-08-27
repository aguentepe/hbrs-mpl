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

#ifndef HBRS_MPL_FN_MINUS_FWD_HBRS_MPL_HPP
#define HBRS_MPL_FN_MINUS_FWD_HBRS_MPL_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/rtsam/fwd.hpp>
#include <hbrs/mpl/dt/submatrix/fwd.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct minus_impl_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const;

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset, Size> const& M1, rtsam<Ring,Order> const& M2) const;

private:
	template<
		typename Ring,
		typename Matrix1,
		typename Matrix2
	>
	decltype(auto)
	impl(Matrix1 const& Mp1, Matrix2 const& Mp2, hana::basic_type<Ring>) const;
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FN_MINUS_IMPLS_HBRS_MPL boost::hana::make_tuple(                                                     \
		hbrs::mpl::detail::minus_impl_rtsam{}                                                                          \
	)

#endif // !HBRS_MPL_FN_MINUS_FWD_HBRS_MPL_HPP
