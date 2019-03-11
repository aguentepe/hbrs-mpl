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

#pragma once

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_SIZE_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_SIZE_HPP

#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/detail/function_object.hpp>
#include <hbrs/mpl/fwd/dt/ctsav.hpp>
#include <hbrs/mpl/fwd/dt/rtsav.hpp>
#include <hbrs/mpl/fwd/dt/submatrix.hpp>
#include <hbrs/mpl/fwd/dt/subsequence.hpp>
#include <hbrs/mpl/fwd/dt/matrix_size.hpp>
#include <hbrs/mpl/fwd/dt/sm.hpp>
#include <hbrs/mpl/fwd/dt/smc.hpp>
#include <hbrs/mpl/fwd/dt/smcs.hpp>
#include <hbrs/mpl/fwd/dt/smr.hpp>
#include <hbrs/mpl/fwd/dt/smrs.hpp>
#include <hbrs/mpl/fwd/dt/sms.hpp>
#include <hbrs/mpl/fwd/dt/srv.hpp>
#include <hbrs/mpl/fwd/dt/scv.hpp>
#include <hbrs/mpl/fwd/dt/ctsam.hpp>
#include <hbrs/mpl/fwd/dt/rtsam.hpp>
#include <hbrs/mpl/fwd/dt/rtsacv.hpp>
#include <hbrs/mpl/fwd/dt/rtsarv.hpp>
#include <hbrs/mpl/fwd/dt/zas.hpp>
#include <boost/hana/tuple.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_ctsav,    ctsav_tag,    length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_rtsav,    rtsav_tag,    length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_submatrix,   submatrix_tag,   size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_subsequence, subsequence_tag, size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_sm,       sm_tag,       size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_smc,      smc_tag,      length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_smcs,     smcs_tag,     length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_smr,      smr_tag,      length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_smrs,     smrs_tag,     length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_sms,      sms_tag,      length)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_srv,      srv_tag,      size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_scv,      scv_tag,      size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_ctsam,    ctsam_tag,    size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_rtsam,    rtsam_tag,    size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_rtsacv,   rtsacv_tag,   size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_rtsarv,   rtsarv_tag,   size)
HBRS_MPL_DEF_FO_TRY_METHOD(size_impl_zas,      zas_tag,      length)

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_SIZE_IMPLS boost::hana::make_tuple(                                                  \
		hbrs::mpl::detail::size_impl_ctsav{},                                                                          \
		hbrs::mpl::detail::size_impl_rtsav{},                                                                          \
		hbrs::mpl::detail::size_impl_submatrix{},                                                                      \
		hbrs::mpl::detail::size_impl_subsequence{},                                                                    \
		hbrs::mpl::detail::size_impl_sm{},                                                                             \
		hbrs::mpl::detail::size_impl_smc{},                                                                            \
		hbrs::mpl::detail::size_impl_smcs{},                                                                           \
		hbrs::mpl::detail::size_impl_smr{},                                                                            \
		hbrs::mpl::detail::size_impl_smrs{},                                                                           \
		hbrs::mpl::detail::size_impl_sms{},                                                                            \
		hbrs::mpl::detail::size_impl_srv{},                                                                            \
		hbrs::mpl::detail::size_impl_scv{},                                                                            \
		hbrs::mpl::detail::size_impl_ctsam{},                                                                          \
		hbrs::mpl::detail::size_impl_rtsam{},                                                                          \
		hbrs::mpl::detail::size_impl_rtsacv{},                                                                         \
		hbrs::mpl::detail::size_impl_rtsarv{},                                                                         \
		hbrs::mpl::detail::size_impl_zas{}                                                                             \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_SIZE_HPP
