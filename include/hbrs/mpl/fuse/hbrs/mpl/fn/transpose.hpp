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

#pragma once

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_TRANSPOSE_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_TRANSPOSE_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/srv.hpp>
#include <hbrs/mpl/dt/scv.hpp>
#include <hbrs/mpl/dt/Matrix.hpp>
#include <hbrs/mpl/dt/CVector.hpp>
#include <hbrs/mpl/dt/RVector.hpp>
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
	operator()(Vector && v) const {
		return make_scv(HBRS_MPL_FWD(v));
	}
};

struct transpose_impl_scv {
	template <
		typename Vector,
		typename std::enable_if_t<
			std::is_same< hana::tag_of_t<Vector>, scv_tag>::value
		>* = nullptr
	>
	constexpr auto 
	operator()(Vector && v) const {
		return make_srv(HBRS_MPL_FWD(v));
	}
};

struct transpose_impl_Matrix {
    decltype(auto)
    operator()(Matrix const& M) const {
        Matrix result {M.n(), M.m()};
        for (std::size_t i {0}; i < result.m(); ++i) {
            for (std::size_t j {0}; j < result.n(); ++j) {
                result.at(i, j) = M.at(j, i);
            }
        }
        return result;
    }
};

struct transpose_impl_CVector {
    decltype(auto)
    operator()(CVector<double> const& v) const {
        return RVector{v};
    }
};

struct transpose_impl_RVector {
    decltype(auto)
    operator()(RVector const& v) const {
        return v.transpose();
    }
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_TRANSPOSE_IMPLS boost::hana::make_tuple(                                             \
		hbrs::mpl::detail::transpose_impl_srv{},                                                                       \
		hbrs::mpl::detail::transpose_impl_scv{},                                                                        \
		hbrs::mpl::detail::transpose_impl_Matrix{},\
		hbrs::mpl::detail::transpose_impl_CVector{},\
		hbrs::mpl::detail::transpose_impl_RVector{}\
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_TRANSPOSE_HPP
