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

#pragma once

#ifndef HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_HPP
#define HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_HPP

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

struct multiply_impl_rtsarv_rtsacv {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) const {
		Ring sum {0};
		for (std::size_t i {0}; i < v1.n(); ++i)
			sum += v1.at(i) * v2.at(i);
		return sum;
	}
};

struct multiply_impl_rtsam_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
		BOOST_ASSERT(M1.n() == M2.m());

		rtsam<Ring,Order> result {M1.m(), M2.n()};
		for (std::size_t i {0}; i < result.m(); ++i) {
			for (std::size_t j {0}; j < result.n(); ++j) {
				result.at(make_matrix_index(i, j)) = M1(i, range<std::size_t,std::size_t>(std::size_t{0}, M1.n() - 1)) * M2(range<std::size_t,std::size_t>(std::size_t{0}, M2.m() - 1), j);
			}
		}
		return result;
	}
};

struct multiply_impl_rtsarv_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v, rtsam<Ring,Order> const& M) const {
		BOOST_ASSERT(v.n() == M.m());

		rtsarv<Ring> result(M.n());
		for (std::size_t i {0}; i < result.n(); ++i) {
			result.at(i) = v * M(range<std::size_t,std::size_t>(std::size_t{0}, M.m() - 1), i);
		}
		return result;
	}
};

struct multiply_impl_rtsam_rtsacv {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M, rtsacv<Ring> const& v) const {
		BOOST_ASSERT(M.n() == v.m());

		rtsacv<Ring> result(M.m());
		for (std::size_t i {0}; i < result.m(); ++i) {
			result.at(i) = M(i, range<std::size_t,std::size_t>(std::size_t{0}, M.n() - 1)) * v;
		}
		return result;
	}
};

struct multiply_impl_rtsacv_rtsarv {
	template<
		typename Ring
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) const {
		rtsam<Ring,storage_order::row_major> result {v1.m(), v2.n()};
		for (std::size_t i {0}; i < result.m(); ++i) {
			for (std::size_t j {0}; j < result.n(); ++j) {
				result.at(make_matrix_index(i, j)) = v1.at(i) * v2.at(j);
			}
		}
		return result;
	}
};

struct multiply_impl_rtsam_ring {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> M, Ring const& d) const {
		for (std::size_t i {0}; i < M.m(); ++i) {
			for (std::size_t j {0}; j < M.n(); ++j) {
				M.at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j)) * d;
			}
		}
		return M;
	}
};

struct multiply_impl_ring_rtsam {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(Ring const& d, rtsam<Ring,Order> M) const {
		return multiply(M,d);
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_IMPLS boost::hana::make_tuple(\
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsacv{},\
		hbrs::mpl::detail::multiply_impl_rtsam_rtsam{},\
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsam{},\
		hbrs::mpl::detail::multiply_impl_rtsam_rtsacv{},\
		hbrs::mpl::detail::multiply_impl_rtsacv_rtsarv{},\
		hbrs::mpl::detail::multiply_impl_rtsam_ring{},\
		hbrs::mpl::detail::multiply_impl_ring_rtsam{}\
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_HPP
