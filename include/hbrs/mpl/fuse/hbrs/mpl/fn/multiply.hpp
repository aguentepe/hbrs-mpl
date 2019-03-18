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
#include <hbrs/mpl/dt/ssm.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>
#include <typeinfo>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator* (rtsam<Ring,Order> M, Ring const& d) {
	return multiply(M,d);
}

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator* (Ring const& d, rtsam<Ring,Order> M) {
	return multiply(d,M);
}

template<typename Ring>
decltype(auto)
operator*(Ring const& s, rtsacv<Ring> v) {
	return multiply(v,s);
}

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag  > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag  > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag  > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag  > && std::is_same_v< hana::tag_of_t<T2>, ssm_tag    > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsacv_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsarv_tag > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag  > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag > && std::is_same_v< hana::tag_of_t<T2>, ssm_tag >
	>* = nullptr
>
decltype(auto)
operator*(T1 && t1, T2 && t2) {
	return multiply(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

namespace detail {

struct multiply_impl_rtsarv_rtsacv {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) const {
		Ring sum {0};
		for (std::size_t i {0}; i < v1.size(); ++i)
			sum += v1.at(i) * v2.at(i);
		return sum;
	}
};

struct multiply_impl_rtsam_rtsam {
	template<
		typename Matrix1,
		typename Matrix2,
		typename std::enable_if_t<
			(std::is_same_v< hana::tag_of_t<Matrix1>, rtsam_tag > || std::is_same_v< hana::tag_of_t<Matrix1>, ssm_tag >) &&
			(std::is_same_v< hana::tag_of_t<Matrix2>, rtsam_tag > || std::is_same_v< hana::tag_of_t<Matrix2>, ssm_tag >)
		>* = nullptr
	>
	/* constexpr */ 
	decltype(auto)
	operator()(Matrix1 const& M1, Matrix2 const& M2) const {
		BOOST_ASSERT(M1.n() == M2.m());

		rtsam<decltype(M1.at(make_matrix_index(0,0))), storage_order::row_major> result {M1.m(), M2.n()};
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
		typename Matrix,
		typename std::enable_if_t< std::is_same_v< hana::tag_of_t<Matrix>, rtsam_tag > || std::is_same_v< hana::tag_of_t<Matrix>, ssm_tag > >* = nullptr
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v, Matrix const& M) const {
		BOOST_ASSERT(v.size() == M.m());

		rtsarv<Ring> result(M.n());
		for (std::size_t i {0}; i < result.size(); ++i) {
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
		BOOST_ASSERT(M.n() == v.size());

		rtsacv<Ring> result(M.m());
		for (std::size_t i {0}; i < result.size(); ++i) {
			result.at(i) = M(i, range<std::size_t,std::size_t>(std::size_t{0}, M.n() - 1)) * v;
		}
		return result;
	}
};

struct multiply_impl_rtsacv_rtsarv {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) const {
		rtsam<Ring,storage_order::row_major> result {v1.size(), v2.size()};
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

struct multiply_impl_rtsacv_ring {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsacv<Ring> v, Ring const& s) const {
		for (std::size_t i {0}; i < v.size(); ++i)
			v.at(i) *= s;
		return v;
	}
};

struct multiply_impl_ring_rtsacv {
	template<typename Ring>
	/* constexpr */ 
	decltype(auto)
	operator()(Ring const& s, rtsacv<Ring> v) const {
		return multiply(v,s);
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_IMPLS boost::hana::make_tuple(                                              \
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsacv{},                                                              \
		hbrs::mpl::detail::multiply_impl_rtsam_rtsam{},                                                                \
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsam{},                                                               \
		hbrs::mpl::detail::multiply_impl_rtsam_rtsacv{},                                                               \
		hbrs::mpl::detail::multiply_impl_rtsacv_rtsarv{},                                                              \
		hbrs::mpl::detail::multiply_impl_rtsam_ring{},                                                                 \
		hbrs::mpl::detail::multiply_impl_ring_rtsam{},                                                                 \
		hbrs::mpl::detail::multiply_impl_rtsacv_ring{},                                                                \
		hbrs::mpl::detail::multiply_impl_ring_rtsacv{}                                                                 \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_HPP
