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
#include <hbrs/mpl/dt/submatrix.hpp>
#include <hbrs/mpl/dt/givens_rotation.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <cmath>
#include <typeinfo>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator* (rtsam<Ring,Order> const& M, Ring const& d) {
	return multiply(M,d);
}

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator* (Ring const& d, rtsam<Ring,Order> const& M) {
	return multiply(d,M);
}

template<typename Ring>
decltype(auto)
operator*(Ring const& s, rtsacv<Ring> const& v) {
	return multiply(v,s);
}

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag       > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag          > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, givens_rotation_tag > ||
		std::is_same_v< hana::tag_of_t<T1>, submatrix_tag       > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           > ||
		std::is_same_v< hana::tag_of_t<T1>, submatrix_tag       > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag       > ||
		std::is_same_v< hana::tag_of_t<T1>, submatrix_tag       > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag          > ||
		std::is_same_v< hana::tag_of_t<T1>, submatrix_tag       > && std::is_same_v< hana::tag_of_t<T2>, givens_rotation_tag > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsacv_tag          > && std::is_same_v< hana::tag_of_t<T2>, rtsarv_tag          > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag          > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag          > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag       > ||
		std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag          > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag          > ||
		std::is_same_v< hana::tag_of_t<T1>, givens_rotation_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           > ||
		std::is_same_v< hana::tag_of_t<T1>, givens_rotation_tag > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag       >
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

struct multiply_impl_matrix_matrix {

	template<
		typename Ring,
		storage_order Order
	>
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
		return impl(M1,M2, hana::type_c<Ring>);
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M1, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M2) const {
		return impl(M1,M2, hana::type_c<Ring>);
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M1, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M2) const {
		return impl(M1,M2, hana::type_c<Ring>);
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M1, rtsam<Ring,Order> const& M2) const {
		return impl(M1,M2, hana::type_c<Ring>);
	}

private:
	template<
		typename Ring,
		typename Matrix1,
		typename Matrix2
	>
	decltype(auto)
	impl(Matrix1 && M1, Matrix2 && M2, hana::basic_type<Ring>) const {
		BOOST_ASSERT(M1.n() == M2.m());

		rtsam<Ring, storage_order::row_major> result {M1.m(), M2.n()};
		for (std::size_t i {0}; i < result.m(); ++i) {
			for (std::size_t j {0}; j < result.n(); ++j) {
				result.at(make_matrix_index(i, j)) = M1(i, range<std::size_t,std::size_t>(std::size_t{0}, M1.n() - 1)) * M2(range<std::size_t,std::size_t>(std::size_t{0}, M2.m() - 1), j);
			}
		}
		return result;
	}
};

struct multiply_impl_rtsarv_matrix {

	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v, rtsam<Ring,Order> const& M) const {
		return impl(v,M, hana::type_c<Ring>);
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsarv<Ring> const& v, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M) const {
		return impl(v,M, hana::type_c<Ring>);
	}

private:
	template<
		typename Ring,
		typename RVector,
		typename Matrix
	>
	/* constexpr */ 
	decltype(auto)
	impl(RVector && vp, Matrix && Mp, hana::basic_type<Ring>) const {
		decltype(auto) v {vp};
		decltype(auto) M {Mp};

		BOOST_ASSERT(v.size() == M.m());

		rtsarv<Ring> result(M.n());
		for (std::size_t i {0}; i < result.size(); ++i) {
			result.at(i) = v * M(range<std::size_t,std::size_t>(std::size_t{0}, M.m() - 1), i);
		}
		return result;
	}
};

struct multiply_impl_matrix_rtsacv {

	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& M, rtsacv<Ring> const& v) const {
		return impl(M,v, hana::type_c<Ring>);
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M, rtsacv<Ring> const& v) const {
		return impl(M,v, hana::type_c<Ring>);
	}

private:
	template<
		typename Ring,
		typename Matrix,
		typename CVector
	>
	decltype(auto)
	impl(Matrix && M, CVector && v, hana::basic_type<Ring>) const {
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
	/* decltype(auto) */
	auto
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

/*
 * Chapter 5.1.9 (Applying Givens Rotations) on page 241
 * A = G(i,k,theta)^T * A
 *              --     --T
 *              |       |
 *              |  c s  |
 * A([i,k],:) = |       | * A([i,k],:)
 *              | -s c  |
 *              |       |
 *              --     --
 *
 * Apply the Givens roation on A and return A.
 */
struct multiply_impl_givens_rotation_matrix {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(givens_rotation<Ring> const& G, rtsam<Ring,Order> const& A) const {
		return impl(G,A, std::integral_constant<storage_order, Order>{});
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(givens_rotation<Ring> const& G, submatrix<rtsam<Ring,Order>&, Offset, Size> const& A) const {
		return impl(G,A, std::integral_constant<storage_order, Order>{});
	}
private:
	template<
		typename Ring,
		storage_order Order,
		typename Matrix
	>
	decltype(auto)
	impl(givens_rotation<Ring> const& G, Matrix const& A, std::integral_constant<storage_order, Order>) const {
		decltype(auto) i {G.i()};
		decltype(auto) k {G.k()};
		decltype(auto) theta {G.theta()};

		BOOST_ASSERT(i < A.m());
		BOOST_ASSERT(k < A.m());

		rtsam<Ring,Order> R{A};
		for (std::size_t j {0}; j <= R.n() - 1; ++j) {
			double const tau1 {R.at(make_matrix_index(i, j))};
			double const tau2 {R.at(make_matrix_index(k, j))};
			R.at(make_matrix_index(i, j)) = theta.at(0) * tau1 - theta.at(1) * tau2;
			R.at(make_matrix_index(k, j)) = theta.at(1) * tau1 + theta.at(0) * tau2;
		}
		return R;
	}
};

/*
 * Chapter 5.1.9 (Applying Givens Rotations) on page 241
 * A = A * G(i,k,theta)
 *                           --     --
 *                           |       |
 *                           |  c s  |
 * A(:,[i,k]) = A(:,[i,k]) * |       |
 *                           | -s c  |
 *                           |       |
 *                           --     --
 *
 * Apply the Givens roation on A and return A.
 */
struct multiply_impl_matrix_givens_rotation {
	template<
		typename Ring,
		storage_order Order
	>
	/* constexpr */ 
	decltype(auto)
	operator()(rtsam<Ring,Order> const& A, givens_rotation<Ring> const& G) const {
		return impl(A,G, std::integral_constant<storage_order, Order>{});
	}

	template<
		typename Ring,
		storage_order Order,
		typename Offset,
		typename Size
	>
	decltype(auto)
	operator()(submatrix<rtsam<Ring,Order>&, Offset, Size> const& A, givens_rotation<Ring> const& G) const {
		return impl(A,G, std::integral_constant<storage_order, Order>{});
	}
private:
	template<
		typename Ring,
		storage_order Order,
		typename Matrix
	>
	decltype(auto)
	impl(Matrix const& A, givens_rotation<Ring> const& G, std::integral_constant<storage_order, Order>) const {
		decltype(auto) i {G.i()};
		decltype(auto) k {G.k()};
		decltype(auto) theta {G.theta()};

		BOOST_ASSERT(i < A.n());
		BOOST_ASSERT(k < A.n());

		rtsam<Ring,Order> R{A};
		for (std::size_t j {0}; j <= R.m() - 1; ++j) {
			double const tau1 {R.at(make_matrix_index(j, i))};
			double const tau2 {R.at(make_matrix_index(j, k))};
			R.at(make_matrix_index(j, i)) = theta.at(0) * tau1 - theta.at(1) * tau2;
			R.at(make_matrix_index(j, k)) = theta.at(1) * tau1 + theta.at(0) * tau2;
		}
		return R;
	}
};

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END

#define HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_IMPLS boost::hana::make_tuple(                                              \
		hbrs::mpl::detail::multiply_impl_rtsarv_rtsacv{},                                                              \
		hbrs::mpl::detail::multiply_impl_matrix_matrix{},                                                                \
		hbrs::mpl::detail::multiply_impl_rtsarv_matrix{},                                                               \
		hbrs::mpl::detail::multiply_impl_matrix_rtsacv{},                                                               \
		hbrs::mpl::detail::multiply_impl_rtsacv_rtsarv{},                                                              \
		hbrs::mpl::detail::multiply_impl_rtsam_ring{},                                                                 \
		hbrs::mpl::detail::multiply_impl_ring_rtsam{},                                                                 \
		hbrs::mpl::detail::multiply_impl_rtsacv_ring{},                                                                \
		hbrs::mpl::detail::multiply_impl_ring_rtsacv{},                                                                \
		hbrs::mpl::detail::multiply_impl_givens_rotation_matrix{},                                                     \
		hbrs::mpl::detail::multiply_impl_matrix_givens_rotation{}                                                      \
	)

#endif // !HBRS_MPL_FUSE_HBRS_MPL_FN_MULTIPLY_HPP
