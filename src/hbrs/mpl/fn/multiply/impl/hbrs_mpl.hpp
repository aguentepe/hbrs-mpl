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

#ifndef HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP
#define HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP

#include "../fwd/hbrs_mpl.hpp"

#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/submatrix.hpp>
#include <hbrs/mpl/dt/givens_rotation.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <cmath>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace detail {

template<typename Ring>
auto
multiply_impl_rtsacv_ring::operator()(rtsacv<Ring> v, Ring const& s) const {
	for (std::size_t i {0}; i < v.length(); ++i) {
		v.at(i) *= s;
	}
	return v;
}


template<typename Ring>
decltype(auto)
multiply_impl_ring_rtsacv::operator()(Ring const& s, rtsacv<Ring> v) const {
	return multiply(v,s);
}

template<typename Ring>
decltype(auto)
multiply_impl_rtsarv_rtsacv::operator()(rtsarv<Ring> const& v1, rtsacv<Ring> const& v2) const {
	typedef std::decay_t<Ring> _Ring_;
	_Ring_ sum {0};
	for (std::size_t i {0}; i < v1.length(); ++i)
		sum += v1.at(i) * v2.at(i);
	return sum;
}


template<typename Ring>
decltype(auto)
multiply_impl_rtsacv_rtsarv::operator()(rtsacv<Ring> const& v1, rtsarv<Ring> const& v2) const {
	typedef std::decay_t<Ring> _Ring_;
	//TODO: Implement code to choose best storage_order
	rtsam<_Ring_,storage_order::row_major> result {v1.length(), v2.length()};
	for (std::size_t i {0}; i < (*m)((*size)(result)); ++i) {
		for (std::size_t j {0}; j < (*n)((*size)(result)); ++j) {
			result.at(make_matrix_index(i, j)) = v1.at(i) * v2.at(j);
		}
	}
	return result;
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_matrix_matrix::operator()(rtsam<Ring,Order> const& M1, rtsam<Ring,Order> const& M2) const {
	return impl(M1,M2, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_matrix_matrix::operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M1, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M2) const {
	return impl(M1,M2, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_matrix_matrix::operator()(rtsam<Ring,Order> const& M1, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M2) const {
	return impl(M1,M2, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_matrix_matrix::operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M1, rtsam<Ring,Order> const& M2) const {
	return impl(M1,M2, hana::type_c<Ring>);
}

template<
	typename Ring,
	typename Matrix1,
	typename Matrix2
>
decltype(auto)
multiply_impl_matrix_matrix::impl(Matrix1 const& M1, Matrix2 const& M2, hana::basic_type<Ring>) const {
	BOOST_ASSERT((*n)((*size)(M1)) == (*m)((*size)(M2)));

	typedef std::decay_t<Ring> _Ring_;
	rtsam<_Ring_, storage_order::row_major> result {(*m)((*size)(M1)), (*n)((*size)(M2))};
	for (std::size_t i {0}; i < (*m)((*size)(result)); ++i) {
		for (std::size_t j {0}; j < (*n)((*size)(result)); ++j) {
			result.at(make_matrix_index(i, j)) = M1(i, range<std::size_t,std::size_t>(std::size_t{0}, (*n)((*size)(M1)) - 1)) * M2(range<std::size_t,std::size_t>(std::size_t{0}, (*m)((*size)(M2)) - 1), j);
		}
	}
	return result;
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_rtsarv_matrix::operator()(rtsarv<Ring> const& v, rtsam<Ring,Order> const& M) const {
	return impl(v,M, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_rtsarv_matrix::operator()(rtsarv<Ring> const& v, submatrix<rtsam<Ring,Order>&, Offset,Size> const& M) const {
	return impl(v,M, hana::type_c<Ring>);
}

template<
	typename Ring,
	typename RVector,
	typename Matrix
>
decltype(auto)
multiply_impl_rtsarv_matrix::impl(RVector const& v, Matrix const& M, hana::basic_type<Ring>) const {
	BOOST_ASSERT(v.length() == (*m)((*size)(M)));

	rtsarv<Ring> result((*n)((*size)(M)));
	for (std::size_t i {0}; i < result.length(); ++i) {
		result.at(i) = v * M(range<std::size_t,std::size_t>(std::size_t{0}, (*m)((*size)(M)) - 1), i);
	}
	return result;
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_matrix_rtsacv::operator()(rtsam<Ring,Order> const& M, rtsacv<Ring> const& v) const {
	return impl(M,v, hana::type_c<Ring>);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_matrix_rtsacv::operator()(submatrix<rtsam<Ring,Order>&, Offset,Size> const& M, rtsacv<Ring> const& v) const {
	return impl(M,v, hana::type_c<Ring>);
}

template<
	typename Ring,
	typename Matrix,
	typename CVector
>
decltype(auto)
multiply_impl_matrix_rtsacv::impl(Matrix const& M, CVector const& v, hana::basic_type<Ring>) const {
	BOOST_ASSERT((*n)((*size)(M)) == v.length());

	rtsacv<Ring> result((*m)((*size)(M)));
	for (std::size_t i {0}; i < result.length(); ++i) {
		result.at(i) = M(i, range<std::size_t,std::size_t>(std::size_t{0}, (*n)((*size)(M)) - 1)) * v;
	}
	return result;
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_rtsam_ring::operator()(rtsam<Ring,Order> M, Ring const& d) const {
	for (std::size_t i {0}; i < (*m)((*size)(M)); ++i) {
		for (std::size_t j {0}; j < (*n)((*size)(M)); ++j) {
			M.at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j)) * d;
		}
	}
	return M;
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_ring_rtsam::operator()(Ring const& d, rtsam<Ring,Order> M) const {
	return multiply(M,d);
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_givens_rotation_matrix::operator()(givens_rotation<Ring> const& G, rtsam<Ring,Order> const& A) const {
	return givens_rotation_expression<givens_rotation<Ring> const&, rtsam<Ring,Order> const&>(G,A);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_givens_rotation_matrix::operator()(givens_rotation<Ring> const& G, submatrix<rtsam<Ring,Order>&, Offset, Size> const& A) const {
	return givens_rotation_expression<givens_rotation<Ring> const&, submatrix<rtsam<Ring,Order>&, Offset, Size> const&>(G,A);
}


template<
	typename Ring,
	storage_order Order
>
decltype(auto)
multiply_impl_matrix_givens_rotation::operator()(rtsam<Ring,Order> const& A, givens_rotation<Ring> const& G) const {
	return givens_rotation_expression<rtsam<Ring,Order> const&, givens_rotation<Ring> const&>(A,G);
}

template<
	typename Ring,
	storage_order Order,
	typename Offset,
	typename Size
>
decltype(auto)
multiply_impl_matrix_givens_rotation::operator()(submatrix<rtsam<Ring,Order>&, Offset, Size> const& A, givens_rotation<Ring> const& G) const {
	return givens_rotation_expression<submatrix<rtsam<Ring,Order>&, Offset, Size> const&, givens_rotation<Ring> const&>(A,G);
}

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END
#endif // !HBRS_MPL_FN_MULTIPLY_IMPL_HBRS_MPL_HPP
