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

#ifndef HBRS_MPL_DT_RTSAM_IMPL_HPP
#define HBRS_MPL_DT_RTSAM_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/smr.hpp>
#include <hbrs/mpl/dt/rtsacv.hpp>
#include <hbrs/mpl/dt/rtsarv.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/givens_rotation.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/minus.hpp>
#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/fn/select.hpp>
#include <hbrs/mpl/detail/translate_index.hpp>
#include <hbrs/mpl/dt/exception.hpp>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>
#include <vector>
#include <type_traits>

HBRS_MPL_NAMESPACE_BEGIN

template<
	typename /* type of matrix entries */ Ring,
	storage_order Order
>
struct rtsam {
	rtsam(matrix_size<std::size_t, std::size_t> sz) : rtsam(sz.m(), sz.n()) {}
	
	rtsam(std::vector<Ring> data, matrix_size<std::size_t, std::size_t> sz) : data_{data}, size_{sz} {
		if (sz.m() * sz.n() != data.size()) {
			BOOST_THROW_EXCEPTION(
				incompatible_matrix_sequence_exception{}
				<< errinfo_matrix_size{sz}
				<< errinfo_sequence_size{data.size()}
			);
		}
	}
	
	rtsam(std::size_t m, std::size_t n) : data_(m * n, Ring{0}), size_{m,n} {}

	explicit rtsam(submatrix<rtsam<Ring,Order>&, matrix_index<std::size_t,std::size_t>, matrix_size<std::size_t,std::size_t>> const& M)
	: rtsam(M.size()) {
		for (std::size_t i{0}; i < m(size()); ++i) {
			for (std::size_t j{0}; j < n(size()); ++j) {
				at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j));
			}
		}
	}

    explicit rtsam(rtsacv<Ring> const& v) : data_ {v.data()}, size_{v.m(), 1} {}

    explicit rtsam(rtsarv<Ring> const& v) : data_ {v.data()}, size_{1, v.n()} {}
	
	rtsam(rtsam const&) = default;
	rtsam(rtsam &&) = default;
	
	rtsam&
	operator=(rtsam const&) = default;
	rtsam&
	operator=(rtsam &&) = default;

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
	rtsam&
	operator=(detail::givens_rotation_expression<givens_rotation<Ring> const&, rtsam<Ring,Order> const&> const& e) {
		if (&(e.t2()) != this) {
			*this = e.t2();
		}

		auto i     {e.t1().i()};
		auto k     {e.t1().k()};
		auto theta {e.t1().theta()};

		BOOST_ASSERT(i < m(size()));
		BOOST_ASSERT(k < m(size()));

		for (std::size_t j {0}; j <= n(size()) - 1; ++j) {
			double const tau1 {at(make_matrix_index(i, j))};
			double const tau2 {at(make_matrix_index(k, j))};
			at(make_matrix_index(i, j)) = theta.at(0) * tau1 - theta.at(1) * tau2;
			at(make_matrix_index(k, j)) = theta.at(1) * tau1 + theta.at(0) * tau2;
		}
		return *this;
	}

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
	rtsam&
	operator=(detail::givens_rotation_expression<rtsam<Ring,Order> const&, givens_rotation<Ring> const&> const& e) {
		if (&(e.t1()) != this) {
			*this = e.t1();
		}

		auto i     {e.t2().i()};
		auto k     {e.t2().k()};
		auto theta {e.t2().theta()};

		BOOST_ASSERT(i < n(size()));
		BOOST_ASSERT(k < n(size()));

		for (std::size_t j {0}; j <= m(size()) - 1; ++j) {
			auto const tau1 {at(make_matrix_index(j, i))};
			auto const tau2 {at(make_matrix_index(j, k))};
			at(make_matrix_index(j, i)) = theta.at(0) * tau1 - theta.at(1) * tau2;
			at(make_matrix_index(j, k)) = theta.at(1) * tau1 + theta.at(0) * tau2;
		}
		return *this;
	}

	auto const&
	size() const { return size_; };
	
	decltype(auto)
	order() const { return storage_order_c<Order>; }
	
	template<typename Index>
	decltype(auto)
	at(Index && i) {
		return data_.at(
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		);
	}
	
	template<typename Index>
	decltype(auto)
	at(Index && i) const {
		return data_.at(
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		);
	}
	
	template<typename Index>
	auto
	operator[](Index && i) & { return smr<rtsam &, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	auto
	operator[](Index && i) const& { return smr<rtsam const&, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	auto
	operator[](Index && i) && { return make_smr(std::move(*this), HBRS_MPL_FWD(i)); }
	
    // Several operator() functions to return submatrices.
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const {
		return select(*this, std::make_pair(rows, column));
	}
	decltype(auto)
	operator()(std::size_t const row, range<std::size_t,std::size_t> const& columns) const {
		return select(*this, std::make_pair(row, columns));
	}
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) & {
		return select(*this, std::make_pair(rows, columns));
	}
private:
	std::vector<Ring> data_;
	matrix_size<std::size_t, std::size_t> size_;
};

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		(std::is_same_v< hana::tag_of_t<T1>, rtsam_tag     > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag >) ||
		(std::is_same_v< hana::tag_of_t<T1>, submatrix_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag >)
	>* = nullptr
>
decltype(auto)
operator-(T1 && t1, T2 && t2) {
	return minus(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator*(rtsam<Ring,Order> const& M, Ring const& d) {
	return multiply(M,d);
}

template<typename Ring, storage_order Order>
rtsam<Ring,Order>
operator*(Ring const& d, rtsam<Ring,Order> const& M) {
	return multiply(d,M);
}

template<
	typename T1,
	typename T2,
	typename std::enable_if_t<
		(std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, submatrix_tag       >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, rtsacv_tag          >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsam_tag           > && std::is_same_v< hana::tag_of_t<T2>, givens_rotation_tag >) ||
		(std::is_same_v< hana::tag_of_t<T1>, submatrix_tag       > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           >) ||
		(std::is_same_v< hana::tag_of_t<T1>, rtsarv_tag          > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           >) ||
		(std::is_same_v< hana::tag_of_t<T1>, givens_rotation_tag > && std::is_same_v< hana::tag_of_t<T2>, rtsam_tag           >)
	>* = nullptr
>
decltype(auto)
operator*(T1 && t1, T2 && t2) {
	return multiply(HBRS_MPL_FWD(t1), HBRS_MPL_FWD(t2));
}

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <
	typename Ring,
	hbrs::mpl::storage_order Order
>
struct tag_of< hbrs::mpl::rtsam<Ring, Order> > {
	using type = hbrs::mpl::rtsam_tag;
};

template <>
struct make_impl<hbrs::mpl::rtsam_tag> {
	template <
		typename Ring,
		hbrs::mpl::storage_order Order
	>
	static hbrs::mpl::rtsam<Ring, Order>
	apply(hana::basic_type<Ring>, hbrs::mpl::matrix_size<std::size_t, std::size_t> sz, hbrs::mpl::storage_order_<Order>) {
		return {sz};
	}
	
	template <
		typename Ring,
		hbrs::mpl::storage_order Order
	>
	static hbrs::mpl::rtsam<std::remove_const_t<Ring>, Order>
	apply(std::vector<Ring> data, hbrs::mpl::matrix_size<std::size_t, std::size_t> sz, hbrs::mpl::storage_order_<Order>) {
		return {data, sz};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_RTSAM_IMPL_HPP
