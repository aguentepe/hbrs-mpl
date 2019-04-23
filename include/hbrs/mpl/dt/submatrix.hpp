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

#ifndef HBRS_MPL_DT_SUBMATRIX_HPP
#define HBRS_MPL_DT_SUBMATRIX_HPP

#include <hbrs/mpl/fwd/dt/submatrix.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

#include <hbrs/mpl/dt/subsequence.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/smr.hpp>
#include <hbrs/mpl/dt/givens_rotation.hpp>

#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/plus.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>

#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename Matrix, typename Offset, typename Size>
struct submatrix {
	
	template<typename Matrix_, typename Offset_, typename Size_>
	constexpr 
	submatrix(Matrix_ && mat, Offset_ && o, Size_ && sz)
	: mat_{HBRS_MPL_FWD(mat)}, o_{HBRS_MPL_FWD(o)}, sz_{HBRS_MPL_FWD(sz)}
	{
		/* BOOST_ASSERT( */
		/* 	(*less_equal)( */
		/* 		/1* plus(m(o_), m(sz_)), *1/ */
		/* 		/1* m(hbrs::mpl::size(mat_)) *1/ */
		/* 	) */
		/* ); */
		BOOST_ASSERT(o_.m() + sz_.m() <= mat_.size().m());
		
		/* BOOST_ASSERT( */
		/* 	(*less_equal)( */
		/* 		/1* plus(n(o_), n(sz_)), *1/ */
		/* 		/1* n(hbrs::mpl::size(mat_)) *1/ */
		/* 	) */
		/* ); */
		BOOST_ASSERT(o_.n() + sz_.n() <= mat_.size().n());
	}
	
	constexpr 
	submatrix(submatrix const&) = default;
	constexpr 
	submatrix(submatrix &&) = default;

	submatrix&
	operator=(submatrix const& M) {
		BOOST_ASSERT(m() == M.m());
		BOOST_ASSERT(n() == M.n());
		for (std::size_t i {0}; i < M.m(); ++i) {
			for (std::size_t j {0}; j < M.n(); ++j) {
				at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j));
			}
		}
		return *this;
	}

	submatrix&
	operator=(std::decay_t<Matrix> const& M) {
		BOOST_ASSERT(m() == M.m());
		BOOST_ASSERT(n() == M.n());
		for (std::size_t i {0}; i < M.m(); ++i) {
			for (std::size_t j {0}; j < M.n(); ++j) {
				at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j));
			}
		}
		return *this;
	}

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
	template<typename Ring>
	submatrix&
	operator=(detail::givens_rotation_expression<givens_rotation<Ring> const&, submatrix<Matrix,Offset,Size> const&> const& e) {
		if (&(e.t2()) != this) {
			*this = e.t2();
		}

		decltype(auto) i     {e.t1().i()};
		decltype(auto) k     {e.t1().k()};
		decltype(auto) theta {e.t1().theta()};

		BOOST_ASSERT(i < m());
		BOOST_ASSERT(k < m());

		for (std::size_t j {0}; j <= n() - 1; ++j) {
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
	template<typename Ring>
	submatrix&
	operator=(detail::givens_rotation_expression<submatrix<Matrix,Offset,Size> const&, givens_rotation<Ring> const&> const& e) {
		if (&(e.t1()) != this) {
			*this = e.t1();
		}

		decltype(auto) i     {e.t2().i()};
		decltype(auto) k     {e.t2().k()};
		decltype(auto) theta {e.t2().theta()};

		BOOST_ASSERT(i < n());
		BOOST_ASSERT(k < n());

		for (std::size_t j {0}; j <= m() - 1; ++j) {
			auto const tau1 {at(make_matrix_index(j, i))};
			auto const tau2 {at(make_matrix_index(j, k))};
			at(make_matrix_index(j, i)) = theta.at(0) * tau1 - theta.at(1) * tau2;
			at(make_matrix_index(j, k)) = theta.at(1) * tau1 + theta.at(0) * tau2;
		}
		return *this;
	}

	
	constexpr decltype(auto)
	size() const {
		return (sz_);
	};
	
	auto
	m() const {
		return sz_.m();
	}

	auto
	n() const {
		return sz_.n();
	}

	template<typename Index>
	constexpr decltype(auto)
	at(Index && i) {
		/* return (*hbrs::mpl::at)( */
		/* 	mat_, */
		/* 	make_matrix_index( */
		/* 		(*plus)(m(o_), m(HBRS_MPL_FWD(i))), */
		/* 		(*plus)(n(o_), n(HBRS_MPL_FWD(i))) */
		/* 	) */
		/* ); */
		return mat_.at(make_matrix_index(o_.m() + HBRS_MPL_FWD(i).m(), o_.n() + HBRS_MPL_FWD(i).n()));
	}
	
	template<typename Index>
	constexpr decltype(auto)
	at(Index && i) const {
		/* return (*hbrs::mpl::at)( */
		/* 	mat_, */
		/* 	make_matrix_index( */
		/* 		(*plus)(m(o_), m(HBRS_MPL_FWD(i))), */
		/* 		(*plus)(n(o_), n(HBRS_MPL_FWD(i))) */
		/* 	) */
		/* ); */
		return mat_.at(make_matrix_index(o_.m() + HBRS_MPL_FWD(i).m(), o_.n() + HBRS_MPL_FWD(i).n()));
	}
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) & {
		return make_subsequence(
			smr<submatrix &, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)},
			/* (*n)(o_), */
			o_.n(),
			/* (*n)(sz_) */
			sz_.n()
		);
	}
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) const& {
		return make_subsequence(
			smr<submatrix const&, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)},
			/* (*n)(o_), */
			o_.n(),
			/* (*n)(sz_) */
			sz_.n()
		);
	}
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) && {
		return make_subsequence(
			make_smr(std::move(*this), HBRS_MPL_FWD(i)),
			/* (*n)(o_), */
			o_.n(),
			/* (*n)(sz_) */
			sz_.n()
		);
	}
	
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const;
	decltype(auto)
	operator()(std::size_t const row, range<std::size_t,std::size_t> const& columns) const;
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns) &;

private:
	Matrix mat_;
	Offset const o_;
	Size const sz_;
};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template<typename Matrix, typename Offset, typename Size>
struct tag_of< hbrs::mpl::submatrix<Matrix, Offset, Size> > {
	using type = hbrs::mpl::submatrix_tag;
};

template <>
struct make_impl<hbrs::mpl::submatrix_tag> {
	
	template<typename Matrix, typename Offset, typename Size>
	static constexpr hbrs::mpl::submatrix< 
		std::decay_t<Matrix>,
		std::decay_t<Offset>,
		std::decay_t<Size>
	>
	apply(Matrix && mat, Offset && o, Size && sz) {
		return {HBRS_MPL_FWD(mat), HBRS_MPL_FWD(o), HBRS_MPL_FWD(sz)};
	}
};

/* namespace hana */ } /* namespace boost */ }


#endif // !HBRS_MPL_DT_CTSSRV_HPP
