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
	
	template<typename Matrix_>
	submatrix&
	operator=(Matrix_ && M) {
		BOOST_ASSERT(m() == HBRS_MPL_FWD(M).m());
		BOOST_ASSERT(n() == HBRS_MPL_FWD(M).n());
		for (std::size_t i {0}; i <= HBRS_MPL_FWD(M).m(); ++i) {
			for (std::size_t j {0}; j <= HBRS_MPL_FWD(M).n(); ++j) {
				at(make_matrix_index(i,j)) = HBRS_MPL_FWD(M).at(make_matrix_index(i,j));
			}
		}
		return *this;
	}

	
	constexpr decltype(auto)
	size() const { return (sz_); };
	
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
	
	/* template<typename Index> */
	/* constexpr auto */
	/* operator[](Index && i) & { */
	/* 	return make_subsequence( */
	/* 		smr<submatrix &, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}, */
	/* 		/1* (*n)(o_), *1/ */
	/* 		o_.n(), */
	/* 		/1* (*n)(sz_) *1/ */
	/* 		sz_.n() */
	/* 	); */
	/* } */
	
	/* template<typename Index> */
	/* constexpr auto */
	/* operator[](Index && i) const& { */
	/* 	return make_subsequence( */
	/* 		smr<submatrix const&, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}, */
	/* 		/1* (*n)(o_), *1/ */
	/* 		o_.n(), */
	/* 		/1* (*n)(sz_) *1/ */
	/* 		sz_.n() */
	/* 	); */
	/* } */
	
	/* template<typename Index> */
	/* constexpr auto */
	/* operator[](Index && i) && { */
	/* 	return make_subsequence( */
	/* 		make_smr(std::move(*this), HBRS_MPL_FWD(i)), */
	/* 		/1* (*n)(o_), *1/ */
	/* 		o_.n(), */
	/* 		/1* (*n)(sz_) *1/ */
	/* 		sz_.n() */
	/* 	); */
	/* } */
	
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const;
	decltype(auto)
	operator()(std::size_t const row, range<std::size_t,std::size_t> const& columns) const;
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns);

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
