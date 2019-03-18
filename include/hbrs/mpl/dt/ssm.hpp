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

#ifndef HBRS_MPL_DT_SSM_HPP
#define HBRS_MPL_DT_SSM_HPP

#include <hbrs/mpl/fwd/dt/ssm.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/detail/translate_index.hpp>

#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/dt/smr.hpp>
#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<
	typename Matrix, 
	typename Index1,
	typename Index2
>
struct ssm {
	
	template<typename Matrix_, typename Index1_, typename Index2_>
	/* constexpr */ 
	ssm(Matrix_ && matrix, Index1_ && index1, Index2_ && index2)
	: matrix_{HBRS_MPL_FWD(matrix)}, index1_{HBRS_MPL_FWD(index1)}, index2_{HBRS_MPL_FWD(index2)} {
		BOOST_ASSERT(index1.m() <= index2.m());
		BOOST_ASSERT(index1.n() <= index2.n());
		BOOST_ASSERT(index2.m() < matrix_.m());
		BOOST_ASSERT(index2.n() < matrix_.n());
	}
	
	constexpr 
	ssm(ssm const&) = default;
	constexpr 
	ssm(ssm &&) = default;
	
	constexpr ssm&
	operator=(ssm const&) = default;
	constexpr ssm&
	operator=(ssm &&) = default;
	
	template<
		typename Ring,
		storage_order Order
	>
	ssm&
	operator=(rtsam<Ring,Order> const& M) {
		BOOST_ASSERT(m() == M.m());
		BOOST_ASSERT(n() == M.n());
		for (std::size_t i {0}; i <= M.m(); ++i) {
			for (std::size_t j {0}; j <= M.n(); ++j) {
				matrix_.at(make_matrix_index(i,j)) = M.at(make_matrix_index(i,j));
			}
		}
		return *this;
	}

	/* constexpr decltype(auto) */
	/* size() const { return (size_); }; */
	
	auto
	m() const {
		return index2_.m() - index1_.m() + 1;
	}

	auto
	n() const {
		return index2_.n() - index1_.n() + 1;
	}

	/* auto */
	/* operator std::remove_reference_t<Matrix>() const { */
	/* 	Matrix result{m(),n()}; */
	/* 	for (std::size_t i {0}; i <= m(); ++i) { */
	/* 		for (std::size_t j {0}; j <= n(); ++j) { */
	/* 			result.at(make_matrix_index(i,j)) = at(make_matrix_index(i,j)); */
	/* 		} */
	/* 	} */
	/* 	return result; */
	/* } */

	template<
		typename Index,
		typename std::enable_if_t<std::is_same_v<matrix_index<std::size_t, std::size_t>, Index>>* = nullptr
	>
	constexpr decltype(auto)
	at(Index && i) {
		return matrix_.at(make_matrix_index(index1_.m() + HBRS_MPL_FWD(i).m(), index1_.n() + HBRS_MPL_FWD(i).n()));
	}
	
	template<
		typename Index,
		typename std::enable_if_t<std::is_same_v<matrix_index<std::size_t, std::size_t>, Index>>* = nullptr
	>
	constexpr decltype(auto)
	at(Index && i) const {
		return matrix_.at(make_matrix_index(index1_.m() + HBRS_MPL_FWD(i).m(), index1_.n() + HBRS_MPL_FWD(i).n()));
	}
	
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) & { return smr<ssm &, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) const& { return smr<ssm const&, std::decay_t<Index>>{*this, HBRS_MPL_FWD(i)}; }
	
	template<typename Index>
	constexpr auto
	operator[](Index && i) && { return make_smr(std::move(*this), HBRS_MPL_FWD(i)); }
	
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, std::size_t const column) const;
	decltype(auto)
	operator()(std::size_t const& row, range<std::size_t,std::size_t> const& columns) const;
	decltype(auto)
	operator()(range<std::size_t,std::size_t> const& rows, range<std::size_t,std::size_t> const& columns);
private:
	Matrix matrix_;
	Index1 const index1_;
	Index2 const index2_;
};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template<
	typename Matrix, 
	typename Index1,
	typename Index2
>
struct tag_of< hbrs::mpl::ssm<Matrix, Index1, Index2> > {
	using type = hbrs::mpl::ssm_tag;
};

template <>
struct make_impl<hbrs::mpl::ssm_tag> {
	
	template<
		typename Matrix, 
		typename Index1,
		typename Index2
	>
	static constexpr hbrs::mpl::ssm< 
		std::decay_t<Matrix>,
		std::decay_t<Index1>,
		std::decay_t<Index2>
	>
	apply(Matrix && seq, Index1 && index1, Index2 && index2) {
		return {HBRS_MPL_FWD(seq), HBRS_MPL_FWD(index1), HBRS_MPL_FWD(index2)};
	}
};

/* namespace hana */ } /* namespace boost */ }


#endif // !HBRS_MPL_DT_SSM_HPP
