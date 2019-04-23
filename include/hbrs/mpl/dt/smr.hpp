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

#ifndef HBRS_MPL_DT_SMR_HPP
#define HBRS_MPL_DT_SMR_HPP

#include <hbrs/mpl/fwd/dt/smr.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/size.hpp>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<
	typename Matrix,
	typename Index
>
struct smr {
	template<typename Matrix_, typename Index_>
	constexpr
	smr(Matrix_ && a, Index_ && row) : a_{HBRS_MPL_FWD(a)}, row_{HBRS_MPL_FWD(row)} {}

	constexpr
	smr(smr const&) = default;
	constexpr
	smr(smr &&) = default;
	
	constexpr smr&
	operator=(smr const&) = default;
	constexpr smr&
	operator=(smr &&) = default;

	
	/* Allow copy construction of a row from a row view */
	template<
		typename Matrix_,
		typename Index_,
		typename std::enable_if_t<
			!std::is_lvalue_reference<Matrix>::value &&
			std::is_lvalue_reference<Matrix_>::value
		>* = nullptr
	>
	constexpr
	smr(smr<Matrix, Index> const& rhs) : a_{rhs.a_}, row_{rhs.row_} {}

	template<typename Index_>
	constexpr decltype(auto)
	at(Index_ && column) {
		using hbrs::mpl::at;
		return (*at)(a_, make_matrix_index(row_, HBRS_MPL_FWD(column)));
	}

	template<typename Index_>
	constexpr decltype(auto)
	at(Index_ && column) const {
		using hbrs::mpl::at;
		return (*at)(a_, make_matrix_index(row_, HBRS_MPL_FWD(column)));
	}
	
	constexpr auto
	length() const {
		using mpl::m;
		return (*m)(size(a_));
	}
	
	template<typename Index_>
	constexpr decltype(auto)
	operator[](Index_ && column) {
		return at(HBRS_MPL_FWD(column));
	}
	
	template<typename Index_>
	constexpr decltype(auto)
	operator[](Index_ && column) const {
		return at(HBRS_MPL_FWD(column));
	}

private:
	Matrix a_;
	Index row_;
};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename Matrix, typename Index>
struct tag_of< hbrs::mpl::smr<Matrix, Index> > {
	using type = hbrs::mpl::smr_tag;
};

/* namespace hana */ } /* namespace boost */ }

namespace boost { namespace hana {

template <>
struct make_impl<hbrs::mpl::smr_tag> {
	template<typename Matrix, typename Index>
	static constexpr hbrs::mpl::smr< std::decay_t<Matrix>, std::decay_t<Index> >
	apply(Matrix && a, Index && row) {
		return {HBRS_MPL_FWD(a), HBRS_MPL_FWD(row)};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_SMR_HPP
