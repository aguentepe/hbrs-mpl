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

#ifndef HBRS_MPL_DT_SUBSEQUENCE_HPP
#define HBRS_MPL_DT_SUBSEQUENCE_HPP

#include <hbrs/mpl/fwd/dt/subsequence.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/plus.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>
#include <boost/assert.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename Sequence, typename Offset, typename Size>
struct subsequence {
	
	template<typename Sequence_, typename Offset_, typename Size_>
	constexpr 
	subsequence(Sequence_ && seq, Offset_ && o, Size_ && sz)
	: seq_{HBRS_MPL_FWD(seq)}, o_{HBRS_MPL_FWD(o)}, sz_{HBRS_MPL_FWD(sz)}
	{
		BOOST_ASSERT(
			(*less_equal)(
				plus(o_, sz_), hbrs::mpl::size(seq_)
			)
		);
	}
	
	constexpr 
	subsequence(subsequence const&) = default;
	constexpr 
	subsequence(subsequence &&) = default;
	
	constexpr subsequence&
	operator=(subsequence const&) = default;
	constexpr subsequence&
	operator=(subsequence &&) = default;
	
	constexpr decltype(auto)
	size() const { return (sz_); };
	
	template<typename Index>
	constexpr decltype(auto)
	at(Index && i) {
		return (*hbrs::mpl::at)(seq_, plus(o_, HBRS_MPL_FWD(i)));
	}
	
	template<typename Index>
	constexpr decltype(auto)
	at(Index && i) const {
		return (*hbrs::mpl::at)(seq_, plus(o_, HBRS_MPL_FWD(i)));
	}
	
	template<typename Index>
	constexpr decltype(auto)
	operator[](Index && i) { return at(HBRS_MPL_FWD(i)); }
	
	template<typename Index>
	constexpr decltype(auto)
	operator[](Index && i) const { return at(HBRS_MPL_FWD(i)); }
	
private:
	Sequence seq_;
	Offset const o_;
	Size const sz_;
};

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template<typename Sequence, typename Offset, typename Size>
struct tag_of< hbrs::mpl::subsequence<Sequence, Offset, Size> > {
	using type = hbrs::mpl::subsequence_tag;
};

template <>
struct make_impl<hbrs::mpl::subsequence_tag> {
	
	template<typename Sequence, typename Offset, typename Size>
	static constexpr hbrs::mpl::subsequence< 
		std::decay_t<Sequence>,
		std::decay_t<Offset>,
		std::decay_t<Size>
	>
	apply(Sequence && seq, Offset && o, Size && sz) {
		return {HBRS_MPL_FWD(seq), HBRS_MPL_FWD(o), HBRS_MPL_FWD(sz)};
	}
};

/* namespace hana */ } /* namespace boost */ }


#endif // !HBRS_MPL_DT_CTSSRV_HPP
