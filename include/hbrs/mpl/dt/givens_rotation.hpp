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

#ifndef HBRS_MPL_DT_GIVENS_ROTATION_HPP
#define HBRS_MPL_DT_GIVENS_ROTATION_HPP

#include <hbrs/mpl/fwd/dt/givens_rotation.hpp>

#include <hbrs/mpl/preprocessor/core.hpp>
#include <array>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>

HBRS_MPL_NAMESPACE_BEGIN

template<typename /* type of matrix entries */ Ring>
struct givens_rotation {
	givens_rotation(std::size_t const i, std::size_t const k, std::array<Ring, 2> const theta) : i_{i}, k_{k}, theta_{theta} {}

	decltype(auto)
	i() const {
		return i_;
	}

	decltype(auto)
	k() const {
		return k_;
	}

	decltype(auto)
	theta() const {
		return (theta_);
	}
private:
	std::size_t i_, k_;
	std::array<Ring, 2> theta_;
};

template<typename Ring>
auto
G(std::size_t const i, std::size_t const k, std::array<Ring, 2> const theta) {
	return givens_rotation{i,k,theta};
}

namespace detail {

template<typename T1, typename T2>
struct givens_rotation_expression {
	template<typename T1_, typename T2_>
	constexpr
	givens_rotation_expression(T1_ && t1, T2_ && t2) : t1_{HBRS_MPL_FWD(t1)}, t2_{HBRS_MPL_FWD(t2)} {}

	constexpr decltype(auto)
	t1() const {
		return t1_;
	}

	constexpr decltype(auto)
	t2() const {
		return t2_;
	}
private:
	T1 t1_;
	T2 t2_;
};

/* namespace detail */ }

HBRS_MPL_NAMESPACE_END

namespace boost { namespace hana {

template <typename Ring>
struct tag_of< hbrs::mpl::givens_rotation<Ring> > {
	using type = hbrs::mpl::givens_rotation_tag;
};

template <>
struct make_impl<hbrs::mpl::givens_rotation_tag> {
	template <typename Ring>
	static hbrs::mpl::givens_rotation<Ring>
	apply(hana::basic_type<Ring>, std::size_t i, std::size_t k, std::array<Ring, 2> theta) {
		return {i,k,theta};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_MPL_DT_GIVENS_ROTATION_HPP
