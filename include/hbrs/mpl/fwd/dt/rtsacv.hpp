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

#ifndef HBRS_MPL_FWD_DT_RTSACV_HPP
#define HBRS_MPL_FWD_DT_RTSACV_HPP

#include <hbrs/mpl/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>
#include <hbrs/mpl/dt/range.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace hana = boost::hana;

/* runtime-size array/continuous/dense/shared-memory column vector */
template<typename /* type of vector entries */ Ring>
struct rtsacv {
	explicit rtsacv(std::vector<Ring> data) : data_{data} {}
	
	explicit rtsacv(std::size_t size) : data_(size, Ring{0}) {
		BOOST_ASSERT(size >= 0);
	}

	explicit rtsacv(std::initializer_list<Ring> const& l);
	
	rtsacv(rtsacv const&) = default;
	rtsacv(rtsacv &&) = default;
	
	rtsacv&
	operator=(rtsacv const&) = default;
	rtsacv&
	operator=(rtsacv &&) = default;
	
	auto const
	size() const {
		return data_.size();
	}

	decltype(auto)
	vector() const {
		return data_;
	}

	auto m() const {
		return data_.size();
	}

	template<
		typename Index,
		typename std::enable_if_t<std::is_integral_v<typename std::decay_t<Index>>, int> = 0
	>
	decltype(auto)
	at(Index && i) {
		return data_[HBRS_MPL_FWD(i)];
	}
	
	template<
		typename Index,
		typename std::enable_if_t<std::is_integral_v<typename std::decay_t<Index>>, int> = 0
	>
	decltype(auto)
	at(Index && i) const {
		return data_[HBRS_MPL_FWD(i)];
	}

    auto const*
	begin() const {
        return &data_[0];
    }

	auto const*
	end() const {
		return &data_[data_.size()];
    }

	auto
	operator() (range<std::size_t,std::size_t> const& r) const ;

private:
	std::vector<Ring> data_;
};

struct rtsacv_tag{};
constexpr auto make_rtsacv = hana::make<rtsacv_tag>;
constexpr auto to_rtsacv = hana::to<rtsacv_tag>;

HBRS_MPL_NAMESPACE_END

#endif // !HBRS_MPL_FWD_DT_RTSACV_HPP
