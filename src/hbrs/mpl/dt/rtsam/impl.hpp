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
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/detail/translate_index.hpp>
#include <hbrs/mpl/dt/exception.hpp>

#include <boost/hana/core/make.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/type.hpp>
#include <boost/assert.hpp>
#include <vector>

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
	
	auto const&
	size() const { return size_; };
	
	decltype(auto)
	order() const { return storage_order_c<Order>; }
	
	template<typename Index>
	decltype(auto)
	at(Index && i) {
		return data_[
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		];
	}
	
	template<typename Index>
	decltype(auto)
	at(Index && i) const {
		return data_[
			detail::translate_index(size_, HBRS_MPL_FWD(i), storage_order_c<Order>)
		];
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
	
private:
	std::vector<Ring> data_;
	matrix_size<std::size_t, std::size_t> size_;
};

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
