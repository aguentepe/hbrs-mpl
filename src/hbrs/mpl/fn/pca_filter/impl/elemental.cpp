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

#include <hbrs/mpl/fn/pca_filter.hpp>

#include <hbrs/mpl/dt/el_matrix.hpp>
#include <hbrs/mpl/dt/el_dist_matrix.hpp>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {

//TODO: Why does El::Complex<...> fail?

template auto pca_filter_impl_matrix::operator()(matrix<float>               const&,  std::function<bool(int)> const&) const;
// template auto pca_filter_impl_matrix::operator()(matrix<El::Complex<float>>  const&,  std::function<bool(int)> const&) const;
template auto pca_filter_impl_matrix::operator()(matrix<double>              const&,  std::function<bool(int)> const&) const;
// template auto pca_filter_impl_matrix::operator()(matrix<El::Complex<double>> const&,  std::function<bool(int)>const& ) const;

template auto pca_filter_impl_matrix::operator()(dist_matrix<float>               const&,  std::function<bool(int)> const&) const;
// template auto pca_filter_impl_matrix::operator()(dist_matrix<El::Complex<float>>  const&,  std::function<bool(int)> const&) const;
template auto pca_filter_impl_matrix::operator()(dist_matrix<double>              const&,  std::function<bool(int)> const&) const;
// template auto pca_filter_impl_matrix::operator()(dist_matrix<El::Complex<double>> const&,  std::function<bool(int)> const&) const;

template auto pca_filter_impl_matrix::operator()(matrix<float>               const&,  std::vector<bool> const&) const;
// template auto pca_filter_impl_matrix::operator()(matrix<El::Complex<float>>  const&,  std::vector<bool> const&) const;
template auto pca_filter_impl_matrix::operator()(matrix<double>              const&,  std::vector<bool> const&) const;
// template auto pca_filter_impl_matrix::operator()(matrix<El::Complex<double>> const&,  std::vector<bool> const&) const;

template auto pca_filter_impl_matrix::operator()(dist_matrix<float>               const&,  std::vector<bool> const&) const;
// template auto pca_filter_impl_matrix::operator()(dist_matrix<El::Complex<float>>  const&,  std::vector<bool> const&) const;
template auto pca_filter_impl_matrix::operator()(dist_matrix<double>              const&,  std::vector<bool> const&) const;
// template auto pca_filter_impl_matrix::operator()(dist_matrix<El::Complex<double>> const&,  std::vector<bool> const&) const;

/* namespace detail */ }
HBRS_MPL_NAMESPACE_END