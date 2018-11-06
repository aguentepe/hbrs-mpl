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

#include <hbrs/mpl/fn/multiply.hpp>
#include <hbrs/mpl/dt/scv.hpp>
#include <vector>

ELEMENTAL_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

template auto multiply_impl_Matrix_Matrix::operator()(El::Matrix<float>               const&, El::Matrix<float>               const&) const;
template auto multiply_impl_Matrix_Matrix::operator()(El::Matrix<El::Complex<float>>  const&, El::Matrix<El::Complex<float>>  const&) const;
template auto multiply_impl_Matrix_Matrix::operator()(El::Matrix<double>              const&, El::Matrix<double>              const&) const;
template auto multiply_impl_Matrix_Matrix::operator()(El::Matrix<El::Complex<double>> const&, El::Matrix<El::Complex<double>> const&) const;

template auto multiply_impl_Matrix_scv_vector::operator()(El::Matrix<float>              const&, mpl::scv<std::vector<int>> const&) const;
//TODO: Solve "error: call to 'operator*' is ambiguous" for El::Complex<...> types
// template auto multiply_impl_Matrix_scv_vector::operator()(El::Matrix<El::Complex<float>> const&, mpl::scv<std::vector<int>> const&) const;
template auto multiply_impl_Matrix_scv_vector::operator()(El::Matrix<double>             const&, mpl::scv<std::vector<int>> const&) const;
// template auto multiply_impl_Matrix_scv_vector::operator()(El::Matrix<El::Complex<double>>const&, mpl::scv<std::vector<int>> const&) const;

template auto multiply_impl_Matrix_Scalar::operator()(El::Matrix<float>              , float               const&) const;
template auto multiply_impl_Matrix_Scalar::operator()(El::Matrix<El::Complex<float>> , El::Complex<float>  const&) const;
template auto multiply_impl_Matrix_Scalar::operator()(El::Matrix<double>             , double              const&) const;
template auto multiply_impl_Matrix_Scalar::operator()(El::Matrix<El::Complex<double>>, El::Complex<double> const&) const;

/* namespace detail */ }
ELEMENTAL_NAMESPACE_END