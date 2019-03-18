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

#define BOOST_TEST_MODULE svd_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/fn/svd.hpp>
#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/mpl/detail/gather.hpp>
#include <hbrs/mpl/dt/ctsav.hpp>
#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/range.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/matrix_index.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/decompose_mode.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/svd.hpp>
#include <hbrs/mpl/fn/transpose.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/select.hpp>
#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
//	#include <elemental/dt/matrix.hpp>
#endif
#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
//	#include <matlab/dt/matrix.hpp>
#endif
#include <boost/hana/tuple.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/range.hpp>
#include <boost/hana/plus.hpp>
#include <boost/hana/length.hpp>
#include <boost/hana/minus.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/drop_back.hpp>
#include <boost/hana/value.hpp>
#include <boost/hana/less_equal.hpp>
#include <boost/hana/integral_constant.hpp>

#include "../data.hpp"
#include "../detail.hpp"

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(svd_test)

//using hbrs::mpl::detail::environment_fixture;
//BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(svd_matrix) {
	using namespace hbrs::mpl;
	rtsam<double, storage_order::row_major> A{
		{1, 2, 3,
		 4, 5, 6,
		 7, 8, 9}, make_matrix_size(3,3)};
	//     Matrix B{{1, 2, 3, 4,
	//               5, 6, 7, 8,
	//               9, 10, 11, 12}, 3};
	rtsam<double, storage_order::row_major> C{
		{1, 2,  3,
		 4, 5,  6,
		 7, 8,  9,
		10, 11, 12}, make_matrix_size(4,3)};
	rtsam<double, storage_order::row_major> D{
		{1, 0,  0,  0,
		 0, 6,  7,  0,
		 0, 0, 11, 12,
		 0, 0,  0,  0}, make_matrix_size(4,4)};
			  
	auto ASVD{svd(A,0)};
	auto CSVD{svd(C,0)};
	auto DSVD{svd(D,0)};

	auto AA   {ASVD.u() * ASVD.s() * transpose(ASVD.v())};
	BOOST_TEST(AA == A);
	BOOST_TEST(CSVD.u() * CSVD.s() * transpose(CSVD.v()) == C);
	auto DD   {DSVD.u() * DSVD.s() * transpose(DSVD.v())};
	BOOST_TEST(DD == D);
}

BOOST_AUTO_TEST_CASE(svd_comparison, * utf::tolerance(0.000000001)) {
	/*using namespace hbrs::mpl;*/
	/*using hbrs::mpl::decompose_mode;*/
	
	/*static constexpr auto datasets = hana::make_tuple(*/
	/*	make_sm(*/
	/*		make_ctsav(test::mat_a), make_matrix_size(hana::size_c<test::mat_a_m>, hana::size_c<test::mat_a_n>), row_major_c*/
	/*	),*/
	/*	make_sm(*/
	/*		make_ctsav(test::mat_g), make_matrix_size(hana::size_c<test::mat_g_m>, hana::size_c<test::mat_g_n>), row_major_c*/
	/*	),*/
	/*	make_sm(*/
	/*		make_ctsav(test::mat_j), make_matrix_size(hana::size_c<test::mat_j_m>, hana::size_c<test::mat_j_n>), row_major_c*/
	/*	)*/
	/*);*/
	
	/*static constexpr auto modes = hana::make_tuple(complete_decompose_c, economy_decompose_c, zero_decompose_c);*/
	
	/*std::size_t dataset_nr = 0;*/
	/*hana::for_each(datasets, [&dataset_nr](auto const& dataset) {*/
	/*	BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);*/
	/*	hana::for_each(modes, [&dataset](auto const& mode) {*/
	/*		 BOOST_TEST_MESSAGE(*/
	/*			 "decompose_mode=" <<*/ 
	/*			 (*/
	/*				(mode == decompose_mode::complete)*/ 
	/*				? "complete"*/ 
	/*				: (*/
	/*					(mode == decompose_mode::economy)*/ 
	/*					? "economy"*/ 
	/*					: "zero"*/
	/*				)*/
	/*			)*/
	/*		);*/
			
	/*		auto funs = hana::drop_back(hana::make_tuple(*/
	/*			#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB*/
	/*			[](auto && a, auto mode) {*/
	/*				return matlab::detail::svd_impl_level0{}(matlab::make_matrix(HBRS_MPL_FWD(a)), mode);*/
	/*			},*/
	/*			[](auto && a, auto mode) {*/
	/*				return matlab::detail::svd_impl_level1{}(matlab::make_matrix(HBRS_MPL_FWD(a)), mode);*/
	/*			},*/
	/*			#endif*/
				
	/*			#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL*/
	/*			[](auto && a, auto mode) {*/
	/*				return elemental::detail::svd_impl_Matrix{}(elemental::make_matrix(HBRS_MPL_FWD(a)), mode);*/
	/*			},*/
	/*			#endif*/
	/*			"SEQUENCE_TERMINATOR___REMOVED_BY_DROP_BACK"*/
	/*		));*/
			
	/*		auto results = hana::transform(funs, [&dataset, &mode](auto f) { return f(dataset, mode); });*/
			
	/*		if constexpr(hana::length(results) >= 2u) {*/
	/*			auto compare = [&dataset, &mode, &results](auto i) {*/
	/*				using hbrs::mpl::select;*/
					
	/*				auto j = i+hana::ushort_c<1>;*/
	/*				auto const& svd_i = hana::at(results, i);*/
	/*				auto const& svd_j = hana::at(results, j);*/
	/*				//TODO: Take ++j if j not supported!*/
					
	/*				if constexpr(*/
	/*					is_supported(svd_i) && is_supported(svd_j)*/
	/*				) {*/
	/*					auto sz_ = (*size)(dataset);*/
	/*					auto m_ = (*m)(sz_);*/
	/*					auto n_ = (*n)(sz_);*/
						
	/*					BOOST_TEST_MESSAGE("comparing svd_u of impl nr " << i << " and " << j);*/
	/*					auto const& svd_u_i = (*at)(svd_i, svd_u{});*/
	/*					auto const& svd_u_j = (*at)(svd_j, svd_u{});*/
	/*					if ((m_ > n_) && (mode == decompose_mode::complete)) {*/
	/*						auto rng = std::make_pair(make_matrix_index(0,0), make_matrix_size((int)m_, (int)n_));*/
	/*						auto const svd_u_i_mxn = (*select)(svd_u_i, rng);*/
	/*						auto const svd_u_j_mxn = (*select)(svd_u_j, rng);*/
	/*						HBRS_MPL_TEST_MMEQ(svd_u_i_mxn, svd_u_j_mxn, true);*/
	/*					} else {*/
	/*						HBRS_MPL_TEST_MMEQ(svd_u_i, svd_u_j, true);*/
	/*					}*/
						
						
	/*					BOOST_TEST_MESSAGE("comparing svd_s of impl nr " << i << " and " << j);*/
	/*					auto const& svd_s_i = (*at)(svd_i, svd_s{});*/
	/*					auto const& svd_s_j = (*at)(svd_j, svd_s{});*/
	/*					HBRS_MPL_TEST_MMEQ(svd_s_i, svd_s_j, false);*/
						
						
	/*					BOOST_TEST_MESSAGE("comparing svd_v of impl nr " << i << " and " << j);*/
	/*					auto const& svd_v_i = (*at)(svd_i, svd_v{});*/
	/*					auto const& svd_v_j = (*at)(svd_j, svd_v{});*/
	/*					if (*/
	/*						(m_ < n_) &&*/ 
	/*						((mode == decompose_mode::complete) || (mode == decompose_mode::zero))*/
	/*					) {*/
	/*						auto rng = std::make_pair(make_matrix_index(0,0), make_matrix_size((int)n_, (int)m_));*/
	/*						auto const svd_v_i_mxn = (*select)(svd_v_i, rng);*/
	/*						auto const svd_v_j_mxn = (*select)(svd_v_j, rng);*/
	/*						HBRS_MPL_TEST_MMEQ(svd_v_i_mxn, svd_v_j_mxn, true);*/
	/*					} else {*/
	/*						HBRS_MPL_TEST_MMEQ(svd_v_i, svd_v_j, true);*/
	/*					}*/
						
	/*					BOOST_TEST_MESSAGE("comparing impl nr " << i << " and " << j << " done.");*/
	/*				}*/
	/*			};*/
				
	/*			hana::for_each(*/
	/*				hana::make_range(*/
	/*					hana::size_c<0>,*/ 
	/*					hana::length(results)-hana::size_c<1>*/
	/*				),*/
	/*				compare*/
	/*			);*/
	/*		}*/
			
	/*		{*/
	/*			std::size_t f_nr = 0;*/
	/*			auto rebuilds = hana::transform(results, [&dataset, &mode, &f_nr](auto && svd_result) {*/
	/*				if constexpr(is_not_supported(svd_result)) {*/
	/*					return hbrs::mpl::detail::not_supported{};*/
	/*				} else {*/
	/*					auto sz_ = (*size)(dataset);*/
	/*					auto m_ = (*m)(sz_);*/
	/*					auto n_ = (*n)(sz_);*/
						
	/*					auto && u = (*at)(svd_result, svd_u{});*/
	/*					auto && s = (*at)(svd_result, svd_s{});*/
	/*					auto && v = (*at)(svd_result, svd_v{});*/
						
	/*					auto u_sz = (*size)(u);*/
	/*					auto u_m = (*m)(u_sz);*/
	/*					auto u_n = (*n)(u_sz);*/
	/*					auto s_sz = (*size)(s);*/
	/*					auto s_m = (*m)(s_sz);*/
	/*					auto s_n = (*n)(s_sz);*/
	/*					auto v_sz = (*size)(v);*/
	/*					auto v_m = (*m)(v_sz);*/
	/*					auto v_n = (*n)(v_sz);*/
						
	/*					if (*/
	/*						(mode == decompose_mode::complete) ||*/
	/*						((mode == decompose_mode::zero) && (m_ <= n_))*/
	/*					) {*/
	/*						BOOST_TEST( (*equal)(u_m,m_) );*/
	/*						BOOST_TEST( (*equal)(u_n,m_) );*/
							
	/*						BOOST_TEST( (*equal)(s_m,m_) );*/
	/*						BOOST_TEST( (*equal)(s_n,n_) );*/
							
	/*						BOOST_TEST( (*equal)(v_m,n_) );*/
	/*						BOOST_TEST( (*equal)(v_n,n_) );*/
	/*					} else if (*/
	/*						(mode == decompose_mode::economy) ||*/
	/*						((mode == decompose_mode::zero) && (m_ > n_))*/
	/*					) {*/
	/*						BOOST_TEST( (*equal)(u_m,m_) );*/
	/*						if (m_ > n_) {*/
	/*							BOOST_TEST( (*equal)(u_n,n_) );*/
	/*						} else {*/
	/*							BOOST_TEST( (*equal)(u_n,m_) );*/
	/*						}*/
							
	/*						if (m_ > n_) {*/
	/*							BOOST_TEST( (*equal)(s_m,n_) );*/
	/*							BOOST_TEST( (*equal)(s_n,n_) );*/
	/*						} else if (m_ == n_) {*/
	/*							BOOST_TEST( (*equal)(s_m,m_) );*/
	/*							BOOST_TEST( (*equal)(s_n,n_) );*/
	/*						} else if (m_ < n_) {*/
	/*							BOOST_TEST( (*equal)(s_m,m_) );*/
	/*							BOOST_TEST( (*equal)(s_n,m_) );*/
	/*						}*/
							
	/*						BOOST_TEST( (*equal)(v_m,n_) );*/
	/*						if (m_ < n_) {*/
	/*							BOOST_TEST( (*equal)(v_n,m_) );*/
	/*						} else {*/
	/*							BOOST_TEST( (*equal)(v_n,n_) );*/
	/*						}*/
	/*					} else {*/
	/*						BOOST_ASSERT(false);*/
	/*					}*/
						
	/*					if (*/
	/*						(m_ > n_) &&*/ 
	/*						((mode==decompose_mode::economy) || (mode==decompose_mode::zero))*/
	/*					) {*/
	/*						BOOST_TEST_MESSAGE("Skipping U*U'==I for f_nr=" << f_nr);*/
	/*					} else {*/
	/*						BOOST_TEST_MESSAGE("Testing U*U'==I for f_nr=" << f_nr);*/
	/*						auto uxut = (*multiply)(u, transpose(u));*/
	/*						BOOST_TEST((*size)(uxut) == make_matrix_size(m_, m_));*/
	/*						HBRS_MPL_TEST_IS_IDENTITY(uxut);*/
	/*					}*/
						
	/*					if ( (m_ < n_) && (mode==decompose_mode::economy) ) {*/
	/*						BOOST_TEST_MESSAGE("Skipping V*V'==I for f_nr=" << f_nr);*/
	/*					} else {*/
	/*						BOOST_TEST_MESSAGE("Testing V*V'==I for f_nr=" << f_nr);*/
	/*						auto vxvt = (*multiply)(v, transpose(v));*/
	/*						BOOST_TEST((*size)(vxvt) == make_matrix_size(n_, n_));*/
	/*						HBRS_MPL_TEST_IS_IDENTITY(vxvt);*/
	/*					}*/
						
	/*					BOOST_TEST_MESSAGE("Testing unitary matrices for f_nr=" << f_nr << " done");*/
	/*					++f_nr;*/
						
	/*					return (*multiply)(u, multiply(s,transpose(v)));*/
	/*				}*/
	/*			});*/
				
	/*			hana::for_each(*/
	/*				hana::make_range(hana::ushort_c<0>, hana::length(results)),*/
	/*				[&dataset, &rebuilds](auto i) {*/
	/*					BOOST_TEST_MESSAGE("comparing original and reconstructed dataset of impl nr " << i);*/
						
	/*					auto const& rebuild = hana::at(rebuilds, i);*/
	/*					if constexpr(is_supported(rebuild)) {*/
	/*						HBRS_MPL_TEST_MMEQ(dataset, rebuild, false);*/
	/*					}*/
	/*				}*/
	/*			);*/
	/*			BOOST_TEST_MESSAGE("comparing original and reconstructed datasets done.");*/
	/*		}*/
	/*	});*/
		
	/*	++dataset_nr;*/
	/*});*/
	
/*/**/

/*complete list of test failures:*/

/*3x3 matrix A*/
  /*decompose_mode=complete / decompose_mode=economy / decompose_mode=zero*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s,v> == el_mat0 <u,s,v>*/
    /*matlab0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab1: U*U'==I*/ 
             /*V*V'==I*/
    /*el_mat0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/


/*13x4 matrix A*/
  /*decompose_mode=complete*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u> != el_mat0 <u> && matlab1 <s,v> == el_mat0 <s,v>*/
    /*matlab0: U*U'==I*/
             /*V*V'==I*/
    /*matlab1: U*U'==I*/
             /*V*V'==I*/
    /*el_mat0: U*U'==I*/
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/

  /*decompose_mode=economy*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s,v> == el_mat0 <u,s,v>*/
    /*matlab0: U*U'!=I*/
             /*V*V'==I*/
    /*matlab1: U*U'!=I*/
             /*V*V'==I*/
    /*el_mat0: U*U'!=I*/
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/

  /*decompose_mode=zero*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s,v> == el_mat0 <u,s,v>*/
    /*matlab0: U*U'!=I*/
             /*V*V'==I*/
    /*matlab1: U*U'!=I*/
             /*V*V'==I*/
    /*el_mat0: U*U'!=I*/
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/


/*4x13 matrix A*/
  /*decompose_mode=complete*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s> == el_mat0 <u,s> && matlab1 <v> != el_mat0 <v>*/
    /*matlab0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab1: U*U'==I*/ 
             /*V*V'==I*/
    /*el_mat0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/

  /*decompose_mode=economy*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s,v> == el_mat0 <u,s,v>*/
    /*matlab0: U*U'==I*/ 
             /*V*V'!=I*/
    /*matlab1: U*U'==I*/ 
             /*V*V'!=I*/
    /*el_mat0: U*U'==I*/ 
             /*V*V'!=I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/

  /*decompose_mode=zero*/
    /*matlab0 <u,s,v> == matlab1 <u,s,v>*/
    /*matlab1 <u,s> == el_mat0 <u,s> && matlab1 <v> != el_mat0 <v>*/
    /*matlab0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab1: U*U'==I*/ 
             /*V*V'==I*/
    /*el_mat0: U*U'==I*/ 
             /*V*V'==I*/
    /*matlab0: U*S*V'== A*/
    /*matlab1: U*S*V'== A*/
    /*el_mat0: U*S*V'== A*/


/*summary of expected test failures:*/
/*13x4 matrix A : decompose_mode=complete : matlab1 <u> != el_mat0 <u> && matlab1 <s,v> == el_mat0 <s,v>*/
/*13x4 matrix A : decompose_mode=economy : matlab0: U*U'!=I*/
/*13x4 matrix A : decompose_mode=economy : matlab1: U*U'!=I*/
/*13x4 matrix A : decompose_mode=economy : el_mat0: U*U'!=I*/
/*13x4 matrix A : decompose_mode=zero : matlab0: U*U'!=I*/
/*13x4 matrix A : decompose_mode=zero : matlab1: U*U'!=I*/
/*13x4 matrix A : decompose_mode=zero : el_mat0: U*U'!=I*/
/*4x13 matrix A : decompose_mode=complete : matlab1 <u,s> == el_mat0 <u,s> && matlab1 <v> != el_mat0 <v>*/
/*4x13 matrix A : decompose_mode=economy : matlab0: V*V'!=I*/
/*4x13 matrix A : decompose_mode=economy : matlab1: V*V'!=I*/
/*4x13 matrix A : decompose_mode=economy : el_mat0: V*V'!=I*/
/*4x13 matrix A : decompose_mode=zero : matlab1 <u,s> == el_mat0 <u,s> && matlab1 <v> != el_mat0 <v>*/

/*reasons for test failures:*/ 
 /*- U*U'!=I and V*V'!=I because in decompose_mode=economy n columns of U / m columns of V are returned thus U/V are not*/ 
   /*necessarly orthogonal anymore.*/
 /*- matlab1 <u> != el_mat0 <u> && matlab1 <s,v> == el_mat0 <s,v>*/ 
   /*and*/
   /*matlab1 <u,s> == el_mat0 <u,s> && matlab1 <v> != el_mat0 <v>*/
   /*because in decompose_mode=complete or decompose_mode=zero with m<n only the first min(m,n) columns of u/v must be equal.*/
   /*There may exist multiple orthogonal matrices where the first min(m,n) columns are equal and the last*/ 
   /*(max(m,n)-min(m,n)) columns might be different.*/


 /**/
}

BOOST_AUTO_TEST_SUITE_END()

