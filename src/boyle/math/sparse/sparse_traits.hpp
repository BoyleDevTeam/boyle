/**
 * @file coo_matrix.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-06-23
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <memory_resource>

#ifndef BOYLE_USE_BOOST_UNORDERED
#include "boost/serialization/unordered_map.hpp"
#include <unordered_map>
#else
#include "boost/unordered/unordered_flat_map.hpp"
#endif

#include "boyle/common/utils/aligned_allocator.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/index_pair.hpp"

namespace boyle::math {

template <
    ScalarArithmetic Scalar, std::integral Index,
    Allocatory Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class DokMatrix;

template <
    ScalarArithmetic Scalar, std::integral Index,
    Allocatory Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class LilMatrix;

template <
    ScalarArithmetic Scalar, std::integral Index,
    Allocatory Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class CooMatrix;

template <
    ScalarArithmetic Scalar, std::integral Index,
    Allocatory Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class CscMatrix;

template <
    ScalarArithmetic Scalar, std::integral Index,
    Allocatory Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class CsrMatrix;

namespace pmr {

template <ScalarArithmetic Scalar, std::integral Index>
using DokMatrix = ::boyle::math::DokMatrix<Scalar, Index, std::pmr::polymorphic_allocator<Scalar>>;

template <ScalarArithmetic Scalar, std::integral Index>
using LilMatrix = ::boyle::math::LilMatrix<Scalar, Index, std::pmr::polymorphic_allocator<Scalar>>;

template <ScalarArithmetic Scalar, std::integral Index>
using CooMatrix = ::boyle::math::CooMatrix<Scalar, Index, std::pmr::polymorphic_allocator<Scalar>>;

template <ScalarArithmetic Scalar, std::integral Index>
using CscMatrix = ::boyle::math::CscMatrix<Scalar, Index, std::pmr::polymorphic_allocator<Scalar>>;

template <ScalarArithmetic Scalar, std::integral Index>
using CsrMatrix = ::boyle::math::CsrMatrix<Scalar, Index, std::pmr::polymorphic_allocator<Scalar>>;

} // namespace pmr

template <typename T>
struct SparseTraits final {
    static_assert(false, "SparseTraits is not implemented for non-Matrix type.");
};

template <ScalarArithmetic Scalar, std::integral Index, Allocatory Alloc>
struct SparseTraits<DokMatrix<Scalar, Index, Alloc>> final {
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
    using dictionary_type =
#ifndef BOYLE_USE_BOOST_UNORDERED
        std::unordered_map<
            IndexPair<index_type>, value_type, IndexPairHash<index_type>,
            IndexPairEqual<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<
                std::pair<const ::boyle::math::IndexPair<index_type>, value_type>>>;
#else
        boost::unordered::unordered_flat_map<
            IndexPair<index_type>, value_type, IndexPairHash<index_type>,
            IndexPairEqual<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<
                std::pair<const ::boyle::math::IndexPair<index_type>, value_type>>>;
#endif
};

template <ScalarArithmetic Scalar, std::integral Index, Allocatory Alloc>
struct SparseTraits<LilMatrix<Scalar, Index, Alloc>> final {
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
    using dictionary_type =
#ifndef BOYLE_USE_BOOST_UNORDERED
        std::unordered_map<
            index_type,
            std::unordered_map<
                index_type, value_type, std::hash<index_type>, std::equal_to<index_type>,
                typename std::allocator_traits<allocator_type>::template rebind_alloc<
                    std::pair<const index_type, value_type>>>,
            std::hash<index_type>, std::equal_to<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<std::pair<
                const index_type,
                std::unordered_map<
                    index_type, value_type, std::hash<index_type>, std::equal_to<index_type>,
                    typename std::allocator_traits<allocator_type>::template rebind_alloc<
                        std::pair<const index_type, value_type>>>>>>;
#else
        boost::unordered::unordered_flat_map<
            index_type,
            boost::unordered::unordered_flat_map<
                index_type, value_type, boost::hash<index_type>, std::equal_to<index_type>,
                typename std::allocator_traits<allocator_type>::template rebind_alloc<
                    std::pair<const index_type, value_type>>>,
            boost::hash<index_type>, std::equal_to<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<std::pair<
                const index_type,
                boost::unordered::unordered_flat_map<
                    index_type, value_type, boost::hash<index_type>, std::equal_to<index_type>,
                    typename std::allocator_traits<allocator_type>::template rebind_alloc<
                        std::pair<const index_type, value_type>>>>>>;
#endif
};

template <ScalarArithmetic Scalar, std::integral Index, Allocatory Alloc>
struct SparseTraits<CooMatrix<Scalar, Index, Alloc>> final {
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
    using dictionary_type =
#ifndef BOYLE_USE_BOOST_UNORDERED
        std::unordered_map<
            IndexPair<index_type>, size_type, IndexPairHash<index_type>, IndexPairEqual<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<
                std::pair<const ::boyle::math::IndexPair<index_type>, size_type>>>;
#else
        boost::unordered::unordered_flat_map<
            IndexPair<index_type>, size_type, IndexPairHash<index_type>, IndexPairEqual<index_type>,
            typename std::allocator_traits<allocator_type>::template rebind_alloc<
                std::pair<const ::boyle::math::IndexPair<index_type>, size_type>>>;
#endif
};

template <ScalarArithmetic Scalar, std::integral Index, Allocatory Alloc>
struct SparseTraits<CscMatrix<Scalar, Index, Alloc>> final {
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
};

template <ScalarArithmetic Scalar, std::integral Index, Allocatory Alloc>
struct SparseTraits<CsrMatrix<Scalar, Index, Alloc>> final {
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
};

} // namespace boyle::math
