/**
 * @file triplet.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-12-16
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <ostream>
#include <tuple>
#include <utility>

#include "boost/serialization/access.hpp"

#include "boyle/math/concepts.hpp"

#define DECLARE_TRIPLET(ClassName, VAR_0, VAR_1, VAR_2)                                         \
    template <typename T0, typename T1 = T0, typename T2 = T0>                                  \
    class ClassName final : public ::boyle::math::Triplet<T0, T1, T2> {                         \
        template <unsigned int, typename T>                                                     \
        friend struct TripletElement;                                                           \
        using ::boyle::math::Triplet<T0, T1, T2>::m_first;                                      \
        using ::boyle::math::Triplet<T0, T1, T2>::m_second;                                     \
        using ::boyle::math::Triplet<T0, T1, T2>::m_third;                                      \
                                                                                                \
      public:                                                                                   \
        [[using gnu: always_inline]]                                                            \
        ClassName() noexcept = default;                                                         \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(T0 cv) noexcept                                                     \
            requires std::same_as<T0, T1> && std::same_as<T0, T2>                               \
            : ::boyle::math::Triplet<T0>{cv} {}                                                 \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(T0 c##VAR_0, T1 c##VAR_1, T2 c##VAR_2) noexcept                     \
            : ::boyle::math::Triplet<T0, T1, T2>{                                               \
                  std::move(c##VAR_0), std::move(c##VAR_1), std::move(c##VAR_2)                 \
              } {}                                                                              \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(const ClassName& other) noexcept = default;                         \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(ClassName&& other) noexcept = default;                              \
        [[using gnu: always_inline]]                                                            \
        constexpr auto operator=(const ClassName& other) noexcept -> ClassName& = default;      \
        [[using gnu: always_inline]]                                                            \
        constexpr auto operator=(ClassName&& other) noexcept -> ClassName& = default;           \
        ~ClassName() noexcept override = default;                                               \
                                                                                                \
        template <                                                                              \
            std::convertible_to<T0> U0, std::convertible_to<T1> U1, std::convertible_to<T2> U2> \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(ClassName<U0, U1, U2> other) noexcept                               \
            : ::boyle::math::Triplet<T0, T1, T2>{                                               \
                  std::move(other.m_first), std::move(other.m_second), std::move(m_third)       \
              } {}                                                                              \
                                                                                                \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(std::tuple<T0, T1, T2> tuple) noexcept                              \
            : ::boyle::math::Triplet<T0, T1, T2>{std::move(tuple)} {}                           \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_0() noexcept -> T0& {                                                \
            return m_first;                                                                     \
        }                                                                                       \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_0() const noexcept -> const T0& {                                    \
            return m_first;                                                                     \
        }                                                                                       \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_1() noexcept -> T1& {                                                \
            return m_second;                                                                    \
        }                                                                                       \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_1() const noexcept -> const T1& {                                    \
            return m_second;                                                                    \
        }                                                                                       \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_2() noexcept -> T2& {                                                \
            return m_third;                                                                     \
        }                                                                                       \
                                                                                                \
        [[using gnu: pure, always_inline, leaf]]                                                \
        constexpr auto VAR_2() const noexcept -> const T2& {                                    \
            return m_third;                                                                     \
        }                                                                                       \
    }

namespace boyle::math {

template <typename T0, typename T1 = T0, typename T2 = T0>
class Triplet {
    template <typename CharT, typename U0, typename U1, typename U2>
    friend inline auto operator<<(
        std::basic_ostream<CharT>& os, const Triplet<U0, U1, U2>& obj
    ) noexcept -> std::basic_ostream<CharT>&;
    friend class boost::serialization::access;
    template <unsigned int, typename T>
    friend struct TripletElement;

  public:
    [[using gnu: always_inline]]
    Triplet() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Triplet(T0 cv) noexcept
        requires std::same_as<T0, T1> && std::same_as<T0, T2>
        : m_first{cv}, m_second{cv}, m_third{cv} {}
    [[using gnu: always_inline]]
    constexpr Triplet(T0 cv0, T1 cv1, T2 cv2) noexcept
        : m_first{std::move(cv0)}, m_second{std::move(cv1)}, m_third{std::move(cv2)} {}
    [[using gnu: always_inline]]
    constexpr Triplet(const Triplet& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Triplet(Triplet&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Triplet& other) noexcept -> Triplet& = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Triplet&& other) noexcept -> Triplet& = default;
    virtual ~Triplet() noexcept = 0;

    [[using gnu: always_inline]]
    constexpr Triplet(std::tuple<T0, T1, T2> tuple) noexcept
        : m_first{std::move(std::get<0>(tuple))}, m_second{std::move(std::get<1>(tuple))},
          m_third{std::move(std::get<2>(tuple))} {}

    [[using gnu: pure, always_inline]]
    constexpr operator std::tuple<T0, T1, T2>() const noexcept {
        return std::make_tuple<T0, T1, T2>(m_first, m_second, m_third);
    }

  protected:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_first;
        archive & m_second;
        archive & m_third;
        return;
    }

    T0 m_first{};
    T1 m_second{};
    T2 m_third{};
};

template <typename T0, typename T1, typename T2>
inline Triplet<T0, T1, T2>::~Triplet() noexcept = default;

template <typename CharT, typename T0, typename T1, typename T2>
[[using gnu: always_inline]]
inline auto operator<<(
    std::basic_ostream<CharT>& os, const boyle::math::Triplet<T0, T1, T2>& obj
) noexcept -> std::basic_ostream<CharT>& {
    os << "(" << obj.m_first << ", " << obj.m_second << ", " << obj.m_third << ")";
    return os;
}

template <typename T>
concept Tripletic = derivedFromTemplate<T, Triplet>;

template <unsigned int Num, typename T>
struct TripletElement;

template <Tripletic T>
struct TripletElement<0, T> {
    using type = decltype(std::declval<T>().m_first);
};

template <Tripletic T>
struct TripletElement<1, T> {
    using type = decltype(std::declval<T>().m_second);
};

template <Tripletic T>
struct TripletElement<2, T> {
    using type = decltype(std::declval<T>().m_third);
};

template <unsigned int Num, Tripletic T>
using TripletElementT = typename TripletElement<Num, T>::type;

} // namespace boyle::math
