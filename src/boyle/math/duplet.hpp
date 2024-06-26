/**
 * @file duplet.hpp
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
#include <utility>

#include "boost/serialization/access.hpp"

#include "boyle/math/concepts.hpp"

#define DECLARE_DUPLET(ClassName, VAR_0, VAR_1)                                                 \
    template <typename T0, typename T1 = T0>                                                    \
    class ClassName final : public ::boyle::math::Duplet<T0, T1> {                              \
        template <unsigned int, typename T>                                                     \
        friend struct DupletElement;                                                            \
        using ::boyle::math::Duplet<T0, T1>::m_first;                                           \
        using ::boyle::math::Duplet<T0, T1>::m_second;                                          \
                                                                                                \
      public:                                                                                   \
        [[using gnu: always_inline]]                                                            \
        ClassName() noexcept = default;                                                         \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(T0 cv) noexcept                                                     \
            requires std::same_as<T0, T1>                                                       \
            : ::boyle::math::Duplet<T0>{cv} {}                                                  \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(T0 c##VAR_0, T1 c##VAR_1) noexcept                                  \
            : ::boyle::math::Duplet<T0, T1>{std::move(c##VAR_0), std::move(c##VAR_1)} {}        \
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
        template <std::convertible_to<T0> U0, std::convertible_to<T1> U1>                       \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(ClassName<U0, U1> other) noexcept                                   \
            : ::boyle::math::Duplet<T0, T1>{std::move(other.first), std::move(other.second)} {} \
                                                                                                \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName(std::pair<T0, T1> pair) noexcept                                    \
            : ::boyle::math::Duplet<T0, T1>{std::move(pair)} {}                                 \
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
    }

namespace boyle::math {

template <typename T0, typename T1 = T0>
class Duplet {
    template <typename CharT, typename U0, typename U1>
    friend inline auto operator<<(std::basic_ostream<CharT>& os, const Duplet<U0, U1>& obj) noexcept
        -> std::basic_ostream<CharT>&;
    friend class boost::serialization::access;
    template <unsigned int, typename T>
    friend struct DupletElement;

  public:
    [[using gnu: always_inline]]
    Duplet() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Duplet(T0 cv) noexcept
        requires std::same_as<T0, T1>
        : m_first{cv}, m_second{cv} {}
    [[using gnu: always_inline]]
    constexpr Duplet(T0 cv0, T1 cv1) noexcept
        : m_first{std::move(cv0)}, m_second{std::move(cv1)} {}
    [[using gnu: always_inline]]
    constexpr Duplet(const Duplet& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Duplet(Duplet&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Duplet& other) noexcept -> Duplet& = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Duplet&& other) noexcept -> Duplet& = default;
    virtual ~Duplet() noexcept = 0;

    [[using gnu: always_inline]]
    constexpr Duplet(std::pair<T0, T1> pair) noexcept
        : m_first{std::move(pair.first)}, m_second{std::move(pair.second)} {}

    [[using gnu: pure, always_inline]]
    constexpr operator std::pair<T0, T1>() const noexcept {
        return std::make_pair<T0, T1>(m_first, m_second);
    }

  protected:
    [[using gnu: always_inline]]
    constexpr auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept
        -> void {
        archive & m_first;
        archive & m_second;
        return;
    }

    T0 m_first{};
    T1 m_second{};
};

template <typename T0, typename T1>
inline Duplet<T0, T1>::~Duplet() noexcept = default;

template <typename CharT, typename T0, typename T1>
[[using gnu: always_inline]]
inline auto operator<<(
    std::basic_ostream<CharT>& os, const boyle::math::Duplet<T0, T1>& obj
) noexcept -> std::basic_ostream<CharT>& {
    os << "(" << obj.m_first << ", " << obj.m_second << ")";
    return os;
}

template <typename T>
concept Dupletic = derivedFromTemplate<T, Duplet>;

template <unsigned int Num, typename T>
struct DupletElement;

template <Dupletic T>
struct DupletElement<0, T> final {
    using type = decltype(std::declval<T>().m_first);
};

template <Dupletic T>
struct DupletElement<1, T> final {
    using type = decltype(std::declval<T>().m_second);
};

template <unsigned int Num, Dupletic T>
using DupletElementT = typename DupletElement<Num, T>::type;

} // namespace boyle::math
