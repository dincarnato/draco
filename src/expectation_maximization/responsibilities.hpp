#pragma once

#include <compare>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <type_traits>

namespace expectation_maximization::responsibilities {

struct Responsibilities;

template <typename R> struct Rows {
  static_assert(std::is_same_v<std::remove_const_t<R>, Responsibilities>);
  using data_type =
      std::conditional_t<std::is_const_v<R>, const double, double>;
  using value_type = std::span<data_type>;
  using difference_type = std::ptrdiff_t;

  constexpr Rows() noexcept = default;
  constexpr Rows(R &responsibilities) noexcept
      : data_(responsibilities.raw_data_pointer()),
        n_clusters_(responsibilities.n_clusters()) {};
  constexpr Rows(R &responsibilities, std::uint32_t row_offset) noexcept
      : data_(responsibilities.raw_data_pointer() +
              static_cast<std::size_t>(row_offset) *
                  static_cast<std::size_t>(responsibilities.n_clusters())),
        n_clusters_(responsibilities.n_clusters()) {};

  constexpr std::span<data_type> operator*() const {
    return std::span(data_, n_clusters_);
  }

  constexpr Rows &operator++() noexcept {
    data_ += n_clusters_;
    return *this;
  }

  constexpr Rows operator++(int) noexcept {
    auto other = *this;
    data_ += n_clusters_;
    return other;
  }

  constexpr Rows &operator+=(difference_type offset) noexcept {
    data_ += n_clusters_ * offset;
    return *this;
  }

  constexpr Rows operator+(difference_type offset) const noexcept {
    auto other = *this;
    other.data_ += n_clusters_ * offset;
    return other;
  }

  friend constexpr Rows operator+(difference_type offset,
                                  Rows const &it) noexcept {
    return it + offset;
  }

  constexpr Rows &operator--() noexcept {
    data_ -= n_clusters_;
    return *this;
  }

  constexpr Rows operator--(int) noexcept {
    auto other = *this;
    data_ -= n_clusters_;
    return other;
  }

  constexpr Rows &operator-=(difference_type offset) noexcept {
    data_ -= n_clusters_ * offset;
    return *this;
  }

  constexpr Rows operator-(difference_type offset) const noexcept {
    auto other = *this;
    other.data_ -= n_clusters_ * offset;
    return other;
  }

  constexpr difference_type operator-(Rows const &other) const noexcept {
    return (data_ - other.data_) / n_clusters_;
  }

  constexpr std::span<data_type>
  operator[](difference_type offset) const noexcept {
    return std::span(data_, offset * n_clusters_, n_clusters_);
  }

  constexpr std::strong_ordering operator<=>(const Rows &other) const noexcept {
    return data_ <=> other.data_;
  }

  constexpr bool operator==(const Rows &other) const noexcept = default;

protected:
  data_type *data_;
  std::uint8_t n_clusters_;
};

static_assert(std::random_access_iterator<Rows<Responsibilities>>);
static_assert(std::random_access_iterator<Rows<const Responsibilities>>);

template <typename R> struct RowsRange {
  static_assert(std::is_same_v<std::remove_const_t<R>, Responsibilities>);
  constexpr RowsRange(R &responsibilities) noexcept
      : responsibilities_(&responsibilities) {}

  constexpr Rows<R> begin() const noexcept {
    return Rows<R>(*responsibilities_);
  }

  constexpr Rows<R> end() const noexcept {
    return Rows<R>(*responsibilities_, responsibilities_->n_rows());
  }

protected:
  R *responsibilities_;
};

template <typename R> struct Cluster {
  static_assert(std::is_same_v<std::remove_const_t<R>, Responsibilities>);
  using data_type =
      std::conditional_t<std::is_const_v<R>, const double, double>;
  using value_type = double;
  using difference_type = std::ptrdiff_t;

  constexpr Cluster() noexcept = default;
  constexpr Cluster(std::span<data_type> data, std::uint8_t n_clusters,
                    std::uint8_t cluster_index) noexcept
      : data_(data), n_clusters_(n_clusters), cluster_index_(cluster_index) {}

  constexpr auto operator*() const noexcept {
    return data_ | std::views::drop(cluster_index_) |
           std::views::stride(n_clusters_);
  }

  constexpr data_type &operator[](std::uint32_t row_index) const noexcept {
    return data_[row_index * n_clusters_ + cluster_index_];
  }

protected:
  std::span<data_type> data_;
  std::uint8_t n_clusters_;
  std::uint8_t cluster_index_{0};
};

template <typename R> struct Clusters {
  static_assert(std::is_same_v<std::remove_const_t<R>, Responsibilities>);
  using data_type =
      std::conditional_t<std::is_const_v<R>, const double, double>;
  using value_type = Cluster<R>;
  using difference_type = std::int16_t;

  constexpr Clusters() noexcept = default;
  constexpr Clusters(R &responsibilities) noexcept
      : data_(responsibilities.raw_data()),
        n_clusters_(responsibilities.n_clusters()) {};
  constexpr Clusters(R &responsibilities, std::uint8_t cluster_offset) noexcept
      : data_(responsibilities.raw_data()),
        n_clusters_(responsibilities.n_clusters()),
        cluster_index_(cluster_offset) {};

  constexpr value_type operator*() const {
    return Cluster<R>(data_, n_clusters_, cluster_index_);
  }

  constexpr Clusters &operator++() noexcept {
    ++cluster_index_;
    return *this;
  }

  constexpr Clusters operator++(int) noexcept {
    auto other = *this;
    ++cluster_index_;
    return other;
  }

  constexpr Clusters &operator+=(difference_type offset) noexcept {
    cluster_index_ = static_cast<std::uint8_t>(
        static_cast<difference_type>(cluster_index_) + offset);
    return *this;
  }

  constexpr Clusters operator+(difference_type offset) const noexcept {
    auto other = *this;
    other.cluster_index_ = static_cast<std::uint8_t>(
        static_cast<difference_type>(other.cluster_index_) + offset);
    return other;
  }

  friend constexpr Clusters operator+(difference_type offset,
                                      Clusters const &it) noexcept {
    return it + offset;
  }

  constexpr Clusters &operator--() noexcept {
    --cluster_index_;
    return *this;
  }

  constexpr Clusters operator--(int) noexcept {
    auto other = *this;
    --cluster_index_;
    return other;
  }

  constexpr Clusters &operator-=(difference_type offset) noexcept {
    cluster_index_ = static_cast<std::uint8_t>(
        static_cast<difference_type>(cluster_index_) - offset);
    return *this;
  }

  constexpr Clusters operator-(difference_type offset) const noexcept {
    auto other = *this;
    other.cluster_index_ = static_cast<std::uint8_t>(
        static_cast<difference_type>(other.cluster_index_) - offset);
    return other;
  }

  constexpr difference_type operator-(Clusters const &other) const noexcept {
    return static_cast<difference_type>(cluster_index_) -
           static_cast<difference_type>(other.cluster_index_);
  }

  constexpr Cluster<R> operator[](difference_type offset) const noexcept {
    return Cluster(data_, n_clusters_, cluster_index_ + offset);
  }

  constexpr std::strong_ordering
  operator<=>(const Clusters &other) const noexcept {
    return cluster_index_ <=> other.cluster_index_;
  }

  constexpr bool operator==(const Clusters &other) const noexcept {
    return cluster_index_ == other.cluster_index_;
  }

protected:
  std::span<data_type> data_;
  std::uint8_t n_clusters_;
  std::uint8_t cluster_index_{0};
};

static_assert(std::random_access_iterator<Clusters<Responsibilities>>);
static_assert(std::random_access_iterator<Clusters<const Responsibilities>>);

template <typename R> struct ClustersRange {
  static_assert(std::is_same_v<std::remove_const_t<R>, Responsibilities>);
  constexpr ClustersRange(R &responsibilities) noexcept
      : responsibilities_(&responsibilities) {}

  constexpr Clusters<R> begin() const noexcept {
    return Clusters<R>(*responsibilities_);
  }

  constexpr Clusters<R> end() const noexcept {
    return Clusters<R>(*responsibilities_, responsibilities_->n_clusters());
  }

protected:
  R *responsibilities_;
};

struct Responsibilities {
  constexpr Responsibilities(std::uint32_t n_rows, std::uint8_t n_clusters)
      : data_(new double[static_cast<std::size_t>(n_rows) *
                         static_cast<std::size_t>(n_clusters)]),
        n_rows_(n_rows), n_clusters_(n_clusters) {}

  constexpr auto rows(this auto &self)
      -> RowsRange<std::remove_reference_t<decltype(self)>> {
    return RowsRange(self);
  }

  constexpr auto clusters(this auto &self)
      -> ClustersRange<std::remove_reference_t<decltype(self)>> {
    return ClustersRange(self);
  }

  constexpr auto row_clusters(this auto &&self,
                              std::uint32_t row_index) noexcept {
    return std::span(self.data_.get(),
                     static_cast<std::size_t>(row_index) *
                         static_cast<std::size_t>(self.n_clusters_),
                     static_cast<std::size_t>(self.n_clusters_));
  }

  constexpr std::uint32_t n_rows() const noexcept { return n_rows_; }
  constexpr std::uint8_t n_clusters() const noexcept { return n_clusters_; }

  constexpr auto raw_data_pointer(this auto &self) noexcept {
    return self.data_.get();
  }
  constexpr auto raw_data(this auto &self) noexcept {
    return std::span(self.data_.get(),
                     static_cast<std::size_t>(self.n_rows_) *
                         static_cast<std::size_t>(self.n_clusters_));
  }

protected:
  std::unique_ptr<double[]> data_;
  std::uint32_t n_rows_;
  std::uint8_t n_clusters_;
};

static_assert(
    std::is_same_v<decltype(std::declval<Responsibilities &>().rows().begin()),
                   Rows<Responsibilities>>);
static_assert(std::is_same_v<
              decltype(std::declval<Responsibilities const &>().rows().begin()),
              Rows<const Responsibilities>>);

} // namespace expectation_maximization::responsibilities

namespace expectation_maximization {
using Responsibilities = responsibilities::Responsibilities;
} // namespace expectation_maximization
