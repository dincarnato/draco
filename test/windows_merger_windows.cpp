#include "windows_merger_windows.hpp"

#include "windows_merger_window_base_impl.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>

#include <range/v3/all.hpp>

using namespace windows_merger;

static void
test_default_construction() {
  WindowsMergerWindows wmw;
  assert(wmw.clusters_size() == 0);
  assert(wmw.bases_capacity() == 0);
  assert(wmw.windows_size() == 0);
  assert(wmw.windows_capacity() == 0);
}

static void
test_construction_with_clusters() {
  WindowsMergerWindows wmw(5);
  assert(wmw.clusters_size() == 5);
  assert(wmw.bases_capacity() == 0);
  assert(wmw.windows_size() == 0);
  assert(wmw.windows_capacity() == 0);
}

static void
test_construction_with_clusters_and_bases() {
  WindowsMergerWindows wmw(5, 70);
  assert(wmw.clusters_size() == 5);
  assert(wmw.bases_capacity() == 70);
  assert(wmw.windows_size() == 0);
  assert(wmw.windows_capacity() == 0);
}

static void
test_construction_with_clusters_bases_and_windows() {
  WindowsMergerWindows wmw(5, 70, 100);
  assert(wmw.clusters_size() == 5);
  assert(wmw.bases_capacity() == 70);
  assert(wmw.windows_size() == 100);
  assert(wmw.windows_capacity() == 100);
}

static void
test_simple_const_accessor() {
  const WindowsMergerWindows wmw(5, 70, 20);
  for (std::size_t index = 0; index < 20; ++index) {
    auto accessor = wmw[index];
    static_assert(std::is_same_v<
                  decltype(accessor),
                  WindowsMergerWindowAccessor<const WindowsMergerWindows>>);

    assert(std::as_const(wmw)[index].index() == index);
    assert(accessor.index() == index);
    assert(accessor.begin_index() == 0);
    assert(accessor.end_index() == 0);
    assert(accessor.size() == 0);
    assert(accessor.coverages().size() == 0);
    assert(accessor.clusters_size() == 5);
    (void)accessor;
  }
}

static std::vector<WindowsMergerWindowBase>
create_random_bases(std::size_t n_clusters, std::size_t n_bases) {
  std::mt19937 random_gen(std::random_device{}());
  std::uniform_real_distribution<double> real_dist(0., 1.);
  std::uniform_int_distribution<unsigned> unsigned_dist(100, 300);

  std::vector<double> weights(n_clusters);
  std::vector<WindowsMergerWindowBase> out_bases(
      n_bases,
      WindowsMergerWindowBase(
          static_cast<WindowsMergerTraits::clusters_size_type>(n_clusters)));
  for (auto& window_base : out_bases) {
    std::generate(ranges::begin(weights), ranges::end(weights),
                  [&] { return real_dist(random_gen); });
    std::transform(ranges::begin(weights), ranges::end(weights),
                   ranges::begin(weights),
                   [normalizer = std::accumulate(ranges::begin(weights),
                                                 ranges::end(weights), 0.)](
                       double weight) { return weight / normalizer; });

    using clusters_size_type =
        typename WindowsMergerWindowBase::clusters_size_type;
    for (clusters_size_type cluster_index = 0; cluster_index < n_clusters;
         ++cluster_index)
      window_base.weight(cluster_index) = weights[cluster_index];
    window_base.coverage() = unsigned_dist(random_gen);
  }

  return out_bases;
}

static std::vector<typename WindowsMergerTraits::bases_size_type>
create_random_starts(WindowsMergerTraits::windows_size_type n_windows,
                     WindowsMergerTraits::bases_size_type n_bases,
                     WindowsMergerTraits::bases_size_type sequence_size) {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  thread_local std::mt19937 random_gen(std::random_device{}());
  std::uniform_int_distribution<bases_size_type> start_dist(0, sequence_size -
                                                                   n_bases);
  std::vector<bases_size_type> all_starts(n_windows);
  std::generate(ranges::begin(all_starts), ranges::end(all_starts),
                [&] { return start_dist(random_gen); });

  return all_starts;
}

static void
test_simple_accessor_emplace_back() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);
  {
    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add) {
      auto new_base_accessor = accessor.emplace_back(window_base);
      assert(new_base_accessor == window_base);
      assert(new_base_accessor == accessor.back());
      assert(std::size(new_base_accessor.weights()) == n_clusters);
    }
  }

  {
    auto accessor = wmw[1];
    auto bases_copied = bases_to_add;
    std::size_t window_base_index = 0;
    for (auto& window_base : bases_copied) {
      auto new_base_accessor = accessor.emplace_back(std::move(window_base));
      assert(new_base_accessor == bases_to_add[window_base_index]);
      assert(new_base_accessor == accessor.back());
      assert(std::size(new_base_accessor.weights()) == n_clusters);
      ++window_base_index;
    }
  }
}

static void
test_simple_accessor_push_back() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);
  {
    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add) {
      accessor.push_back(window_base);
      assert(accessor.back() == window_base);
      assert(std::size(accessor.back().weights()) == n_clusters);
    }
  }

  {
    auto accessor = wmw[1];
    auto bases_copied = bases_to_add;
    std::size_t window_base_index = 0;
    for (auto& window_base : bases_copied) {
      accessor.push_back(std::move(window_base));
      assert(accessor.back() == bases_to_add[window_base_index]);
      assert(std::size(accessor.back().weights()) == n_clusters);
      ++window_base_index;
    }
  }
}

static void
test_front_back_accessors() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);

  {
    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add)
      accessor.push_back(window_base);
  }

  {
    auto accessor = wmw[19];
    for (const auto& window_base : bases_to_add)
      accessor.push_back(window_base);
  }

  assert(&wmw.front()[0].weight(0) == &wmw[0][0].weight(0));
  assert(&wmw.back()[0].weight(0) == &wmw[19][0].weight(0));
  assert(&std::as_const(wmw).front()[0].weight(0) ==
         &std::as_const(wmw)[0][0].weight(0));
  assert(&std::as_const(wmw).back()[0].weight(0) ==
         &std::as_const(wmw)[19][0].weight(0));
}

static void
test_base_accessor_const_iterator() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);

  auto accessor = wmw[0];
  for (const auto& window_base : bases_to_add)
    accessor.emplace_back(window_base);

  for (typename WindowsMergerTraits::bases_size_type window_base_index = 0;
       window_base_index < n_bases_to_add; ++window_base_index) {
    const auto& window_base = accessor[window_base_index];
    auto weights_accessor = window_base.weights();
    assert(std::size(weights_accessor) == n_clusters);
    auto&& base_added = bases_to_add[window_base_index];

    {
      typename WindowsMergerTraits::clusters_size_type cluster_index = 0;
      std::size_t parsed_weights = 0;
      for (auto&& weight : weights_accessor) {
        assert(weight == base_added.weight(cluster_index++));
        ++parsed_weights;
      }

      assert(parsed_weights == n_clusters);
    }

    {
      typename WindowsMergerTraits::clusters_size_type cluster_index = 0;
      std::size_t parsed_weights = 0;
      auto weights_accessor = std::as_const(base_added).weights();
      for (auto&& weight : weights_accessor) {
        assert(weight == base_added.weight(cluster_index++));
        ++parsed_weights;
      }

      assert(parsed_weights == n_clusters);
    }
  }
}

static void
test_base_accessor_iterator() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);

  auto accessor = wmw[0];
  for (const auto& window_base : bases_to_add)
    accessor.emplace_back(window_base);

  for (typename WindowsMergerTraits::bases_size_type window_base_index = 0;
       window_base_index < n_bases_to_add; ++window_base_index) {
    const auto& window_base = accessor[window_base_index];
    auto weights_accessor = window_base.weights();

    ranges::fill(ranges::begin(weights_accessor), ranges::end(weights_accessor),
                 TinyFraction(0.2));
    assert(ranges::all_of(weights_accessor, [](auto&& weight) {
      return weight == TinyFraction(0.2);
    }));

    {
      auto base_added = std::as_const(bases_to_add)[window_base_index];
      auto weights_accessor = base_added.weights();
      ranges::fill(ranges::begin(weights_accessor),
                   ranges::end(weights_accessor), TinyFraction(0.3));
      assert(ranges::all_of(weights_accessor, [](auto&& weight) {
        return weight == TinyFraction(0.3);
      }));
    }
  }
}

static void
test_base_accessor_weights_accessor() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);

  auto accessor = wmw[0];
  for (const auto& window_base : bases_to_add)
    accessor.emplace_back(window_base);

  for (typename WindowsMergerTraits::bases_size_type window_base_index = 0;
       window_base_index < n_bases_to_add; ++window_base_index) {
    const auto& window_base = accessor[window_base_index];
    auto weights_accessor = window_base.weights();

    for (WindowsMergerTraits::clusters_size_type cluster_index = 0;
         cluster_index < n_clusters; ++cluster_index) {
      assert(weights_accessor[cluster_index] ==
             window_base.weight(cluster_index));

      constexpr TinyFraction new_weight(0.2);
      weights_accessor[cluster_index] = new_weight;
      assert(weights_accessor[cluster_index] == new_weight);
      assert(window_base.weight(cluster_index) == new_weight);
    }
  }
}

static void
test_base_comparisons() {
  WindowsMergerWindowBase base(3);
  auto&& weights = base.weights();
  weights[0] = TinyFraction(0.5);
  weights[1] = TinyFraction(0.3);
  weights[2] = TinyFraction(0.2);
  base.coverage() = 125;

  {
    WindowsMergerWindowBase test_base(2);
    test_base.coverage() = 125;

    assert(not(test_base == base));
    assert(test_base != base);
  }

  {
    WindowsMergerWindowBase test_base = base;
    assert(test_base == base);
    assert(not(test_base != base));

    test_base.weight(2) = TinyFraction(0.3);
    assert(not(test_base == base));
    assert(test_base != base);
  }

  {
    WindowsMergerWindowBase test_base = base;
    test_base.coverage() = 128;
    assert(not(test_base == base));
    assert(test_base != base);
  }

  WindowsMergerWindows wmw(3, 70, 20);
  {
    auto accessor = wmw[0];
    accessor.push_back(base);
    assert(accessor[0] == base);
    assert(!(accessor[0] != base));
    assert(accessor[0] == accessor[0]);
    assert(!(accessor[0] != accessor[0]));

    {
      auto test_base = base;
      test_base.coverage() = 128;
      accessor.push_back(test_base);

      assert(!(accessor[1] == base));
      assert(accessor[1] != base);
      assert(accessor[1] == test_base);
      assert(!(accessor[1] != test_base));
      assert(accessor[1] == accessor[1]);
      assert(!(accessor[1] != accessor[1]));
    }

    {
      auto test_base = base;
      test_base.weight(2) = TinyFraction(0.3);
      accessor.push_back(test_base);

      assert(!(accessor[2] == base));
      assert(accessor[2] != base);
      assert(accessor[2] == test_base);
      assert(!(accessor[2] != test_base));
      assert(accessor[2] == accessor[2]);
      assert(!(accessor[2] != accessor[2]));
    }

    assert(accessor[0] != accessor[1]);
    assert(accessor[0] != accessor[2]);
    assert(accessor[1] != accessor[2]);
  }

  {
    auto accessor = std::as_const(wmw)[0];
    assert(accessor[0] == base);
    assert(!(accessor[0] != base));
    assert(accessor[0] == accessor[0]);
    assert(!(accessor[0] != accessor[0]));

    {
      auto test_base = base;
      test_base.coverage() = 128;

      assert(!(accessor[1] == base));
      assert(accessor[1] != base);
      assert(accessor[1] == test_base);
      assert(!(accessor[1] != test_base));
      assert(accessor[1] == accessor[1]);
      assert(!(accessor[1] != accessor[1]));
    }

    {
      auto test_base = base;
      test_base.weight(2) = TinyFraction(0.3);

      assert(!(accessor[2] == base));
      assert(accessor[2] != base);
      assert(accessor[2] == test_base);
      assert(!(accessor[2] != test_base));
      assert(accessor[2] == accessor[2]);
      assert(!(accessor[2] != accessor[2]));
    }

    assert(accessor[0] != accessor[1]);
    assert(accessor[0] != accessor[2]);
    assert(accessor[1] != accessor[2]);
  }

  WindowsMergerWindows new_wmw(5, 70, 20);

  WindowsMergerWindowBase test_base(5);
  {
    auto&& weights = test_base.weights();
    weights[0] = TinyFraction(0.5);
    weights[1] = TinyFraction(0.3);
    weights[2] = TinyFraction(0.2);
    weights[3] = TinyFraction(0.);
    weights[4] = TinyFraction(0.);
    test_base.coverage() = 125;
  }

  {
    auto&& accessor = new_wmw[0];
    accessor.push_back(test_base);
    assert(accessor[0] == test_base);
    assert(!(accessor[0] != test_base));
    assert(!(accessor[0] == base));
    assert(accessor[0] != base);
    assert(accessor[0] == accessor[0]);
    assert(!(accessor[0] != accessor[0]));

    {
      auto&& wmw_accessor = wmw[0];
      assert(accessor[0] != wmw_accessor[0]);
      assert(accessor[0] != wmw_accessor[1]);
      assert(accessor[0] != wmw_accessor[2]);
    }
    {
      auto&& wmw_accessor = std::as_const(wmw)[0];
      assert(accessor[0] != wmw_accessor[0]);
      assert(accessor[0] != wmw_accessor[1]);
      assert(accessor[0] != wmw_accessor[2]);
    }
  }

  {
    auto&& accessor = std::as_const(new_wmw)[0];
    assert(accessor[0] == test_base);
    assert(!(accessor[0] != test_base));
    assert(!(accessor[0] == base));
    assert(accessor[0] != base);
    assert(accessor[0] == accessor[0]);
    assert(!(accessor[0] != accessor[0]));

    {
      auto&& wmw_accessor = wmw[0];
      assert(accessor[0] != wmw_accessor[0]);
      assert(accessor[0] != wmw_accessor[1]);
      assert(accessor[0] != wmw_accessor[2]);
    }
    {
      auto&& wmw_accessor = std::as_const(wmw)[0];
      assert(accessor[0] != wmw_accessor[0]);
      assert(accessor[0] != wmw_accessor[1]);
      assert(accessor[0] != wmw_accessor[2]);
    }
  }

  {
    WindowsMergerWindows new_wmw(wmw.clusters_size(), wmw.bases_capacity(),
                                 wmw.windows_capacity());
    {
      auto&& accessor = std::as_const(wmw)[0];
      auto&& new_accessor = new_wmw[0];
      for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 3;
           ++base_index) {
        WindowsMergerWindowBase base(wmw.clusters_size());
        base.coverage() = accessor[base_index].coverage();
        for (WindowsMergerTraits::clusters_size_type cluster_index = 0;
             cluster_index < wmw.clusters_size(); ++cluster_index)
          base.weight(cluster_index) =
              accessor[base_index].weight(cluster_index);

        new_accessor.push_back(std::move(base));
      }
    }

    {
      auto&& accessor_a = wmw[0];
      auto&& accessor_b = new_wmw[0];

      assert(accessor_a[0] == accessor_b[0]);
      assert(accessor_a[0] != accessor_b[1]);
      assert(accessor_a[0] != accessor_b[2]);
      assert(accessor_a[1] == accessor_b[1]);
      assert(accessor_a[1] != accessor_b[2]);
      assert(accessor_a[2] == accessor_b[2]);
    }

    {
      auto&& accessor_a = std::as_const(wmw)[0];
      auto&& accessor_b = new_wmw[0];

      assert(accessor_a[0] == accessor_b[0]);
      assert(accessor_a[0] != accessor_b[1]);
      assert(accessor_a[0] != accessor_b[2]);
      assert(accessor_a[1] == accessor_b[1]);
      assert(accessor_a[1] != accessor_b[2]);
      assert(accessor_a[2] == accessor_b[2]);
    }

    {
      auto&& accessor_a = wmw[0];
      auto&& accessor_b = std::as_const(new_wmw)[0];

      assert(accessor_a[0] == accessor_b[0]);
      assert(accessor_a[0] != accessor_b[1]);
      assert(accessor_a[0] != accessor_b[2]);
      assert(accessor_a[1] == accessor_b[1]);
      assert(accessor_a[1] != accessor_b[2]);
      assert(accessor_a[2] == accessor_b[2]);
    }

    {
      auto&& accessor_a = std::as_const(wmw)[0];
      auto&& accessor_b = std::as_const(new_wmw)[0];

      assert(accessor_a[0] == accessor_b[0]);
      assert(accessor_a[0] != accessor_b[1]);
      assert(accessor_a[0] != accessor_b[2]);
      assert(accessor_a[1] == accessor_b[1]);
      assert(accessor_a[1] != accessor_b[2]);
      assert(accessor_a[2] == accessor_b[2]);
    }
  }
}

static void
test_accessor_resizing_push() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 1, 20);
  auto accessor = wmw[0];
  for (const auto& window_base : bases_to_add) {
    auto new_base_accessor = accessor.emplace_back(window_base);
    assert(new_base_accessor == window_base);
    assert(new_base_accessor == accessor.back());
  }

  assert(wmw.bases_capacity() >= std::size(bases_to_add));
  for (WindowsMergerTraits::bases_size_type base_index = 0;
       base_index < std::size(bases_to_add); ++base_index)
    assert(accessor[base_index] == bases_to_add[base_index]);
}

static void
test_accessor_resizing_push_on_empty() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  {
    WindowsMergerWindows wmw(5, 0, 20);
    {
      auto accessor = wmw[0];
      assert(accessor.index() == 0);
      assert(accessor.begin_index() == 0);
      assert(accessor.end_index() == 0);
      assert(accessor.size() == 0);
      assert(accessor.coverages().size() == 0);
    }

    {
      auto accessor = std::as_const(wmw)[0];
      assert(accessor.index() == 0);
      assert(accessor.begin_index() == 0);
      assert(accessor.end_index() == 0);
      assert(accessor.size() == 0);
      assert(accessor.coverages().size() == 0);
    }

    {
      auto accessor = std::move(wmw)[0];
      assert(accessor.index() == 0);
      assert(accessor.begin_index() == 0);
      assert(accessor.end_index() == 0);
      assert(accessor.size() == 0);
      assert(accessor.coverages().size() == 0);
    }
  }

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  {
    WindowsMergerWindows wmw(5, 0, 20);
    auto accessor = wmw[0];
    assert(accessor.index() == 0);
    assert(accessor.begin_index() == 0);
    assert(accessor.end_index() == 0);
    assert(accessor.size() == 0);
    assert(accessor.coverages().size() == 0);

    for (const auto& window_base : bases_to_add) {
      auto new_base_accessor = accessor.emplace_back(window_base);
      assert(new_base_accessor == window_base);
      assert(new_base_accessor == accessor.back());
    }

    assert(wmw.bases_capacity() >= std::size(bases_to_add));
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < std::size(bases_to_add); ++base_index)
      assert(accessor[base_index] == bases_to_add[base_index]);
  }

  {
    WindowsMergerWindows wmw(5, 0, 20);
    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add) {
      accessor.push_back(window_base);
      assert(accessor.back() == window_base);
    }

    assert(wmw.bases_capacity() >= std::size(bases_to_add));
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < std::size(bases_to_add); ++base_index)
      assert(accessor[base_index] == bases_to_add[base_index]);
  }

  {
    WindowsMergerWindows wmw(5, 0, 20);
    auto accessor = wmw[0];
    for (auto& window_base : bases_to_add) {
      accessor.push_back(window_base);
      assert(accessor.back() == window_base);
    }

    assert(wmw.bases_capacity() >= std::size(bases_to_add));
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < std::size(bases_to_add); ++base_index)
      assert(accessor[base_index] == bases_to_add[base_index]);
  }
}

static void
test_windows_reshape_noreshape() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  {
    WindowsMergerWindows wmw(5);
    wmw.reshape(0, typename WindowsMergerWindows::reserver_type(0));
  }

  {
    WindowsMergerWindows wmw(5);
    wmw.reshape(0, typename WindowsMergerWindows::resizer_type(0));
  }

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, 20);

  {
    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto first_element_address =
        std::addressof(static_cast<TinyFraction&>(wmw[0].front().weight(0)));
    wmw.reshape(0, typename WindowsMergerWindows::reserver_type(0));
    assert(first_element_address == std::addressof(static_cast<TinyFraction&>(
                                        wmw[0].front().weight(0))));
  }

  {
    auto first_element_address =
        std::addressof(static_cast<TinyFraction&>(wmw[0].front().weight(0)));
    wmw.reshape(70, typename WindowsMergerWindows::reserver_type(20));
    assert(first_element_address == std::addressof(static_cast<TinyFraction&>(
                                        wmw[0].front().weight(0))));
  }

  {
    auto first_element_address =
        std::addressof(static_cast<TinyFraction&>(wmw[0].front().weight(0)));
    wmw.reshape(70, typename WindowsMergerWindows::resizer_type(20));
    assert(first_element_address == std::addressof(static_cast<TinyFraction&>(
                                        wmw[0].front().weight(0))));
  }
}

static void
test_windows_reshape_from_empty() {
  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 20;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  auto add_to_wmw = [&](WindowsMergerWindows& wmw) {
    if (wmw.windows_size() == 0)
      wmw.reshape(0, WindowsMergerWindows::resizer_type(1));

    auto accessor = wmw[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);

    assert(wmw.bases_capacity() == 70);
    assert(wmw.windows_capacity() == 20);
  };

  auto check_wmw = [&](const WindowsMergerWindows& wmw) {
    auto accessor = wmw[0];
    assert(std::size(accessor) == n_bases_to_add);
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < std::size(accessor); ++base_index) {
      assert(accessor[base_index] == bases_to_add[base_index]);
    }

    assert(wmw.bases_capacity() == 70);
    assert(wmw.windows_capacity() == 30);
  };

  {
    WindowsMergerWindows wmw(5, 0, 20);
    assert(wmw.windows_size() == 20);
    wmw.reshape(70, typename WindowsMergerWindows::reserver_type(20));
    add_to_wmw(wmw);
    assert(wmw.windows_size() == 20);

    wmw.reshape(70, typename WindowsMergerWindows::reserver_type(30));
    check_wmw(wmw);
    assert(wmw.windows_capacity() == 30);
    assert(wmw.windows_size() == 20);
  }

  {
    WindowsMergerWindows wmw(5, 0, 20);
    assert(wmw.windows_size() == 20);
    wmw.reshape(70, typename WindowsMergerWindows::resizer_type(20));
    add_to_wmw(wmw);
    assert(wmw.windows_size() == 20);

    wmw.reshape(70, typename WindowsMergerWindows::resizer_type(30));
    check_wmw(wmw);
    assert(wmw.windows_size() == 30);
  }

  {
    WindowsMergerWindows wmw(5, 70, 0);
    assert(wmw.windows_size() == 0);
    wmw.reshape(70, typename WindowsMergerWindows::reserver_type(20));
    add_to_wmw(wmw);
    assert(wmw.windows_size() == 1);

    wmw.reshape(70, typename WindowsMergerWindows::reserver_type(30));
    check_wmw(wmw);
    assert(wmw.windows_size() == 1);
  }

  {
    WindowsMergerWindows wmw(5, 70, 0);
    assert(wmw.windows_size() == 0);
    wmw.reshape(70, typename WindowsMergerWindows::resizer_type(20));
    add_to_wmw(wmw);
    assert(wmw.windows_size() == 20);

    wmw.reshape(70, typename WindowsMergerWindows::resizer_type(30));
    check_wmw(wmw);
    assert(wmw.windows_size() == 30);
  }
}

static void
test_windows_reshape_shrink() {
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
  WindowsMergerWindows wmw(5, 70, n_windows);

  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < n_windows; ++window_index) {
    auto accessor = wmw[window_index];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  wmw.reshape(50, typename WindowsMergerWindows::resizer_type(10));
  assert(wmw.bases_capacity() == 70);
  assert(wmw.windows_capacity() == 20);
  assert(wmw.windows_size() == 10);
  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < 10; ++window_index) {
    auto accessor = wmw[window_index];
    const auto accessor_size = std::size(accessor);
    assert(accessor_size == n_bases_to_add);
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < accessor_size; ++base_index)
      assert(accessor[base_index] == bases_to_add[base_index]);
  }
}

static void
test_windows_emplace_back() {
  constexpr WindowsMergerTraits::clusters_size_type n_clusters = 5;
  {
    WindowsMergerWindows wmw(n_clusters);
    for (std::size_t index = 0; index < 10; ++index) {
      auto new_window = wmw.emplace_back();
      assert(wmw.windows_size() == index + 1);
      assert(wmw.windows_capacity() >= index + 1);

      auto new_bases = create_random_bases(n_clusters, 30);
      for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
           ++base_index)
        new_window.emplace_back(new_bases[base_index]);

      assert(std::size(new_window) == 30);
      assert(wmw.bases_capacity() >= 30);
      for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
           ++base_index)
        assert(new_window[base_index] == new_bases[base_index]);
    }
  }

  {
    WindowsMergerWindows wmw(n_clusters);
    for (std::size_t index = 0; index < 10; ++index) {
      wmw.emplace_back();
      assert(wmw.windows_size() == index + 1);
      assert(wmw.windows_capacity() >= index + 1);
    }

    for (std::size_t index = 0; index < 10; ++index) {
      auto new_window = wmw[index];
      auto new_bases = create_random_bases(n_clusters, 30);
      for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
           ++base_index)
        new_window.emplace_back(new_bases[base_index]);

      assert(std::size(new_window) == 30);
      assert(wmw.bases_capacity() >= 30);
      for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
           ++base_index)
        assert(new_window[base_index] == new_bases[base_index]);
    }
  }
}

static void
test_window_set_begin_index() {
  constexpr WindowsMergerTraits::clusters_size_type n_clusters = 5;
  WindowsMergerWindows wmw(n_clusters);

  auto new_window = wmw.emplace_back();
  assert(wmw.windows_size() == 1);
  assert(wmw.windows_capacity() >= 1);
  assert(wmw.bases_capacity() == 0);
  assert(new_window.begin_index() == 0);
  assert(new_window.end_index() == 0);

  new_window.set_begin_index(300);
  assert(wmw.windows_size() == 1);
  assert(wmw.windows_capacity() >= 1);
  assert(wmw.bases_capacity() > 0);
  assert(new_window.begin_index() == 300);
  assert(new_window.end_index() == 300);

  auto new_bases = create_random_bases(n_clusters, 30);
  for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
       ++base_index)
    new_window.emplace_back(new_bases[base_index]);

  assert(std::size(new_window) == 30);
  assert(wmw.bases_capacity() >= 30);
  assert(new_window.begin_index() == 300);
  assert(new_window.end_index() == 330);
  for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
       ++base_index)
    assert(new_window[base_index] == new_bases[base_index]);
}

static void
test_window_coverage_accessor() {
  constexpr std::size_t n_clusters = 5;
  const auto new_bases = create_random_bases(n_clusters, 30);

  WindowsMergerWindows wmw(n_clusters);
  auto new_window = wmw.emplace_back();

  {
    auto coverages = new_window.coverages();
    assert(coverages.index() == 0);
    assert(coverages.size() == 0);
    assert(ranges::begin(coverages) == nullptr);
    assert(ranges::end(coverages) == nullptr);
  }

  {
    auto coverages = std::as_const(wmw)[0].coverages();
    assert(coverages.size() == 0);
    assert(ranges::begin(coverages) == nullptr);
    assert(ranges::end(coverages) == nullptr);
  }

  for (WindowsMergerTraits::bases_size_type base_index = 0; base_index < 30;
       ++base_index)
    new_window.emplace_back(new_bases[base_index]);

  auto perform_test = [&new_bases](auto&& new_window) {
    auto new_window_coverages = new_window.coverages();
    assert(ranges::equal(
        new_bases, new_window_coverages, {},
        [](auto&& base_accessor) { return base_accessor.coverage(); }));
    assert(ranges::equal(
        std::as_const(new_bases), new_window_coverages, {},
        [](auto&& base_accessor) { return base_accessor.coverage(); }));
    for (typename WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < 30; ++base_index) {
      assert(new_window_coverages[base_index] ==
             new_bases[base_index].coverage());
    }
    assert(new_window_coverages.front() == new_bases.front().coverage());
    assert(new_window_coverages.back() == new_bases.back().coverage());
  };

  perform_test(new_window);
  perform_test(std::as_const(wmw).back());
}

static void
test_window_accessor_iterator() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr windows_size_type n_windows = 20;
  constexpr bases_size_type n_bases = 30;
  std::vector<decltype(create_random_bases(n_clusters, n_bases))> raw_windows(
      n_windows);
  std::generate(ranges::begin(raw_windows), ranges::end(raw_windows),
                [] { return create_random_bases(n_clusters, n_bases); });

  WindowsMergerWindows wmw(n_clusters);
  for (auto&& window : std::as_const(raw_windows)) {
    auto&& new_window = wmw.emplace_back();
    // Cannot use std::back_inserter :(
    for (auto&& base : std::as_const(window))
      new_window.push_back(base);
  }

  auto perform_test = [&raw_windows](auto&& wmw) {
    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      auto accessor = wmw[window_index];
      const auto& raw_window = raw_windows[window_index];

      {
        auto accessor_iter = ranges::begin(accessor);
        const auto accessor_end = ranges::end(accessor);
        auto raw_window_iter = ranges::begin(raw_window);
        typename decltype(accessor_iter)::difference_type position = 0;

        {
          auto iter = accessor_iter;
          ++iter;
          ++iter;
          iter++;
          const auto inc_iter = accessor_iter + 3;
          {
            auto new_iter = accessor_iter;
            new_iter += 3;
            assert(new_iter == inc_iter);
          }

          assert(inc_iter == iter);
          assert(inc_iter >= iter);
          assert(inc_iter <= iter);
          assert(3 + accessor_iter == iter);
          assert(3 + accessor_iter >= accessor_iter);
          assert(accessor_iter <= accessor_iter + 3);
          assert(accessor_iter != iter);
          iter--;
          --iter;
          --iter;
          assert(inc_iter - 3 == iter);
          assert(inc_iter != iter);

          {
            auto new_iter = inc_iter;
            new_iter -= 3;
            assert(new_iter == iter);
          }
        }

        for (; accessor_iter < accessor_end;
             ++accessor_iter, ++raw_window_iter, ++position) {
          assert(*accessor_iter == *raw_window_iter);
          assert(ranges::distance(ranges::begin(accessor), accessor_iter) ==
                 position);
          assert(ranges::begin(accessor)[position] == *accessor_iter);
        }

        assert(raw_window_iter == ranges::end(raw_window));
      }

      {
        auto accessor_iter = std::rbegin(accessor);
        const auto accessor_end = std::rend(accessor);
        auto raw_window_iter = std::rbegin(raw_window);
        typename decltype(accessor_iter)::difference_type position = 0;
        for (; accessor_iter < accessor_end;
             ++accessor_iter, ++raw_window_iter, ++position) {
          assert(*accessor_iter == *raw_window_iter);
          assert(ranges::distance(std::rbegin(accessor), accessor_iter) ==
                 position);
          assert(std::rbegin(accessor)[position] == *accessor_iter);
        }

        assert(raw_window_iter == std::rend(raw_window));
      }
    }
  };

  perform_test(wmw);
  perform_test(std::as_const(wmw));
}

static void
test_accessor_assignment() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr bases_size_type n_bases_to_add = 30;

  auto perform_test = [&](bases_size_type initial_capacity, auto get_rhs_wmw,
                          auto get_lhs_wmw, auto assign_accessor) {
    WindowsMergerWindows wmw_lhs(5, initial_capacity, 1);
    {
      auto starts = create_random_starts(1, initial_capacity, 300);
      auto bases_to_add = create_random_bases(n_clusters, initial_capacity);
      auto accessor = wmw_lhs[0];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      accessor.set_begin_index(starts[0]);
    }

    WindowsMergerWindows wmw_rhs(5, n_bases_to_add, 1);
    {
      auto starts = create_random_starts(1, n_bases_to_add, 300);
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw_rhs[0];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      accessor.set_begin_index(starts[0]);
    }

    assert(not ranges::equal(wmw_lhs[0], wmw_rhs[0]));
    {
      auto accessor = get_rhs_wmw(wmw_rhs)[0];
      assign_accessor(get_lhs_wmw(wmw_lhs)[0], accessor);
    }
    // I can do this even in case of move because I know that the source state
    // is unchanged
    assert(wmw_lhs[0].begin_index() == wmw_rhs[0].begin_index());
    assert(ranges::equal(wmw_lhs[0], wmw_rhs[0]));
  };

  for (bases_size_type bases_initial_capacity :
       std::array{bases_size_type(70), bases_size_type(13)}) {

    perform_test(
        bases_initial_capacity, [](auto& wmw) -> decltype(auto) { return wmw; },
        [](auto& wmw) -> decltype(auto) { return wmw; },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return wmw; },
                 [](auto& wmw) -> decltype(auto) { return wmw; },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });

    perform_test(
        bases_initial_capacity,
        [](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
        [](auto& wmw) -> decltype(auto) { return wmw; },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
                 [](auto& wmw) -> decltype(auto) { return wmw; },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });

    perform_test(
        bases_initial_capacity, [](auto& wmw) -> decltype(auto) { return wmw; },
        [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return wmw; },
                 [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });

    perform_test(
        bases_initial_capacity,
        [](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
        [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
                 [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });

    perform_test(
        bases_initial_capacity,
        [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
        [](auto& wmw) -> decltype(auto) { return wmw; },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
                 [](auto& wmw) -> decltype(auto) { return wmw; },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });

    perform_test(
        bases_initial_capacity,
        [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
        [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
        [](auto&& accessor_a, auto&& accessor_b) { accessor_a = accessor_b; });

    perform_test(bases_initial_capacity,
                 [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
                 [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
                 [](auto&& accessor_a, auto&& accessor_b) {
                   accessor_a = std::move(accessor_b);
                 });
  }
}

static void
test_accessor_swap() {
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
      added_windows;
  auto added_windows_starts =
      create_random_starts(n_windows, n_bases_to_add, 300);

  WindowsMergerWindows wmw(5, 70, n_windows);
  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < n_windows; ++window_index) {
    auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
    auto accessor = wmw[window_index];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    accessor.set_begin_index(added_windows_starts[window_index]);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  {
    auto wmw_0 = wmw[0];
    auto wmw_1 = wmw[1];
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];
    assert(
        not ranges::equal(wmw[0], wmw_1 | ranges::view::take(wmw[0].size())));
    assert(ranges::equal(wmw_0,
                         added_windows_0 | ranges::view::take(wmw_0.size())));
    assert(ranges::equal(wmw_1,
                         added_windows_1 | ranges::view::take(wmw_1.size())));
    assert(wmw_0.begin_index() == added_windows_starts[0]);
    assert(wmw_1.begin_index() == added_windows_starts[1]);
  }

  ranges::swap(wmw[0], wmw[1]);

  {
    auto&& wmw_0 = wmw[0];
    auto&& wmw_1 = wmw[1];
    assert(wmw_0.index() == 0);
    assert(wmw_1.index() == 1);
    assert(wmw_0.begin_index() == added_windows_starts[1]);
    assert(wmw_1.begin_index() == added_windows_starts[0]);
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];
    assert(ranges::equal(wmw_0,
                         added_windows_1 | ranges::view::take(wmw_0.size())));
    assert(ranges::equal(wmw_1,
                         added_windows_0 | ranges::view::take(wmw_1.size())));
  }
}

static void
test_accessor_assignment_reshape() {
  constexpr std::size_t n_clusters = 5;

  WindowsMergerWindows wmw_a(5, 30, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  WindowsMergerWindows wmw_b(5, 70, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
  }
  wmw_a[0] = wmw_b[0];
  assert(wmw_a.bases_capacity() == 70);
  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    assert(ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
  }
}

static void
test_accessor_assignment_destroy() {
  constexpr std::size_t n_clusters = 5;

  WindowsMergerWindows wmw_a(5, 70, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  WindowsMergerWindows wmw_b(5, 30, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
    assert(wmw_a0.size() == 70);
    assert(wmw_a0.coverages().size() == 70);

    wmw_a0 = wmw_b0;

    assert(wmw_a0.size() == 30);
    assert(wmw_a0.coverages().size() == 30);
    assert(ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
  }
}

static void
test_accessor_move_assignment_reshape() {
  constexpr std::size_t n_clusters = 5;

  WindowsMergerWindows wmw_a(5, 30, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  WindowsMergerWindows wmw_b(5, 70, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];

    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
    wmw_a0 = std::move(wmw_b)[0];
    assert(wmw_a.bases_capacity() == 70);
    assert(ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
  }
}

static void
test_accessor_move_assignment_destroy() {
  constexpr std::size_t n_clusters = 5;

  WindowsMergerWindows wmw_a(5, 70, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  WindowsMergerWindows wmw_b(5, 30, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];

    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
    assert(wmw_a0.size() == 70);
    assert(wmw_a0.coverages().size() == 70);

    wmw_a0 = std::move(wmw_b)[0];

    assert(wmw_a0.size() == 30);
    assert(wmw_a0.coverages().size() == 30);
    assert(ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
  }
}

static void
test_accessor_swap_reshape_lhs() {
  constexpr std::size_t n_clusters = 5;

  std::vector<std::vector<WindowsMergerWindowBase>> added_windows;

  WindowsMergerWindows wmw_a(5, 30, 1);
  WindowsMergerWindows wmw_b(5, 70, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];

    assert(wmw_a.bases_capacity() == 30);
    assert(wmw_b.bases_capacity() == 70);
    assert(std::size(wmw_a0) == 30);
    assert(std::size(wmw_b0) == 70);
    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_a0,
                         added_windows_0 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_b0,
                         added_windows_1 | ranges::view::take(wmw_b0.size())));
  }

  ranges::swap(wmw_a[0], wmw_b[0]);

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];

    assert(wmw_a.bases_capacity() == 70);
    assert(wmw_b.bases_capacity() == 70);
    assert(std::size(wmw_a0) == 70);
    assert(std::size(wmw_b0) == 30);
    assert(ranges::equal(wmw_a0,
                         added_windows_1 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_b0,
                         added_windows_0 | ranges::view::take(wmw_b0.size())));
  }
}

static void
test_accessor_swap_reshape_rhs() {
  constexpr std::size_t n_clusters = 5;

  std::vector<std::vector<WindowsMergerWindowBase>> added_windows;

  WindowsMergerWindows wmw_a(5, 70, 1);
  WindowsMergerWindows wmw_b(5, 30, 1);
  {
    auto bases_to_add = create_random_bases(n_clusters, 70);
    auto accessor = wmw_a[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  {
    auto bases_to_add = create_random_bases(n_clusters, 30);
    auto accessor = wmw_b[0];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];

    assert(wmw_a.bases_capacity() == 70);
    assert(wmw_b.bases_capacity() == 30);
    assert(std::size(wmw_a0) == 70);
    assert(std::size(wmw_b0) == 30);
    assert(
        not ranges::equal(wmw_a0, wmw_b0 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_a0,
                         added_windows_0 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_b0,
                         added_windows_1 | ranges::view::take(wmw_b0.size())));
  }

  ranges::swap(wmw_a[0], wmw_b[0]);

  {
    auto&& wmw_a0 = wmw_a[0];
    auto&& wmw_b0 = wmw_b[0];
    auto&& added_windows_0 = added_windows[0];
    auto&& added_windows_1 = added_windows[1];

    assert(wmw_a.bases_capacity() == 70);
    assert(wmw_b.bases_capacity() == 70);
    assert(std::size(wmw_a0) == 30);
    assert(std::size(wmw_b0) == 70);
    assert(ranges::equal(wmw_a0,
                         added_windows_1 | ranges::view::take(wmw_a0.size())));
    assert(ranges::equal(wmw_b0,
                         added_windows_0 | ranges::view::take(wmw_b0.size())));
  }
}

static void
test_base_accessor_assignment() {
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;
  auto starts = create_random_starts(n_windows, n_bases_to_add, 300);

  auto perform_test = [&](auto get_rhs_wmw, auto get_lhs_wmw, auto assign_op) {
    WindowsMergerWindows wmw(5, 70, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw[window_index];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      accessor.set_begin_index(starts[window_index]);
    }

    assert(not ranges::equal(wmw[0], wmw[1]));
    {
      auto lhs_accessor = get_lhs_wmw(wmw)[0];
      auto rhs_accessor = get_rhs_wmw(wmw)[1];

      auto rhs_iter = ranges::begin(rhs_accessor);
      const auto rhs_end_iter = ranges::end(rhs_accessor);
      auto lhs_iter = ranges::begin(lhs_accessor);
      for (; rhs_iter < rhs_end_iter; ++rhs_iter, ++lhs_iter)
        assign_op(*lhs_iter, *rhs_iter);
      wmw[0].set_begin_index(wmw[1].begin_index());
    }
    // I can do this even in case of move because I know that the source state
    // is unchanged
    assert(ranges::equal(wmw[0], wmw[1]));
  };

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = rhs; });

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::as_const(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return wmw; },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });

  perform_test([](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& lhs, auto&& rhs) { lhs = std::move(rhs); });
}

static void
test_window_accessor_to_window_constructor() {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [](auto wmw_get) {
    std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
        added_windows;

    WindowsMergerWindows wmw(n_clusters, 70, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw[window_index];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      added_windows.emplace_back(bases_to_add);
    }

    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      WindowsMergerWindow window(wmw_get(wmw)[window_index]);
      assert(ranges::equal(
          added_windows[window_index],
          window | ranges::view::take(added_windows[window_index].size())));
    }
  };

  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); });
}

static void
test_window_accessor_to_window_assignment() {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [](auto wmw_get, auto get_accessor) {
    std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
        added_windows;

    WindowsMergerWindows wmw(n_clusters, 70, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw[window_index];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      added_windows.emplace_back(std::move(bases_to_add));
    }

    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      WindowsMergerWindow window(n_clusters);
      auto accessor = wmw_get(wmw)[window_index];
      window = get_accessor(accessor);
      assert(ranges::equal(
          added_windows[window_index],
          window | ranges::view::take(added_windows[window_index].size())));
    }
  };

  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return wmw; },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
}

static void
test_window_base_accessor_to_base_constructor() {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  using bases_size_type = WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [](auto wmw_get) {
    std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
        added_windows;

    WindowsMergerWindows wmw(n_clusters, 70, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw[window_index];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      added_windows.emplace_back(bases_to_add);
    }

    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      auto window_accessor = wmw_get(wmw)[window_index];
      const auto& added_window = added_windows[window_index];
      for (bases_size_type base_index = 0; base_index < n_bases_to_add;
           ++base_index) {
        WindowsMergerWindowBase base(window_accessor[base_index]);
        assert(base == added_window[base_index]);
      }
    }
  };

  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); });
}

static void
test_window_base_accessor_to_base_assignment() {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  using bases_size_type = WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [](auto wmw_get, auto get_accessor) {
    std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
        added_windows;

    WindowsMergerWindows wmw(n_clusters, 70, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
      auto accessor = wmw[window_index];
      for (const auto& window_base : bases_to_add)
        accessor.emplace_back(window_base);
      added_windows.emplace_back(bases_to_add);
    }

    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      auto window_accessor = wmw_get(wmw)[window_index];
      const auto& added_window = added_windows[window_index];
      for (bases_size_type base_index = 0; base_index < n_bases_to_add;
           ++base_index) {
        WindowsMergerWindowBase base;
        auto accessor = window_accessor[base_index];
        base = get_accessor(accessor);
        assert(base == added_window[base_index]);
      }
    }
  };

  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& accessor) -> decltype(auto) { return accessor; });
  perform_test([](auto&& wmw) -> decltype(auto) { return wmw; },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test([](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
               [](auto&& accessor) -> decltype(auto) {
                 return std::as_const(accessor);
               });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return wmw; },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return std::as_const(wmw); },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
  perform_test(
      [](auto&& wmw) -> decltype(auto) { return std::move(wmw); },
      [](auto&& accessor) -> decltype(auto) { return std::move(accessor); });
}

static void
test_windows_iterators() {
  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_bases_to_add = 30;

  std::vector<decltype(create_random_bases(n_clusters, n_bases_to_add))>
      added_windows;

  WindowsMergerWindows wmw(5, 70, n_windows);
  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < n_windows; ++window_index) {
    auto bases_to_add = create_random_bases(n_clusters, n_bases_to_add);
    auto accessor = wmw[window_index];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  auto perform_test = [&added_windows](auto&& wmw) {
    auto test_with_iterators = [](auto wmw_begin, auto wmw_end,
                                  auto added_windows_begin,
                                  auto added_windows_end) {
      if (ranges::distance(wmw_begin, wmw_end) !=
          ranges::distance(added_windows_begin, added_windows_end))
        return false;

      {
        typename WindowsMergerTraits::bases_size_type index = 0;
        auto wmw_iter = wmw_begin;

        assert(ranges::next(wmw_iter) == ranges::next(wmw_begin));
        assert(ranges::next(wmw_iter) != wmw_begin);
        assert(ranges::next(wmw_iter) > wmw_begin);
        assert(ranges::next(wmw_iter) >= wmw_begin);
        assert(wmw_begin < ranges::next(wmw_iter));
        assert(wmw_begin <= ranges::next(wmw_iter));
        {
          auto iter = wmw_begin;
          ++iter;
          ++iter;
          assert(wmw_begin + 2 == iter);
        }

        {
          auto iter = wmw_begin;
          ++iter;
          ++iter;
          assert(2 + wmw_begin == iter);
        }
        assert(ranges::next(wmw_iter, 2) - 2 == wmw_begin);

        auto added_windows_iter = added_windows_begin;
        for (; wmw_iter < wmw_end; ++wmw_iter, ++added_windows_iter, ++index) {
          auto&& wmw_window = *wmw_iter;
          assert(wmw_begin[index].index() == (*wmw_iter).index() and
                 wmw_begin[index].begin_index() == (*wmw_iter).begin_index() and
                 wmw_begin[index].end_index() == (*wmw_iter).end_index() and
                 ranges::equal(
                     wmw_begin[index],
                     wmw_window | ranges::view::take(wmw_begin[index].size())));

          if (added_windows_iter == added_windows_end)
            return false;

          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_end)
          return false;
      }

      {
        auto wmw_iter = wmw_end;
        auto added_windows_iter = added_windows_end;
        while (wmw_iter > wmw_begin) {
          --wmw_iter;
          --added_windows_iter;

          if (added_windows_iter == added_windows_end)
            return false;

          auto&& wmw_window = *wmw_iter;
          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_begin)
          return false;
      }

      {
        auto wmw_iter = wmw_begin;
        auto added_windows_iter = added_windows_begin;
        for (; wmw_iter < wmw_end; wmw_iter++, added_windows_iter++) {
          if (added_windows_iter == added_windows_end)
            return false;

          auto&& wmw_window = *wmw_iter;
          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_end)
          return false;
      }

      {
        auto wmw_iter = wmw_end;
        auto added_windows_iter = added_windows_end;
        while (wmw_iter > wmw_begin) {
          wmw_iter--;
          added_windows_iter--;

          if (added_windows_iter == added_windows_end)
            return false;

          auto&& wmw_window = *wmw_iter;
          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_begin)
          return false;
      }

      {
        auto wmw_iter = wmw_begin;
        auto added_windows_iter = added_windows_begin;
        for (; wmw_iter < wmw_end; wmw_iter += 2, added_windows_iter += 2) {
          if (added_windows_iter == added_windows_end)
            return false;

          auto&& wmw_window = *wmw_iter;
          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_end)
          return false;
      }

      {
        auto wmw_iter = wmw_end;
        auto added_windows_iter = added_windows_end;
        while (wmw_iter > ranges::next(wmw_begin)) {
          wmw_iter -= 2;
          added_windows_iter -= 2;

          if (added_windows_iter == added_windows_end)
            return false;

          auto&& wmw_window = *wmw_iter;
          auto&& added_window = *added_windows_iter;
          if (not ranges::equal(wmw_window, added_window))
            return false;
        }

        if (added_windows_iter != added_windows_begin)
          return false;
      }

      return true;
    };

    using wmw_type = decltype(wmw);
    using wmw_iter_type = decltype(std::forward<wmw_type>(wmw).begin());
    using wmw_iter_expected_type = std::conditional_t<
        std::is_lvalue_reference_v<wmw_type>,
        std::conditional_t<
            std::is_const_v<std::remove_reference_t<wmw_type>>,
            WindowsMergerWindowsIterator<const WindowsMergerWindows>,
            WindowsMergerWindowsIterator<WindowsMergerWindows>>,
        WindowsMergerWindowsIterator<WindowsMergerWindows&&>>;

    using wmw_deref_iter_type = decltype(*std::forward<wmw_type>(wmw).begin());
    using wmw_deref_iter_expected_type = std::conditional_t<
        std::is_lvalue_reference_v<wmw_type>,
        std::conditional_t<
            std::is_const_v<std::remove_reference_t<wmw_type>>,
            WindowsMergerWindowAccessor<const WindowsMergerWindows>,
            WindowsMergerWindowAccessor<WindowsMergerWindows>>,
        WindowsMergerWindowAccessor<WindowsMergerWindows&&>>;

    static_assert(std::is_same_v<wmw_iter_type, wmw_iter_expected_type>);
    static_assert(
        std::is_same_v<wmw_deref_iter_type, wmw_deref_iter_expected_type>);

    return test_with_iterators(std::forward<wmw_type>(wmw).begin(),
                               std::forward<wmw_type>(wmw).end(),
                               ranges::begin(added_windows),
                               ranges::end(added_windows)) and
           test_with_iterators(std::forward<wmw_type>(wmw).rbegin(),
                               std::forward<wmw_type>(wmw).rend(),
                               std::rbegin(added_windows),
                               std::rend(added_windows));
  };

  assert(perform_test(std::as_const(wmw)));
  assert(perform_test(wmw));
  assert(perform_test(std::move(wmw)));
}

static void
test_concrete_base() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr bases_size_type n_bases = 30;

  WindowsMergerWindows wmw(n_clusters, n_bases, 1);
  const auto accessor = wmw[0];

  auto bases = create_random_bases(n_clusters, n_bases);
  for (auto&& base : bases)
    accessor.emplace_back(std::move(base));

  bases = create_random_bases(n_clusters, n_bases);
  ranges::copy(bases, ranges::begin(accessor));
  assert(ranges::equal(accessor, bases));

  bases = create_random_bases(n_clusters, n_bases);
  ranges::move(bases, ranges::begin(accessor));
  assert(ranges::equal(accessor, bases));

  bases = create_random_bases(n_clusters, n_bases);
  ranges::copy(accessor, ranges::begin(bases));
  assert(ranges::equal(accessor, bases));

  bases = create_random_bases(n_clusters, n_bases);
  ranges::move(accessor, ranges::begin(bases));
  assert(ranges::equal(accessor, bases));
}

static void
test_concrete_window() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr bases_size_type n_bases_to_add = 30;

  auto perform_test = [&](bases_size_type initial_bases_capacity) {
    std::vector<WindowsMergerWindow> added_windows(n_windows);
    auto generate_added_windows = [&](bool use_push_back = false) {
      std::generate(
          ranges::begin(added_windows), ranges::end(added_windows),
          [window_index = windows_size_type(0), use_push_back]() mutable {
            WindowsMergerWindow window(n_clusters);
            auto random_bases = create_random_bases(n_clusters, n_bases_to_add);
            for (auto&& base : random_bases) {
              if (use_push_back)
                window.push_back(base);
              else
                window.emplace_back(base);
            }

            assert(window.front() == random_bases.front());
            assert(window.back() == random_bases.back());
            assert(std::as_const(window).front() == random_bases.front());
            assert(std::as_const(window).back() == random_bases.back());
            assert(window.size() == random_bases.size());
            assert(window.coverages().size() == random_bases.size());
            window.set_begin_index(
                static_cast<bases_size_type>(window_index++));
            return window;
          });
    };

// For now, I explicitly don't want copy and move constructors for
// WindowsMergerWindows, because I need to decide if they have to be shallow or
// deep.
//
// Copy or copy not; there is no shallow.
//  - Master Yoda
#define CREATE_WMW()                                                           \
  WindowsMergerWindows wmw(n_clusters, initial_bases_capacity, n_windows);     \
                                                                               \
  generate_added_windows();                                                    \
  for (WindowsMergerTraits::windows_size_type window_index = 0;                \
       window_index < n_windows; ++window_index) {                             \
    assert(added_windows[window_index].begin_index() == window_index);         \
    auto accessor = wmw[window_index];                                         \
    accessor.set_begin_index(static_cast<bases_size_type>(window_index));      \
                                                                               \
    auto random_bases =                                                        \
        create_random_bases(n_clusters, initial_bases_capacity);               \
    for (auto& window_base : random_bases) {                                   \
      accessor.emplace_back(window_base);                                      \
      assert(accessor.back() == window_base);                                  \
    }                                                                          \
    assert(accessor.front() == std::as_const(random_bases).front());           \
    assert(accessor.back() == std::as_const(random_bases).back());             \
    assert(accessor.size() == std::as_const(random_bases).size());             \
    assert(accessor.coverages().size() == std::as_const(random_bases).size()); \
  }                                                                            \
  for (windows_size_type window_index = 0; window_index < n_windows;           \
       ++window_index) {                                                       \
    assert(wmw[window_index].begin_index() == window_index);                   \
  }

    {
      generate_added_windows(true);
      CREATE_WMW()
      ranges::copy(added_windows, ranges::begin(wmw));
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows();
      CREATE_WMW()
      ranges::move(added_windows, ranges::begin(wmw));
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows(true);
      CREATE_WMW()
      ranges::copy(added_windows, std::move(wmw).begin());
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows();
      CREATE_WMW()
      ranges::move(added_windows, std::move(wmw).begin());
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows(true);
      CREATE_WMW()
      ranges::copy(wmw, ranges::begin(added_windows));
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows();
      CREATE_WMW()
      ranges::move(wmw, ranges::begin(added_windows));
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows(true);
      CREATE_WMW()
      ranges::copy(std::move(wmw), ranges::begin(added_windows));
      assert(ranges::equal(wmw, added_windows));
    }

    {
      generate_added_windows();
      CREATE_WMW()
      ranges::move(std::move(wmw), ranges::begin(added_windows));
      assert(ranges::equal(wmw, added_windows));
    }
  };

  perform_test(70);
  perform_test(13);
}

static void
test_concrete_window_iterator() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [](auto get_begin, auto get_end) {
    WindowsMergerWindow window(n_clusters);
    auto random_bases = create_random_bases(n_clusters, n_bases_to_add);
    for (auto&& base : random_bases)
      window.emplace_back(base);

    assert(ranges::equal(get_begin(window), get_end(window),
                         get_begin(random_bases), get_end(random_bases)));

    for (bases_size_type base_index = 0; base_index < n_bases_to_add;
         ++base_index) {
      assert(std::as_const(window)[base_index] == random_bases[base_index]);
      assert(window[base_index] == random_bases[base_index]);
    }
  };

  perform_test([](auto&& data) -> decltype(auto) { return data.begin(); },
               [](auto&& data) -> decltype(auto) { return data.end(); });

  perform_test([](auto&& data) -> decltype(auto) { return data.rbegin(); },
               [](auto&& data) -> decltype(auto) { return data.rend(); });

  perform_test(
      [](auto&& data) -> decltype(auto) { return std::as_const(data).begin(); },
      [](auto&& data) -> decltype(auto) { return std::as_const(data).end(); });

  perform_test(
      [](auto&& data) -> decltype(auto) {
        return std::as_const(data).rbegin();
      },
      [](auto&& data) -> decltype(auto) { return std::as_const(data).rend(); });

  perform_test(
      [](auto&& data) -> decltype(auto) { return std::move(data).begin(); },
      [](auto&& data) -> decltype(auto) { return std::move(data).end(); });

  perform_test(
      [](auto&& data) -> decltype(auto) { return std::move(data).rbegin(); },
      [](auto&& data) -> decltype(auto) { return std::move(data).rend(); });
}

static void
test_concrete_window_coverage_accessor() {
  using bases_size_type = WindowsMergerWindow::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  const auto added_bases = create_random_bases(n_clusters, n_bases_to_add);

  WindowsMergerWindow window(n_clusters);
  for (auto&& base : added_bases)
    window.emplace_back(std::move(base));

  auto perform_test = [&](auto window_get) {
    auto added_bases_iter = ranges::begin(added_bases);
    auto coverages = window_get(window).coverages();
    assert(coverages.size() == n_bases_to_add);
    assert(coverages.front() == added_bases.front().coverage());
    assert(coverages.back() == added_bases.back().coverage());

    {
      auto coverages_iter = ranges::begin(coverages);
      const auto coverages_end = ranges::end(coverages);

      for (; coverages_iter < coverages_end;
           ++coverages_iter, ++added_bases_iter) {
        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == ranges::end(added_bases));
    }

    for (bases_size_type base_index = 0; base_index < n_bases_to_add;
         ++base_index)
      assert(coverages[base_index] == added_bases[base_index].coverage());
  };

  perform_test([](auto&& window) -> decltype(auto) { return window; });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); });
}

static void
test_concrete_window_coverage_iterator() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  const auto added_bases = create_random_bases(n_clusters, n_bases_to_add);

  auto perform_test = [&](auto window_get, auto get_begin, auto get_end) {
    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_begin(added_bases);
      auto coverages_iter = get_begin(coverages);
      const auto coverages_end = get_end(coverages);

      assert(*(coverages_iter + 2) == (added_bases_iter + 2)->coverage());
      assert(*(2 + coverages_iter) == (added_bases_iter + 2)->coverage());
      {
        auto iter = coverages_iter;
        ++iter;
        ++iter;
        assert(iter >= coverages_iter);
        assert(coverages_iter <= iter);
        assert(iter - 2 == coverages_iter);
      }
      assert(coverages_iter <= coverages_iter);
      assert(coverages_iter >= coverages_iter);
      assert(ranges::distance(coverages_iter, coverages_end) == n_bases_to_add);

      for (; coverages_iter < coverages_end;
           ++coverages_iter, ++added_bases_iter) {
        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_end(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_end(added_bases);
      auto coverages_iter = get_end(coverages);
      const auto coverages_begin = get_begin(coverages);

      for (; coverages_iter > coverages_begin;) {
        --coverages_iter;
        --added_bases_iter;

        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_begin(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_begin(added_bases);
      auto coverages_iter = get_begin(coverages);
      const auto coverages_end = get_end(coverages);

      for (; coverages_iter < coverages_end;
           coverages_iter++, added_bases_iter++) {
        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_end(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_end(added_bases);
      auto coverages_iter = get_end(coverages);
      const auto coverages_begin = get_begin(coverages);

      for (; coverages_iter > coverages_begin;) {
        coverages_iter--;
        added_bases_iter--;

        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_begin(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_begin(added_bases);
      auto coverages_iter = get_begin(coverages);
      const auto coverages_end = get_end(coverages);

      for (; coverages_iter < coverages_end;
           coverages_iter += 2, added_bases_iter += 2) {
        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_end(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();

      auto added_bases_iter = get_end(added_bases);
      auto coverages_iter = get_end(coverages);
      const auto coverages_begin = get_begin(coverages);

      for (; coverages_iter > coverages_begin;) {
        coverages_iter -= 2;
        added_bases_iter -= 2;

        assert(*coverages_iter == added_bases_iter->coverage());
      }
      assert(added_bases_iter == get_begin(added_bases));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto coverages = window_get(window).coverages();
      const auto coverages_begin = get_begin(coverages);
      const auto added_bases_begin = get_begin(added_bases);
      for (bases_size_type base_index = 0; base_index < n_bases_to_add;
           ++base_index) {
        assert(coverages_begin[base_index] ==
               added_bases_begin[base_index].coverage());
      }
    }
  };

  perform_test([](auto&& window) -> decltype(auto) { return window; },
               [](auto&& data) -> decltype(auto) { return data.begin(); },
               [](auto&& data) -> decltype(auto) { return data.end(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); },
      [](auto&& data) -> decltype(auto) { return data.begin(); },
      [](auto&& data) -> decltype(auto) { return data.end(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); },
      [](auto&& data) -> decltype(auto) { return data.begin(); },
      [](auto&& data) -> decltype(auto) { return data.end(); });

  perform_test([](auto&& window) -> decltype(auto) { return window; },
               [](auto&& data) -> decltype(auto) { return data.rbegin(); },
               [](auto&& data) -> decltype(auto) { return data.rend(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); },
      [](auto&& data) -> decltype(auto) { return data.rbegin(); },
      [](auto&& data) -> decltype(auto) { return data.rend(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); },
      [](auto&& data) -> decltype(auto) { return data.rbegin(); },
      [](auto&& data) -> decltype(auto) { return data.rend(); });
}

static void
test_concrete_window_weights_accessor() {
  using clusters_size_type = WindowsMergerWindow::clusters_size_type;
  using bases_size_type = WindowsMergerWindow::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  auto perform_test = [&](auto window_get) {
    const auto added_bases = create_random_bases(n_clusters, n_bases_to_add);
    WindowsMergerWindow window(n_clusters);
    for (auto&& base : added_bases)
      window.emplace_back(base);

    auto linear_weights = window_get(window).linear_weights();
    assert(linear_weights.size() == n_bases_to_add * n_clusters);
    assert(linear_weights.front() == added_bases.front().weight(0));
    assert(linear_weights.back() == added_bases.back().weight(n_clusters - 1));

    {
      auto linear_weights_iter = ranges::begin(linear_weights);
      const auto linear_weights_end = ranges::end(linear_weights);
      assert(ranges::distance(linear_weights_iter, linear_weights_end) ==
             n_clusters * n_bases_to_add);
      auto added_bases_iter = ranges::begin(added_bases);
      clusters_size_type weight_index = 0;

      for (; linear_weights_iter < linear_weights_end;
           ++linear_weights_iter, ++weight_index) {
        if (weight_index == n_clusters) {
          weight_index = 0;
          ++added_bases_iter;
        }
        assert(*linear_weights_iter == added_bases_iter->weight(weight_index));
      }
      assert(added_bases_iter == std::prev(ranges::end(added_bases)));
      assert(weight_index == n_clusters);
    }

    {
      bases_size_type linear_weight_index = 0;
      for (bases_size_type base_index = 0; base_index < n_bases_to_add;
           ++base_index) {
        const auto& base = added_bases[base_index];
        for (clusters_size_type cluster_index = 0; cluster_index < n_clusters;
             ++cluster_index) {
          assert(linear_weights[linear_weight_index++] ==
                 base.weight(cluster_index));
        }
      }
    }
  };

  perform_test([](auto&& window) -> decltype(auto) { return window; });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); });
}

static void
test_concrete_window_weights_iterator() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  using weight_type = typename WindowsMergerTraits::weight_type;

  constexpr std::size_t n_clusters = 5;
  constexpr std::size_t n_bases_to_add = 30;

  const auto added_bases = create_random_bases(n_clusters, n_bases_to_add);
  const auto added_linear_weights = [&added_bases] {
    std::vector<weight_type> linear_weights;
    linear_weights.reserve(n_clusters * n_bases_to_add);

    for (const auto& base : added_bases) {
      for (auto&& weight : base.weights()) {
        linear_weights.emplace_back(weight);
      }
    }

    return linear_weights;
  }();

  auto perform_test = [&](auto window_get, auto get_begin, auto get_end) {
    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_begin(added_linear_weights);
      auto linear_weights_iter = get_begin(linear_weights);
      const auto linear_weights_end = get_end(linear_weights);

      assert(*(linear_weights_iter + 2) == *(added_linear_weights_iter + 2));
      assert(*(2 + linear_weights_iter) == *(added_linear_weights_iter + 2));
      {
        auto iter = linear_weights_iter;
        ++iter;
        ++iter;
        assert(iter >= linear_weights_iter);
        assert(linear_weights_iter <= iter);
        assert(iter - 2 == linear_weights_iter);
      }
      assert(linear_weights_iter <= linear_weights_iter);
      assert(linear_weights_iter >= linear_weights_iter);
      assert(ranges::distance(linear_weights_iter, linear_weights_end) ==
             n_bases_to_add * n_clusters);

      for (; linear_weights_iter < linear_weights_end;
           ++linear_weights_iter, ++added_linear_weights_iter) {
        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_end(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_end(added_linear_weights);
      auto linear_weights_iter = get_end(linear_weights);
      const auto linear_weights_begin = get_begin(linear_weights);

      for (; linear_weights_iter > linear_weights_begin;) {
        --linear_weights_iter;
        --added_linear_weights_iter;

        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_begin(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_begin(added_linear_weights);
      auto linear_weights_iter = get_begin(linear_weights);
      const auto linear_weights_end = get_end(linear_weights);

      for (; linear_weights_iter < linear_weights_end;
           linear_weights_iter++, added_linear_weights_iter++) {
        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_end(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_end(added_linear_weights);
      auto linear_weights_iter = get_end(linear_weights);
      const auto linear_weights_begin = get_begin(linear_weights);

      for (; linear_weights_iter > linear_weights_begin;) {
        linear_weights_iter--;
        added_linear_weights_iter--;

        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_begin(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_begin(added_linear_weights);
      auto linear_weights_iter = get_begin(linear_weights);
      const auto linear_weights_end = get_end(linear_weights);

      for (; linear_weights_iter < linear_weights_end;
           linear_weights_iter += 2, added_linear_weights_iter += 2) {
        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_end(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();

      auto added_linear_weights_iter = get_end(added_linear_weights);
      auto linear_weights_iter = get_end(linear_weights);
      const auto linear_weights_begin = get_begin(linear_weights);

      for (; linear_weights_iter > linear_weights_begin;) {
        linear_weights_iter -= 2;
        added_linear_weights_iter -= 2;

        assert(*linear_weights_iter == *added_linear_weights_iter);
      }
      assert(added_linear_weights_iter == get_begin(added_linear_weights));
    }

    {
      WindowsMergerWindow window(n_clusters);
      for (auto&& base : added_bases)
        window.emplace_back(base);

      auto linear_weights = window_get(window).linear_weights();
      const auto linear_weights_begin = get_begin(linear_weights);
      const auto added_linear_weights_begin = get_begin(added_linear_weights);
      for (bases_size_type weight_index = 0;
           weight_index < n_bases_to_add * n_clusters; ++weight_index) {
        assert(linear_weights_begin[weight_index] ==
               added_linear_weights_begin[weight_index]);
      }
    }
  };

  perform_test([](auto&& window) -> decltype(auto) { return window; },
               [](auto&& data) -> decltype(auto) { return data.begin(); },
               [](auto&& data) -> decltype(auto) { return data.end(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); },
      [](auto&& data) -> decltype(auto) { return data.begin(); },
      [](auto&& data) -> decltype(auto) { return data.end(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); },
      [](auto&& data) -> decltype(auto) { return data.begin(); },
      [](auto&& data) -> decltype(auto) { return data.end(); });

  perform_test([](auto&& window) -> decltype(auto) { return window; },
               [](auto&& data) -> decltype(auto) { return data.rbegin(); },
               [](auto&& data) -> decltype(auto) { return data.rend(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::as_const(window); },
      [](auto&& data) -> decltype(auto) { return data.rbegin(); },
      [](auto&& data) -> decltype(auto) { return data.rend(); });
  perform_test(
      [](auto&& window) -> decltype(auto) { return std::move(window); },
      [](auto&& data) -> decltype(auto) { return data.rbegin(); },
      [](auto&& data) -> decltype(auto) { return data.rend(); });
}

static void
test_windows_sort() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr windows_size_type n_windows = 80;
  constexpr bases_size_type n_bases = 30;
  constexpr bases_size_type max_sequence_size = 300;
  std::vector<decltype(create_random_bases(n_clusters, n_bases))> raw_windows(
      n_windows);
  std::generate(ranges::begin(raw_windows), ranges::end(raw_windows),
                [] { return create_random_bases(n_clusters, n_bases); });
  WindowsMergerWindows wmw(n_clusters, n_bases, n_windows);
  auto all_starts = create_random_starts(n_windows, n_bases, max_sequence_size);

  {
    windows_size_type window_index = 0;
    for (const auto& window : raw_windows) {
      auto window_accessor = wmw[window_index];
      window_accessor.set_begin_index(all_starts[window_index]);
      for (const auto& base : window)
        window_accessor.emplace_back(base);

      ++window_index;
    }
  }

  {
    std::vector<std::size_t> indices(n_windows);
    std::iota(ranges::begin(indices), ranges::end(indices), std::size_t(0));
    ranges::stable_sort(indices, ranges::less{},
                        [&](std::size_t index) { return all_starts[index]; });

    decltype(raw_windows) new_raw_windows;
    new_raw_windows.reserve(n_windows);
    std::transform(ranges::begin(indices), ranges::end(indices),
                   std::back_inserter(new_raw_windows),
                   [&](std::size_t index) -> decltype(auto) {
                     return std::move(raw_windows[index]);
                   });

    raw_windows = std::move(new_raw_windows);
    ranges::stable_sort(all_starts);
  }

  ranges::stable_sort(wmw, ranges::less{},
                      [](auto&& window) { return window.begin_index(); });
  assert(ranges::is_sorted(wmw, ranges::less{},
                           [](auto&& window) { return window.begin_index(); }));

  for (windows_size_type window_index = 0; window_index < n_windows;
       ++window_index) {
    auto&& window_accessor = wmw[window_index];
    auto&& raw_window = raw_windows[window_index];

    assert(window_accessor.begin_index() == all_starts[window_index]);
    assert(window_accessor.end_index() == all_starts[window_index] + n_bases);
    assert(ranges::equal(window_accessor, raw_window));
  }
}

static void
test_window_accessor_resize() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr bases_size_type n_bases = 30;

  WindowsMergerWindows wmw(n_clusters, n_bases, 1);
  const auto accessor = wmw[0];

  auto bases = create_random_bases(n_clusters, n_bases);
  for (auto&& base : bases)
    accessor.emplace_back(base);
  assert(ranges::equal(bases, accessor));

  accessor.resize(60);
  assert(wmw.bases_capacity() == 60);
  assert(accessor.size() == 60);
  assert(accessor.coverages().size() == 60);
  assert(ranges::equal(bases, accessor | ranges::view::take(bases.size())));

  WindowsMergerWindowBase default_base(n_clusters);
  ranges::fill(default_base.weights(), TinyFraction(0.));
  default_base.coverage() = 0;

  assert(ranges::all_of(
      ranges::next(ranges::begin(accessor), n_bases), ranges::end(accessor),
      [&default_base](auto&& base) { return base == default_base; }));

  accessor.resize(20);
  assert(wmw.bases_capacity() == 60);
  assert(accessor.size() == 20);
  assert(accessor.coverages().size() == 20);
  assert(ranges::equal(ranges::begin(bases),
                       ranges::next(ranges::begin(bases), 20),
                       ranges::begin(accessor), ranges::end(accessor)));
}

static void
test_move_constructor() {
  {
    WindowsMergerWindows wmw;
    WindowsMergerWindows new_wmw(std::move(wmw));

    assert(new_wmw.clusters_size() == 0);
    assert(new_wmw.bases_capacity() == 0);
    assert(new_wmw.windows_size() == 0);
    assert(new_wmw.windows_capacity() == 0);
  }

  {
    WindowsMergerWindows wmw(5);
    WindowsMergerWindows new_wmw(std::move(wmw));

    assert(new_wmw.clusters_size() == 5);
    assert(new_wmw.bases_capacity() == 0);
    assert(new_wmw.windows_size() == 0);
    assert(new_wmw.windows_capacity() == 0);
  }

  {
    WindowsMergerWindows wmw(5, 70);
    WindowsMergerWindows new_wmw(std::move(wmw));

    assert(new_wmw.clusters_size() == 5);
    assert(new_wmw.bases_capacity() == 70);
    assert(new_wmw.windows_size() == 0);
    assert(new_wmw.windows_capacity() == 0);

    assert(wmw.bases_capacity() == 0);
    assert(wmw.windows_size() == 0);
    assert(wmw.windows_capacity() == 0);
  }

  {
    WindowsMergerWindows wmw(5, 70, 100);
    WindowsMergerWindows new_wmw(std::move(wmw));

    assert(new_wmw.clusters_size() == 5);
    assert(new_wmw.bases_capacity() == 70);
    assert(new_wmw.windows_size() == 100);
    assert(new_wmw.windows_capacity() == 100);

    assert(wmw.bases_capacity() == 0);
    assert(wmw.windows_size() == 0);
    assert(wmw.windows_capacity() == 0);
  }

  {
    using bases_size_type = typename WindowsMergerTraits::bases_size_type;
    using windows_size_type = typename WindowsMergerTraits::windows_size_type;

    constexpr std::size_t n_clusters = 5;
    constexpr windows_size_type n_windows = 20;
    constexpr bases_size_type n_bases = 30;
    std::vector<decltype(create_random_bases(n_clusters, n_bases))> raw_windows(
        n_windows);
    std::generate(ranges::begin(raw_windows), ranges::end(raw_windows),
                  [] { return create_random_bases(n_clusters, n_bases); });

    WindowsMergerWindows wmw(n_clusters);
    for (auto&& window : std::as_const(raw_windows)) {
      auto&& new_window = wmw.emplace_back();
      for (auto&& base : std::as_const(window))
        new_window.push_back(base);
    }

    WindowsMergerWindows new_wmw(std::move(wmw));
    assert(new_wmw.clusters_size() == n_clusters);
    assert(new_wmw.windows_size() == n_windows);
    assert(new_wmw.windows_capacity() >= n_windows);
    assert(new_wmw.bases_capacity() >= n_bases);

    assert(wmw.bases_capacity() == 0);
    assert(wmw.windows_size() == 0);
    assert(wmw.windows_capacity() == 0);

    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      const auto& raw_window = raw_windows[window_index];
      auto&& wmw_window = std::as_const(new_wmw)[window_index];

      for (bases_size_type base_index = 0; base_index < n_bases; ++base_index)
        assert(wmw_window[base_index] == raw_window[base_index]);
    }
  }
}

static void
test_move_assignment() {
  auto perform_test = [](auto&& create_wmw) {
    {
      WindowsMergerWindows wmw;
      WindowsMergerWindows new_wmw = create_wmw();
      new_wmw = std::move(wmw);

      assert(new_wmw.clusters_size() == 0);
      assert(new_wmw.bases_capacity() == 0);
      assert(new_wmw.windows_size() == 0);
      assert(new_wmw.windows_capacity() == 0);
    }

    {
      WindowsMergerWindows wmw(5);
      WindowsMergerWindows new_wmw = create_wmw();
      new_wmw = std::move(wmw);

      assert(new_wmw.clusters_size() == 5);
      assert(new_wmw.bases_capacity() == 0);
      assert(new_wmw.windows_size() == 0);
      assert(new_wmw.windows_capacity() == 0);
    }

    {
      WindowsMergerWindows wmw(5, 70);
      WindowsMergerWindows new_wmw = create_wmw();
      new_wmw = std::move(wmw);

      assert(new_wmw.clusters_size() == 5);
      assert(new_wmw.bases_capacity() == 70);
      assert(new_wmw.windows_size() == 0);
      assert(new_wmw.windows_capacity() == 0);

      assert(wmw.bases_capacity() == 0);
      assert(wmw.windows_size() == 0);
      assert(wmw.windows_capacity() == 0);
    }

    {
      WindowsMergerWindows wmw(5, 70, 100);
      WindowsMergerWindows new_wmw = create_wmw();
      new_wmw = std::move(wmw);

      assert(new_wmw.clusters_size() == 5);
      assert(new_wmw.bases_capacity() == 70);
      assert(new_wmw.windows_size() == 100);
      assert(new_wmw.windows_capacity() == 100);

      assert(wmw.bases_capacity() == 0);
      assert(wmw.windows_size() == 0);
      assert(wmw.windows_capacity() == 0);
    }

    {
      using bases_size_type = typename WindowsMergerTraits::bases_size_type;
      using windows_size_type = typename WindowsMergerTraits::windows_size_type;

      constexpr std::size_t n_clusters = 5;
      constexpr windows_size_type n_windows = 20;
      constexpr bases_size_type n_bases = 30;
      std::vector<decltype(create_random_bases(n_clusters, n_bases))>
          raw_windows(n_windows);
      std::generate(ranges::begin(raw_windows), ranges::end(raw_windows),
                    [] { return create_random_bases(n_clusters, n_bases); });

      WindowsMergerWindows wmw(n_clusters);
      for (auto&& window : std::as_const(raw_windows)) {
        auto&& new_window = wmw.emplace_back();
        for (auto&& base : std::as_const(window))
          new_window.push_back(base);
      }

      WindowsMergerWindows new_wmw = create_wmw();
      new_wmw = std::move(wmw);

      assert(new_wmw.clusters_size() == n_clusters);
      assert(new_wmw.windows_size() == n_windows);
      assert(new_wmw.windows_capacity() >= n_windows);
      assert(new_wmw.bases_capacity() >= n_bases);

      assert(wmw.bases_capacity() == 0);
      assert(wmw.windows_size() == 0);
      assert(wmw.windows_capacity() == 0);

      for (windows_size_type window_index = 0; window_index < n_windows;
           ++window_index) {
        const auto& raw_window = raw_windows[window_index];
        auto&& wmw_window = std::as_const(new_wmw)[window_index];

        for (bases_size_type base_index = 0; base_index < n_bases; ++base_index)
          assert(wmw_window[base_index] == raw_window[base_index]);
      }
    }
  };

  perform_test([] { return WindowsMergerWindows{}; });
  perform_test([] { return WindowsMergerWindows(5); });
  perform_test([] { return WindowsMergerWindows(5, 70); });
  perform_test([] { return WindowsMergerWindows(5, 70, 100); });
  perform_test([] {
    using bases_size_type = typename WindowsMergerTraits::bases_size_type;
    using windows_size_type = typename WindowsMergerTraits::windows_size_type;

    constexpr std::size_t n_clusters = 5;
    constexpr windows_size_type n_windows = 20;
    constexpr bases_size_type n_bases = 30;

    WindowsMergerWindows wmw(n_clusters);
    for (windows_size_type window_index = 0; window_index < n_windows;
         ++window_index) {
      auto bases = create_random_bases(n_clusters, n_bases);
      auto&& new_window = wmw.emplace_back();
      for (auto&& base : bases)
        new_window.push_back(std::move(base));
    }

    return wmw;
  });
}

static void
test_window_accessor_clear() {
  using bases_size_type = typename WindowsMergerTraits::bases_size_type;
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr std::size_t n_clusters = 5;
  constexpr windows_size_type n_windows = 20;
  constexpr bases_size_type n_bases = 30;
  std::vector<decltype(create_random_bases(n_clusters, n_bases))> raw_windows(
      n_windows);
  std::generate(ranges::begin(raw_windows), ranges::end(raw_windows),
                [] { return create_random_bases(n_clusters, n_bases); });

  WindowsMergerWindows wmw(n_clusters);
  for (auto&& window : std::as_const(raw_windows)) {
    auto&& new_window = wmw.emplace_back();
    for (auto&& base : std::as_const(window))
      new_window.push_back(base);
  }

  const auto window = wmw[n_windows / 2];
  window.clear();
  assert(window.size() == 0);
  assert(window.coverages().size() == 0);
}

static void
test_reorder_clusters() {
  using clusters_size_type = WindowsMergerWindow::clusters_size_type;
  using bases_size_type = WindowsMergerWindow::bases_size_type;

  constexpr clusters_size_type n_clusters = 8;
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr bases_size_type n_bases = 70;

  std::vector<decltype(create_random_bases(n_clusters, n_bases))> added_windows;
  auto added_windows_starts = create_random_starts(n_windows, n_bases, 300);

  WindowsMergerWindows wmw(n_clusters, n_bases, n_windows);
  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < n_windows; ++window_index) {
    auto bases_to_add = create_random_bases(n_clusters, n_bases);
    auto accessor = wmw[window_index];
    for (const auto& window_base : bases_to_add)
      accessor.emplace_back(window_base);
    accessor.set_begin_index(added_windows_starts[window_index]);
    added_windows.emplace_back(std::move(bases_to_add));
  }

  std::vector<clusters_size_type> new_clusters_order(n_clusters);
  ranges::iota(new_clusters_order, clusters_size_type(0));
  {
    std::mt19937 random_gen(std::random_device{}());
    ranges::shuffle(new_clusters_order, random_gen);
  }

  for (auto&& window : wmw) {
    window.reorder_clusters(new_clusters_order);
  }

  {
    auto wmw_window_iter = ranges::cbegin(wmw);
    auto const wmw_end = ranges::cend(wmw);
    auto added_windows_iter = ranges::cbegin(added_windows);

    for (; wmw_window_iter < wmw_end; ++wmw_window_iter, ++added_windows_iter) {
      auto&& window = *wmw_window_iter;
      auto&& added_window = *added_windows_iter;

      auto window_iter = ranges::cbegin(window);
      auto const window_end = ranges::cend(window);
      auto added_window_iter = ranges::cbegin(added_window);
      for (; window_iter < window_end; ++window_iter, ++added_window_iter) {
        auto&& base = *window_iter;
        auto&& added_window_base = *added_window_iter;

        for (auto cluster_index = clusters_size_type(0);
             cluster_index < n_clusters; ++cluster_index) {
          assert(base.weight(cluster_index) ==
                 added_window_base.weight(new_clusters_order[cluster_index]));
        }
      }
    }
  }
}

int
main() {
  static_assert(std::is_nothrow_default_constructible_v<WindowsMergerWindows>);
  static_assert(std::is_move_constructible_v<WindowsMergerWindows>);
  static_assert(std::is_move_assignable_v<WindowsMergerWindows>);
  /*
  static_assert(std::is_copy_constructible_v<WindowsMergerWindows>);
  static_assert(std::is_copy_assignable_v<WindowsMergerWindows>);
  */

  test_default_construction();
  test_construction_with_clusters();
  test_construction_with_clusters_and_bases();
  test_construction_with_clusters_bases_and_windows();
  test_simple_const_accessor();
  test_simple_accessor_emplace_back();
  test_simple_accessor_push_back();
  test_front_back_accessors();
  test_base_accessor_const_iterator();
  test_base_accessor_iterator();
  test_base_accessor_weights_accessor();
  test_base_comparisons();
  test_accessor_resizing_push();
  test_accessor_resizing_push_on_empty();
  test_windows_reshape_noreshape();
  test_windows_reshape_from_empty();
  test_windows_reshape_shrink();
  test_windows_emplace_back();
  test_window_set_begin_index();
  test_window_coverage_accessor();
  test_window_accessor_iterator();
  test_accessor_assignment();
  test_accessor_swap();
  test_accessor_assignment_reshape();
  test_accessor_move_assignment_destroy();
  test_accessor_move_assignment_reshape();
  test_accessor_assignment_destroy();
  test_accessor_swap_reshape_lhs();
  test_accessor_swap_reshape_rhs();
  test_windows_iterators();
  test_concrete_base();
  test_concrete_window();
  test_concrete_window_iterator();
  test_base_accessor_assignment();
  test_window_accessor_to_window_constructor();
  test_window_accessor_to_window_assignment();
  test_window_base_accessor_to_base_constructor();
  test_window_base_accessor_to_base_assignment();
  test_concrete_window_coverage_accessor();
  test_concrete_window_coverage_iterator();
  test_concrete_window_weights_accessor();
  test_concrete_window_weights_iterator();
  test_windows_sort();
  test_window_accessor_resize();
  test_window_accessor_clear();
  test_move_constructor();
  test_move_assignment();
  test_reorder_clusters();
}
