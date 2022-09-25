#include "windows_merger_cache_indices.hpp"

#include <algorithm>
#include <cassert>
#include <random>

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

using namespace windows_merger;

static void test_default_construction() {
  WindowsMergerCacheIndices wmci;
  assert(wmci.line_capacity() == 0);
  assert(wmci.size() == 0);
  assert(wmci.capacity() == 0);
}

static void test_construction_with_clusters_and_bases() {
  WindowsMergerCacheIndices wmci(70);
  assert(wmci.line_capacity() == 70);
  assert(wmci.size() == 0);
  assert(wmci.capacity() == 0);
}

static void test_construction_with_clusters_bases_and_windows() {
  WindowsMergerCacheIndices wmci(70, 100);
  assert(wmci.line_capacity() == 70);
  assert(wmci.size() == 100);
  assert(wmci.capacity() == 100);
}

static void test_simple_const_accessor() {
  const WindowsMergerCacheIndices wmci(70, 20);
  for (std::size_t index = 0; index < 20; ++index) {
    auto accessor = wmci[index];
    static_assert(std::is_same_v<decltype(accessor),
                                 WindowsMergerCacheIndicesAccessor<
                                     const WindowsMergerCacheIndices>>);

    assert(accessor.index() == index);
    assert(accessor.size() == 0);
    assert(accessor.empty());
    (void)accessor;
  }
}

static std::vector<WindowsMergerTraits::windows_size_type>
create_random_indices(std::size_t n_indices) {
  using indices_type = std::vector<WindowsMergerTraits::windows_size_type>;
  std::mt19937 random_gen(std::random_device{}());
  std::uniform_real_distribution<double> real_dist(0., 1.);
  std::uniform_int_distribution<std::size_t> indices_size_dist(1, 20);
  std::uniform_int_distribution<unsigned> indices_dist(0, 300);

  indices_type out_indices(n_indices);
  auto start_iter = ranges::begin(out_indices);
  do {
    std::generate(start_iter, ranges::end(out_indices),
                  [&] { return indices_dist(random_gen); });
    std::sort(ranges::begin(out_indices), ranges::end(out_indices));
    start_iter =
        std::unique(ranges::begin(out_indices), ranges::end(out_indices));
  } while (start_iter != ranges::end(out_indices));

  return out_indices;
}

static void test_simple_accessor_emplace_back() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, 20);
  {
    auto accessor = wmci[0];
    for (const auto &index : indices_to_add) {
      auto &&new_index = accessor.emplace_back(index);
      assert(new_index == index);
    }

    assert(std::size(accessor) == n_indices_to_add);
    assert(not accessor.empty());
  }
}

static void test_simple_accessor_push_back() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, 20);
  {
    auto accessor = wmci[0];
    for (const auto &index : indices_to_add) {
      accessor.push_back(index);
      assert(accessor.back() == index);
    }

    assert(std::size(accessor) == n_indices_to_add);
    assert(not accessor.empty());
  }
}

static void test_index_accessor_const_iterator() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, 20);

  auto accessor = wmci[0];
  for (const auto &index : indices_to_add)
    accessor.emplace_back(index);

  assert(not accessor.empty());
  for (typename WindowsMergerTraits::windows_size_type index_index = 0;
       index_index < n_indices_to_add; ++index_index) {
    const auto &index = accessor[index_index];
    auto &&index_added = indices_to_add[index_index];

    assert(index == index_added);
  }
}

static void test_index_accessor_iterator() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, 20);

  auto accessor = wmci[0];
  for (const auto &index : indices_to_add)
    accessor.emplace_back(index);

  for (typename WindowsMergerTraits::windows_size_type index_index = 0;
       index_index < n_indices_to_add; ++index_index) {
    auto &index = accessor[index_index];
    auto &&index_added = indices_to_add[index_index];

    ++index;
    assert(index != index_added);
    assert(index == index_added + 1);
  }
}

static void test_accessor_resizing_push() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(1, 20);
  auto accessor = wmci[0];
  for (const auto &index : indices_to_add) {
    auto const &new_index = accessor.emplace_back(index);
    assert(new_index == index);
    assert(new_index == accessor.back());
  }

  assert(wmci.line_capacity() >= std::size(indices_to_add));
  for (WindowsMergerTraits::windows_size_type index_index = 0;
       index_index < std::size(indices_to_add); ++index_index)
    assert(accessor[index_index] == indices_to_add[index_index]);
}

static void test_accessor_resizing_push_on_empty() {
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  {
    WindowsMergerCacheIndices wmci(0, 20);
    auto accessor = wmci[0];
    for (const auto &index : indices_to_add) {
      auto const &new_index = accessor.emplace_back(index);
      assert(new_index == index);
      assert(new_index == accessor.back());
    }

    assert(wmci.line_capacity() >= std::size(indices_to_add));
    for (WindowsMergerTraits::windows_size_type index_index = 0;
         index_index < std::size(indices_to_add); ++index_index)
      assert(accessor[index_index] == indices_to_add[index_index]);
  }

  {
    WindowsMergerCacheIndices wmci(0, 20);
    auto accessor = wmci[0];
    for (const auto &indices : indices_to_add) {
      accessor.push_back(indices);
      assert(accessor.back() == indices);
    }

    assert(wmci.line_capacity() >= std::size(indices_to_add));
    for (WindowsMergerTraits::windows_size_type index_index = 0;
         index_index < std::size(indices_to_add); ++index_index)
      assert(accessor[index_index] == indices_to_add[index_index]);
  }

  {
    WindowsMergerCacheIndices wmci(0, 20);
    auto accessor = wmci[0];
    for (auto &index : indices_to_add) {
      accessor.push_back(index);
      assert(accessor.back() == index);
    }

    assert(wmci.line_capacity() >= std::size(indices_to_add));
    for (WindowsMergerTraits::windows_size_type index_index = 0;
         index_index < std::size(indices_to_add); ++index_index)
      assert(accessor[index_index] == indices_to_add[index_index]);
  }
}

static void test_windows_reshape_noreshape() {
  constexpr std::size_t n_indices_to_add = 30;

  {
    WindowsMergerCacheIndices wmci;
    wmci.reshape(0, typename WindowsMergerCacheIndices::reserver_type(0));
  }

  {
    WindowsMergerCacheIndices wmci;
    wmci.reshape(0, typename WindowsMergerCacheIndices::resizer_type(0));
  }

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, 20);

  {
    auto accessor = wmci[0];
    for (const auto &window_base : indices_to_add)
      accessor.emplace_back(window_base);
  }

  {
    auto first_element_address = &wmci[0].front();
    wmci.reshape(0, typename WindowsMergerCacheIndices::reserver_type(0));
    assert(first_element_address == &wmci[0].front());
  }

  {
    auto first_element_address = &wmci[0].front();
    wmci.reshape(70, typename WindowsMergerCacheIndices::reserver_type(20));
    assert(first_element_address == &wmci[0].front());
  }

  {
    auto first_element_address = &wmci[0].front();
    wmci.reshape(70, typename WindowsMergerCacheIndices::resizer_type(20));
    assert(first_element_address == &wmci[0].front());
  }
}

static void test_windows_reshape_from_empty() {
  constexpr std::size_t n_indices_to_add = 20;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  auto add_to_wmci = [&](WindowsMergerCacheIndices &wmci) {
    if (wmci.size() == 0)
      wmci.reshape(0, WindowsMergerCacheIndices::resizer_type(1));

    auto accessor = wmci[0];
    for (const auto &window_base : indices_to_add)
      accessor.emplace_back(window_base);

    assert(wmci.line_capacity() == 70);
    assert(wmci.capacity() == 20);
  };

  auto check_wmci = [&](const WindowsMergerCacheIndices &wmci) {
    auto accessor = wmci[0];
    assert(std::size(accessor) == n_indices_to_add);
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < std::size(accessor); ++base_index) {
      assert(accessor[base_index] == indices_to_add[base_index]);
    }

    assert(wmci.line_capacity() == 70);
    assert(wmci.capacity() == 30);
  };

  {
    WindowsMergerCacheIndices wmci(0, 20);
    assert(wmci.size() == 20);
    wmci.reshape(70, typename WindowsMergerCacheIndices::reserver_type(20));
    add_to_wmci(wmci);
    assert(wmci.size() == 20);

    wmci.reshape(70, typename WindowsMergerCacheIndices::reserver_type(30));
    check_wmci(wmci);
    assert(wmci.capacity() == 30);
    assert(wmci.size() == 20);
  }

  {
    WindowsMergerCacheIndices wmci(0, 20);
    assert(wmci.size() == 20);
    wmci.reshape(70, typename WindowsMergerCacheIndices::resizer_type(20));
    add_to_wmci(wmci);
    assert(wmci.size() == 20);

    wmci.reshape(70, typename WindowsMergerCacheIndices::resizer_type(30));
    check_wmci(wmci);
    assert(wmci.size() == 30);
  }

  {
    WindowsMergerCacheIndices wmci(70, 0);
    assert(wmci.size() == 0);
    wmci.reshape(70, typename WindowsMergerCacheIndices::reserver_type(20));
    add_to_wmci(wmci);
    assert(wmci.size() == 1);

    wmci.reshape(70, typename WindowsMergerCacheIndices::reserver_type(30));
    check_wmci(wmci);
    assert(wmci.size() == 1);
  }

  {
    WindowsMergerCacheIndices wmci(70, 0);
    assert(wmci.size() == 0);
    wmci.reshape(70, typename WindowsMergerCacheIndices::resizer_type(20));
    add_to_wmci(wmci);
    assert(wmci.size() == 20);

    wmci.reshape(70, typename WindowsMergerCacheIndices::resizer_type(30));
    check_wmci(wmci);
    assert(wmci.size() == 30);
  }
}

static void test_windows_reshape_shrink() {
  constexpr WindowsMergerTraits::windows_size_type n_windows = 20;
  constexpr std::size_t n_indices_to_add = 30;

  auto indices_to_add = create_random_indices(n_indices_to_add);
  WindowsMergerCacheIndices wmci(70, n_windows);

  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < n_windows; ++window_index) {
    auto accessor = wmci[window_index];
    for (const auto &window_base : indices_to_add)
      accessor.emplace_back(window_base);
  }

  wmci.reshape(50, typename WindowsMergerCacheIndices::resizer_type(10));
  assert(wmci.line_capacity() == 70);
  assert(wmci.capacity() == 20);
  assert(wmci.size() == 10);
  for (WindowsMergerTraits::windows_size_type window_index = 0;
       window_index < 10; ++window_index) {
    auto accessor = wmci[window_index];
    const auto accessor_size = std::size(accessor);
    assert(accessor_size == 30);
    for (WindowsMergerTraits::bases_size_type base_index = 0;
         base_index < accessor_size; ++base_index)
      assert(accessor[base_index] == indices_to_add[base_index]);
  }
}

static void test_windows_iterators() {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  constexpr windows_size_type n_windows = 20;
  constexpr std::size_t n_indices_to_add = 30;

  auto perform_test = [&](auto get_begin, auto get_end) {
    std::vector<std::vector<WindowsMergerTraits::windows_size_type>>
        indices_to_add(n_windows);
    WindowsMergerCacheIndices wmci(70, n_windows);

    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto accessor = wmci[window_index];

      auto &current_indices_to_add = indices_to_add[window_index];
      current_indices_to_add = create_random_indices(n_indices_to_add);

      for (const auto &window_base : current_indices_to_add)
        accessor.emplace_back(window_base);
    }

    {
      const auto added_indices_iter_begin = get_begin(indices_to_add);
      const auto indices_iter_begin = get_begin(wmci);
      const auto indices_iter_end = get_end(wmci);
      auto added_indices_iter = added_indices_iter_begin;
      auto indices_iter = indices_iter_begin;

      {
        auto iter = indices_iter;
        ++iter;
        ++iter;
        assert(iter == indices_iter + 2);
      }
      assert(2 + indices_iter == indices_iter + 2);
      assert(indices_iter_end - indices_iter == n_windows);
      assert((indices_iter + 2) - 2 == indices_iter);
      assert(indices_iter != indices_iter + 1);
      assert(indices_iter <= indices_iter);
      assert(indices_iter >= indices_iter);
      assert(indices_iter <= indices_iter + 1);
      assert(indices_iter + 1 >= indices_iter);
      assert(indices_iter < indices_iter + 1);
      assert(indices_iter + 1 > indices_iter);

      for (; indices_iter < indices_iter_end;
           ++added_indices_iter, ++indices_iter) {
        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
        assert(ranges::equal(indices_line | ranges::view::reverse,
                             added_indices_line | ranges::view::reverse));
      }

      for (std::ptrdiff_t index = 0;
           index < static_cast<std::ptrdiff_t>(n_windows); ++index) {
        auto &&indices_line = indices_iter_begin[index];
        auto &&added_indices_line = added_indices_iter_begin[index];

        assert(ranges::equal(indices_line, added_indices_line));
      }
    }

    {
      auto added_indices_iter = get_begin(indices_to_add);
      auto indices_iter = get_begin(wmci);
      const auto indices_iter_end = get_end(wmci);

      for (; indices_iter < indices_iter_end;
           added_indices_iter++, indices_iter++) {
        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
      }
    }

    {
      auto added_indices_iter = get_begin(indices_to_add);
      auto indices_iter = get_begin(wmci);
      const auto indices_iter_end = get_end(wmci);

      for (; indices_iter < indices_iter_end;
           added_indices_iter += 2, indices_iter += 2) {
        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
      }
    }

    {
      auto added_indices_iter = get_end(indices_to_add);
      auto indices_iter = get_end(wmci);
      const auto indices_iter_begin = get_begin(wmci);

      while (indices_iter > indices_iter_begin) {
        --indices_iter;
        --added_indices_iter;

        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
      }
    }

    {
      auto added_indices_iter = get_end(indices_to_add);
      auto indices_iter = get_end(wmci);
      const auto indices_iter_begin = get_begin(wmci);

      while (indices_iter > indices_iter_begin) {
        indices_iter--;
        added_indices_iter--;

        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
      }
    }

    {
      auto added_indices_iter = get_end(indices_to_add);
      auto indices_iter = get_end(wmci);
      const auto indices_iter_begin = get_begin(wmci);

      while (indices_iter > indices_iter_begin) {
        indices_iter -= 2;
        added_indices_iter -= 2;

        auto &&indices_line = *indices_iter;
        auto &&added_indices_line = *added_indices_iter;
        assert(ranges::equal(indices_line, added_indices_line));
      }
    }
  };

  perform_test([](auto &&window) { return ranges::begin(window); },
               [](auto &&window) { return ranges::end(window); });

  perform_test(
      [](auto &&window) { return ranges::begin(std::as_const(window)); },
      [](auto &&window) { return ranges::end(std::as_const(window)); });

  perform_test([](auto &&window) { return std::rbegin(window); },
               [](auto &&window) { return std::rend(window); });

  perform_test([](auto &&window) { return std::rbegin(std::as_const(window)); },
               [](auto &&window) { return std::rend(std::as_const(window)); });
}

static void test_move_constructor() {
  {
    WindowsMergerCacheIndices wmci;
    WindowsMergerCacheIndices new_wmci(std::move(wmci));

    assert(new_wmci.line_capacity() == 0);
    assert(new_wmci.size() == 0);
    assert(new_wmci.capacity() == 0);
  }

  {
    WindowsMergerCacheIndices wmci(70);
    WindowsMergerCacheIndices new_wmci(std::move(wmci));

    assert(new_wmci.line_capacity() == 70);
    assert(new_wmci.size() == 0);
    assert(new_wmci.capacity() == 0);

    assert(wmci.line_capacity() == 0);
  }

  {
    WindowsMergerCacheIndices wmci(70, 100);
    WindowsMergerCacheIndices new_wmci(std::move(wmci));

    assert(new_wmci.line_capacity() == 70);
    assert(new_wmci.size() == 100);
    assert(new_wmci.capacity() == 100);

    assert(wmci.line_capacity() == 0);
    assert(wmci.size() == 0);
    assert(wmci.capacity() == 0);
  }

  {
    using windows_size_type = WindowsMergerTraits::windows_size_type;
    constexpr windows_size_type n_windows = 20;
    constexpr std::size_t n_indices_to_add = 30;

    WindowsMergerCacheIndices wmci(n_indices_to_add, n_windows);
    std::vector<std::vector<WindowsMergerTraits::windows_size_type>>
        indices_to_add(n_windows);

    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto accessor = wmci[window_index];

      auto &current_indices_to_add = indices_to_add[window_index];
      current_indices_to_add = create_random_indices(n_indices_to_add);

      for (const auto &window_base : current_indices_to_add)
        accessor.emplace_back(window_base);
    }

    WindowsMergerCacheIndices new_wmci(std::move(wmci));
    assert(new_wmci.line_capacity() >= n_indices_to_add);
    assert(new_wmci.size() == n_windows);
    assert(new_wmci.capacity() >= n_windows);

    assert(wmci.line_capacity() == 0);
    assert(wmci.size() == 0);
    assert(wmci.capacity() == 0);

    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto &&new_window = std::as_const(new_wmci)[window_index];
      const auto &indices = indices_to_add[window_index];

      for (windows_size_type index_index = 0; index_index < n_indices_to_add;
           ++index_index)
        assert(new_window[index_index] == indices[index_index]);
    }
  }
}

static void test_move_assignment() {

  auto perform_test = [](auto &&create_wmci) {
    {
      WindowsMergerCacheIndices wmci;
      WindowsMergerCacheIndices new_wmci = create_wmci();
      new_wmci = std::move(wmci);

      assert(new_wmci.line_capacity() == 0);
      assert(new_wmci.size() == 0);
      assert(new_wmci.capacity() == 0);
    }

    {
      WindowsMergerCacheIndices wmci(70);
      WindowsMergerCacheIndices new_wmci = create_wmci();
      new_wmci = std::move(wmci);

      assert(new_wmci.line_capacity() == 70);
      assert(new_wmci.size() == 0);
      assert(new_wmci.capacity() == 0);

      assert(wmci.line_capacity() == 0);
    }

    {
      WindowsMergerCacheIndices wmci(70, 100);
      WindowsMergerCacheIndices new_wmci = create_wmci();
      new_wmci = std::move(wmci);

      assert(new_wmci.line_capacity() == 70);
      assert(new_wmci.size() == 100);
      assert(new_wmci.capacity() == 100);

      assert(wmci.line_capacity() == 0);
      assert(wmci.size() == 0);
      assert(wmci.capacity() == 0);
    }

    {
      using windows_size_type = WindowsMergerTraits::windows_size_type;
      constexpr windows_size_type n_windows = 20;
      constexpr std::size_t n_indices_to_add = 30;

      WindowsMergerCacheIndices wmci(n_indices_to_add, n_windows);
      std::vector<std::vector<WindowsMergerTraits::windows_size_type>>
          indices_to_add(n_windows);

      for (WindowsMergerTraits::windows_size_type window_index = 0;
           window_index < n_windows; ++window_index) {
        auto accessor = wmci[window_index];

        auto &current_indices_to_add = indices_to_add[window_index];
        current_indices_to_add = create_random_indices(n_indices_to_add);

        for (const auto &window_base : current_indices_to_add)
          accessor.emplace_back(window_base);
      }

      WindowsMergerCacheIndices new_wmci = create_wmci();
      new_wmci = std::move(wmci);

      assert(new_wmci.line_capacity() >= n_indices_to_add);
      assert(new_wmci.size() == n_windows);
      assert(new_wmci.capacity() >= n_windows);

      assert(wmci.line_capacity() == 0);
      assert(wmci.size() == 0);
      assert(wmci.capacity() == 0);

      for (WindowsMergerTraits::windows_size_type window_index = 0;
           window_index < n_windows; ++window_index) {
        auto &&new_window = std::as_const(new_wmci)[window_index];
        const auto &indices = indices_to_add[window_index];

        for (windows_size_type index_index = 0; index_index < n_indices_to_add;
             ++index_index)
          assert(new_window[index_index] == indices[index_index]);
      }
    }
  };

  perform_test([] { return WindowsMergerCacheIndices{}; });
  perform_test([] { return WindowsMergerCacheIndices(70); });
  perform_test([] { return WindowsMergerCacheIndices(70, 100); });
  perform_test([] {
    using windows_size_type = WindowsMergerTraits::windows_size_type;
    constexpr windows_size_type n_windows = 20;
    constexpr std::size_t n_indices_to_add = 30;

    WindowsMergerCacheIndices wmci(n_indices_to_add, n_windows);
    for (WindowsMergerTraits::windows_size_type window_index = 0;
         window_index < n_windows; ++window_index) {
      auto &&window = wmci[window_index];
      auto &&indices = create_random_indices(n_indices_to_add);
      assert(window.size() == 0);

      for (auto &&index : indices)
        window.emplace_back(std::move(index));
    }

    return wmci;
  });
}

static void test_window_accessor_resize() {
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr windows_size_type n_indices = 30;

  WindowsMergerCacheIndices wmci(n_indices, 1);
  const auto accessor = wmci[0];

  auto bases = create_random_indices(n_indices);
  for (auto &&base : bases)
    accessor.emplace_back(base);
  assert(ranges::equal(bases, accessor));

  accessor.resize(60);
  assert(wmci.line_capacity() == 60);
  assert(accessor.size() == 60);
  assert(ranges::equal(bases, accessor | ranges::view::take(bases.size())));

  accessor.resize(20);
  assert(wmci.line_capacity() == 60);
  assert(accessor.size() == 20);
  assert(ranges::equal(ranges::begin(bases),
                       ranges::next(ranges::begin(bases), 20),
                       ranges::begin(accessor), ranges::end(accessor)));
}

static void test_window_accessor_clear() {
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr windows_size_type n_indices = 30;

  WindowsMergerCacheIndices wmci(n_indices, 1);
  const auto accessor = wmci[0];

  auto bases = create_random_indices(n_indices);
  for (auto &&base : bases)
    accessor.emplace_back(base);
  assert(ranges::equal(bases, accessor));
  accessor.clear();
  assert(accessor.size() == 0);
}

static void test_window_accessor_assignment() {
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr windows_size_type n_indices_to_add = 30;

  auto perform_test = [&](windows_size_type initial_capacity, auto get_rhs_wmci,
                          auto get_lhs_wmci, auto assign_accessor) {
    WindowsMergerCacheIndices wmci_lhs(initial_capacity, 1);
    {
      auto bases_to_add = create_random_indices(initial_capacity);
      auto accessor = wmci_lhs[0];
      for (const auto &window_base : bases_to_add)
        accessor.emplace_back(window_base);
    }

    WindowsMergerCacheIndices wmci_rhs(n_indices_to_add, 1);
    {
      auto bases_to_add = create_random_indices(n_indices_to_add);
      auto accessor = wmci_rhs[0];
      for (const auto &window_base : bases_to_add)
        accessor.emplace_back(window_base);
    }

    assert(not ranges::equal(wmci_lhs[0], wmci_rhs[0]));
    {
      auto accessor = get_rhs_wmci(wmci_rhs)[0];
      assign_accessor(get_lhs_wmci(wmci_lhs)[0], accessor);
    }
    // I can do this even in case of move because I know that the source state
    // is unchanged
    assert(ranges::equal(wmci_lhs[0], wmci_rhs[0]));
  };

  for (windows_size_type indices_initial_capacity :
       std::array{windows_size_type(70), windows_size_type(13)}) {

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::as_const(wmci); },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::as_const(wmci); },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::as_const(wmci); },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::as_const(wmci); },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &wmci) -> decltype(auto) { return wmci; },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) { accessor_a = accessor_b; });

    perform_test(
        indices_initial_capacity,
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &wmci) -> decltype(auto) { return std::move(wmci); },
        [](auto &&accessor_a, auto &&accessor_b) {
          accessor_a = std::move(accessor_b);
        });
  }
}

static void test_line() {
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;

  constexpr windows_size_type n_indices = 30;

  {
    WindowsMergerCacheIndices wmci(n_indices, 1);
    const auto accessor = wmci[0];

    auto bases = create_random_indices(n_indices);
    for (auto &&base : bases)
      accessor.emplace_back(base);
    assert(ranges::equal(bases, accessor));

    WindowsMergerCacheIndicesLine line(accessor);
    assert(ranges::equal(accessor, line));
    {
      WindowsMergerCacheIndices wmci_new(n_indices, 1);
      const auto accessor = wmci_new[0];
      accessor = line;
      assert(ranges::equal(accessor, line));

      assert(wmci == wmci_new);
    }
  }

  {
    WindowsMergerCacheIndices wmci(n_indices, 1);
    const auto accessor = wmci[0];

    auto bases = create_random_indices(n_indices);
    for (auto &&base : bases)
      accessor.emplace_back(base);
    assert(ranges::equal(bases, accessor));

    WindowsMergerCacheIndicesLine line(std::move(accessor));
    assert(ranges::equal(accessor, line));
    {
      WindowsMergerCacheIndices wmci_new(n_indices, 1);
      const auto accessor = wmci_new[0];
      accessor = std::move(line);

      assert(wmci == wmci_new);
    }
  }
}

int main() {
  static_assert(
      std::is_nothrow_default_constructible_v<WindowsMergerCacheIndices>);
  static_assert(std::is_move_constructible_v<WindowsMergerCacheIndices>);
  static_assert(std::is_move_assignable_v<WindowsMergerCacheIndices>);
  /*
  static_assert(std::is_copy_constructible_v<WindowsMergerCacheIndices>);
  static_assert(std::is_copy_assignable_v<WindowsMergerCacheIndices>);
  */

  test_default_construction();
  test_construction_with_clusters_and_bases();
  test_construction_with_clusters_bases_and_windows();
  test_simple_const_accessor();
  test_simple_accessor_emplace_back();
  test_simple_accessor_push_back();
  test_index_accessor_const_iterator();
  test_index_accessor_iterator();
  test_accessor_resizing_push();
  test_accessor_resizing_push_on_empty();
  test_windows_reshape_noreshape();
  test_windows_reshape_from_empty();
  test_windows_reshape_shrink();
  test_windows_iterators();
  test_move_constructor();
  test_move_assignment();
  test_window_accessor_resize();
  test_window_accessor_clear();
  test_window_accessor_assignment();
  test_line();
}
