#include "math.hpp"

#include <algorithm>
#include <array>

static void test1() {
  constexpr std::array<int, 5> data1{5, 30, 100, 1000, 10000};
  constexpr std::array<int, 5> data2{125, 12345, 7, 35, 980};
  constexpr std::array<int, 5> expectedResult{2, 3, 0, 4, 1};

  auto result = getBestMatchingIndices(
      std::begin(data1), std::end(data1), std::begin(data2), std::end(data2),
      [](auto a, auto b) { return std::abs(a - b); }, std::greater<int>());
  if (not std::equal(std::begin(expectedResult), std::end(expectedResult),
                     std::begin(result), std::end(result)))
    std::exit(-1);
}

static void test2() {
  const std::vector data1{5, 30, 100, 1000, 10000};
  const std::vector data2{125, 12345, 7, 35, 980};
  constexpr std::array<int, 5> expectedResult{2, 3, 0, 4, 1};

  auto result = getBestMatchingIndices(
      std::begin(data1), std::end(data1), std::begin(data2), std::end(data2),
      [](auto a, auto b) { return std::abs(a - b); }, std::greater<int>());
  if (not std::equal(std::begin(expectedResult), std::end(expectedResult),
                     std::begin(result), std::end(result)))
    std::exit(-1);
}

static void test_asymmetric() {
  constexpr std::array<int, 3> data1{5, 100, 10000};
  constexpr std::array<int, 5> data2{125, 35, 7, 980, 12345};
  constexpr std::array<int, 3> expectedResult{2, 0, 4};

  auto result = getBestMatchingIndices(
      std::begin(data1), std::end(data1), std::begin(data2), std::end(data2),
      [](auto a, auto b) { return std::abs(a - b); }, std::greater<int>());
  if (not std::equal(std::begin(expectedResult), std::end(expectedResult),
                     std::begin(result), std::end(result)))
    std::exit(-1);
}

static void test_less_matching() {
  constexpr std::array<int, 5> data1{5, 30, 100, 1000, 10000};
  constexpr std::array<int, 2> data2{32, 2};
  constexpr std::array<int, 2> expectedResult{1, 0};

  auto result = getBestMatchingIndices(
      std::begin(data1), std::end(data1), std::begin(data2), std::end(data2),
      [](auto a, auto b) { return std::abs(a - b); }, std::greater<int>());
  if (not std::equal(std::begin(expectedResult), std::end(expectedResult),
                     std::begin(result), std::end(result)))
    std::exit(-1);
}

int main() {
  test1();
  test2();
  test_asymmetric();
  test_less_matching();
}
