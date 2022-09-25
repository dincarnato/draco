#pragma once

#include <string>
#include <vector>

enum class BaseSecondaryStructure {
  single_strand,
  double_strand_open,
  double_strand_close
};

struct Helix;

class RnaSecondaryStructure : public std::vector<BaseSecondaryStructure> {
public:
  using base_type = std::vector<BaseSecondaryStructure>;

  RnaSecondaryStructure() = default;
  explicit RnaSecondaryStructure(const std::string &rawStructure);
  explicit RnaSecondaryStructure(const char *const rawStructure);
  RnaSecondaryStructure(std::size_t length);

  RnaSecondaryStructure(const RnaSecondaryStructure &) noexcept;
  RnaSecondaryStructure(RnaSecondaryStructure &&) noexcept;
  RnaSecondaryStructure &operator=(const RnaSecondaryStructure &) noexcept;
  RnaSecondaryStructure &operator=(RnaSecondaryStructure &&) noexcept;

  RnaSecondaryStructure &operator+=(const RnaSecondaryStructure &);
  RnaSecondaryStructure &operator+=(RnaSecondaryStructure &&);
  RnaSecondaryStructure operator+(const RnaSecondaryStructure &) const;
  RnaSecondaryStructure operator+(RnaSecondaryStructure &&) const;

  unsigned distance(const RnaSecondaryStructure &other) const;
  const std::string &raw() const;
  std::string str() const;
  void regenRaw();

  bool isBalanced() const noexcept;
  std::vector<Helix> getHelices(bool keepLonelyPairs = true) const;

private:
  std::string rawStructure;
};
