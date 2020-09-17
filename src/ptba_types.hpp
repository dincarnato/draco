#pragma once

#include <vector>

struct PerturbedEigengap : std::vector<double> {
  using base_type = std::vector<double>;
  using base_type::base_type;
  using base_type::operator=;
};

struct PerturbedEigengaps : std::vector<PerturbedEigengap> {
  using base_type = std::vector<PerturbedEigengap>;
  using base_type::base_type;
  using base_type::operator=;
};
