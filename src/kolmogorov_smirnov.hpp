#pragma once

#include <cstdint>

double KScdf(int n, double x);
double KSfbar(int n, double x);

constexpr std::size_t kolmogorov_smirnov_table_size = 10'000;
double kolmogorov_smirnov_critical_value(unsigned n, double alpha);
