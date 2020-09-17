#include "kolmogorov_smirnov.hpp"
#include "kolmogorov_smirnov_general.hpp"
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <vector>

constexpr int nExact = 500;
constexpr int nKolmo = 100'000;
constexpr int mFact = 30;

static double DurbinMatrix(int n, double d);

static constexpr double LnFactorial[mFact + 1] = {0.,
                                                  0.,
                                                  0.6931471805599453,
                                                  1.791759469228055,
                                                  3.178053830347946,
                                                  4.787491742782046,
                                                  6.579251212010101,
                                                  8.525161361065415,
                                                  10.60460290274525,
                                                  12.80182748008147,
                                                  15.10441257307552,
                                                  17.50230784587389,
                                                  19.98721449566188,
                                                  22.55216385312342,
                                                  25.19122118273868,
                                                  27.89927138384088,
                                                  30.67186010608066,
                                                  33.50507345013688,
                                                  36.39544520803305,
                                                  39.33988418719949,
                                                  42.33561646075348,
                                                  45.3801388984769,
                                                  48.47118135183522,
                                                  51.60667556776437,
                                                  54.7847293981123,
                                                  58.00360522298051,
                                                  61.26170176100199,
                                                  64.55753862700632,
                                                  67.88974313718154,
                                                  71.257038967168,
                                                  74.65823634883016};

static double
getLogFactorial(int n) {
  if (n <= mFact) {
    return LnFactorial[n];

  } else {
    double x = static_cast<double>(n + 1);
    double y = 1.0 / (x * x);
    double z = ((-(5.95238095238E-4 * y) + 7.936500793651E-4) * y -
                2.7777777777778E-3) *
                   y +
               8.3333333333333E-2;
    z = ((x - 0.5) * log(x) - x) + 9.1893853320467E-1 + z / x;
    return z;
  }
}

static double
rapfac(int n) {
  int i;
  double res = 1.0 / n;
  for (i = 2; i <= n; i++) {
    res *= static_cast<double>(i) / n;
  }
  return res;
}

static double
KSPlusbarAsymp(int n, double x) {
  double t = (6.0 * n * x + 1);
  double z = t * t / (18.0 * n);
  double v = 1.0 - (2.0 * z * z - 4.0 * z - 1.0) / (18.0 * n);
  if (v <= 0.0)
    return 0.0;
  v = v * exp(-z);
  if (v >= 1.0)
    return 1.0;
  return v;
}

static double
KSPlusbarUpper(int n, double x) {
  const double EPSILON = 1.0E-12;
  double q;
  double Sum = 0.0;
  double term;
  double t;
  double LogCom;
  double LOGJMAX;
  int j;
  int jdiv;
  int jmax = static_cast<int>(n * (1.0 - x));

  if (n > 200000)
    return KSPlusbarAsymp(n, x);

  if ((1.0 - x - static_cast<double>(jmax) / n) <= 0.0)
    jmax--;

  if (n > 3000)
    jdiv = 2;
  else
    jdiv = 3;

  j = jmax / jdiv + 1;
  LogCom = getLogFactorial(n) - getLogFactorial(j) - getLogFactorial(n - j);
  LOGJMAX = LogCom;

  while (j <= jmax) {
    q = static_cast<double>(j) / n + x;
    term = LogCom + (j - 1) * log(q) + (n - j) * log1p(-q);
    t = exp(term);
    Sum += t;
    LogCom += log(static_cast<double>(n - j) / (j + 1));
    if (t <= Sum * EPSILON)
      break;
    j++;
  }

  j = jmax / jdiv;
  LogCom = LOGJMAX + log(static_cast<double>(j + 1) / (n - j));

  while (j > 0) {
    q = static_cast<double>(j) / n + x;
    term = LogCom + (j - 1) * log(q) + (n - j) * log1p(-q);
    t = exp(term);
    Sum += t;
    LogCom += log(static_cast<double>(j) / (n - j + 1));
    if (t <= Sum * EPSILON)
      break;
    j--;
  }

  Sum *= x;
  Sum += exp(n * log1p(-x));
  return Sum;
}

static double
Pelz(int n, double x) {
  constexpr int JMAX = 20;
  constexpr double EPS = 1.0e-10;
  constexpr double C = 2.506628274631001;   /* sqrt(2*Pi) */
  constexpr double C2 = 1.2533141373155001; /* sqrt(Pi/2) */
  constexpr double PI2 = M_PI * M_PI;
  constexpr double PI4 = PI2 * PI2;

  const double RACN = sqrt(static_cast<double>(n));
  const double z = RACN * x;
  const double z2 = z * z;
  const double z4 = z2 * z2;
  const double z6 = z4 * z2;
  const double w = PI2 / (2.0 * z * z);
  double ti, term, tom;
  double sum;
  int j;

  term = 1;
  j = 0;
  sum = 0;
  while (j <= JMAX && term > EPS * sum) {
    ti = j + 0.5;
    term = exp(-ti * ti * w);
    sum += term;
    j++;
  }
  sum *= C / z;

  term = 1;
  tom = 0;
  j = 0;
  while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
    ti = j + 0.5;
    term = (PI2 * ti * ti - z2) * exp(-ti * ti * w);
    tom += term;
    j++;
  }
  sum += tom * C2 / (RACN * 3.0 * z4);

  term = 1;
  tom = 0;
  j = 0;
  while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
    ti = j + 0.5;
    term = 6 * z6 + 2 * z4 + PI2 * (2 * z4 - 5 * z2) * ti * ti +
           PI4 * (1 - 2 * z2) * ti * ti * ti * ti;
    term *= exp(-ti * ti * w);
    tom += term;
    j++;
  }
  sum += tom * C2 / (n * 36.0 * z * z6);

  term = 1;
  tom = 0;
  j = 1;
  while (j <= JMAX && term > EPS * tom) {
    ti = j;
    term = PI2 * ti * ti * exp(-ti * ti * w);
    tom += term;
    j++;
  }
  sum -= tom * C2 / (n * 18.0 * z * z2);

  term = 1;
  tom = 0;
  j = 0;
  while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
    ti = j + 0.5;
    ti = ti * ti;
    term = -30 * z6 - 90 * z6 * z2 + PI2 * (135 * z4 - 96 * z6) * ti +
           PI4 * (212 * z4 - 60 * z2) * ti * ti +
           PI2 * PI4 * ti * ti * ti * (5 - 30 * z2);
    term *= exp(-ti * w);
    tom += term;
    j++;
  }
  sum += tom * C2 / (RACN * n * 3240.0 * z4 * z6);

  term = 1;
  tom = 0;
  j = 1;
  while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
    ti = j * j;
    term = (3 * PI2 * ti * z2 - PI4 * ti * ti) * exp(-ti * w);
    tom += term;
    j++;
  }
  sum += tom * C2 / (RACN * n * 108.0 * z6);

  return sum;
}

static void
CalcFloorCeil(int n, double t, std::vector<double>& A,
              std::vector<double>& Atflo, std::vector<double>& Atcei) {
  int i;
  auto ell = static_cast<int>(t);
  double z = t - ell;
  double w = ceil(t) - t;

  if (z > 0.5) {
    for (i = 2; i <= 2 * n + 2; i += 2)
      Atflo[i] = i / 2 - 2 - ell;
    for (i = 1; i <= 2 * n + 2; i += 2)
      Atflo[i] = i / 2 - 1 - ell;

    for (i = 2; i <= 2 * n + 2; i += 2)
      Atcei[i] = i / 2 + ell;
    for (i = 1; i <= 2 * n + 2; i += 2)
      Atcei[i] = i / 2 + 1 + ell;

  } else if (z > 0.0) {
    for (i = 1; i <= 2 * n + 2; i++)
      Atflo[i] = i / 2 - 1 - ell;

    for (i = 2; i <= 2 * n + 2; i++)
      Atcei[i] = i / 2 + ell;
    Atcei[1] = 1 + ell;

  } else { /* z == 0 */
    for (i = 2; i <= 2 * n + 2; i += 2)
      Atflo[i] = i / 2 - 1 - ell;
    for (i = 1; i <= 2 * n + 2; i += 2)
      Atflo[i] = i / 2 - ell;

    for (i = 2; i <= 2 * n + 2; i += 2)
      Atcei[i] = i / 2 - 1 + ell;
    for (i = 1; i <= 2 * n + 2; i += 2)
      Atcei[i] = i / 2 + ell;
  }

  if (w < z)
    z = w;
  A[0] = A[1] = 0;
  A[2] = z;
  A[3] = 1 - A[2];
  for (i = 4; i <= 2 * n + 1; i++)
    A[i] = A[i - 2] + 1;
  A[2 * n + 2] = n;
}

static double
Pomeranz(int n, double x) {
  const double EPS = 1.0e-15;
  const int ENO = 350;
  const double RENO = ldexp(1.0, ENO);
  int coreno;
  const double t = n * x;
  double w, sum, minsum;
  int i, j, k, s;
  int r1, r2;
  int jlow, jup, klow, kup, kup0;

  std::vector<double> A(2 * n + 3);
  std::vector<double> Atflo(2 * n + 3);
  std::vector<double> Atcei(2 * n + 3);
  arma::mat V(2, n + 2);
  arma::mat H(4, n + 2);

  CalcFloorCeil(n, t, A, Atflo, Atcei);

  for (j = 1; j <= n + 1; j++)
    V(0, j) = 0;
  for (j = 2; j <= n + 1; j++)
    V(1, j) = 0;
  V(1, 1) = RENO;
  coreno = 1;

  H(0, 0) = 1;
  w = 2.0 * A[2] / n;
  for (j = 1; j <= n + 1; j++)
    H(0, j) = w * H(0, j - 1) / j;

  H(1, 0) = 1;
  w = (1.0 - 2.0 * A[2]) / n;
  for (j = 1; j <= n + 1; j++)
    H(1, j) = w * H(1, j - 1) / j;

  H(2, 0) = 1;
  w = A[2] / n;
  for (j = 1; j <= n + 1; j++)
    H(2, j) = w * H(2, j - 1) / j;

  H(3, 0) = 1;
  for (j = 1; j <= n + 1; j++)
    H(3, j) = 0;

  r1 = 0;
  r2 = 1;
  for (i = 2; i <= 2 * n + 2; i++) {
    jlow = 2 + static_cast<int>(Atflo[i]);
    if (jlow < 1)
      jlow = 1;
    jup = static_cast<int>(Atcei[i]);
    if (jup > n + 1)
      jup = n + 1;

    klow = 2 + static_cast<int>(Atflo[i - 1]);
    if (klow < 1)
      klow = 1;
    kup0 = static_cast<int>(Atcei[i - 1]);

    w = (A[i] - A[i - 1]) / n;
    s = -1;
    for (j = 0; j < 4; j++) {
      if (fabs(w - H(j, 1)) <= EPS) {
        s = j;
        break;
      }
    }

    minsum = RENO;
    r1 = (r1 + 1) & 1;
    r2 = (r2 + 1) & 1;

    for (j = jlow; j <= jup; j++) {
      kup = kup0;
      if (kup > j)
        kup = j;
      sum = 0;
      for (k = kup; k >= klow; k--)
        sum += V(r1, k) * H(s, j - k);
      V(r2, j) = sum;
      if (sum < minsum)
        minsum = sum;
    }

    if (minsum < 1.0e-280) {
      for (j = jlow; j <= jup; j++)
        V(r2, j) *= RENO;
      coreno++;
    }
  }

  sum = V(r2, n + 1);
  w = getLogFactorial(n) - coreno * ENO * std::log(2) + log(sum);
  if (w >= 0.)
    return 1.;
  return exp(w);
}

static double
cdfSpecial(int n, double x) {
  if ((n * x * x >= 18.0) || (x >= 1.0))
    return 1.0;

  if (x <= 0.5 / n)
    return 0.0;

  if (n == 1)
    return 2.0 * x - 1.0;

  if (x <= 1.0 / n) {
    double t = 2.0 * x * n - 1.0;
    double w;
    if (n <= nExact) {
      w = rapfac(n);
      return w * pow(t, static_cast<double>(n));
    }
    w = getLogFactorial(n) + n * log(t / n);
    return exp(w);
  }

  if (x >= 1.0 - 1.0 / n) {
    return 1.0 - 2.0 * pow(1.0 - x, static_cast<double>(n));
  }

  return -1.0;
}

double
KScdf(int n, double x) {
  const double w = n * x * x;
  double u = cdfSpecial(n, x);
  if (u >= 0.0)
    return u;

  if (n <= nExact) {
    if (w < 0.754693)
      return DurbinMatrix(n, x);
    if (w < 4.0)
      return Pomeranz(n, x);
    return 1.0 - KSfbar(n, x);
  }

  if ((w * x * n <= 7.0) && (n <= nKolmo))
    return DurbinMatrix(n, x);

  return Pelz(n, x);
}

static double
fbarSpecial(int n, double x) {
  const double w = n * x * x;

  if ((w >= 370.0) || (x >= 1.0))
    return 0.0;
  if ((w <= 0.0274) || (x <= 0.5 / n))
    return 1.0;
  if (n == 1)
    return 2.0 - 2.0 * x;

  if (x <= 1.0 / n) {
    double z;
    double t = 2.0 * x * n - 1.0;
    if (n <= nExact) {
      z = rapfac(n);
      return 1.0 - z * pow(t, static_cast<double>(n));
    }
    z = getLogFactorial(n) + n * log(t / n);
    return 1.0 - exp(z);
  }

  if (x >= 1.0 - 1.0 / n) {
    return 2.0 * pow(1.0 - x, static_cast<double>(n));
  }
  return -1.0;
}

double
KSfbar(int n, double x) {
  const double w = n * x * x;
  double v = fbarSpecial(n, x);
  if (v >= 0.0)
    return v;

  if (n <= nExact) {
    if (w < 4.0)
      return 1.0 - KScdf(n, x);
    else
      return 2.0 * KSPlusbarUpper(n, x);
  }

  if (w >= 2.65)
    return 2.0 * KSPlusbarUpper(n, x);

  return 1.0 - KScdf(n, x);
}

constexpr double norm = 1e140;
constexpr double inorm = 1. / 1e140;
constexpr double logNorm = 140.;

static void mMultiply(std::vector<double> const& A,
                      std::vector<double> const& B, std::vector<double>& C,
                      int m);

static void mPower(std::vector<double> const& A, int eA, std::vector<double>& V,
                   int* eV, int m, int n);

static double
DurbinMatrix(int n, double d) {
  int k, m, i, j, g, eH, eQ;
  double h, s;
  k = static_cast<int>(n * d) + 1;
  m = 2 * k - 1;
  h = k - n * d;
  std::vector<double> H(m * m);
  std::vector<double> Q(m * m);

  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i - j + 1 < 0)
        H[i * m + j] = 0;
      else
        H[i * m + j] = 1;
  for (i = 0; i < m; i++) {
    H[i * m] -= pow(h, static_cast<double>(i + 1));
    H[(m - 1) * m + i] -= pow(h, static_cast<double>(m - i));
  }
  H[(m - 1) * m] +=
      (2 * h - 1 > 0 ? pow(2 * h - 1, static_cast<double>(m)) : 0);
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i - j + 1 > 0)
        for (g = 1; g <= i - j + 1; g++)
          H[i * m + j] /= g;
  eH = 0;
  mPower(H, eH, Q, &eQ, m, n);
  s = Q[(k - 1) * m + k - 1];

  for (i = 1; i <= n; i++) {
    s = s * static_cast<double>(i) / n;
    if (s < inorm) {
      s *= norm;
      eQ -= logNorm;
    }
  }
  s *= pow(10., static_cast<double>(eQ));
  return s;
}

static void
mMultiply(std::vector<double> const& A, std::vector<double> const& B,
          std::vector<double>& C, int m) {
  int i, j, k;
  double s;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      s = 0.;
      for (k = 0; k < m; k++)
        s += A[i * m + k] * B[k * m + j];
      C[i * m + j] = s;
    }
}

static void
renormalize(std::vector<double>& V, int m, int* p) {
  int i;
  for (i = 0; i < m * m; i++)
    V[i] *= inorm;
  *p += logNorm;
}

static void
mPower(std::vector<double> const& A, int eA, std::vector<double>& V, int* eV,
       int m, int n) {
  int eB, i;
  if (n == 1) {
    for (i = 0; i < m * m; i++)
      V[i] = A[i];
    *eV = eA;
    return;
  }
  mPower(A, eA, V, eV, m, n / 2);
  std::vector<double> B(m * m);
  mMultiply(V, V, B, m);
  eB = 2 * (*eV);
  if (B[(m / 2) * m + (m / 2)] > norm)
    renormalize(B, m, &eB);

  if (n % 2 == 0) {
    for (i = 0; i < m * m; i++)
      V[i] = B[i];
    *eV = eB;
  } else {
    mMultiply(A, B, V, m);
    *eV = eA + eB;
  }

  if (V[(m / 2) * m + (m / 2)] > norm)
    renormalize(V, m, eV);
}

double
kolmogorov_smirnov_critical_value(unsigned n, double alpha) {
  thread_local std::map<int, std::vector<double>> cached_probabilities;

  if (alpha <= 0. or alpha > 1.)
    return std::numeric_limits<double>::quiet_NaN();

  if (n > 50) {
    auto const value_iter =
        std::upper_bound(std::begin(kolmogorov_smirnov_general_table),
                         std::end(kolmogorov_smirnov_general_table), alpha,
                         std::greater<double>());
    return static_cast<double>(std::distance(
               std::begin(kolmogorov_smirnov_general_table), value_iter)) /
           static_cast<double>(kolmogorov_smirnov_general_table.size()) *
           std::sqrt(kolmogorog_smirnov_general_n) /
           std::sqrt(static_cast<double>(n));
  }

  auto const& table = [&]() -> decltype(auto) {
    if (auto table_iter = cached_probabilities.find(n);
        table_iter == std::end(cached_probabilities)) {
      std::vector<double> fbars(kolmogorov_smirnov_table_size);
      for (std::size_t index = 0; index < kolmogorov_smirnov_table_size;
           ++index) {
        double x = static_cast<double>(index + 1) /
                   static_cast<double>(kolmogorov_smirnov_table_size);
        fbars[index] = KSfbar(n, x);
      }
      return cached_probabilities.emplace(n, std::move(fbars)).first->second;
    } else {
      return table_iter->second;
    }
  }();

  auto const value_iter = std::upper_bound(std::begin(table), std::end(table),
                                           alpha, std::greater<double>());
  return static_cast<double>(std::distance(std::begin(table), value_iter)) /
         static_cast<double>(kolmogorov_smirnov_table_size);
}
