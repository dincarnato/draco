#include "parallel/blocking_queue.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <thread>

static constexpr std::size_t nLoops = 10'000;
static constexpr std::size_t minVecSize = 10;
static constexpr std::size_t maxVecSize = 100;
static constexpr std::size_t nConsumers = 10;

static std::atomic<size_t> mapReads = 0;
static std::atomic_bool mapWriting = false;
static std::map<std::vector<double>, std::atomic<std::ptrdiff_t>>
    availableElements;

static void queue_filler(parallel::blocking_queue<std::vector<double>> &queue) {
  std::mt19937 randomGen(std::random_device{}());
  std::uniform_real_distribution<> randomDist;
  std::uniform_int_distribution<std::size_t> randomVecSize(minVecSize,
                                                           maxVecSize);

  for (std::size_t loopIndex = 0; loopIndex < nLoops; ++loopIndex) {
    std::vector<double> vec(randomVecSize(randomGen));
    std::generate(std::begin(vec), std::end(vec),
                  [&] { return randomDist(randomGen); });
    {
      mapWriting.store(true, std::memory_order_release);
      while (mapReads.load(std::memory_order_acquire) != 0)
        ;

      if (auto mappedElement = availableElements.find(vec);
          mappedElement == std::end(availableElements))
        availableElements.emplace(vec, 1);
      else
        mappedElement->second.fetch_add(1, std::memory_order_acq_rel);

      mapWriting.store(false, std::memory_order_release);
    }

    queue.push(std::move(vec));
  }

  queue.finish();
}

std::atomic<std::size_t> allVecsPopped = 0;

static void
queue_consumer(parallel::blocking_queue<std::vector<double>> &queue) {
  std::mt19937 randomGen(std::random_device{}());

  std::size_t vecsPopped = 0;
  for (;;) {
    auto popped = queue.pop();
    if (not popped) {
      assert(queue.finished());
      break;
    }

    ++vecsPopped;
    const std::vector<double> &vec = *popped;
    {
      for (;;) {
        if (mapWriting.load(std::memory_order_acquire))
          continue;

        mapReads.fetch_add(1, std::memory_order_acq_rel);
        if (not mapWriting.load(std::memory_order_acquire))
          break;

        mapReads.fetch_sub(1, std::memory_order_acq_rel);
      }

      auto mappedElement = availableElements.find(vec);
      assert(mappedElement != std::end(availableElements));

      mapReads.fetch_sub(1, std::memory_order_acq_rel);

      assert(mappedElement->second.fetch_sub(1, std::memory_order_acq_rel) >=
             0);
    }

    double result = 10.;
    for (unsigned loop = 0; loop < 50; ++loop) {
      result =
          std::accumulate(std::begin(vec), std::end(vec), result,
                          [](double accumulator, double value) {
                            double newValue =
                                accumulator * std::log(value > 0 ? value : 1.);
                            return std::log(newValue > 0 ? newValue : 10.);
                          });
    }
    (void)result;
  }

  allVecsPopped += vecsPopped;
}

int main() {
  using namespace std::chrono_literals;
  parallel::blocking_queue<std::vector<double>> queue(
      static_cast<std::size_t>(1.5 * nConsumers));

  std::thread filler([&] { queue_filler(queue); });

  std::this_thread::sleep_for(10ms);
  std::vector<std::thread> consumers;
  consumers.reserve(nConsumers);

  for (std::size_t index = 0; index < nConsumers; ++index)
    consumers.emplace_back([&] { queue_consumer(queue); });

  filler.join();

  for (auto &consumer : consumers)
    consumer.join();

  assert(allVecsPopped = nLoops);
  assert(std::all_of(std::begin(availableElements), std::end(availableElements),
                     [](const auto &vecAndCount) {
                       return vecAndCount.second.load(
                                  std::memory_order_acquire) == 0;
                     }));
}
