#pragma once
#include <condition_variable>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

// Runs work(item, result) on worker threads while the calling thread produces items in order and
// consumes their results in the same order. produce(item) returns false when no items remain.
// With one thread everything runs inline. An exception from work is rethrown on the calling thread
// when its result is due, after all earlier results were consumed; an exception from produce is
// rethrown after the results already produced were consumed. Workers are joined before returning.
template <class Item, class Result, class Produce, class Work, class Consume>
void runOrdered(int threads, Produce produce, Work work, Consume consume) {
	if (threads <= 1) {
		while (true) {
			Item item;
			if (!produce(item)) { return; }
			Result result;
			work(item, result);
			consume(result);
		}
	}

	struct Task {
		Item item;
		Result result;
		std::exception_ptr error;
		bool done = false;
	};
	std::mutex mutex;
	std::condition_variable workAvailable;
	std::condition_variable resultReady;
	std::deque<Task*> pending;
	bool stopping = false;
	std::deque<std::unique_ptr<Task> > inFlight; // in production order; outlives the workers

	auto worker = [&]() {
		std::unique_lock<std::mutex> lock(mutex);
		while (true) {
			workAvailable.wait(lock, [&] { return stopping || !pending.empty(); });
			if (pending.empty()) { return; }
			Task* task = pending.front();
			pending.pop_front();
			lock.unlock();
			try {
				work(task->item, task->result);
			} catch (...) {
				task->error = std::current_exception();
			}
			lock.lock();
			task->done = true;
			resultReady.notify_all();
		}
	};

	struct Joiner {
		std::vector<std::thread> pool;
		std::mutex& mutex;
		std::condition_variable& workAvailable;
		std::deque<Task*>& pending;
		bool& stopping;
		~Joiner() {
			{
				std::lock_guard<std::mutex> lock(mutex);
				stopping = true;
				pending.clear();
			}
			workAvailable.notify_all();
			for (auto& thread : pool) { thread.join(); }
		}
	} joiner{{}, mutex, workAvailable, pending, stopping};
	for (int i = 0; i < threads; ++i) { joiner.pool.emplace_back(worker); }

	// Enough queued items to keep every worker busy; more would only hold memory.
	const size_t window = static_cast<size_t>(threads) + 2;
	bool more = true;
	std::exception_ptr produceError;
	while (true) {
		while (more && inFlight.size() < window) {
			auto task = std::make_unique<Task>();
			try {
				more = produce(task->item);
			} catch (...) {
				produceError = std::current_exception();
				more = false;
			}
			if (!more) { break; }
			{
				std::lock_guard<std::mutex> lock(mutex);
				pending.push_back(task.get());
			}
			inFlight.push_back(std::move(task));
			workAvailable.notify_one();
		}
		if (inFlight.empty()) { break; }
		Task* next = inFlight.front().get();
		{
			std::unique_lock<std::mutex> lock(mutex);
			resultReady.wait(lock, [next] { return next->done; });
		}
		if (next->error) { std::rethrow_exception(next->error); }
		consume(next->result);
		inFlight.pop_front();
	}
	if (produceError) { std::rethrow_exception(produceError); }
}
