#include "threadpool.h"
#include <iostream>

void ThreadPool::Start(int numthreads) {
  const int maxnumthreads = std::thread::hardware_concurrency();
  if (numthreads <= 0) numthreads = 1;
  if (numthreads >= maxnumthreads) numthreads = maxnumthreads-1;
  std::cout << "Number of threads = " << numthreads << std::endl;
    
  threads.resize(numthreads);
  for (int i = 0; i < numthreads; i++) {
    threads.at(i) = std::thread(&ThreadPool::ThreadLoop, this);
  }
}

void ThreadPool::QueueJob(const std::function<void()>& job) {
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    jobs.push(job);
  }
  mutex_condition.notify_one();
}

bool ThreadPool::Busy() {
  bool poolbusy;
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    poolbusy = !jobs.empty();
  }
  return poolbusy;
}

void ThreadPool::Stop() {
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    should_terminate = true;
  }
  mutex_condition.notify_all();
  for (std::thread& active_thread : threads) {
    active_thread.join();
  }
  threads.clear();
}

void ThreadPool::ThreadLoop(ThreadPool* pool) {
  while (true) {
    std::function<void()> job;
    {
      std::unique_lock<std::mutex> lock(pool->queue_mutex);
      pool->mutex_condition.wait(lock, [pool] {
	return pool->should_terminate || !pool->jobs.empty();
      });
      if (pool->should_terminate) {
	return;
      }
      job = pool->jobs.front();
      pool->jobs.pop();
    }
    job();
  }
}


