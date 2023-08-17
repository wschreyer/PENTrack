#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>

class ThreadPool {
public:
    void Start(int numthread);
    void QueueJob(const std::function<void()>& job);
    bool Busy();
    void Stop();

private:
    static void ThreadLoop(ThreadPool* pool);

    std::vector<std::thread> threads;
    std::queue<std::function<void()>> jobs;
    std::mutex queue_mutex;
    std::condition_variable mutex_condition;
    bool should_terminate = false;
};

#endif // THREADPOOL_H
