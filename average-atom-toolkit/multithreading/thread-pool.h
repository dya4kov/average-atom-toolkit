#pragma once

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

namespace aatk { namespace multithreading { 

#define THREAD_POOL_RES_IT typename std::result_of<F(iterator, Args...)>::type
#define THREAD_POOL_RES    typename std::result_of<F(std::size_t, Args...)>::type

class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type>;

    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    auto for_each(iterator&& a, iterator&& b, F&& f, Args&&... args)
        -> typename std::enable_if<std::is_void<THREAD_POOL_RES_IT>::value, void>::type;

    template<std::size_t ntasks_per_thread = 16, class F, class... Args>
    auto for_each(std::size_t&& a, std::size_t&& b, F&& f, Args&&... args)
        -> typename std::enable_if<std::is_void<THREAD_POOL_RES>::value, void>::type;

    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    auto for_each(iterator&& a, iterator&& b, F&& f, Args&&... args)
        -> typename std::enable_if<!std::is_void<THREAD_POOL_RES_IT>::value, 
            std::vector<THREAD_POOL_RES_IT>>::type;

    ~ThreadPool();
    std::size_t size() const { return workers.size(); }
    void resize(size_t nthreads = 0);
    void stop();
    void start();
    ThreadPool(const ThreadPool& tp) {
        started = false;
        resize(tp.size());
    }
    ThreadPool& operator=(const ThreadPool& tp) {
        started = false;
        resize(tp.size());
        return *this;
    }
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void()> > tasks;
        
    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool started;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads) : started(false) {
    resize(threads);
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type> {
    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(!started)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task](){ (*task)(); });
    }
    condition.notify_one();
    return res;
}

// change amount of workers
inline void ThreadPool::resize(size_t threads) {
    if (started) stop();
    if (threads > 0) start();
    for(size_t i = 0; i < threads; ++i)
        workers.emplace_back(
            [this]
            {
                for(;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return !this->started || !this->tasks.empty(); });
                        if(!this->started && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            }
        );
}

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto ThreadPool::for_each(iterator&& a, iterator&& b, F&& f, Args&&... args) 
    -> typename std::enable_if<std::is_void<THREAD_POOL_RES_IT>::value, void>::type {
    auto size = (b - a); if (size == 0) return;
    std::size_t num_tasks = ntasks_per_thread * workers.size();
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<void>> results(num_tasks);

    auto pos = a;
    auto end = b;
    for (auto &result : results) {
        auto to = std::min(pos + bin, end);
        result = enqueue( 
            [pos, to, f, &args...](){ 
                for (auto it = pos; it != to; ++it) 
                    f(it, std::forward<Args>(args)...);
            } 
        );
        pos = pos + bin;
    }

    for (auto& result : results) 
        result.wait();
};

template<std::size_t ntasks_per_thread, class F, class... Args>
auto ThreadPool::for_each(std::size_t&& a, std::size_t&& b, F&& f, Args&&... args) 
    -> typename std::enable_if<std::is_void<THREAD_POOL_RES>::value, void>::type {
    auto size = (b - a); if (size == 0) return;
    std::size_t num_tasks = ntasks_per_thread * workers.size();
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<void>> results(num_tasks);

    auto pos = a;
    auto end = b;
    for (auto &result : results) {
        auto to = std::min(pos + bin, end);
        result = enqueue( 
            [pos, to, f, &args...](){ 
                for (auto i = pos; i != to; ++i) 
                    f(i, std::forward<Args>(args)...);
            } 
        );
        pos = pos + bin;
    }

    for (auto& result : results) 
        result.wait();
};

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto ThreadPool::for_each(iterator&& a, iterator&& b, F&& f, Args&&... args) 
    -> typename std::enable_if<!std::is_void<THREAD_POOL_RES_IT>::value, 
            std::vector<THREAD_POOL_RES_IT>>::type {
    auto size = (b - a); 
    if (size == 0) return std::vector<THREAD_POOL_RES_IT>();
    std::size_t num_tasks = ntasks_per_thread * workers.size();
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<THREAD_POOL_RES_IT>> results_future(num_tasks);
    std::vector<THREAD_POOL_RES_IT> results(num_tasks);

    auto pos = a;
    auto end = b;
    for (auto &result : results_future) {
        auto to = std::min(pos + bin, end);
        result = enqueue( 
            [pos, to, f, &args...](){ 
                for (auto it = pos; it != to; ++it) 
                    f(it, std::forward<Args>(args)...);
            } 
        );
        pos = pos + bin;
    }

    for (auto i = 0; i < results.size(); i++) 
        results[i] = results_future[i].get();
    return results;
};

// stop joins all threads
inline void ThreadPool::stop() {
    if (started || size() > 0) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            started = false;
        }

        condition.notify_all();
        
        for(std::thread &worker: workers)
            worker.join();

        workers.clear();
    }
}

inline void ThreadPool::start() {
    started = true;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool() {
    stop();
}

static ThreadPool dummy_pool(0);

} /*multithreading*/ } /*aatk*/
