/* Copyright 2020 by Mahrud Sayrafi */

/*
  https://en.cppreference.com/w/cpp/thread/thread
  https://en.cppreference.com/w/cpp/atomic/atomic
  https://en.cppreference.com/w/cpp/thread/mutex
  https://www.boost.org/doc/libs/1_74_0/doc/html/thread.html#thread.overview
  https://software.intel.com/content/www/us/en/develop/articles/choosing-the-right-threading-framework.html
  https://latedev.wordpress.com/2013/02/24/investigating-c11-threads/
*/

#include "strings-exports.h"

#include <atomic>
#include <fstream>
#include <iostream>
#include <sstream>
#include <mutex>
#include <chrono>
#include <thread>

thread_local std::stringbuf outbuf, errbuf;
thread_local std::ostream outstream(&outbuf), errstream(&errbuf);
std::mutex stdout_mutex, stderr_mutex;

void threads_initialize_stream(std::stringbuf* buf, std::ostream* stream)
{
  // instead of reallocating memory, skip to the beginning
  buf->str("");
  stream->clear();
  stream->seekp(0);
  stream->rdbuf(buf);
}

extern "C" {

void threads_initialize(bool live)
{
  threads_initialize_stream(&outbuf, &outstream);
  if (live)
    threads_initialize_stream(&errbuf, &errstream);
  else
    errstream.rdbuf(&outbuf);
}

void threads_flush(int n)
{
  std::lock_guard<std::mutex> stdout_guard(stdout_mutex), stderr_guard(stderr_mutex);
  std::cout << n << " out:\n" << outbuf.str() << std::flush;
  std::cerr << n << " err:\n" << errbuf.str() << std::flush;
}

void threads_compute(const M2_string str, int n)
{
  errstream << "\t";
  for(int i = 0; i < 100; i++)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      errstream << M2_tocharstar(str);
    }
  errstream << std::endl;
  outstream << "\tstr = '" << errbuf.str() << "'" << std::endl;
  threads_flush(n);
}

} /* extern "C" */
