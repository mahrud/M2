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
#include <chrono>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>

extern thread_local std::stringbuf outbuf, errbuf;
extern thread_local std::ostream outstream, errstream;

extern "C" {

void streams_Sinitialize(int, bool);
size_t streams_Sflush(int);

void threads_Tinitialize(bool live)
{
  streams_Sinitialize(1, live); // stdout
  streams_Sinitialize(2, live); // stderr
}

void threads_Tcompute(const M2_string str, int n)
{
  errstream << "\t";
  for(int i = 0; i < 100; i++)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      errstream << M2_tocharstar(str);
    }
  errstream << std::endl;
  outstream << "\tstr = '" << errbuf.str() << "'" << std::endl;
  streams_Sflush(1);
}

} /* extern "C" */
