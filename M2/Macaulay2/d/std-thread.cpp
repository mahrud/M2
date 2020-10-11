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
#include <mutex>
#include <thread>

extern "C" {

int threads_print(const M2_string buffer, int n)
{
  std::cout << n << ":\t" << M2_tocharstar(buffer) << std::endl;
  return n;
}

} /* extern "C" */
