/* Copyright 2020 by Mahrud Sayrafi */

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

extern "C" {

void streams_initialize(std::stringbuf* buf, std::ostream* stream)
{
  // instead of reallocating memory, skip to the beginning
  buf->str("");
  stream->clear();
  stream->seekp(0);
  stream->rdbuf(buf);
}

void streams_flush(int n)
{
  std::lock_guard<std::mutex> stdout_guard(stdout_mutex), stderr_guard(stderr_mutex);
  std::cout << n << " out:\n" << outbuf.str() << std::flush;
  std::cerr << n << " err:\n" << errbuf.str() << std::flush;
}

} /* extern "C" */
