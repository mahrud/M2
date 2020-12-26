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

void initialize(std::stringbuf* buf, std::ostream* stream)
{
  // instead of reallocating memory, skip to the beginning
  buf->str("");
  stream->clear();
  stream->seekp(0);
  stream->rdbuf(buf);
}

extern "C" {

// TODO: ability to use stdout or stderr instead of n
void streams_Sinitialize(int n, bool live)
{
  switch (n)
    {
      case 1:
        initialize(&outbuf, &outstream);
	break;
      case 2:
        if (live)
          initialize(&errbuf, &errstream);
        else
          errstream.rdbuf(&outbuf);
	break;
    }
}

void streams_Sflush()
{
  std::lock_guard<std::mutex> stdout_guard(stdout_mutex), stderr_guard(stderr_mutex);
  std::cout << " out:\n" << outbuf.str() << std::flush;
  std::cerr << " err:\n" << errbuf.str() << std::flush;
  // TODO: should we empty the buffer here?
  outbuf.str("");
  errbuf.str("");
}

M2_string streams_StoString(int n)
{
  M2_string ret;
  switch (n)
    {
      case 1:
	ret = M2_tostring(outbuf.str().c_str());
	outbuf.str("");
	break;
      case 2:
	ret = M2_tostring(errbuf.str().c_str());
	errbuf.str("");
	break;
    }
  return ret;
}

} /* extern "C" */
