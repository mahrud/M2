/* Copyright 2020 by Mahrud Sayrafi */

#include "strings-exports.h"

#include <assert.h>

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

size_t streams_Sflush(int fd)
{
  if (fd == 1)
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << outbuf.str() << std::flush;
      outbuf.str("");
      return std::cout.tellp();
    }
  else if (fd == 2)
    {
      std::lock_guard<std::mutex> stderr_guard(stderr_mutex);
      std::cerr << errbuf.str() << std::flush;
      errbuf.str("");
      return std::cerr.tellp();
    }
  else
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << " SHOULD HAVE BEEN TO FILE" << std::endl;
    }
  return -1;
}

size_t streams_Swrite(int fd, std::stringbuf* buf, size_t offset, size_t size)
{
  // std::ios::sync_with_stdio(false); // TODO: try this
  if (size == 0) return 0;
  assert(offset + size <= buf->str().size());
  std::string str(buf->str().substr(offset, size));
  if (fd == 1)
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << str << std::flush;
    }
  else if (fd == 2)
    {
      std::lock_guard<std::mutex> stderr_guard(stderr_mutex);
      std::cerr << str << std::flush;
    }
  else
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << " SHOULD HAVE BEEN TO FILE:\n" << str << std::flush;
    }
  buf->pubseekoff(0, std::ios_base::beg, std::ios_base::in | std::ios_base::out);
  return size;
}

M2_string streams_Stostring(std::stringbuf* buf, size_t offset, size_t size)
{
  if (size == 0) return M2_tostring("");
  assert(offset + size <= buf->str().size());
  return M2_tostring(buf->str().substr(offset, size).c_str());
}

int streams_Sgetc(std::stringbuf* buf) { return buf->sgetc(); }
size_t streams_Sgetn(M2_string s, size_t n, std::stringbuf* buf)
{
  return buf->sgetn(M2_tocharstar(s), n);
}

int streams_Sputc(char c, std::stringbuf* buf) { return buf->sputc(c); }
size_t streams_Sputn(M2_string s, size_t n, std::stringbuf* buf)
{
  return buf->sputn(M2_tocharstar(s), n);
}

int streams_Srewindc(std::stringbuf* buf) { buf->sungetc(); return buf->sbumpc(); }

void streams_Sempty(std::stringbuf* buf) { return buf->str(""); }

size_t streams_Slength(std::stringbuf* buf) { return buf->str().size(); }

std::stringbuf* streams_Snew(M2_string str)
{
  return new std::stringbuf(M2_tocharstar(str));
}

} /* extern "C" */
