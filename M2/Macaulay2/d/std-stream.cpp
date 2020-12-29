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

std::mutex stdout_mutex, stderr_mutex;
thread_local std::stringstream *outstream, *errstream;
thread_local std::stringbuf *outbuf = outstream->rdbuf(),
                            *errbuf = errstream->rdbuf();

extern "C" {

// TODO: use std::make_unique for these two?
// New stringstream from M2_string
std::stringstream* streams_Sstringstream(M2_string str)
{
  return new std::stringstream(
      M2_tocharstar(str),
      std::ios_base::in | std::ios_base::out | std::ios_base::ate);
}

// New stringbuf from stringstream
std::stringbuf* streams_Sstringbuf(std::stringstream stream)
{
  return stream.rdbuf();
}

// Get the number of characters not yet printed
size_t streams_Slength(std::stringstream* stream) { return stream->rdbuf()->in_avail(); }

// Empty the stream, deallocating
void streams_Sempty(std::stringstream* stream) { stream->clear(); return stream->str(""); }

// Flush the buffered content to the appropriate file descriptor
size_t streams_Sflush(int fd, std::stringstream* stream)
{
  if (streams_Slength(stream) == 0) return stream->tellg();
  if (fd == 1)
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << stream->rdbuf() << std::flush;
      return stream->tellg();
    }
  else if (fd == 2)
    {
      std::lock_guard<std::mutex> stderr_guard(stderr_mutex);
      std::cerr << stream->rdbuf() << std::flush;
      return stream->tellg();
    }
  else
    {
      std::lock_guard<std::mutex> stdout_guard(stdout_mutex);
      std::cout << " SHOULD HAVE BEEN TO FILE" << std::endl;
    }
  return -1;
}

/*
size_t streams_Swrite(int fd, std::stringstream* stream, size_t offset, size_t size)
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
*/

M2_string streams_Stostring(std::stringstream* stream, size_t offset, size_t size)
{
  if (size == 0) return M2_tostring("");
  auto str = stream->rdbuf()->str();
  assert(offset + size <= str.size());
  return M2_tostring(str.substr(offset, size).c_str());
}

int streams_Sgetc(std::stringstream* stream) { return stream->get(); }
size_t streams_Sgetn(M2_string s, size_t n, std::stringstream* stream)
{
  stream->read(M2_tocharstar(s), n);
  return n;
}

// Get the last character printed
int streams_Srewindc(std::stringstream* stream) { return stream->unget().get(); }

int streams_Sputc(char c, std::stringstream* stream) { stream->put(c); return c; }
size_t streams_Sputn(M2_string s, size_t n, std::stringstream* stream)
{
  *stream << std::string(M2_tocharstar(s)) << std::flush;
  //  stream->rdbuf()->sputn(M2_tocharstar(s), n);
  return n;
}

} /* extern "C" */
