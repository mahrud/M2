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

// TODO: use std::make_unique
// New stringbuf from M2_string
std::stringbuf* streams_newstringbuf(M2_string str)
{
  return new std::stringbuf(
      M2_tocharstar(str),
      std::ios_base::in | std::ios_base::out | std::ios_base::ate);
}

// Get the number of characters not yet printed
size_t streams_Slength(std::stringbuf* buf) { return buf->in_avail(); }

// Returns the character as an integer
int streams_Sgetc(std::stringbuf* buf) { return buf->sgetc(); }
size_t streams_Sgetn(M2_string s, size_t n, std::stringbuf* buf)
{
  return buf->sgetn(M2_tocharstar(s), n);
}

// Returns the character as an integer and advance put pointer
int streams_Sputc(char c, std::stringbuf* buf) { return buf->sputc(c); }
// Note: does not advance put pointer
size_t streams_Sputn(M2_string s, size_t n, std::stringbuf* buf)
{
  return buf->sputn(M2_tocharstar(s), n);
}

M2_string streams_buftostring(std::stringbuf* buf, size_t offset, size_t size)
{
  if (size == 0) return M2_tostring("");
  auto str = buf->str();
  assert(offset + size <= str.size());
  return M2_tostring(str.substr(offset, size).c_str());
}

// Empty the stream, deallocating
void streams_cleanbuf(std::stringbuf* buf) { buf->str(""); }

/***************************************************************************/

// TODO: use std::make_unique
// New stringstream from M2_string
std::stringstream* streams_newstringstream(M2_string str)
{
  return new std::stringstream(
      M2_tocharstar(str),
      std::ios_base::in | std::ios_base::out | std::ios_base::ate);
}

// Get the internal stringbuf of a stringstream
std::stringbuf* streams_getstringbuf(std::stringstream stream)
{
  return stream.rdbuf();
}

// Empty the stream, deallocating
void streams_cleanstream(std::stringstream* stream) { stream->clear(); return stream->str(""); }

// Flush the buffered content to the appropriate file descriptor
size_t streams_Sflush(int fd, std::stringstream* stream)
{
  if (streams_Slength(stream->rdbuf()) == 0) return stream->tellg();
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

// Note: these advance the read pointer
int streams_Sget(std::stringstream* stream) { return stream->get(); }
size_t streams_Sread(M2_string s, size_t n, std::stringstream* stream)
{
  stream->read(M2_tocharstar(s), n);
  return n;
}

// Get the last character printed
int streams_Srewindc(std::stringstream* stream) { return stream->unget().get(); }

int streams_Sput(char c, std::stringstream* stream) { stream->put(c); return c; }
size_t streams_Swrite(M2_string s, size_t n, std::stringstream* stream)
{
  *stream << std::string(M2_tocharstar(s)) << std::flush;
  // TODO: if we ever need to write binary data, we should use this instead:
  //  stream->write(M2_tocharstar(s), n);
  return n;
}

} /* extern "C" */
