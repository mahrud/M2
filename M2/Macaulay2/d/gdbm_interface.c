/**********************************************
 *                  dbm stuff                 *
 **********************************************/

#include <stdio.h>
#include <string.h>
#include <gdbm.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include "M2mem.h"

#define TRUE 1
#define FALSE 0
#define ERROR (-1)
#define OPEN_RETRY_USECONDS 100000

typedef char bool;
static int numfiles = 0;
static GDBM_FILE *gdbm_files = NULL;
static datum *lastkeys = NULL;
static bool *hadlastkeys = NULL;
static pthread_mutex_t gdbm_mutex = PTHREAD_MUTEX_INITIALIZER;
static __thread gdbm_error last_gdbm_errno = GDBM_NO_ERROR;

static void lock_dbms(void) {
     if (pthread_mutex_lock(&gdbm_mutex) != 0) abort();
     }

static void unlock_dbms(void) {
     if (pthread_mutex_unlock(&gdbm_mutex) != 0) abort();
     }

static void clear_lastkey(int handle) {
     if (hadlastkeys[handle]) {
	  free(lastkeys[handle].dptr);
	  lastkeys[handle].dptr = NULL;
	  lastkeys[handle].dsize = 0;
	  hadlastkeys[handle] = FALSE;
	  }
     }

static void initialize_file_slots(int first, int last) {
     int i;
     for (i=first; i<last; i++) {
	  gdbm_files[i] = NULL;
	  lastkeys[i].dptr = NULL;
	  lastkeys[i].dsize = 0;
	  hadlastkeys[i] = FALSE;
	  }
     }

static bool is_lock_error(gdbm_error error) {
     return error == GDBM_CANT_BE_READER || error == GDBM_CANT_BE_WRITER;
     }

void close_all_dbms(void) {
     int i;
     lock_dbms();
     for (i=0; i<numfiles; i++) {
	  clear_lastkey(i);
	  if (gdbm_files[i] != NULL) {
	       gdbm_close(gdbm_files[i]);
	       gdbm_files[i] = NULL;
	       }
	  }
     unlock_dbms();
     }

#include "M2-exports.h"

int system_dbmopen(M2_string filename, M2_bool mutable) {
     int gdbm_handle;
     int flags = mutable ? GDBM_WRCREAT : GDBM_READER;
     int mode = 0666;
     char *FileName = M2_tocharstar(filename);
     GDBM_FILE f;
     lock_dbms();
     for (;;) {
	  f = gdbm_open(FileName, 4096, flags, mode, NULL);
	  last_gdbm_errno = gdbm_errno;
	  if (f != NULL || !is_lock_error(last_gdbm_errno)) break;
	  unlock_dbms();
	  usleep(OPEN_RETRY_USECONDS);
	  lock_dbms();
	  }
     freemem(FileName);
     if (f == NULL) {
	  unlock_dbms();
	  return ERROR;
	  }
     if (numfiles == 0) {
	  numfiles = 10;
	  gdbm_files = (GDBM_FILE *) getmem(numfiles * sizeof(GDBM_FILE));
	  lastkeys = (datum *) getmem(numfiles * sizeof(datum));
	  hadlastkeys = (bool *) getmem(numfiles * sizeof(bool));
	  initialize_file_slots(0, numfiles);
	  gdbm_handle = 0;
	  }
     else {
	  for (gdbm_handle=0; TRUE ; gdbm_handle++) {
	       if (gdbm_handle==numfiles) {
		    GDBM_FILE *p;
		    datum *lk;
		    bool *hlk;
		    int j;
		    int oldnumfiles = numfiles;
		    numfiles *= 2;
		    p = (GDBM_FILE *) getmem(numfiles * sizeof(GDBM_FILE));
		    lk = (datum *) getmem(numfiles * sizeof(datum));
		    hlk = (bool *) getmem(numfiles * sizeof(bool));
		    for (j=0; j<gdbm_handle; j++) p[j] = gdbm_files[j];
		    for (j=0; j<gdbm_handle; j++) lk[j] = lastkeys[j];
		    for (j=0; j<gdbm_handle; j++) hlk[j] = hadlastkeys[j];
		    gdbm_files = p;
		    lastkeys = lk;
		    hadlastkeys = hlk;
	  	    initialize_file_slots(oldnumfiles, numfiles);
		    break;
		    }
	       else if (gdbm_files[gdbm_handle] == NULL) break;
	       }
	  }
     gdbm_files[gdbm_handle] = f;
     unlock_dbms();
     return gdbm_handle;
     }

int system_dbmclose(int handle) {
     lock_dbms();
     clear_lastkey(handle);
     gdbm_close(gdbm_files[handle]);
     last_gdbm_errno = gdbm_errno;
     gdbm_files[handle] = NULL;
     unlock_dbms();
     return 0;
     }

static datum todatum(M2_string x) {
     datum y;
     y.dptr = (char *)x->array;
     y.dsize = x->len;
     return y;
     }

static M2_string fromdatum(datum y) {
     M2_string x;
     if (y.dptr == NULL) return NULL;
     x = (M2_string)getmem(sizeofarray(x,y.dsize));
     x->len = y.dsize;
     memcpy(x->array, y.dptr, y.dsize);
     return x;
     }

static M2_string fromdatum_and_free(datum y) {
     M2_string x = fromdatum(y);
     free(y.dptr);
     y.dptr = NULL;
     y.dsize = 0;
     return x;
     }

int system_dbmstore(int handle, M2_string key, M2_string content) {
     int ret;
     lock_dbms();
     clear_lastkey(handle);
     ret = gdbm_store(gdbm_files[handle],todatum(key),todatum(content),GDBM_REPLACE);
     last_gdbm_errno = gdbm_errno;
     unlock_dbms();
     return ret;
     }

M2_string /* or NULL */ system_dbmfetch(int handle, M2_string key) {
     M2_string ret;
     lock_dbms();
     ret = fromdatum_and_free(gdbm_fetch(gdbm_files[handle],todatum(key)));
     last_gdbm_errno = gdbm_errno;
     unlock_dbms();
     return ret;
     }

int system_dbmdelete(int handle, M2_string key) {
     int ret;
     lock_dbms();
     clear_lastkey(handle);
     ret = gdbm_delete(gdbm_files[handle],todatum(key));
     last_gdbm_errno = gdbm_errno;
     unlock_dbms();
     return ret;
     }

M2_string /* or NULL */ system_dbmfirst(int handle) {
     M2_string ret;
     lock_dbms();
     clear_lastkey(handle);
     lastkeys[handle] = gdbm_firstkey(gdbm_files[handle]);
     last_gdbm_errno = gdbm_errno;
     hadlastkeys[handle] = TRUE;
     ret = fromdatum(lastkeys[handle]);
     unlock_dbms();
     return ret;
     }

M2_string /* or NULL */ system_dbmnext(int handle) {
     M2_string ret;
     lock_dbms();
     if (hadlastkeys[handle]) {
	  datum x = gdbm_nextkey(gdbm_files[handle], lastkeys[handle]);
	  last_gdbm_errno = gdbm_errno;
	  free(lastkeys[handle].dptr);
	  lastkeys[handle] = x;
	  ret = fromdatum(lastkeys[handle]);
	  }
     else {
	  lastkeys[handle] = gdbm_firstkey(gdbm_files[handle]);
	  last_gdbm_errno = gdbm_errno;
	  hadlastkeys[handle] = TRUE;
	  ret = fromdatum(lastkeys[handle]);
	  }
     unlock_dbms();
     return ret;
     }

int system_dbmreorganize(int handle) {
     int ret;
     lock_dbms();
     clear_lastkey(handle);
     ret = gdbm_reorganize(gdbm_files[handle]);
     last_gdbm_errno = gdbm_errno;
     unlock_dbms();
     return ret;
     }

M2_string system_dbmstrerror(void) {
     M2_string ret;
     lock_dbms();
     ret = M2_tostring(gdbm_strerror(last_gdbm_errno));
     unlock_dbms();
     return ret;
     }

/*
 Local Variables:
 compile-command: "echo \"make: Entering directory \\`$M2BUILDDIR/Macaulay2/d'\" && make -C $M2BUILDDIR/Macaulay2/d gdbm_interface.o "
 End:
*/
