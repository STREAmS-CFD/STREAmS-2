#ifdef __linux__

#define _GNU_SOURCE
#include <stdio.h>
#include <sys/sysinfo.h>
#include <unistd.h>

void getmemory(long int* totmem, long int* avmem) {
    /* PAGESIZE is POSIX: http://pubs.opengroup.org/onlinepubs/9699919799/
     * but PHYS_PAGES and AVPHYS_PAGES are glibc extensions. I bet those are
     * parsed from /proc/meminfo. */
    //sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGESIZE)
    //sysconf(_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE)

    /* glibc extensions. man says they are parsed from /proc/meminfo. */
    *totmem = get_phys_pages()   * sysconf(_SC_PAGESIZE);
    *avmem  = get_avphys_pages() * sysconf(_SC_PAGESIZE);
}

#else

void getmemory(long int* totmem, long int* avmem) {
    *totmem = 0L;
    *avmem  = 0L;
}

#endif
