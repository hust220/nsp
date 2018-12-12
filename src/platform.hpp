#pragma once

#if defined(__APPLE__) && defined(__GNUC__)  
#  define JN_OS_MACX  
#elif defined(__MACOSX__)  
#  define JN_OS_MACX  
#elif defined(macintosh)  
#  define JN_OS_MAC9  
#elif defined(__CYGWIN__)  
#  define JN_OS_CYGWIN  
#elif defined(MSDOS) || defined(_MSDOS)  
#  define JN_OS_MSDOS  
#elif defined(__OS2__)  
#  if defined(__EMX__)  
#    define JN_OS_OS2EMX  
#  else  
#    define JN_OS_OS2  
#  endif  
#elif !defined(SAG_COM) && (defined(WIN64) || defined(_WIN64) || defined(__WIN64__))  
#  define JN_OS_WIN32  
#  define JN_OS_WIN64  
#  define JN_OS_WIN
#elif !defined(SAG_COM) && (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))  
#  define JN_OS_WIN32  
#  define JN_OS_WIN
#elif defined(__MWERKS__) && defined(__INTEL__)  
#  define JN_OS_WIN32  
#  define JN_OS_WIN
#elif defined(__sun) || defined(sun)  
#  define JN_OS_SOLARIS  
#elif defined(hpux) || defined(__hpux)  
#  define JN_OS_HPUX  
#elif defined(__ultrix) || defined(ultrix)  
#  define JN_OS_ULTRIX  
#elif defined(sinix)  
#  define JN_OS_RELIANT  
#elif defined(__linux__) || defined(__linux)  
#  define JN_OS_LINUX  
#elif defined(__FreeBSD__)  
#  define JN_OS_FREEBSD  
#  define JN_OS_BSD4  
#elif defined(__NetBSD__)  
#  define JN_OS_NETBSD  
#  define JN_OS_BSD4  
#elif defined(__OpenBSD__)  
#  define JN_OS_OPENBSD  
#  define JN_OS_BSD4  
#elif defined(__bsdi__)  
#  define JN_OS_BSDI  
#  define JN_OS_BSD4  
#elif defined(__sgi)  
#  define JN_OS_IRIX  
#elif defined(__osf__)  
#  define JN_OS_OSF  
#elif defined(_AIX)  
#  define JN_OS_AIX  
#elif defined(__Lynx__)  
#  define JN_OS_LYNX  
#elif defined(__GNU_HURD__)  
#  define JN_OS_HURD  
#elif defined(__DGUX__)  
#  define JN_OS_DGUX  
#elif defined(__QNXNTO__)  
#  define JN_OS_QNX6  
#elif defined(__QNX__)  
#  define JN_OS_QNX  
#elif defined(_SEQUENT_)  
#  define JN_OS_DYNIX  
#elif defined(_SCO_DS)                   /* SCO OpenServer 5 + GCC */  
#  define JN_OS_SCO  
#elif defined(__USLC__)                  /* all SCO platforms + UDK or OUDK */  
#  define JN_OS_UNIXWARE  
#  define JN_OS_UNIXWARE7  
#elif defined(__svr4__) && defined(i386) /* Open UNIX 8 + GCC */  
#  define JN_OS_UNIXWARE  
#  define JN_OS_UNIXWARE7  
#else  
#  error "Qt has not been ported to this OS - talk to qt-bugs@trolltech.com"  
#endif  

#if defined(JN_OS_MAC9) || defined(JN_OS_MACX)  
#  define JN_OS_MAC  
#endif  

#if defined(JN_OS_MAC9) || defined(JN_OS_MSDOS) || defined(JN_OS_OS2) || defined(JN_OS_WIN32) || defined(JN_OS_WIN64)  
#  undef JN_OS_UNIX  
#elif !defined(JN_OS_UNIX)  
#  define JN_OS_UNIX  
#endif
