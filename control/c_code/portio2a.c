# 1 "portio2a.c"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "portio2a.c"
# 11 "portio2a.c"
# 1 "/usr/include/stdio.h" 1 3 4
# 28 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/features.h" 1 3 4
# 361 "/usr/include/features.h" 3 4
# 1 "/usr/include/sys/cdefs.h" 1 3 4
# 373 "/usr/include/sys/cdefs.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 374 "/usr/include/sys/cdefs.h" 2 3 4
# 362 "/usr/include/features.h" 2 3 4
# 385 "/usr/include/features.h" 3 4
# 1 "/usr/include/gnu/stubs.h" 1 3 4



# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 5 "/usr/include/gnu/stubs.h" 2 3 4




# 1 "/usr/include/gnu/stubs-64.h" 1 3 4
# 10 "/usr/include/gnu/stubs.h" 2 3 4
# 386 "/usr/include/features.h" 2 3 4
# 29 "/usr/include/stdio.h" 2 3 4





# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 211 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 3 4
typedef long unsigned int size_t;
# 35 "/usr/include/stdio.h" 2 3 4

# 1 "/usr/include/bits/types.h" 1 3 4
# 28 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/bits/types.h" 2 3 4


typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;







typedef long int __quad_t;
typedef unsigned long int __u_quad_t;
# 131 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/typesizes.h" 1 3 4
# 132 "/usr/include/bits/types.h" 2 3 4


typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned long int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned long int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;

typedef int __daddr_t;
typedef long int __swblk_t;
typedef int __key_t;


typedef int __clockid_t;


typedef void * __timer_t;


typedef long int __blksize_t;




typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;


typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;


typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;

typedef long int __ssize_t;



typedef __off64_t __loff_t;
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;


typedef long int __intptr_t;


typedef unsigned int __socklen_t;
# 37 "/usr/include/stdio.h" 2 3 4
# 45 "/usr/include/stdio.h" 3 4
struct _IO_FILE;



typedef struct _IO_FILE FILE;





# 65 "/usr/include/stdio.h" 3 4
typedef struct _IO_FILE __FILE;
# 75 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/libio.h" 1 3 4
# 32 "/usr/include/libio.h" 3 4
# 1 "/usr/include/_G_config.h" 1 3 4
# 15 "/usr/include/_G_config.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 16 "/usr/include/_G_config.h" 2 3 4




# 1 "/usr/include/wchar.h" 1 3 4
# 83 "/usr/include/wchar.h" 3 4
typedef struct
{
  int __count;
  union
  {

    unsigned int __wch;



    char __wchb[4];
  } __value;
} __mbstate_t;
# 21 "/usr/include/_G_config.h" 2 3 4

typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;
# 53 "/usr/include/_G_config.h" 3 4
typedef int _G_int16_t __attribute__ ((__mode__ (__HI__)));
typedef int _G_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int _G_uint16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int _G_uint32_t __attribute__ ((__mode__ (__SI__)));
# 33 "/usr/include/libio.h" 2 3 4
# 53 "/usr/include/libio.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h" 1 3 4
# 40 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 54 "/usr/include/libio.h" 2 3 4
# 170 "/usr/include/libio.h" 3 4
struct _IO_jump_t; struct _IO_FILE;
# 180 "/usr/include/libio.h" 3 4
typedef void _IO_lock_t;





struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;



  int _pos;
# 203 "/usr/include/libio.h" 3 4
};


enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};
# 271 "/usr/include/libio.h" 3 4
struct _IO_FILE {
  int _flags;




  char* _IO_read_ptr;
  char* _IO_read_end;
  char* _IO_read_base;
  char* _IO_write_base;
  char* _IO_write_ptr;
  char* _IO_write_end;
  char* _IO_buf_base;
  char* _IO_buf_end;

  char *_IO_save_base;
  char *_IO_backup_base;
  char *_IO_save_end;

  struct _IO_marker *_markers;

  struct _IO_FILE *_chain;

  int _fileno;



  int _flags2;

  __off_t _old_offset;



  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];



  _IO_lock_t *_lock;
# 319 "/usr/include/libio.h" 3 4
  __off64_t _offset;
# 328 "/usr/include/libio.h" 3 4
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;

  int _mode;

  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];

};


typedef struct _IO_FILE _IO_FILE;


struct _IO_FILE_plus;

extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;
# 364 "/usr/include/libio.h" 3 4
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);







typedef __ssize_t __io_write_fn (void *__cookie, __const char *__buf,
     size_t __n);







typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);


typedef int __io_close_fn (void *__cookie);
# 416 "/usr/include/libio.h" 3 4
extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);
# 460 "/usr/include/libio.h" 3 4
extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) __attribute__ ((__nothrow__));
extern int _IO_ferror (_IO_FILE *__fp) __attribute__ ((__nothrow__));

extern int _IO_peekc_locked (_IO_FILE *__fp);





extern void _IO_flockfile (_IO_FILE *) __attribute__ ((__nothrow__));
extern void _IO_funlockfile (_IO_FILE *) __attribute__ ((__nothrow__));
extern int _IO_ftrylockfile (_IO_FILE *) __attribute__ ((__nothrow__));
# 490 "/usr/include/libio.h" 3 4
extern int _IO_vfscanf (_IO_FILE * __restrict, const char * __restrict,
   __gnuc_va_list, int *__restrict);
extern int _IO_vfprintf (_IO_FILE *__restrict, const char *__restrict,
    __gnuc_va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);

extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);

extern void _IO_free_backup_area (_IO_FILE *) __attribute__ ((__nothrow__));
# 76 "/usr/include/stdio.h" 2 3 4




typedef __gnuc_va_list va_list;
# 93 "/usr/include/stdio.h" 3 4
typedef __off64_t off_t;
# 103 "/usr/include/stdio.h" 3 4
typedef __ssize_t ssize_t;









typedef _G_fpos64_t fpos_t;


# 161 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/bits/stdio_lim.h" 1 3 4
# 162 "/usr/include/stdio.h" 2 3 4



extern struct _IO_FILE *stdin;
extern struct _IO_FILE *stdout;
extern struct _IO_FILE *stderr;









extern int remove (__const char *__filename) __attribute__ ((__nothrow__));

extern int rename (__const char *__old, __const char *__new) __attribute__ ((__nothrow__));




extern int renameat (int __oldfd, __const char *__old, int __newfd,
       __const char *__new) __attribute__ ((__nothrow__));



# 197 "/usr/include/stdio.h" 3 4
extern FILE *tmpfile (void) __asm__ ("" "tmpfile64") ;
# 208 "/usr/include/stdio.h" 3 4
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__)) ;





extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__)) ;
# 226 "/usr/include/stdio.h" 3 4
extern char *tempnam (__const char *__dir, __const char *__pfx)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;








extern int fclose (FILE *__stream);




extern int fflush (FILE *__stream);

# 251 "/usr/include/stdio.h" 3 4
extern int fflush_unlocked (FILE *__stream);
# 265 "/usr/include/stdio.h" 3 4

# 282 "/usr/include/stdio.h" 3 4
extern FILE *fopen (__const char *__restrict __filename, __const char *__restrict __modes) __asm__ ("" "fopen64")

  ;
extern FILE *freopen (__const char *__restrict __filename, __const char *__restrict __modes, FILE *__restrict __stream) __asm__ ("" "freopen64")


  ;






# 305 "/usr/include/stdio.h" 3 4
extern FILE *fdopen (int __fd, __const char *__modes) __attribute__ ((__nothrow__)) ;
# 318 "/usr/include/stdio.h" 3 4
extern FILE *fmemopen (void *__s, size_t __len, __const char *__modes)
  __attribute__ ((__nothrow__)) ;




extern FILE *open_memstream (char **__bufloc, size_t *__sizeloc) __attribute__ ((__nothrow__)) ;






extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__));



extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
      int __modes, size_t __n) __attribute__ ((__nothrow__));





extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
         size_t __size) __attribute__ ((__nothrow__));


extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__));








extern int fprintf (FILE *__restrict __stream,
      __const char *__restrict __format, ...);




extern int printf (__const char *__restrict __format, ...);

extern int sprintf (char *__restrict __s,
      __const char *__restrict __format, ...) __attribute__ ((__nothrow__));





extern int vfprintf (FILE *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg);




extern int vprintf (__const char *__restrict __format, __gnuc_va_list __arg);

extern int vsprintf (char *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg) __attribute__ ((__nothrow__));





extern int snprintf (char *__restrict __s, size_t __maxlen,
       __const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
        __const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));

# 416 "/usr/include/stdio.h" 3 4
extern int vdprintf (int __fd, __const char *__restrict __fmt,
       __gnuc_va_list __arg)
     __attribute__ ((__format__ (__printf__, 2, 0)));
extern int dprintf (int __fd, __const char *__restrict __fmt, ...)
     __attribute__ ((__format__ (__printf__, 2, 3)));








extern int fscanf (FILE *__restrict __stream,
     __const char *__restrict __format, ...) ;




extern int scanf (__const char *__restrict __format, ...) ;

extern int sscanf (__const char *__restrict __s,
     __const char *__restrict __format, ...) __attribute__ ((__nothrow__));
# 447 "/usr/include/stdio.h" 3 4
extern int fscanf (FILE *__restrict __stream, __const char *__restrict __format, ...) __asm__ ("" "__isoc99_fscanf")

                               ;
extern int scanf (__const char *__restrict __format, ...) __asm__ ("" "__isoc99_scanf")
                              ;
extern int sscanf (__const char *__restrict __s, __const char *__restrict __format, ...) __asm__ ("" "__isoc99_sscanf")

                          __attribute__ ((__nothrow__));
# 467 "/usr/include/stdio.h" 3 4








extern int vfscanf (FILE *__restrict __s, __const char *__restrict __format,
      __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;





extern int vscanf (__const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;


extern int vsscanf (__const char *__restrict __s,
      __const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
# 498 "/usr/include/stdio.h" 3 4
extern int vfscanf (FILE *__restrict __s, __const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vfscanf")



     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (__const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vscanf")

     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (__const char *__restrict __s, __const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vsscanf")



     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
# 526 "/usr/include/stdio.h" 3 4









extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);





extern int getchar (void);

# 554 "/usr/include/stdio.h" 3 4
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);
# 565 "/usr/include/stdio.h" 3 4
extern int fgetc_unlocked (FILE *__stream);











extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);





extern int putchar (int __c);

# 598 "/usr/include/stdio.h" 3 4
extern int fputc_unlocked (int __c, FILE *__stream);







extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);






extern int getw (FILE *__stream);


extern int putw (int __w, FILE *__stream);








extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;






extern char *gets (char *__s) ;

# 660 "/usr/include/stdio.h" 3 4
extern __ssize_t __getdelim (char **__restrict __lineptr,
          size_t *__restrict __n, int __delimiter,
          FILE *__restrict __stream) ;
extern __ssize_t getdelim (char **__restrict __lineptr,
        size_t *__restrict __n, int __delimiter,
        FILE *__restrict __stream) ;







extern __ssize_t getline (char **__restrict __lineptr,
       size_t *__restrict __n,
       FILE *__restrict __stream) ;








extern int fputs (__const char *__restrict __s, FILE *__restrict __stream);





extern int puts (__const char *__s);






extern int ungetc (int __c, FILE *__stream);






extern size_t fread (void *__restrict __ptr, size_t __size,
       size_t __n, FILE *__restrict __stream) ;




extern size_t fwrite (__const void *__restrict __ptr, size_t __size,
        size_t __n, FILE *__restrict __s) ;

# 732 "/usr/include/stdio.h" 3 4
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
         size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (__const void *__restrict __ptr, size_t __size,
          size_t __n, FILE *__restrict __stream) ;








extern int fseek (FILE *__stream, long int __off, int __whence);




extern long int ftell (FILE *__stream) ;




extern void rewind (FILE *__stream);

# 776 "/usr/include/stdio.h" 3 4
extern int fseeko (FILE *__stream, __off64_t __off, int __whence) __asm__ ("" "fseeko64")

                  ;
extern __off64_t ftello (FILE *__stream) __asm__ ("" "ftello64");








# 801 "/usr/include/stdio.h" 3 4
extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos) __asm__ ("" "fgetpos64")
                                          ;
extern int fsetpos (FILE *__stream, __const fpos_t *__pos) __asm__ ("" "fsetpos64")
                                                            ;






# 819 "/usr/include/stdio.h" 3 4


extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__));

extern int feof (FILE *__stream) __attribute__ ((__nothrow__)) ;

extern int ferror (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;








extern void perror (__const char *__s);






# 1 "/usr/include/bits/sys_errlist.h" 1 3 4
# 27 "/usr/include/bits/sys_errlist.h" 3 4
extern int sys_nerr;
extern __const char *__const sys_errlist[];
# 849 "/usr/include/stdio.h" 2 3 4




extern int fileno (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
# 868 "/usr/include/stdio.h" 3 4
extern FILE *popen (__const char *__command, __const char *__modes) ;





extern int pclose (FILE *__stream);





extern char *ctermid (char *__s) __attribute__ ((__nothrow__));
# 908 "/usr/include/stdio.h" 3 4
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__));



extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__)) ;


extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__));
# 938 "/usr/include/stdio.h" 3 4

# 12 "portio2a.c" 2
# 1 "/usr/include/stdlib.h" 1 3 4
# 33 "/usr/include/stdlib.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 323 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 3 4
typedef int wchar_t;
# 34 "/usr/include/stdlib.h" 2 3 4








# 1 "/usr/include/bits/waitflags.h" 1 3 4
# 43 "/usr/include/stdlib.h" 2 3 4
# 1 "/usr/include/bits/waitstatus.h" 1 3 4
# 65 "/usr/include/bits/waitstatus.h" 3 4
# 1 "/usr/include/endian.h" 1 3 4
# 37 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/endian.h" 1 3 4
# 38 "/usr/include/endian.h" 2 3 4
# 61 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/byteswap.h" 1 3 4
# 28 "/usr/include/bits/byteswap.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/bits/byteswap.h" 2 3 4
# 62 "/usr/include/endian.h" 2 3 4
# 66 "/usr/include/bits/waitstatus.h" 2 3 4

union wait
  {
    int w_status;
    struct
      {

 unsigned int __w_termsig:7;
 unsigned int __w_coredump:1;
 unsigned int __w_retcode:8;
 unsigned int:16;







      } __wait_terminated;
    struct
      {

 unsigned int __w_stopval:8;
 unsigned int __w_stopsig:8;
 unsigned int:16;






      } __wait_stopped;
  };
# 44 "/usr/include/stdlib.h" 2 3 4
# 68 "/usr/include/stdlib.h" 3 4
typedef union
  {
    union wait *__uptr;
    int *__iptr;
  } __WAIT_STATUS __attribute__ ((__transparent_union__));
# 96 "/usr/include/stdlib.h" 3 4


typedef struct
  {
    int quot;
    int rem;
  } div_t;



typedef struct
  {
    long int quot;
    long int rem;
  } ldiv_t;







__extension__ typedef struct
  {
    long long int quot;
    long long int rem;
  } lldiv_t;


# 140 "/usr/include/stdlib.h" 3 4
extern size_t __ctype_get_mb_cur_max (void) __attribute__ ((__nothrow__)) ;




extern double atof (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern int atoi (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern long int atol (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





__extension__ extern long long int atoll (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





extern double strtod (__const char *__restrict __nptr,
        char **__restrict __endptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern float strtof (__const char *__restrict __nptr,
       char **__restrict __endptr) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

extern long double strtold (__const char *__restrict __nptr,
       char **__restrict __endptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern long int strtol (__const char *__restrict __nptr,
   char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

extern unsigned long int strtoul (__const char *__restrict __nptr,
      char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




__extension__
extern long long int strtoq (__const char *__restrict __nptr,
        char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

__extension__
extern unsigned long long int strtouq (__const char *__restrict __nptr,
           char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





__extension__
extern long long int strtoll (__const char *__restrict __nptr,
         char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

__extension__
extern unsigned long long int strtoull (__const char *__restrict __nptr,
     char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

# 311 "/usr/include/stdlib.h" 3 4
extern char *l64a (long int __n) __attribute__ ((__nothrow__)) ;


extern long int a64l (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;




# 1 "/usr/include/sys/types.h" 1 3 4
# 28 "/usr/include/sys/types.h" 3 4






typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;




typedef __loff_t loff_t;





typedef __ino64_t ino_t;
# 61 "/usr/include/sys/types.h" 3 4
typedef __dev_t dev_t;




typedef __gid_t gid_t;




typedef __mode_t mode_t;




typedef __nlink_t nlink_t;




typedef __uid_t uid_t;
# 99 "/usr/include/sys/types.h" 3 4
typedef __pid_t pid_t;





typedef __id_t id_t;
# 116 "/usr/include/sys/types.h" 3 4
typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;





typedef __key_t key_t;
# 133 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/time.h" 1 3 4
# 58 "/usr/include/time.h" 3 4


typedef __clock_t clock_t;



# 74 "/usr/include/time.h" 3 4


typedef __time_t time_t;



# 92 "/usr/include/time.h" 3 4
typedef __clockid_t clockid_t;
# 104 "/usr/include/time.h" 3 4
typedef __timer_t timer_t;
# 134 "/usr/include/sys/types.h" 2 3 4
# 147 "/usr/include/sys/types.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 148 "/usr/include/sys/types.h" 2 3 4



typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
# 195 "/usr/include/sys/types.h" 3 4
typedef int int8_t __attribute__ ((__mode__ (__QI__)));
typedef int int16_t __attribute__ ((__mode__ (__HI__)));
typedef int int32_t __attribute__ ((__mode__ (__SI__)));
typedef int int64_t __attribute__ ((__mode__ (__DI__)));


typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));

typedef int register_t __attribute__ ((__mode__ (__word__)));
# 220 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/sys/select.h" 1 3 4
# 31 "/usr/include/sys/select.h" 3 4
# 1 "/usr/include/bits/select.h" 1 3 4
# 23 "/usr/include/bits/select.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 24 "/usr/include/bits/select.h" 2 3 4
# 32 "/usr/include/sys/select.h" 2 3 4


# 1 "/usr/include/bits/sigset.h" 1 3 4
# 24 "/usr/include/bits/sigset.h" 3 4
typedef int __sig_atomic_t;




typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;
# 35 "/usr/include/sys/select.h" 2 3 4



typedef __sigset_t sigset_t;





# 1 "/usr/include/time.h" 1 3 4
# 120 "/usr/include/time.h" 3 4
struct timespec
  {
    __time_t tv_sec;
    long int tv_nsec;
  };
# 45 "/usr/include/sys/select.h" 2 3 4

# 1 "/usr/include/bits/time.h" 1 3 4
# 75 "/usr/include/bits/time.h" 3 4
struct timeval
  {
    __time_t tv_sec;
    __suseconds_t tv_usec;
  };
# 47 "/usr/include/sys/select.h" 2 3 4


typedef __suseconds_t suseconds_t;





typedef long int __fd_mask;
# 67 "/usr/include/sys/select.h" 3 4
typedef struct
  {






    __fd_mask __fds_bits[1024 / (8 * (int) sizeof (__fd_mask))];


  } fd_set;






typedef __fd_mask fd_mask;
# 99 "/usr/include/sys/select.h" 3 4

# 109 "/usr/include/sys/select.h" 3 4
extern int select (int __nfds, fd_set *__restrict __readfds,
     fd_set *__restrict __writefds,
     fd_set *__restrict __exceptfds,
     struct timeval *__restrict __timeout);
# 121 "/usr/include/sys/select.h" 3 4
extern int pselect (int __nfds, fd_set *__restrict __readfds,
      fd_set *__restrict __writefds,
      fd_set *__restrict __exceptfds,
      const struct timespec *__restrict __timeout,
      const __sigset_t *__restrict __sigmask);



# 221 "/usr/include/sys/types.h" 2 3 4


# 1 "/usr/include/sys/sysmacros.h" 1 3 4
# 30 "/usr/include/sys/sysmacros.h" 3 4
__extension__
extern unsigned int gnu_dev_major (unsigned long long int __dev)
     __attribute__ ((__nothrow__));
__extension__
extern unsigned int gnu_dev_minor (unsigned long long int __dev)
     __attribute__ ((__nothrow__));
__extension__
extern unsigned long long int gnu_dev_makedev (unsigned int __major,
            unsigned int __minor)
     __attribute__ ((__nothrow__));
# 224 "/usr/include/sys/types.h" 2 3 4





typedef __blksize_t blksize_t;
# 249 "/usr/include/sys/types.h" 3 4
typedef __blkcnt64_t blkcnt_t;



typedef __fsblkcnt64_t fsblkcnt_t;



typedef __fsfilcnt64_t fsfilcnt_t;
# 271 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/bits/pthreadtypes.h" 1 3 4
# 23 "/usr/include/bits/pthreadtypes.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 24 "/usr/include/bits/pthreadtypes.h" 2 3 4
# 50 "/usr/include/bits/pthreadtypes.h" 3 4
typedef unsigned long int pthread_t;


typedef union
{
  char __size[56];
  long int __align;
} pthread_attr_t;



typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;
# 76 "/usr/include/bits/pthreadtypes.h" 3 4
typedef union
{
  struct __pthread_mutex_s
  {
    int __lock;
    unsigned int __count;
    int __owner;

    unsigned int __nusers;



    int __kind;

    int __spins;
    __pthread_list_t __list;
# 101 "/usr/include/bits/pthreadtypes.h" 3 4
  } __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;




typedef union
{
  struct
  {
    int __lock;
    unsigned int __futex;
    __extension__ unsigned long long int __total_seq;
    __extension__ unsigned long long int __wakeup_seq;
    __extension__ unsigned long long int __woken_seq;
    void *__mutex;
    unsigned int __nwaiters;
    unsigned int __broadcast_seq;
  } __data;
  char __size[48];
  __extension__ long long int __align;
} pthread_cond_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;



typedef unsigned int pthread_key_t;



typedef int pthread_once_t;





typedef union
{

  struct
  {
    int __lock;
    unsigned int __nr_readers;
    unsigned int __readers_wakeup;
    unsigned int __writer_wakeup;
    unsigned int __nr_readers_queued;
    unsigned int __nr_writers_queued;
    int __writer;
    int __shared;
    unsigned long int __pad1;
    unsigned long int __pad2;


    unsigned int __flags;
  } __data;
# 187 "/usr/include/bits/pthreadtypes.h" 3 4
  char __size[56];
  long int __align;
} pthread_rwlock_t;

typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;





typedef volatile int pthread_spinlock_t;




typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;
# 272 "/usr/include/sys/types.h" 2 3 4



# 321 "/usr/include/stdlib.h" 2 3 4






extern long int random (void) __attribute__ ((__nothrow__));


extern void srandom (unsigned int __seed) __attribute__ ((__nothrow__));





extern char *initstate (unsigned int __seed, char *__statebuf,
   size_t __statelen) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));



extern char *setstate (char *__statebuf) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));







struct random_data
  {
    int32_t *fptr;
    int32_t *rptr;
    int32_t *state;
    int rand_type;
    int rand_deg;
    int rand_sep;
    int32_t *end_ptr;
  };

extern int random_r (struct random_data *__restrict __buf,
       int32_t *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern int srandom_r (unsigned int __seed, struct random_data *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));

extern int initstate_r (unsigned int __seed, char *__restrict __statebuf,
   size_t __statelen,
   struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4)));

extern int setstate_r (char *__restrict __statebuf,
         struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));






extern int rand (void) __attribute__ ((__nothrow__));

extern void srand (unsigned int __seed) __attribute__ ((__nothrow__));




extern int rand_r (unsigned int *__seed) __attribute__ ((__nothrow__));







extern double drand48 (void) __attribute__ ((__nothrow__));
extern double erand48 (unsigned short int __xsubi[3]) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern long int lrand48 (void) __attribute__ ((__nothrow__));
extern long int nrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern long int mrand48 (void) __attribute__ ((__nothrow__));
extern long int jrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern void srand48 (long int __seedval) __attribute__ ((__nothrow__));
extern unsigned short int *seed48 (unsigned short int __seed16v[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
extern void lcong48 (unsigned short int __param[7]) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





struct drand48_data
  {
    unsigned short int __x[3];
    unsigned short int __old_x[3];
    unsigned short int __c;
    unsigned short int __init;
    unsigned long long int __a;
  };


extern int drand48_r (struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int erand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int lrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int nrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int mrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int jrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int srand48_r (long int __seedval, struct drand48_data *__buffer)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));

extern int seed48_r (unsigned short int __seed16v[3],
       struct drand48_data *__buffer) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern int lcong48_r (unsigned short int __param[7],
        struct drand48_data *__buffer)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));









extern void *malloc (size_t __size) __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;

extern void *calloc (size_t __nmemb, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;










extern void *realloc (void *__ptr, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__warn_unused_result__));

extern void free (void *__ptr) __attribute__ ((__nothrow__));




extern void cfree (void *__ptr) __attribute__ ((__nothrow__));



# 1 "/usr/include/alloca.h" 1 3 4
# 25 "/usr/include/alloca.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 26 "/usr/include/alloca.h" 2 3 4







extern void *alloca (size_t __size) __attribute__ ((__nothrow__));






# 498 "/usr/include/stdlib.h" 2 3 4





extern void *valloc (size_t __size) __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;




extern int posix_memalign (void **__memptr, size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




extern void abort (void) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));



extern int atexit (void (*__func) (void)) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 531 "/usr/include/stdlib.h" 3 4





extern int on_exit (void (*__func) (int __status, void *__arg), void *__arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern void exit (int __status) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));
# 554 "/usr/include/stdlib.h" 3 4






extern void _Exit (int __status) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));






extern char *getenv (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




extern char *__secure_getenv (__const char *__name)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern int putenv (char *__string) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int setenv (__const char *__name, __const char *__value, int __replace)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));


extern int unsetenv (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern int clearenv (void) __attribute__ ((__nothrow__));
# 606 "/usr/include/stdlib.h" 3 4
extern char *mktemp (char *__template) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 623 "/usr/include/stdlib.h" 3 4
extern int mkstemp (char *__template) __asm__ ("" "mkstemp64")
     __attribute__ ((__nonnull__ (1))) ;
# 645 "/usr/include/stdlib.h" 3 4
extern int mkstemps (char *__template, int __suffixlen) __asm__ ("" "mkstemps64")
                     __attribute__ ((__nonnull__ (1))) ;
# 663 "/usr/include/stdlib.h" 3 4
extern char *mkdtemp (char *__template) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 712 "/usr/include/stdlib.h" 3 4





extern int system (__const char *__command) ;

# 734 "/usr/include/stdlib.h" 3 4
extern char *realpath (__const char *__restrict __name,
         char *__restrict __resolved) __attribute__ ((__nothrow__)) ;






typedef int (*__compar_fn_t) (__const void *, __const void *);
# 752 "/usr/include/stdlib.h" 3 4



extern void *bsearch (__const void *__key, __const void *__base,
        size_t __nmemb, size_t __size, __compar_fn_t __compar)
     __attribute__ ((__nonnull__ (1, 2, 5))) ;



extern void qsort (void *__base, size_t __nmemb, size_t __size,
     __compar_fn_t __compar) __attribute__ ((__nonnull__ (1, 4)));
# 771 "/usr/include/stdlib.h" 3 4
extern int abs (int __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;
extern long int labs (long int __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;



__extension__ extern long long int llabs (long long int __x)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;







extern div_t div (int __numer, int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;
extern ldiv_t ldiv (long int __numer, long int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;




__extension__ extern lldiv_t lldiv (long long int __numer,
        long long int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;

# 808 "/usr/include/stdlib.h" 3 4
extern char *ecvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *fcvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *gcvt (double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3))) ;




extern char *qecvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qfcvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qgcvt (long double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3))) ;




extern int ecvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int fcvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));

extern int qecvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int qfcvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));







extern int mblen (__const char *__s, size_t __n) __attribute__ ((__nothrow__)) ;


extern int mbtowc (wchar_t *__restrict __pwc,
     __const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__)) ;


extern int wctomb (char *__s, wchar_t __wchar) __attribute__ ((__nothrow__)) ;



extern size_t mbstowcs (wchar_t *__restrict __pwcs,
   __const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__));

extern size_t wcstombs (char *__restrict __s,
   __const wchar_t *__restrict __pwcs, size_t __n)
     __attribute__ ((__nothrow__));








extern int rpmatch (__const char *__response) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 896 "/usr/include/stdlib.h" 3 4
extern int getsubopt (char **__restrict __optionp,
        char *__const *__restrict __tokens,
        char **__restrict __valuep)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2, 3))) ;
# 948 "/usr/include/stdlib.h" 3 4
extern int getloadavg (double __loadavg[], int __nelem)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 964 "/usr/include/stdlib.h" 3 4

# 13 "portio2a.c" 2
# 1 "/usr/include/string.h" 1 3 4
# 29 "/usr/include/string.h" 3 4





# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 35 "/usr/include/string.h" 2 3 4









extern void *memcpy (void *__restrict __dest,
       __const void *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern void *memmove (void *__dest, __const void *__src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));






extern void *memccpy (void *__restrict __dest, __const void *__restrict __src,
        int __c, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));





extern void *memset (void *__s, int __c, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int memcmp (__const void *__s1, __const void *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 95 "/usr/include/string.h" 3 4
extern void *memchr (__const void *__s, int __c, size_t __n)
      __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 126 "/usr/include/string.h" 3 4


extern char *strcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncpy (char *__restrict __dest,
        __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern char *strcat (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncat (char *__restrict __dest, __const char *__restrict __src,
        size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcmp (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern int strncmp (__const char *__s1, __const char *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcoll (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern size_t strxfrm (char *__restrict __dest,
         __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));






# 1 "/usr/include/xlocale.h" 1 3 4
# 28 "/usr/include/xlocale.h" 3 4
typedef struct __locale_struct
{

  struct __locale_data *__locales[13];


  const unsigned short int *__ctype_b;
  const int *__ctype_tolower;
  const int *__ctype_toupper;


  const char *__names[13];
} *__locale_t;


typedef __locale_t locale_t;
# 163 "/usr/include/string.h" 2 3 4


extern int strcoll_l (__const char *__s1, __const char *__s2, __locale_t __l)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2, 3)));

extern size_t strxfrm_l (char *__dest, __const char *__src, size_t __n,
    __locale_t __l) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4)));





extern char *strdup (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));






extern char *strndup (__const char *__string, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));
# 210 "/usr/include/string.h" 3 4

# 235 "/usr/include/string.h" 3 4
extern char *strchr (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 262 "/usr/include/string.h" 3 4
extern char *strrchr (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 281 "/usr/include/string.h" 3 4



extern size_t strcspn (__const char *__s, __const char *__reject)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern size_t strspn (__const char *__s, __const char *__accept)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 314 "/usr/include/string.h" 3 4
extern char *strpbrk (__const char *__s, __const char *__accept)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 342 "/usr/include/string.h" 3 4
extern char *strstr (__const char *__haystack, __const char *__needle)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strtok (char *__restrict __s, __const char *__restrict __delim)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));




extern char *__strtok_r (char *__restrict __s,
    __const char *__restrict __delim,
    char **__restrict __save_ptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 3)));

extern char *strtok_r (char *__restrict __s, __const char *__restrict __delim,
         char **__restrict __save_ptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 3)));
# 397 "/usr/include/string.h" 3 4


extern size_t strlen (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern size_t strnlen (__const char *__string, size_t __maxlen)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern char *strerror (int __errnum) __attribute__ ((__nothrow__));

# 427 "/usr/include/string.h" 3 4
extern int strerror_r (int __errnum, char *__buf, size_t __buflen) __asm__ ("" "__xpg_strerror_r") __attribute__ ((__nothrow__))

                        __attribute__ ((__nonnull__ (2)));
# 445 "/usr/include/string.h" 3 4
extern char *strerror_l (int __errnum, __locale_t __l) __attribute__ ((__nothrow__));





extern void __bzero (void *__s, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));



extern void bcopy (__const void *__src, void *__dest, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern void bzero (void *__s, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int bcmp (__const void *__s1, __const void *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 489 "/usr/include/string.h" 3 4
extern char *index (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 517 "/usr/include/string.h" 3 4
extern char *rindex (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));




extern int ffs (int __i) __attribute__ ((__nothrow__)) __attribute__ ((__const__));
# 536 "/usr/include/string.h" 3 4
extern int strcasecmp (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strncasecmp (__const char *__s1, __const char *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 559 "/usr/include/string.h" 3 4
extern char *strsep (char **__restrict __stringp,
       __const char *__restrict __delim)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strsignal (int __sig) __attribute__ ((__nothrow__));


extern char *__stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));



extern char *__stpncpy (char *__restrict __dest,
   __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpncpy (char *__restrict __dest,
        __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
# 646 "/usr/include/string.h" 3 4

# 14 "portio2a.c" 2
# 1 "/usr/include/sys/stat.h" 1 3 4
# 39 "/usr/include/sys/stat.h" 3 4
# 1 "/usr/include/time.h" 1 3 4
# 40 "/usr/include/sys/stat.h" 2 3 4
# 105 "/usr/include/sys/stat.h" 3 4


# 1 "/usr/include/bits/stat.h" 1 3 4
# 46 "/usr/include/bits/stat.h" 3 4
struct stat
  {
    __dev_t st_dev;




    __ino_t st_ino;







    __nlink_t st_nlink;
    __mode_t st_mode;

    __uid_t st_uid;
    __gid_t st_gid;

    int __pad0;

    __dev_t st_rdev;




    __off_t st_size;



    __blksize_t st_blksize;

    __blkcnt_t st_blocks;
# 91 "/usr/include/bits/stat.h" 3 4
    struct timespec st_atim;
    struct timespec st_mtim;
    struct timespec st_ctim;
# 106 "/usr/include/bits/stat.h" 3 4
    long int __unused[3];
# 115 "/usr/include/bits/stat.h" 3 4
  };
# 108 "/usr/include/sys/stat.h" 2 3 4
# 219 "/usr/include/sys/stat.h" 3 4
extern int stat (__const char *__restrict __file, struct stat *__restrict __buf) __asm__ ("" "stat64") __attribute__ ((__nothrow__))

     __attribute__ ((__nonnull__ (1, 2)));
extern int fstat (int __fd, struct stat *__buf) __asm__ ("" "fstat64") __attribute__ ((__nothrow__))
     __attribute__ ((__nonnull__ (2)));
# 245 "/usr/include/sys/stat.h" 3 4
extern int fstatat (int __fd, __const char *__restrict __file, struct stat *__restrict __buf, int __flag) __asm__ ("" "fstatat64") __attribute__ ((__nothrow__))


                 __attribute__ ((__nonnull__ (2, 3)));
# 269 "/usr/include/sys/stat.h" 3 4
extern int lstat (__const char *__restrict __file, struct stat *__restrict __buf) __asm__ ("" "lstat64") __attribute__ ((__nothrow__))


     __attribute__ ((__nonnull__ (1, 2)));
# 286 "/usr/include/sys/stat.h" 3 4
extern int chmod (__const char *__file, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int lchmod (__const char *__file, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));




extern int fchmod (int __fd, __mode_t __mode) __attribute__ ((__nothrow__));





extern int fchmodat (int __fd, __const char *__file, __mode_t __mode,
       int __flag)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2))) ;






extern __mode_t umask (__mode_t __mask) __attribute__ ((__nothrow__));
# 323 "/usr/include/sys/stat.h" 3 4
extern int mkdir (__const char *__path, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int mkdirat (int __fd, __const char *__path, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));






extern int mknod (__const char *__path, __mode_t __mode, __dev_t __dev)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int mknodat (int __fd, __const char *__path, __mode_t __mode,
      __dev_t __dev) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));





extern int mkfifo (__const char *__path, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int mkfifoat (int __fd, __const char *__path, __mode_t __mode)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));





extern int utimensat (int __fd, __const char *__path,
        __const struct timespec __times[2],
        int __flags)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));




extern int futimens (int __fd, __const struct timespec __times[2]) __attribute__ ((__nothrow__));
# 412 "/usr/include/sys/stat.h" 3 4
extern int __fxstat (int __ver, int __fildes, struct stat *__stat_buf) __asm__ ("" "__fxstat64") __attribute__ ((__nothrow__))

     __attribute__ ((__nonnull__ (3)));
extern int __xstat (int __ver, __const char *__filename, struct stat *__stat_buf) __asm__ ("" "__xstat64") __attribute__ ((__nothrow__))

     __attribute__ ((__nonnull__ (2, 3)));
extern int __lxstat (int __ver, __const char *__filename, struct stat *__stat_buf) __asm__ ("" "__lxstat64") __attribute__ ((__nothrow__))

     __attribute__ ((__nonnull__ (2, 3)));
extern int __fxstatat (int __ver, int __fildes, __const char *__filename, struct stat *__stat_buf, int __flag) __asm__ ("" "__fxstatat64") __attribute__ ((__nothrow__))


                    __attribute__ ((__nonnull__ (3, 4)));
# 444 "/usr/include/sys/stat.h" 3 4
extern int __xmknod (int __ver, __const char *__path, __mode_t __mode,
       __dev_t *__dev) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4)));

extern int __xmknodat (int __ver, int __fd, __const char *__path,
         __mode_t __mode, __dev_t *__dev)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 5)));
# 536 "/usr/include/sys/stat.h" 3 4

# 15 "portio2a.c" 2
# 1 "/usr/include/errno.h" 1 3 4
# 32 "/usr/include/errno.h" 3 4




# 1 "/usr/include/bits/errno.h" 1 3 4
# 25 "/usr/include/bits/errno.h" 3 4
# 1 "/usr/include/linux/errno.h" 1 3 4



# 1 "/usr/include/asm/errno.h" 1 3 4
# 1 "/usr/include/asm-generic/errno.h" 1 3 4



# 1 "/usr/include/asm-generic/errno-base.h" 1 3 4
# 5 "/usr/include/asm-generic/errno.h" 2 3 4
# 1 "/usr/include/asm/errno.h" 2 3 4
# 5 "/usr/include/linux/errno.h" 2 3 4
# 26 "/usr/include/bits/errno.h" 2 3 4
# 47 "/usr/include/bits/errno.h" 3 4
extern int *__errno_location (void) __attribute__ ((__nothrow__)) __attribute__ ((__const__));
# 37 "/usr/include/errno.h" 2 3 4
# 59 "/usr/include/errno.h" 3 4

# 16 "portio2a.c" 2
# 1 "/usr/include/time.h" 1 3 4
# 30 "/usr/include/time.h" 3 4








# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 39 "/usr/include/time.h" 2 3 4



# 1 "/usr/include/bits/time.h" 1 3 4
# 43 "/usr/include/time.h" 2 3 4
# 131 "/usr/include/time.h" 3 4


struct tm
{
  int tm_sec;
  int tm_min;
  int tm_hour;
  int tm_mday;
  int tm_mon;
  int tm_year;
  int tm_wday;
  int tm_yday;
  int tm_isdst;


  long int tm_gmtoff;
  __const char *tm_zone;




};








struct itimerspec
  {
    struct timespec it_interval;
    struct timespec it_value;
  };


struct sigevent;
# 180 "/usr/include/time.h" 3 4



extern clock_t clock (void) __attribute__ ((__nothrow__));


extern time_t time (time_t *__timer) __attribute__ ((__nothrow__));


extern double difftime (time_t __time1, time_t __time0)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern time_t mktime (struct tm *__tp) __attribute__ ((__nothrow__));





extern size_t strftime (char *__restrict __s, size_t __maxsize,
   __const char *__restrict __format,
   __const struct tm *__restrict __tp) __attribute__ ((__nothrow__));

# 217 "/usr/include/time.h" 3 4
extern size_t strftime_l (char *__restrict __s, size_t __maxsize,
     __const char *__restrict __format,
     __const struct tm *__restrict __tp,
     __locale_t __loc) __attribute__ ((__nothrow__));
# 230 "/usr/include/time.h" 3 4



extern struct tm *gmtime (__const time_t *__timer) __attribute__ ((__nothrow__));



extern struct tm *localtime (__const time_t *__timer) __attribute__ ((__nothrow__));





extern struct tm *gmtime_r (__const time_t *__restrict __timer,
       struct tm *__restrict __tp) __attribute__ ((__nothrow__));



extern struct tm *localtime_r (__const time_t *__restrict __timer,
          struct tm *__restrict __tp) __attribute__ ((__nothrow__));





extern char *asctime (__const struct tm *__tp) __attribute__ ((__nothrow__));


extern char *ctime (__const time_t *__timer) __attribute__ ((__nothrow__));







extern char *asctime_r (__const struct tm *__restrict __tp,
   char *__restrict __buf) __attribute__ ((__nothrow__));


extern char *ctime_r (__const time_t *__restrict __timer,
        char *__restrict __buf) __attribute__ ((__nothrow__));




extern char *__tzname[2];
extern int __daylight;
extern long int __timezone;




extern char *tzname[2];



extern void tzset (void) __attribute__ ((__nothrow__));



extern int daylight;
extern long int timezone;





extern int stime (__const time_t *__when) __attribute__ ((__nothrow__));
# 313 "/usr/include/time.h" 3 4
extern time_t timegm (struct tm *__tp) __attribute__ ((__nothrow__));


extern time_t timelocal (struct tm *__tp) __attribute__ ((__nothrow__));


extern int dysize (int __year) __attribute__ ((__nothrow__)) __attribute__ ((__const__));
# 328 "/usr/include/time.h" 3 4
extern int nanosleep (__const struct timespec *__requested_time,
        struct timespec *__remaining);



extern int clock_getres (clockid_t __clock_id, struct timespec *__res) __attribute__ ((__nothrow__));


extern int clock_gettime (clockid_t __clock_id, struct timespec *__tp) __attribute__ ((__nothrow__));


extern int clock_settime (clockid_t __clock_id, __const struct timespec *__tp)
     __attribute__ ((__nothrow__));






extern int clock_nanosleep (clockid_t __clock_id, int __flags,
       __const struct timespec *__req,
       struct timespec *__rem);


extern int clock_getcpuclockid (pid_t __pid, clockid_t *__clock_id) __attribute__ ((__nothrow__));




extern int timer_create (clockid_t __clock_id,
    struct sigevent *__restrict __evp,
    timer_t *__restrict __timerid) __attribute__ ((__nothrow__));


extern int timer_delete (timer_t __timerid) __attribute__ ((__nothrow__));


extern int timer_settime (timer_t __timerid, int __flags,
     __const struct itimerspec *__restrict __value,
     struct itimerspec *__restrict __ovalue) __attribute__ ((__nothrow__));


extern int timer_gettime (timer_t __timerid, struct itimerspec *__value)
     __attribute__ ((__nothrow__));


extern int timer_getoverrun (timer_t __timerid) __attribute__ ((__nothrow__));
# 417 "/usr/include/time.h" 3 4

# 17 "portio2a.c" 2



# 1 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h" 1
# 28 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h"
typedef unsigned long long int u_integer;
typedef long long int integer;
# 43 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h"
typedef double real;
typedef integer um_data_t;
# 21 "portio2a.c" 2


# 1 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_portio.h" 1
# 24 "portio2a.c" 2


# 1 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_pio_timer.h" 1
# 15 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_pio_timer.h"
# 1 "/usr/include/sys/time.h" 1 3 4
# 29 "/usr/include/sys/time.h" 3 4
# 1 "/usr/include/bits/time.h" 1 3 4
# 30 "/usr/include/sys/time.h" 2 3 4
# 39 "/usr/include/sys/time.h" 3 4

# 57 "/usr/include/sys/time.h" 3 4
struct timezone
  {
    int tz_minuteswest;
    int tz_dsttime;
  };

typedef struct timezone *__restrict __timezone_ptr_t;
# 73 "/usr/include/sys/time.h" 3 4
extern int gettimeofday (struct timeval *__restrict __tv,
    __timezone_ptr_t __tz) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));




extern int settimeofday (__const struct timeval *__tv,
    __const struct timezone *__tz)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int adjtime (__const struct timeval *__delta,
      struct timeval *__olddelta) __attribute__ ((__nothrow__));




enum __itimer_which
  {

    ITIMER_REAL = 0,


    ITIMER_VIRTUAL = 1,



    ITIMER_PROF = 2

  };



struct itimerval
  {

    struct timeval it_interval;

    struct timeval it_value;
  };






typedef int __itimer_which_t;




extern int getitimer (__itimer_which_t __which,
        struct itimerval *__value) __attribute__ ((__nothrow__));




extern int setitimer (__itimer_which_t __which,
        __const struct itimerval *__restrict __new,
        struct itimerval *__restrict __old) __attribute__ ((__nothrow__));




extern int utimes (__const char *__file, __const struct timeval __tvp[2])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));



extern int lutimes (__const char *__file, __const struct timeval __tvp[2])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int futimes (int __fd, __const struct timeval __tvp[2]) __attribute__ ((__nothrow__));
# 191 "/usr/include/sys/time.h" 3 4

# 16 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_pio_timer.h" 2
# 1 "/usr/include/unistd.h" 1 3 4
# 28 "/usr/include/unistd.h" 3 4

# 203 "/usr/include/unistd.h" 3 4
# 1 "/usr/include/bits/posix_opt.h" 1 3 4
# 204 "/usr/include/unistd.h" 2 3 4



# 1 "/usr/include/bits/environments.h" 1 3 4
# 23 "/usr/include/bits/environments.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 24 "/usr/include/bits/environments.h" 2 3 4
# 208 "/usr/include/unistd.h" 2 3 4
# 227 "/usr/include/unistd.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 228 "/usr/include/unistd.h" 2 3 4
# 256 "/usr/include/unistd.h" 3 4
typedef __useconds_t useconds_t;
# 268 "/usr/include/unistd.h" 3 4
typedef __intptr_t intptr_t;






typedef __socklen_t socklen_t;
# 288 "/usr/include/unistd.h" 3 4
extern int access (__const char *__name, int __type) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 305 "/usr/include/unistd.h" 3 4
extern int faccessat (int __fd, __const char *__file, int __type, int __flag)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2))) ;
# 334 "/usr/include/unistd.h" 3 4
extern __off64_t lseek (int __fd, __off64_t __offset, int __whence) __asm__ ("" "lseek64") __attribute__ ((__nothrow__))

             ;
# 350 "/usr/include/unistd.h" 3 4
extern int close (int __fd);






extern ssize_t read (int __fd, void *__buf, size_t __nbytes) ;





extern ssize_t write (int __fd, __const void *__buf, size_t __n) ;
# 385 "/usr/include/unistd.h" 3 4
extern ssize_t pread (int __fd, void *__buf, size_t __nbytes, __off64_t __offset) __asm__ ("" "pread64")

                    ;
extern ssize_t pwrite (int __fd, __const void *__buf, size_t __nbytes, __off64_t __offset) __asm__ ("" "pwrite64")

                     ;
# 414 "/usr/include/unistd.h" 3 4
extern int pipe (int __pipedes[2]) __attribute__ ((__nothrow__)) ;
# 429 "/usr/include/unistd.h" 3 4
extern unsigned int alarm (unsigned int __seconds) __attribute__ ((__nothrow__));
# 441 "/usr/include/unistd.h" 3 4
extern unsigned int sleep (unsigned int __seconds);







extern __useconds_t ualarm (__useconds_t __value, __useconds_t __interval)
     __attribute__ ((__nothrow__));






extern int usleep (__useconds_t __useconds);
# 466 "/usr/include/unistd.h" 3 4
extern int pause (void);



extern int chown (__const char *__file, __uid_t __owner, __gid_t __group)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;



extern int fchown (int __fd, __uid_t __owner, __gid_t __group) __attribute__ ((__nothrow__)) ;




extern int lchown (__const char *__file, __uid_t __owner, __gid_t __group)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;






extern int fchownat (int __fd, __const char *__file, __uid_t __owner,
       __gid_t __group, int __flag)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2))) ;



extern int chdir (__const char *__path) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;



extern int fchdir (int __fd) __attribute__ ((__nothrow__)) ;
# 508 "/usr/include/unistd.h" 3 4
extern char *getcwd (char *__buf, size_t __size) __attribute__ ((__nothrow__)) ;
# 522 "/usr/include/unistd.h" 3 4
extern char *getwd (char *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) __attribute__ ((__deprecated__)) ;




extern int dup (int __fd) __attribute__ ((__nothrow__)) ;


extern int dup2 (int __fd, int __fd2) __attribute__ ((__nothrow__));
# 540 "/usr/include/unistd.h" 3 4
extern char **__environ;







extern int execve (__const char *__path, char *__const __argv[],
     char *__const __envp[]) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));




extern int fexecve (int __fd, char *__const __argv[], char *__const __envp[])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));




extern int execv (__const char *__path, char *__const __argv[])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execle (__const char *__path, __const char *__arg, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execl (__const char *__path, __const char *__arg, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execvp (__const char *__file, char *__const __argv[])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));




extern int execlp (__const char *__file, __const char *__arg, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
# 595 "/usr/include/unistd.h" 3 4
extern int nice (int __inc) __attribute__ ((__nothrow__)) ;




extern void _exit (int __status) __attribute__ ((__noreturn__));





# 1 "/usr/include/bits/confname.h" 1 3 4
# 26 "/usr/include/bits/confname.h" 3 4
enum
  {
    _PC_LINK_MAX,

    _PC_MAX_CANON,

    _PC_MAX_INPUT,

    _PC_NAME_MAX,

    _PC_PATH_MAX,

    _PC_PIPE_BUF,

    _PC_CHOWN_RESTRICTED,

    _PC_NO_TRUNC,

    _PC_VDISABLE,

    _PC_SYNC_IO,

    _PC_ASYNC_IO,

    _PC_PRIO_IO,

    _PC_SOCK_MAXBUF,

    _PC_FILESIZEBITS,

    _PC_REC_INCR_XFER_SIZE,

    _PC_REC_MAX_XFER_SIZE,

    _PC_REC_MIN_XFER_SIZE,

    _PC_REC_XFER_ALIGN,

    _PC_ALLOC_SIZE_MIN,

    _PC_SYMLINK_MAX,

    _PC_2_SYMLINKS

  };


enum
  {
    _SC_ARG_MAX,

    _SC_CHILD_MAX,

    _SC_CLK_TCK,

    _SC_NGROUPS_MAX,

    _SC_OPEN_MAX,

    _SC_STREAM_MAX,

    _SC_TZNAME_MAX,

    _SC_JOB_CONTROL,

    _SC_SAVED_IDS,

    _SC_REALTIME_SIGNALS,

    _SC_PRIORITY_SCHEDULING,

    _SC_TIMERS,

    _SC_ASYNCHRONOUS_IO,

    _SC_PRIORITIZED_IO,

    _SC_SYNCHRONIZED_IO,

    _SC_FSYNC,

    _SC_MAPPED_FILES,

    _SC_MEMLOCK,

    _SC_MEMLOCK_RANGE,

    _SC_MEMORY_PROTECTION,

    _SC_MESSAGE_PASSING,

    _SC_SEMAPHORES,

    _SC_SHARED_MEMORY_OBJECTS,

    _SC_AIO_LISTIO_MAX,

    _SC_AIO_MAX,

    _SC_AIO_PRIO_DELTA_MAX,

    _SC_DELAYTIMER_MAX,

    _SC_MQ_OPEN_MAX,

    _SC_MQ_PRIO_MAX,

    _SC_VERSION,

    _SC_PAGESIZE,


    _SC_RTSIG_MAX,

    _SC_SEM_NSEMS_MAX,

    _SC_SEM_VALUE_MAX,

    _SC_SIGQUEUE_MAX,

    _SC_TIMER_MAX,




    _SC_BC_BASE_MAX,

    _SC_BC_DIM_MAX,

    _SC_BC_SCALE_MAX,

    _SC_BC_STRING_MAX,

    _SC_COLL_WEIGHTS_MAX,

    _SC_EQUIV_CLASS_MAX,

    _SC_EXPR_NEST_MAX,

    _SC_LINE_MAX,

    _SC_RE_DUP_MAX,

    _SC_CHARCLASS_NAME_MAX,


    _SC_2_VERSION,

    _SC_2_C_BIND,

    _SC_2_C_DEV,

    _SC_2_FORT_DEV,

    _SC_2_FORT_RUN,

    _SC_2_SW_DEV,

    _SC_2_LOCALEDEF,


    _SC_PII,

    _SC_PII_XTI,

    _SC_PII_SOCKET,

    _SC_PII_INTERNET,

    _SC_PII_OSI,

    _SC_POLL,

    _SC_SELECT,

    _SC_UIO_MAXIOV,

    _SC_IOV_MAX = _SC_UIO_MAXIOV,

    _SC_PII_INTERNET_STREAM,

    _SC_PII_INTERNET_DGRAM,

    _SC_PII_OSI_COTS,

    _SC_PII_OSI_CLTS,

    _SC_PII_OSI_M,

    _SC_T_IOV_MAX,



    _SC_THREADS,

    _SC_THREAD_SAFE_FUNCTIONS,

    _SC_GETGR_R_SIZE_MAX,

    _SC_GETPW_R_SIZE_MAX,

    _SC_LOGIN_NAME_MAX,

    _SC_TTY_NAME_MAX,

    _SC_THREAD_DESTRUCTOR_ITERATIONS,

    _SC_THREAD_KEYS_MAX,

    _SC_THREAD_STACK_MIN,

    _SC_THREAD_THREADS_MAX,

    _SC_THREAD_ATTR_STACKADDR,

    _SC_THREAD_ATTR_STACKSIZE,

    _SC_THREAD_PRIORITY_SCHEDULING,

    _SC_THREAD_PRIO_INHERIT,

    _SC_THREAD_PRIO_PROTECT,

    _SC_THREAD_PROCESS_SHARED,


    _SC_NPROCESSORS_CONF,

    _SC_NPROCESSORS_ONLN,

    _SC_PHYS_PAGES,

    _SC_AVPHYS_PAGES,

    _SC_ATEXIT_MAX,

    _SC_PASS_MAX,


    _SC_XOPEN_VERSION,

    _SC_XOPEN_XCU_VERSION,

    _SC_XOPEN_UNIX,

    _SC_XOPEN_CRYPT,

    _SC_XOPEN_ENH_I18N,

    _SC_XOPEN_SHM,


    _SC_2_CHAR_TERM,

    _SC_2_C_VERSION,

    _SC_2_UPE,


    _SC_XOPEN_XPG2,

    _SC_XOPEN_XPG3,

    _SC_XOPEN_XPG4,


    _SC_CHAR_BIT,

    _SC_CHAR_MAX,

    _SC_CHAR_MIN,

    _SC_INT_MAX,

    _SC_INT_MIN,

    _SC_LONG_BIT,

    _SC_WORD_BIT,

    _SC_MB_LEN_MAX,

    _SC_NZERO,

    _SC_SSIZE_MAX,

    _SC_SCHAR_MAX,

    _SC_SCHAR_MIN,

    _SC_SHRT_MAX,

    _SC_SHRT_MIN,

    _SC_UCHAR_MAX,

    _SC_UINT_MAX,

    _SC_ULONG_MAX,

    _SC_USHRT_MAX,


    _SC_NL_ARGMAX,

    _SC_NL_LANGMAX,

    _SC_NL_MSGMAX,

    _SC_NL_NMAX,

    _SC_NL_SETMAX,

    _SC_NL_TEXTMAX,


    _SC_XBS5_ILP32_OFF32,

    _SC_XBS5_ILP32_OFFBIG,

    _SC_XBS5_LP64_OFF64,

    _SC_XBS5_LPBIG_OFFBIG,


    _SC_XOPEN_LEGACY,

    _SC_XOPEN_REALTIME,

    _SC_XOPEN_REALTIME_THREADS,


    _SC_ADVISORY_INFO,

    _SC_BARRIERS,

    _SC_BASE,

    _SC_C_LANG_SUPPORT,

    _SC_C_LANG_SUPPORT_R,

    _SC_CLOCK_SELECTION,

    _SC_CPUTIME,

    _SC_THREAD_CPUTIME,

    _SC_DEVICE_IO,

    _SC_DEVICE_SPECIFIC,

    _SC_DEVICE_SPECIFIC_R,

    _SC_FD_MGMT,

    _SC_FIFO,

    _SC_PIPE,

    _SC_FILE_ATTRIBUTES,

    _SC_FILE_LOCKING,

    _SC_FILE_SYSTEM,

    _SC_MONOTONIC_CLOCK,

    _SC_MULTI_PROCESS,

    _SC_SINGLE_PROCESS,

    _SC_NETWORKING,

    _SC_READER_WRITER_LOCKS,

    _SC_SPIN_LOCKS,

    _SC_REGEXP,

    _SC_REGEX_VERSION,

    _SC_SHELL,

    _SC_SIGNALS,

    _SC_SPAWN,

    _SC_SPORADIC_SERVER,

    _SC_THREAD_SPORADIC_SERVER,

    _SC_SYSTEM_DATABASE,

    _SC_SYSTEM_DATABASE_R,

    _SC_TIMEOUTS,

    _SC_TYPED_MEMORY_OBJECTS,

    _SC_USER_GROUPS,

    _SC_USER_GROUPS_R,

    _SC_2_PBS,

    _SC_2_PBS_ACCOUNTING,

    _SC_2_PBS_LOCATE,

    _SC_2_PBS_MESSAGE,

    _SC_2_PBS_TRACK,

    _SC_SYMLOOP_MAX,

    _SC_STREAMS,

    _SC_2_PBS_CHECKPOINT,


    _SC_V6_ILP32_OFF32,

    _SC_V6_ILP32_OFFBIG,

    _SC_V6_LP64_OFF64,

    _SC_V6_LPBIG_OFFBIG,


    _SC_HOST_NAME_MAX,

    _SC_TRACE,

    _SC_TRACE_EVENT_FILTER,

    _SC_TRACE_INHERIT,

    _SC_TRACE_LOG,


    _SC_LEVEL1_ICACHE_SIZE,

    _SC_LEVEL1_ICACHE_ASSOC,

    _SC_LEVEL1_ICACHE_LINESIZE,

    _SC_LEVEL1_DCACHE_SIZE,

    _SC_LEVEL1_DCACHE_ASSOC,

    _SC_LEVEL1_DCACHE_LINESIZE,

    _SC_LEVEL2_CACHE_SIZE,

    _SC_LEVEL2_CACHE_ASSOC,

    _SC_LEVEL2_CACHE_LINESIZE,

    _SC_LEVEL3_CACHE_SIZE,

    _SC_LEVEL3_CACHE_ASSOC,

    _SC_LEVEL3_CACHE_LINESIZE,

    _SC_LEVEL4_CACHE_SIZE,

    _SC_LEVEL4_CACHE_ASSOC,

    _SC_LEVEL4_CACHE_LINESIZE,



    _SC_IPV6 = _SC_LEVEL1_ICACHE_SIZE + 50,

    _SC_RAW_SOCKETS,


    _SC_V7_ILP32_OFF32,

    _SC_V7_ILP32_OFFBIG,

    _SC_V7_LP64_OFF64,

    _SC_V7_LPBIG_OFFBIG,


    _SC_SS_REPL_MAX,


    _SC_TRACE_EVENT_NAME_MAX,

    _SC_TRACE_NAME_MAX,

    _SC_TRACE_SYS_MAX,

    _SC_TRACE_USER_EVENT_MAX,


    _SC_XOPEN_STREAMS,


    _SC_THREAD_ROBUST_PRIO_INHERIT,

    _SC_THREAD_ROBUST_PRIO_PROTECT

  };


enum
  {
    _CS_PATH,


    _CS_V6_WIDTH_RESTRICTED_ENVS,



    _CS_GNU_LIBC_VERSION,

    _CS_GNU_LIBPTHREAD_VERSION,


    _CS_V5_WIDTH_RESTRICTED_ENVS,



    _CS_V7_WIDTH_RESTRICTED_ENVS,



    _CS_LFS_CFLAGS = 1000,

    _CS_LFS_LDFLAGS,

    _CS_LFS_LIBS,

    _CS_LFS_LINTFLAGS,

    _CS_LFS64_CFLAGS,

    _CS_LFS64_LDFLAGS,

    _CS_LFS64_LIBS,

    _CS_LFS64_LINTFLAGS,


    _CS_XBS5_ILP32_OFF32_CFLAGS = 1100,

    _CS_XBS5_ILP32_OFF32_LDFLAGS,

    _CS_XBS5_ILP32_OFF32_LIBS,

    _CS_XBS5_ILP32_OFF32_LINTFLAGS,

    _CS_XBS5_ILP32_OFFBIG_CFLAGS,

    _CS_XBS5_ILP32_OFFBIG_LDFLAGS,

    _CS_XBS5_ILP32_OFFBIG_LIBS,

    _CS_XBS5_ILP32_OFFBIG_LINTFLAGS,

    _CS_XBS5_LP64_OFF64_CFLAGS,

    _CS_XBS5_LP64_OFF64_LDFLAGS,

    _CS_XBS5_LP64_OFF64_LIBS,

    _CS_XBS5_LP64_OFF64_LINTFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_CFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_LDFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_LIBS,

    _CS_XBS5_LPBIG_OFFBIG_LINTFLAGS,


    _CS_POSIX_V6_ILP32_OFF32_CFLAGS,

    _CS_POSIX_V6_ILP32_OFF32_LDFLAGS,

    _CS_POSIX_V6_ILP32_OFF32_LIBS,

    _CS_POSIX_V6_ILP32_OFF32_LINTFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_CFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_LDFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_LIBS,

    _CS_POSIX_V6_ILP32_OFFBIG_LINTFLAGS,

    _CS_POSIX_V6_LP64_OFF64_CFLAGS,

    _CS_POSIX_V6_LP64_OFF64_LDFLAGS,

    _CS_POSIX_V6_LP64_OFF64_LIBS,

    _CS_POSIX_V6_LP64_OFF64_LINTFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_CFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LDFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LIBS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LINTFLAGS,


    _CS_POSIX_V7_ILP32_OFF32_CFLAGS,

    _CS_POSIX_V7_ILP32_OFF32_LDFLAGS,

    _CS_POSIX_V7_ILP32_OFF32_LIBS,

    _CS_POSIX_V7_ILP32_OFF32_LINTFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_CFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_LDFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_LIBS,

    _CS_POSIX_V7_ILP32_OFFBIG_LINTFLAGS,

    _CS_POSIX_V7_LP64_OFF64_CFLAGS,

    _CS_POSIX_V7_LP64_OFF64_LDFLAGS,

    _CS_POSIX_V7_LP64_OFF64_LIBS,

    _CS_POSIX_V7_LP64_OFF64_LINTFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_CFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LDFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LIBS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LINTFLAGS,


    _CS_V6_ENV,

    _CS_V7_ENV

  };
# 607 "/usr/include/unistd.h" 2 3 4


extern long int pathconf (__const char *__path, int __name)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern long int fpathconf (int __fd, int __name) __attribute__ ((__nothrow__));


extern long int sysconf (int __name) __attribute__ ((__nothrow__));



extern size_t confstr (int __name, char *__buf, size_t __len) __attribute__ ((__nothrow__));




extern __pid_t getpid (void) __attribute__ ((__nothrow__));


extern __pid_t getppid (void) __attribute__ ((__nothrow__));




extern __pid_t getpgrp (void) __attribute__ ((__nothrow__));
# 643 "/usr/include/unistd.h" 3 4
extern __pid_t __getpgid (__pid_t __pid) __attribute__ ((__nothrow__));

extern __pid_t getpgid (__pid_t __pid) __attribute__ ((__nothrow__));






extern int setpgid (__pid_t __pid, __pid_t __pgid) __attribute__ ((__nothrow__));
# 669 "/usr/include/unistd.h" 3 4
extern int setpgrp (void) __attribute__ ((__nothrow__));
# 686 "/usr/include/unistd.h" 3 4
extern __pid_t setsid (void) __attribute__ ((__nothrow__));



extern __pid_t getsid (__pid_t __pid) __attribute__ ((__nothrow__));



extern __uid_t getuid (void) __attribute__ ((__nothrow__));


extern __uid_t geteuid (void) __attribute__ ((__nothrow__));


extern __gid_t getgid (void) __attribute__ ((__nothrow__));


extern __gid_t getegid (void) __attribute__ ((__nothrow__));




extern int getgroups (int __size, __gid_t __list[]) __attribute__ ((__nothrow__)) ;
# 719 "/usr/include/unistd.h" 3 4
extern int setuid (__uid_t __uid) __attribute__ ((__nothrow__));




extern int setreuid (__uid_t __ruid, __uid_t __euid) __attribute__ ((__nothrow__));




extern int seteuid (__uid_t __uid) __attribute__ ((__nothrow__));






extern int setgid (__gid_t __gid) __attribute__ ((__nothrow__));




extern int setregid (__gid_t __rgid, __gid_t __egid) __attribute__ ((__nothrow__));




extern int setegid (__gid_t __gid) __attribute__ ((__nothrow__));
# 775 "/usr/include/unistd.h" 3 4
extern __pid_t fork (void) __attribute__ ((__nothrow__));







extern __pid_t vfork (void) __attribute__ ((__nothrow__));





extern char *ttyname (int __fd) __attribute__ ((__nothrow__));



extern int ttyname_r (int __fd, char *__buf, size_t __buflen)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2))) ;



extern int isatty (int __fd) __attribute__ ((__nothrow__));





extern int ttyslot (void) __attribute__ ((__nothrow__));




extern int link (__const char *__from, __const char *__to)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern int linkat (int __fromfd, __const char *__from, int __tofd,
     __const char *__to, int __flags)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4))) ;




extern int symlink (__const char *__from, __const char *__to)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern ssize_t readlink (__const char *__restrict __path,
    char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern int symlinkat (__const char *__from, int __tofd,
        __const char *__to) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 3))) ;


extern ssize_t readlinkat (int __fd, __const char *__restrict __path,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 3))) ;



extern int unlink (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));



extern int unlinkat (int __fd, __const char *__name, int __flag)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));



extern int rmdir (__const char *__path) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));



extern __pid_t tcgetpgrp (int __fd) __attribute__ ((__nothrow__));


extern int tcsetpgrp (int __fd, __pid_t __pgrp_id) __attribute__ ((__nothrow__));






extern char *getlogin (void);







extern int getlogin_r (char *__name, size_t __name_len) __attribute__ ((__nonnull__ (1)));




extern int setlogin (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 890 "/usr/include/unistd.h" 3 4
# 1 "/usr/include/getopt.h" 1 3 4
# 59 "/usr/include/getopt.h" 3 4
extern char *optarg;
# 73 "/usr/include/getopt.h" 3 4
extern int optind;




extern int opterr;



extern int optopt;
# 152 "/usr/include/getopt.h" 3 4
extern int getopt (int ___argc, char *const *___argv, const char *__shortopts)
       __attribute__ ((__nothrow__));
# 891 "/usr/include/unistd.h" 2 3 4







extern int gethostname (char *__name, size_t __len) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern int sethostname (__const char *__name, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;



extern int sethostid (long int __id) __attribute__ ((__nothrow__)) ;





extern int getdomainname (char *__name, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
extern int setdomainname (__const char *__name, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern int vhangup (void) __attribute__ ((__nothrow__));


extern int revoke (__const char *__file) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;







extern int profil (unsigned short int *__sample_buffer, size_t __size,
     size_t __offset, unsigned int __scale)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int acct (__const char *__name) __attribute__ ((__nothrow__));



extern char *getusershell (void) __attribute__ ((__nothrow__));
extern void endusershell (void) __attribute__ ((__nothrow__));
extern void setusershell (void) __attribute__ ((__nothrow__));





extern int daemon (int __nochdir, int __noclose) __attribute__ ((__nothrow__)) ;






extern int chroot (__const char *__path) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;



extern char *getpass (__const char *__prompt) __attribute__ ((__nonnull__ (1)));
# 976 "/usr/include/unistd.h" 3 4
extern int fsync (int __fd);






extern long int gethostid (void);


extern void sync (void) __attribute__ ((__nothrow__));





extern int getpagesize (void) __attribute__ ((__nothrow__)) __attribute__ ((__const__));




extern int getdtablesize (void) __attribute__ ((__nothrow__));
# 1011 "/usr/include/unistd.h" 3 4
extern int truncate (__const char *__file, __off64_t __length) __asm__ ("" "truncate64") __attribute__ ((__nothrow__))

                  __attribute__ ((__nonnull__ (1))) ;
# 1029 "/usr/include/unistd.h" 3 4
extern int ftruncate (int __fd, __off64_t __length) __asm__ ("" "ftruncate64") __attribute__ ((__nothrow__))
                        ;
# 1047 "/usr/include/unistd.h" 3 4
extern int brk (void *__addr) __attribute__ ((__nothrow__)) ;





extern void *sbrk (intptr_t __delta) __attribute__ ((__nothrow__));
# 1068 "/usr/include/unistd.h" 3 4
extern long int syscall (long int __sysno, ...) __attribute__ ((__nothrow__));
# 1094 "/usr/include/unistd.h" 3 4
extern int lockf (int __fd, int __cmd, __off64_t __len) __asm__ ("" "lockf64")
                       ;
# 1122 "/usr/include/unistd.h" 3 4
extern int fdatasync (int __fildes);
# 1151 "/usr/include/unistd.h" 3 4
extern char *ctermid (char *__s) __attribute__ ((__nothrow__));
# 1160 "/usr/include/unistd.h" 3 4

# 17 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_pio_timer.h" 2
# 27 "portio2a.c" 2





# 1 "/usr/include/byteswap.h" 1 3 4
# 23 "/usr/include/byteswap.h" 3 4
# 1 "/usr/include/bits/byteswap.h" 1 3 4
# 24 "/usr/include/byteswap.h" 2 3 4
# 33 "portio2a.c" 2
# 53 "portio2a.c"
extern u_integer io_timer_active ;
void io_update_timer( struct timeval, struct timeval,
      u_integer, integer );
void change_endian(void *, int, int);
void output_buffer( integer *, um_data_t *, integer *, integer *,
                    integer *, real *);
void sync_to_cache( integer *);



void flush_unit_buffer( integer * );
void use_buffer( integer *, um_data_t *, integer *, integer *, real *);


static char message[ 384 ];



FILE *pf[300]=
                     {((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                     };

enum fileop { UM_READ, UM_WRITE, UM_UNSET } prev_fileop [300] ;


static integer io_position[300];


# 1 "/usr/include/malloc.h" 1 3 4
# 25 "/usr/include/malloc.h" 3 4
# 1 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 1 3 4
# 149 "/plush/dugong/usr/bin/../lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 26 "/usr/include/malloc.h" 2 3 4
# 48 "/usr/include/malloc.h" 3 4



extern void *malloc (size_t __size) __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;


extern void *calloc (size_t __nmemb, size_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__malloc__)) ;






extern void *realloc (void *__ptr, size_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__warn_unused_result__));


extern void free (void *__ptr) __attribute__ ((__nothrow__));


extern void cfree (void *__ptr) __attribute__ ((__nothrow__));


extern void *memalign (size_t __alignment, size_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__malloc__)) ;


extern void *valloc (size_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__malloc__)) ;



extern void * pvalloc (size_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__malloc__)) ;



extern void *(*__morecore) (ptrdiff_t __size);


extern void *__default_morecore (ptrdiff_t __size) __attribute__ ((__nothrow__))
       __attribute__ ((__malloc__));



struct mallinfo {
  int arena;
  int ordblks;
  int smblks;
  int hblks;
  int hblkhd;
  int usmblks;
  int fsmblks;
  int uordblks;
  int fordblks;
  int keepcost;
};


extern struct mallinfo mallinfo (void) __attribute__ ((__nothrow__));
# 135 "/usr/include/malloc.h" 3 4
extern int mallopt (int __param, int __val) __attribute__ ((__nothrow__));



extern int malloc_trim (size_t __pad) __attribute__ ((__nothrow__));



extern size_t malloc_usable_size (void *__ptr) __attribute__ ((__nothrow__));


extern void malloc_stats (void) __attribute__ ((__nothrow__));


extern int malloc_info (int __options, FILE *__fp);


extern void *malloc_get_state (void) __attribute__ ((__nothrow__));



extern int malloc_set_state (void *__ptr) __attribute__ ((__nothrow__));




extern void (*__malloc_initialize_hook) (void);

extern void (*__free_hook) (void *__ptr, __const void *)
                             ;
extern void *(*__malloc_hook) (size_t __size, __const void *)
                                  ;
extern void *(*__realloc_hook) (void *__ptr, size_t __size, __const void *)
                                   ;
extern void *(*__memalign_hook) (size_t __alignment, size_t __size, __const void *)

                                    ;
extern void (*__after_morecore_hook) (void);


extern void __malloc_check_init (void) __attribute__ ((__nothrow__));



# 111 "portio2a.c" 2



static um_data_t *unit_buffer[300]=
                     {((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                      ((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),((void *)0),
                     };



static integer unit_intent[300], unit_offset[300];



static char *real_add;





integer buffer_size=524288;

static integer readonly = 0;

int open_flag[300] =
     {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
# 187 "portio2a.c"
integer *the_unit;



integer printstatus = 0;



void



set_printstatus_



(printstatus_f)
integer *printstatus_f;
{
printstatus = *printstatus_f;
}


void
# 223 "portio2a.c"
buffin64_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
struct timeval first;
struct timeval second;
struct timezone tzp;
int k;
the_unit=unit;

if(open_flag[*unit]== 0){



    if(prev_fileop[*unit] == UM_WRITE)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_READ;

    if (io_timer_active){
      gettimeofday( &first, &tzp );
    }

    *length = fread(array,sizeof(um_data_t),*maxlen,pf[*unit]);

    if (io_timer_active){
      gettimeofday( &second, &tzp );
      io_update_timer( first, second, sizeof(um_data_t) * (*length),
                       0 );
    }


change_endian(array, sizeof(integer), *maxlen);


*status=-1.0;
k=feof(pf[*unit]);

    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      snprintf(message, 384 ,
       "BUFFIN: C I/O Error - Return code = %d", k);
      fprintf(stdout, "%s\n", message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      snprintf(message, 384 ,
       "BUFFIN: C I/O Error - Return code = %d", k);
      fprintf(stdout, "%s\n", message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

    io_position[*unit]=io_position[*unit]+(*length * sizeof(um_data_t));
}


void
# 308 "portio2a.c"
buffout64_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
  integer elsize;
  the_unit=unit;

  elsize = sizeof(um_data_t);

  if(open_flag[*unit]== 0){


    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_WRITE;



    change_endian(array, sizeof(integer), *maxlen);



    use_buffer(unit, array, maxlen, length, status);





    change_endian(array, sizeof(integer), *maxlen);

  }
  else
    *status=3.0;
}


void use_buffer
(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;

{
  integer record_length;
  integer record_offset;
  integer i;

  um_data_t *buffer_address;



  record_length=*maxlen;
  record_offset=0;
  buffer_address=unit_buffer[*unit];



  while(unit_offset[*unit]+record_length > buffer_size) {





    for(i=unit_offset[*unit]; i<buffer_size; i++) {
      buffer_address[i]=array[record_offset];
      record_length=record_length-1;
      record_offset=record_offset+1;
    }



    unit_offset[*unit]=buffer_size;
    flush_unit_buffer(unit);
  }






  for(i=0; i<record_length; i++) {
    buffer_address[i+unit_offset[*unit]]=array[record_offset];
    record_offset=record_offset+1;
  }




  unit_offset[*unit]=unit_offset[*unit]+record_length;

  *length=*maxlen;
  *status=-1.0;

}


void output_buffer(unit, array, maxlen, length, elsize, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
integer *elsize;
real *status;
{
int k;
struct timeval first;
struct timeval second;
struct timezone tzp;

the_unit=unit;


    if (io_timer_active){
      gettimeofday( &first, &tzp );
    }

    *length = fwrite(array,*elsize,*maxlen,pf[*unit]);

    if (io_timer_active){
      gettimeofday( &second, &tzp);
      io_update_timer( first, second, (*elsize) * (*length),
                       1);
    }

    *status=-1.0;
    k=feof(pf[*unit]);

    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      snprintf(message, 384 ,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      fprintf(stdout, "%s\n", message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      snprintf(message, 384 ,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      fprintf(stdout, "%s\n", message);
      *status=1.0;
    }


    unit_offset[*unit]=0;


    io_position[*unit]=io_position[*unit]+((*length) * (*elsize));
}


void
# 486 "portio2a.c"
setpos_single_




(unit, word_address,err)
integer *unit;
integer *word_address;
integer *err;

{
  int k;

  off_t byte_address;




  the_unit=unit;
  *err = 0;




  if(unit_intent[*unit] != readonly) {



    if(io_position[*unit]+unit_offset[*unit]*sizeof(um_data_t)
        == *word_address*sizeof(um_data_t)) {


      return;

    }


      flush_unit_buffer(unit);
}



    if (io_position[*unit] != *word_address*sizeof(um_data_t)) {
      byte_address=(*word_address)*sizeof(um_data_t);

      k = fseeko(pf[*unit],byte_address,0);



      if(k!=0){
        perror("\nSETPOS: Seek Failed");
        fprintf(stderr,
         "SETPOS: Unit %d to Word Address %d Failed with Error Code %d",
         (int) *unit, (int) *word_address, k);
         *err = 1;
         abort();
      }
    io_position[*unit]=byte_address;
    }

}


void flush_unit_buffer(integer *unit)
{
  integer k;
  integer buffer_length;
  integer len;
  integer elsize;

  integer *buffer_length_address;
  integer *length;

  real stat;
  real minus_1=-1.0;

  um_data_t *buffer_address;
  real *status;

  elsize = sizeof(um_data_t);


  if(unit_offset[*unit] > 0) {

    buffer_address=unit_buffer[*unit];
    buffer_length=unit_offset[*unit];
    k=buffer_length;
    buffer_length_address=&buffer_length;
    length=&len;
    status=&stat;
    output_buffer(unit, buffer_address,
                  buffer_length_address, length, &elsize, status);



    if(*status != minus_1 || *length != k) {

      fprintf(stderr,
       "\nFLUSH_UNIT_BUFFER: Error Flushing Buffered Data on PE %d\n",
       0);
      fprintf(stderr, "FLUSH_UNIT_BUFFER: Status is %7.1f\n", *status);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length Requested was %lld\n", (long long) k);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length written   was %lld\n\n", (long long) *length);
      abort();
    }
    else {
      unit_offset[*unit]=0;
    }
  }
}



void
# 614 "portio2a.c"
open_single_




(unit,file_name, char_len,intent,environ_var_flag,err)
char file_name[];
integer *unit;

integer *char_len;

integer *intent ;

integer *environ_var_flag;


integer *err;


{
   char *fname;
   struct timeval first;
   struct timeval second;
   struct timezone tzp;
   char *gname;
   integer k;
   int malloc_error;


   enum filestat { old, new };
   enum filestat filestatus;


   fname = calloc(*char_len + 1,1);
   the_unit=unit;
   pf[*unit] = ((void *)0);


   strncpy( fname, file_name, *char_len );
   fname[ *char_len ] = '\0';
   sscanf( fname, "%s", fname );

   if (io_timer_active){
     gettimeofday( &first, &tzp );
   }

   if ( *environ_var_flag == 0 )
   { gname = getenv( fname );
      if ( gname == ((void *)0) ) {
        snprintf(message, 384 ,
         "OPEN:  WARNING: Environment variable %s not set",
         fname);
        fprintf(stdout, "%s\n", message);
        open_flag[*unit]=1;
        *err=1;
        free (fname);
        return;
      }
   }
   else
      gname = fname;




   if ( access( gname, 0 ) == 0 ) {

      if (printstatus >= 2) {
      snprintf(message, 384 ,
       "OPEN:  File %s to be Opened on Unit %d Exists",
       gname, (int) *unit );
      fprintf(stdout, "%s\n", message);
      }
      filestatus = old;
   }
   else {

      if (printstatus >= 2) {
      snprintf(message, 384 ,
       "OPEN:  File %s to be Opened on Unit %d does not Exist",
       gname, (int) *unit);
      fprintf(stdout, "%s\n", message);
      }
      filestatus = new;
   }


   if ( filestatus == old ) {

      if ( *intent == readonly ) {

         if ( ( pf[*unit] = fopen( gname, "rb" ) ) == ((void *)0) ) {
            perror("OPEN:  File Open Failed");
            snprintf(message, 384 ,
              "OPEN:  Unable to Open File %s for Reading", gname );
            fprintf(stdout, "%s\n", message);
         }
      }
      else {

         if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == ((void *)0) ) {
            perror("OPEN:  File Open Failed");
            snprintf(message, 384 ,
              "OPEN:  Unable to Open File %s for Read/Write", gname );
            fprintf(stdout, "%s\n", message);
         }
      }
   }



   if ( filestatus == new ) {


      pf[*unit] = ((void *)0);

      if ( *intent == readonly ) {
         snprintf(message, 384 , "OPEN:  **WARNING: FILE NOT FOUND" );
         fprintf(stdout, "%s\n", message);
         snprintf(message, 384 ,
          "OPEN:  Ignored Request to Open File %s for Reading",
           gname );
         fprintf(stdout, "%s\n", message);
      }
      else {


          if ( ( pf[*unit] = fopen( gname, "w+b" ) ) == ((void *)0) ) {
            perror("OPEN:  File Creation Failed");
            snprintf(message, 384 ,
             "OPEN:  Unable to Open File %s for Read/Write", gname );
            fprintf(stdout, "%s\n", message);
          }
          else {

            snprintf(message, 384 , "OPEN:  File %s Created on Unit %d",
                     gname, (int) *unit );
            fprintf(stdout, "%s\n", message);

          }
      }
   }




   if( pf[*unit] == ((void *)0) ) {
      *err = 1;
      open_flag[*unit]=1;
   }
   else {
      if ( setvbuf( pf[*unit], ((void *)0), 0, 8192 ) != 0 ) {
         perror("\n**Warning: setvbuf failed");
         *err=1;
         open_flag[*unit]=1;
       }
       else
       {
         *err = 0;
         open_flag[*unit]=0;



        }
   }
   io_position[*unit]=0;
   prev_fileop[*unit]=UM_UNSET;





  if(*intent != readonly) {




    if(unit_buffer[*unit] == ((void *)0)) {




      k=buffer_size*sizeof(um_data_t);



      real_add=malloc(k);



      if( real_add == ((void *)0)) {
        malloc_error=(*__errno_location ());
        fprintf(stderr,
         "OPEN: PE %d had C I/O Error: failed in MALLOC\n",
         0);
        fprintf(stderr,
         "OPEN: Return code = %d, while claiming %lld bytes\n",
         malloc_error, (long long) k);
        abort();
      }
      else {
        snprintf(message, 384 ,
         "OPEN:  Claimed %lld Bytes (%lld Words) for Buffering",
         (long long) k, (long long) k/sizeof(um_data_t));
        fprintf(stdout, "%s\n", message);
        snprintf(message, 384 ,
         "OPEN:  Buffer Address is %16lX", (unsigned long) real_add);
        fprintf(stdout, "%s\n", message);
      }



      unit_buffer[*unit]=(um_data_t *)real_add;



      unit_offset[*unit]=0;
      unit_intent[*unit]=1;
    }
  }
  else {
    unit_intent[*unit]=0;
  }



   if (io_timer_active){
     gettimeofday( &second, &tzp );
     io_update_timer( first, second, 0, 2 );
   }

free (fname);
}



void
# 863 "portio2a.c"
close_single_




(unit,file_name,char_len,environ_var_flag,delete,err)
char file_name[];
integer *unit;
integer *char_len;
integer *delete;
integer *environ_var_flag;

integer *err;
{
char *fname;
char *gname;
struct timeval first;
struct timeval second;
struct timezone tzp;
int i;
integer k;

the_unit=unit;
fname = calloc(*char_len + 1,1);

if(open_flag[*unit]== 0){


      if (io_timer_active){
        gettimeofday( &first, &tzp );
      }


      flush_unit_buffer(unit);


      k=fclose(pf[*unit]);


        strncpy(fname,file_name,*char_len);
        fname[*char_len] = '\0';
        for (i=0; i<*char_len; i++){

            if (fname[i] == ' '){
               fname[i] = '\0';
               break;
            }
         }

        if(*environ_var_flag == 0)
          { gname = getenv( fname );
            if ( gname == ((void *)0) ) {
              snprintf(message, 384 ,
               "CLOSE: WARNING: Environment variable %s not set",
               fname);
              fprintf(stdout, "%s\n", message);
            open_flag[*unit]=1;
            *err=1;
            free( fname );
            return;
            }
          }
        else
          gname=fname;

      if(k==0){


        if(*delete != 0){
          k=remove(gname);
          if( k != 0){
            fprintf(stderr,
             "CLOSE: Cannot Delete File %s",gname);
            *err = 1;
            abort();
          }
          else{
            open_flag[*unit]=1;
            *err = 0;
            if (printstatus >= 2) {
              snprintf(message, 384 ,
               "CLOSE: File %s Deleted", gname);
              fprintf(stdout, "%s\n", message);
            }
          }

        }
        else{

           open_flag[*unit]=1;
          if (printstatus >= 2) {
            snprintf(message, 384 ,
             "CLOSE: File %s Closed on Unit %d",
             gname, (int) *unit);
            fprintf(stdout, "%s\n", message);
          }
        }
      }

    else {
          snprintf(message, 384 ,
           "CLOSE: Cannot Close File %s on Unit %d",
           gname, (int) *unit);
          fprintf(stdout, "%s\n", message);
    }

    if (io_timer_active){
      gettimeofday( &second, &tzp );
      io_update_timer( first, second, 0, 3 );
    }

}

else {

      if (printstatus >= 2) {
  snprintf(message, 384 ,
   "CLOSE: WARNING: Unit %d Not Opened", (int) *unit);
  fprintf(stdout, "%s\n", message);
      }
}

free( fname );
}



void
# 1003 "portio2a.c"
buffo32_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
  integer nsize;
  integer elsize;


  elsize = 4;

  if (open_flag[*unit]==0){

    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_WRITE;


    nsize = (*maxlen+1)/2;





    change_endian(array, 4, nsize);




    if (*maxlen % (sizeof(um_data_t)/4) == 0)
    {


      use_buffer(unit, array, &nsize, length, status);
    }
    else
    {

      flush_unit_buffer(unit);
      output_buffer(unit, array, maxlen, length, &elsize, status);
    }





    change_endian(array, 4, nsize);




    if (*maxlen % (sizeof(um_data_t)/4) == 0)
    {
      *length = *length*(sizeof(um_data_t)/4);
    }

  }
}

void
# 1084 "portio2a.c"
buffin32_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
  struct timeval first;
  struct timeval second;
  struct timezone tzp;
  int nsize;
  int k;

    if (open_flag[*unit]==0){

      if(prev_fileop[*unit] == UM_WRITE)
      {
        sync_to_cache(unit);
      }
      prev_fileop[*unit] = UM_READ;

      if (io_timer_active){
         gettimeofday( &first, &tzp );
      }

      *length = fread(array,4,*maxlen,pf [*unit]);

      if (io_timer_active){
        gettimeofday( &second, &tzp );
        io_update_timer( first, second, 4 * (*maxlen),
                         0 );
      }



      nsize = (*maxlen+1)/2;



      change_endian(array, 4, nsize);


      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        snprintf(message, 384 ,"C I/O Error: failed in BUFFIN32\n");
        fprintf(stdout, "%s\n", message);
        snprintf(message, 384 ,"Return code = %d\n",k);
        fprintf(stdout, "%s\n", message);
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        snprintf(message, 384 , "C I/O Error: failed in BUFFIN32\n");
        fprintf(stdout, "%s\n", message);
        snprintf(message, 384 , "Return code = %d\n",k);
        fprintf(stdout, "%s\n", message);
        *status=1.0;
      }
    }
    else
      *status=3.0;

    io_position[*unit]=io_position[*unit]+(*length*4);
}

void
# 1170 "portio2a.c"
buffin16_single_




(unit, array, maxlen, length, status)
integer *unit;



short array[];

integer *maxlen;
integer *length;
real *status;
{
  struct timeval first;
  struct timeval second;
  struct timezone tzp;
  integer k;

    if (open_flag[*unit]==0){

      if(prev_fileop[*unit] == UM_WRITE)
      {
        sync_to_cache(unit);
      }
      prev_fileop[*unit] = UM_READ;
# 1232 "portio2a.c"
      if (io_timer_active){
         gettimeofday( &first, &tzp );
      }

      *length = fread(array,2,*maxlen,pf [*unit]);

      if (io_timer_active){
        gettimeofday( &second, &tzp );
        io_update_timer( first, second, 2 * (*maxlen),
                         0 );
      }


      change_endian(array, 2, *maxlen);


      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        snprintf(message, 384 , "C I/O Error: failed in BUFFIN16\n");
        fprintf(stdout, "%s\n", message);
        snprintf(message, 384 , "Return code = %lld\n", (long long) k);
        fprintf(stdout, "%s\n", message);
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        snprintf(message, 384 , "C I/O Error: failed in BUFFIN16\n");
        fprintf(stdout, "%s\n", message);
        snprintf(message, 384 , "Return code = %lld\n", (long long) k);
        fprintf(stdout, "%s\n", message);
        *status=1.0;

      }
    }
    else
      *status=3.0;

    io_position[*unit]=io_position[*unit]+(*length*2);
}

void
# 1288 "portio2a.c"
buffin8_single_




(unit, array, maxlen, length, status)
integer *unit;
char array[];
integer *maxlen;
integer *length;
real *status;
{
  int k;

  if(open_flag[*unit]== 0){
    *length = fread(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      snprintf(message, 384 , "C I/O Error: failed in BUFFIN8\n");
      fprintf(stdout, "%s\n", message);
      snprintf(message, 384 , "Return code = %d\n",k);
      fprintf(stdout, "%s\n", message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      snprintf(message, 384 ,"C I/O Error: failed in BUFFIN8\n");
      fprintf(stdout, "%s\n", message);
      snprintf(message, 384 ,"Return code = %d\n",k);
      fprintf(stdout, "%s\n", message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

}

void



buffou8_



(unit, array, maxlen, length, status)
integer *unit;
char array[];
integer *maxlen;
integer *length;
real *status;
{
  int k;

  if(open_flag[*unit]== 0){


    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_WRITE;

    *length = fwrite(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      snprintf(message, 384 , "C I/O Error: failed in BUFFOU8\n");
      fprintf(stdout, "%s\n", message);
      snprintf(message, 384 , "Return code = %d\n",k);
      fprintf(stdout, "%s\n", message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      snprintf(message, 384 , "C I/O Error: failed in BUFFOU8\n");
      fprintf(stdout, "%s\n", message);
      snprintf(message, 384 , "Return code = %d\n",k);
      fprintf(stdout, "%s\n", message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

}
void
# 1395 "portio2a.c"
setpos8_single_




(unit, byte_address)
integer *unit;
integer *byte_address;
{

off_t k;





    k = fseeko(pf[*unit],*byte_address,0);



    if(k!=0){
         snprintf(message, 384 , "ERROR detected in SETPOS\n");
         fprintf(stdout, "%s\n", message);
         snprintf(message, 384 , "word_address = %lld\n", (long long) *byte_address);
         fprintf(stdout, "%s\n", message);
         snprintf(message, 384 , "Return code from fseek = %d\n",(int) k);
         fprintf(stdout, "%s\n", message);
    }
}

void



getpos8_



(unit, byte_address)
integer *unit;
integer *byte_address;
{

off_t k;





    k = ftello(pf[*unit]);




    *byte_address = k;
}

void
# 1465 "portio2a.c"
setpos32_single_




(unit,word32_address,err)

integer *unit;
integer *word32_address;
integer *err;


{
  int k;

  off_t byte_address;




  *err=0;
  if (open_flag[*unit]==0) {



  if(unit_intent[*unit] != readonly) {



    if(io_position[*unit]+unit_offset[*unit]*sizeof(um_data_t) == *word32_address*4) {



      return;

    }


      flush_unit_buffer(unit);
}



    byte_address=(*word32_address)*4;


    k=fseeko(pf[*unit],byte_address,0);




    if (k != 0) {
      k=(*__errno_location ());
      perror("\nSETPOS32: Seek failed");
      fprintf(stderr,
        "SETPOS32: Unit %d to 32bit Word Address %d failed. Error: %d",
        (int) *unit, (int) *word32_address, k);
      the_unit=unit;
      *err=1;
      abort();
    }

  }
}

void



getpos32_



(unit,word32_address)

integer *unit;
integer *word32_address;

{

  off_t byte_address;




  if (open_flag[*unit]==0) {

    byte_address=ftello(pf[*unit]);



    *word32_address = byte_address/4;

    if (byte_address != *word32_address * 4){
      snprintf(message, 384 ,
       "GETPOS32: byte address (%16lX) is not aligned on a %d byte boundary!",
       (unsigned long) byte_address, 4);
      fprintf(stdout, "%s\n", message);
      abort();

    }

  }
}

void



getpos_



(unit, word_address)
integer *unit;
integer *word_address;
{

    off_t byte_address;




    the_unit=unit;

     byte_address = ftello(pf[*unit]);



    *word_address = byte_address/sizeof(um_data_t);






    if (byte_address != *word_address * sizeof(um_data_t)){
      fprintf(stderr,
       "GETPOS: byte address (%16lX) is not aligned on a %lu byte boundary!",
       (unsigned long) byte_address, (unsigned long) sizeof(um_data_t));
      abort();
    }
# 1618 "portio2a.c"
    if (unit_offset[*unit]!=0)
      *word_address+=unit_offset[*unit];

      if(byte_address != io_position[*unit]) {
        fprintf(stderr,
         "GETPOS: IO_POSITION is %lld, but FTELL gives %lld",
         (long long) io_position[*unit], (long long) byte_address);
        the_unit=unit;
        abort();
      }

}





void



flush_buffer_



(unit, icode)
integer *unit ;
integer *icode ;
{
int i ;
  if(open_flag[*unit]== 0){


      flush_unit_buffer(unit);



      i = fflush(pf[*unit]);
      *icode = i;
  }
  else {
    if(pf[*unit] != ((void *)0)) {
      fprintf(stderr,
       "FLUSH_BUFFER: File Pointer for Unopened Unit %d is %16lX",
       (int) *unit, (unsigned long) pf[*unit]);
      the_unit=unit;
      abort();
    }
  }
}
# 1680 "portio2a.c"
void change_endian (void *ptr_array, int nbytes, int nsize) {



  integer *ptr_array_integer;
  int *ptr_array_int;
  short *ptr_array_short;
  int i;

  if (nbytes == 8)
  {


     ptr_array_integer = (integer *) ptr_array;
#pragma omp parallel for private(i) shared(nsize, ptr_array_integer)
     for (i=0 ; i<nsize ; i++)
       ptr_array_integer[i] = (__extension__ ({ register unsigned long __v, __x = (ptr_array_integer[i]); if (__builtin_constant_p (__x)) __v = ((((__x) & 0xff00000000000000ull) >> 56) | (((__x) & 0x00ff000000000000ull) >> 40) | (((__x) & 0x0000ff0000000000ull) >> 24) | (((__x) & 0x000000ff00000000ull) >> 8) | (((__x) & 0x00000000ff000000ull) << 8) | (((__x) & 0x0000000000ff0000ull) << 24) | (((__x) & 0x000000000000ff00ull) << 40) | (((__x) & 0x00000000000000ffull) << 56)); else __asm__ ("bswap %q0" : "=r" (__v) : "0" (__x)); __v; }));
  }
  else if (nbytes == 4)
  {

     ptr_array_int = (int *) ptr_array;
#pragma omp parallel for private(i) shared(nsize, ptr_array_int)
     for (i=0 ; i<nsize*sizeof(integer)/4 ; i++)
       ptr_array_int[i] = (__extension__ ({ register unsigned int __v, __x = (ptr_array_int[i]); if (__builtin_constant_p (__x)) __v = ((((__x) & 0xff000000) >> 24) | (((__x) & 0x00ff0000) >> 8) | (((__x) & 0x0000ff00) << 8) | (((__x) & 0x000000ff) << 24)); else __asm__ ("bswap %0" : "=r" (__v) : "0" (__x)); __v; }));
  }
  else
  {


     ptr_array_short = (short *) ptr_array;
#pragma omp parallel for private(i) shared(nsize, ptr_array_short)
     for (i=0 ; i<nsize ; i++)
       ptr_array_short[i] = (__extension__ ({ register unsigned short int __v, __x = (unsigned short int) (ptr_array_short[i]); if (__builtin_constant_p (__x)) __v = ((unsigned short int) ((((__x) >> 8) & 0xff) | (((__x) & 0xff) << 8))); else __asm__ ("rorw $8, %w0" : "=r" (__v) : "0" (__x) : "cc"); __v; }));

  }

  return;
}
# 1729 "portio2a.c"
void



is_unit_open_local_



(unit, ret_code)
integer *unit;
integer *ret_code;
{

  *ret_code = open_flag[*unit];
  return;
}

void
sync_to_cache
(integer * unit)
{
  if(open_flag[*unit]== 0){


    flush_unit_buffer(unit);


    fflush(pf[*unit]);
  }
}


void



sync_single_



(integer * unit,integer * icode)
{
  *icode=0;
  sync_to_cache(unit);
  fsync(fileno(pf[*unit]));
}






void



getextent_



(integer * unit,integer * extent,integer * icode)
{

  off_t c_pos;
  off_t e_pos;





  *icode=0;


  flush_unit_buffer(unit);




  c_pos=ftello(pf[*unit]);
  fseeko(pf[*unit],0,2);
  e_pos=ftello(pf[*unit]);
  fseeko(pf[*unit],c_pos,0);







  *extent=e_pos;
}


void



portioshutdown_



()
{

}


void
# 1849 "portio2a.c"
buffout_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
# 1878 "portio2a.c"
  buffout64_single_
  (unit, array, maxlen, length, status);
# 1910 "portio2a.c"
}


void
# 1926 "portio2a.c"
buffin_single_




(unit, array, maxlen, length, status)
integer *unit;
um_data_t array[];
integer *maxlen;
integer *length;
real *status;
{
# 1955 "portio2a.c"
  buffin64_single_
  (unit, array, maxlen, length, status);
# 1987 "portio2a.c"
}
# 1998 "portio2a.c"
void portioinit_



(integer* b_size,
 integer* dummyA,
 integer* dummyB,
 integer* dummyC,
 integer* dummyD,
 integer* io_timing)
{
buffer_size = *b_size;



snprintf(message, 384 ,
 "Buffered I/O active. Buffer size set to %lld, %lu byte words ", (long long) *b_size, (unsigned long) sizeof(um_data_t));
fprintf(stdout, "%s\n", message);

}




void portioaddhelper_



(integer * role)
{



}




void portiodettachallhelpers_



()
{



}
