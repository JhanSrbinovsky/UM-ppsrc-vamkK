# 1 "pio_data_conv.c"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "pio_data_conv.c"
# 58 "pio_data_conv.c"
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

# 59 "pio_data_conv.c" 2
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

# 60 "pio_data_conv.c" 2
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

# 61 "pio_data_conv.c" 2

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

# 63 "pio_data_conv.c" 2



# 1 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h" 1
# 28 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h"
typedef unsigned long long int u_integer;
typedef long long int integer;
# 43 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_fort2c_2a.h"
typedef double real;
typedef integer um_data_t;
# 67 "pio_data_conv.c" 2







integer *the_unit;
static char message[256];


# 1 "/short/p66/jxs599/UM_ROUTDIR/jxs599/vamkk/umatmos/inc/c_data_conv.h" 1
# 79 "pio_data_conv.c" 2






integer



ieee2ieg_



(type, num, ieg_num_out, offset_out, cri_num_in, stride,
    size_num_in, size_num_out)
integer *type, *num, *offset_out, *stride, *size_num_in, *size_num_out;
integer ieg_num_out[], cri_num_in[];
{

void read_number();
void write_number();

integer errorcode=0 ;
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer);
int out_offset ;


integer sign=0, expo=0, mant=0, mant_bits_read=0;
integer unit;

unit = -1;
the_unit = &unit;



if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid num = %d. Return code = %d\n",
     (int) *num, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid stride = %d. Return code = %d\n",
     (int) *stride, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



 offset = word_length - (int) *offset_out - (int) *size_num_out;




if (offset < 0 || offset%( (int) *size_num_out ) != 0
               || offset > word_length
               || (int) *offset_out < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEEE2IEG: Error - Invalid bitoff = %d. Return code = %d\n",
     (int) *offset_out, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}




if ((int) *type == 2){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 3){



  switch ((int) *size_num_in){
    case 32: type_num_in = 3; break;

    case 64: type_num_in = 4; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 32: type_num_out = 3; break;

    case 64: type_num_out = 4; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid forlen %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 5){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IEG: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
          "IEEE2IEG: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else {

  errorcode = -2 ;
  sprintf(message,
   "IEEE2IEG: Error - unsupported data type = %d. Return code = %d\n",
     (int) *type, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



out_offset = offset ;

for (i=0; i < *num; i++){
  read_number(type_num_in, (int) *size_num_in, 0, *cri_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, out_offset,
               ieg_num_out, sign, expo, mant, mant_bits_read);

  if ( out_offset == 0 ) {
    out_offset = word_length - (int) *size_num_out ;
    ieg_num_out++ ;
  }
  else {
    out_offset = out_offset - (int) *size_num_out ;
  }

  cri_num_in++ ;

}
 return errorcode ;
}






integer



ieg2ieee_



(type, num, ieg_num_in, offset_in, cri_num_out, stride,
    size_num_out, size_num_in)
integer *type, *num, *offset_in, *stride, *size_num_in, *size_num_out;
integer ieg_num_in[], cri_num_out[];
{

void read_number();
void write_number();

integer errorcode=0 ;
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer);
int in_offset ;


integer sign=0, expo=0, mant=0, mant_bits_read=0;

integer unit;

unit = -1;
the_unit = &unit;



if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid num = %d. Return code = %d\n",
     (int) *num, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid stride = %d. Return code = %d\n",
     (int) *stride, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



 offset = word_length - (int) *offset_in - (int) *size_num_in;




if (offset < 0 || offset%( (int) *size_num_in ) != 0
               || offset > word_length
               || (int) *offset_in < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEG2IEEE: Error - Invalid bitoff = %d. Return code = %d\n",
     (int) *offset_in, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}





if ((int) *type == 2){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 3){



  switch ((int) *size_num_in){
    case 32: type_num_in = 3; break;

    case 64: type_num_in = 4; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 32: type_num_out = 3; break;

    case 64: type_num_out = 4; break;

    default: errorcode = -6 ;
        sprintf(message,
           "IEG2IEEE: Error - Invalid forlen %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 5){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
          "IEG2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEG2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else {



  errorcode = -2 ;
  sprintf(message,
   "IEG2IEEE: Error - unsupported data type = %d. Return code = %d\n",
     (int) *type, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



in_offset = offset ;

for (i=0; i<*num; i++){
  read_number(type_num_in, (int) *size_num_in, in_offset, *ieg_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, 0, cri_num_out,
                      sign, expo, mant, mant_bits_read);

  if ( in_offset == 0 ) {
    in_offset = word_length - (int) *size_num_in ;
    ieg_num_in++ ;
  }
  else {
    in_offset = in_offset - (int) *size_num_in ;
  }

  cri_num_out++ ;

}
 return errorcode ;
}






integer



ieee2ibm_



(type, num, ibm_num_out, offset_out, cri_num_in, stride,
    size_num_in, size_num_out)
integer *type, *num, *offset_out, *stride, *size_num_in, *size_num_out;
integer ibm_num_out[], cri_num_in[];
{

void read_number();
void write_number();

integer errorcode=0 ;
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer);
int out_offset ;


integer sign=0, expo=0, mant=0, mant_bits_read=0;

integer unit;

unit = -1;
the_unit = &unit;



if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid num = %d. Return code = %d\n",
     (int) *num, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid stride = %d. Return code = %d\n",
     (int) *stride, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}


 offset = word_length - (int) *offset_out - (int) *size_num_out;




if (offset < 0 || offset%( (int) *size_num_out ) != 0
               || offset > word_length
               || (int) *offset_out < 0){
  errorcode = -4 ;
  sprintf(message,
     "IEEE2IBM: Error - Invalid bitoff = %d. Return code = %d\n",
     (int) *offset_out, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}




if ((int) *type == 2){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 3){



  switch ((int) *size_num_in){
    case 32: type_num_in = 3; break;

    case 64: type_num_in = 4; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 32: type_num_out = 1; break;

    case 64: type_num_out = 2; break;

    default: errorcode = -6 ;
        sprintf(message,
           "IEEE2IBM: Error - Invalid forlen %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 5){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid natlen = %d . Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IEEE2IBM: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else {



  errorcode = -2 ;
  sprintf(message,
   "IEEE2IBM: Error - unsupported data type = %d. Return code = %d\n",
     (int) *type, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



out_offset = offset ;

for (i=0; i < *num; i++){

  read_number(type_num_in, (int) *size_num_in, 0, *cri_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, out_offset,
               ibm_num_out, sign, expo, mant, mant_bits_read);

  if ( out_offset == 0 ) {
    out_offset = word_length - (int) *size_num_out ;
    ibm_num_out++ ;
  }
  else {
    out_offset = out_offset - (int) *size_num_out ;
  }

  cri_num_in++ ;

}
 return errorcode ;
}






integer



ibm2ieee_



(type, num, ibm_num_in, offset_in, cri_num_out, stride,
size_num_out, size_num_in)
integer *type, *num, *offset_in, *stride, *size_num_in, *size_num_out;
integer ibm_num_in[], cri_num_out[];
{

void read_number();
void write_number();

integer errorcode=0 ;
int i;
int type_num_in, type_num_out, offset;
int word_length=8*sizeof(integer);
int in_offset ;


integer sign=0, expo=0, mant=0, mant_bits_read=0;

integer unit;

unit = -1;
the_unit = &unit;



if ((int) *num <= 0) {
  errorcode = -3 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid num = %d. Return code = %d\n",
     (int) *num, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}

if ((int) *stride != 1) {
  errorcode = -7 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid stride = %d. Return code = %d\n",
     (int) *stride, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}


 offset = word_length - (int) *offset_in - (int) *size_num_in;




if (offset < 0 || offset%( (int) *size_num_in ) != 0
               || offset > word_length
               || (int) *offset_in < 0){
  errorcode = -4 ;
  sprintf(message,
     "IBM2IEEE: Error - Invalid bitoff = %d. Return code = %d\n",
     (int) *offset_in, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}




if ((int) *type == 2){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 3){



  switch ((int) *size_num_in){
    case 32: type_num_in = 1; break;

    case 64: type_num_in = 2; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 32: type_num_out = 3; break;

    case 64: type_num_out = 4; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else if ((int) *type == 5){



  switch ((int) *size_num_in){
    case 16: type_num_in = 6; break;
    case 32: type_num_in = 7; break;

    case 64: type_num_in = 8; break;

    default: errorcode = -5 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid natlen = %d. Return code = %d\n",
           (int) *size_num_in, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
  switch ((int) *size_num_out){
    case 16: type_num_out = 6; break;
    case 32: type_num_out = 7; break;

    case 64: type_num_out = 8; break;

    default: errorcode = -6 ;
        sprintf(message,
         "IBM2IEEE: Error - Invalid forlen = %d. Return code = %d\n",
           (int) *size_num_out, (int) errorcode);
        fprintf(stdout, "%s\n", message);
        return errorcode ;
  }
}
else {

  errorcode = -2 ;
  sprintf(message,
   "IBM2IEEE: Error - unsupported data type = %d. Return code = %d\n",
     (int) *type, (int) errorcode);
  fprintf(stdout, "%s\n", message);
  return errorcode ;
}



in_offset = offset ;

for (i=0; i<*num; i++){
  read_number(type_num_in, (int) *size_num_in, in_offset, *ibm_num_in,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (int) *size_num_out, 0, cri_num_out,
                      sign, expo, mant, mant_bits_read);

  if ( in_offset == 0 ) {
    in_offset = word_length - (int) *size_num_in ;
    ibm_num_in++ ;
  }
  else {
    in_offset = in_offset - (int) *size_num_in ;
  }

  cri_num_out++ ;

}
 return errorcode ;
}







void read_number
(format,size,offset,in_number,out_sign,out_expo,out_mant,mant_bits_read)
int format, size, offset;
integer in_number;
integer *out_sign, *out_expo, *out_mant, *mant_bits_read;
{

integer expo_bits, expo_bias;
integer sign_mask, expo_mask, mant_mask, in_number_mask;
integer one ;
int i, first_bit;
int word_length=8*sizeof(integer) ;

switch (format){
  case 1:

    sign_mask = 0x80000000;
    expo_mask = 0x7F000000;
    expo_bits = 7;
    expo_bias = 64;
    mant_mask = 0x00FFFFFF;
    *mant_bits_read = 24;
  break;


  case 2:

    sign_mask = 0x8000000000000000LL;
    expo_mask = 0x7F00000000000000LL;
    expo_bits = 7;
    expo_bias = 64;
    mant_mask = 0x00FFFFFFFFFFFFFFLL;
    *mant_bits_read = 56;
  break;


  case 3:

    sign_mask = 0x80000000;
    expo_mask = 0x7F800000;
    expo_bits = 8;
    expo_bias = 127;
    mant_mask = 0x007FFFFF;
    *mant_bits_read = 23;
  break;


  case 4:

    sign_mask = 0x8000000000000000LL;
    expo_mask = 0x7FF0000000000000LL;
    expo_bits = 11;
    expo_bias = 1023;
    mant_mask = 0x000FFFFFFFFFFFFFLL;
    *mant_bits_read = 52;
  break;

  case 6:

    sign_mask = 0x8000;
    mant_mask = 0x7FFF;
    *mant_bits_read = 15;
    expo_bits = 0;
    expo_bias = 0;
  break;

  case 7:

    sign_mask = 0x80000000;
    mant_mask = 0x7FFFFFFF;
    *mant_bits_read = 31;
    expo_bits = 0;
    expo_bias = 0;
  break;


  case 8:

    sign_mask = 0x8000000000000000LL;
    mant_mask = 0x7FFFFFFFFFFFFFFFLL;
    *mant_bits_read = 63;
    expo_bits = 0;
    expo_bias = 0;
  break;

}

if ( offset != 0 || size != (*mant_bits_read + expo_bits + 1)
     || (word_length - offset -size != 0) ){



  in_number = in_number >> offset;
  one = 1 ;
  in_number_mask = ( (one << size) - 1 );
  in_number = in_number & in_number_mask;


  one = 1 ;
  sign_mask = one << (size - 1);
}





  *out_sign= (in_number & sign_mask) >> (*mant_bits_read + expo_bits);



if ((in_number & mant_mask) == 0 && (in_number & expo_mask) == 0) {
  *out_mant = 0 ;
  *out_expo = 0 ;
}
else{
if( format == 1 || format ==2 ){





  first_bit = 0;
  for (i = 1; first_bit == 0; i++){
    if((((in_number & mant_mask) >> (*mant_bits_read - i)) & 1 ) == 1){
      first_bit = i;
    }
    if(i == *mant_bits_read) first_bit = *mant_bits_read;
  }



  *out_mant = (in_number & mant_mask) << first_bit;




  *out_expo =((((in_number & expo_mask) >> *mant_bits_read)
                - expo_bias) * 4) - first_bit;


  *mant_bits_read = *mant_bits_read + 1;

}
else if(format == 3 || format == 4){




  one = 1 ;
  one = one << *mant_bits_read ;

  *out_mant= (in_number & mant_mask) | one ;


  *out_expo= ((in_number & expo_mask) >> *mant_bits_read ) - expo_bias;


  *mant_bits_read = *mant_bits_read + 1;
}
else {

  *out_mant = mant_mask & in_number;
  *out_expo = 0;
}
}
# 1070 "pio_data_conv.c"
}







void write_number
(format,size,offset,out_number,in_sign,in_expo,in_mant,mant_bits_read)
int format, size, offset;
integer *out_number;
integer in_sign, in_expo, in_mant, mant_bits_read;
{
integer mant_bits_write, expo_bits_write;
integer sign_mask, expo_mask, mant_mask;
integer sign, mant, expo, conv_number, conv_number_mask;
integer first_bit, temp;
integer one, add_one ;
integer mant_bits_diff, expo_bias;
int i, expo_diff;
int word_length=8*sizeof(integer) ;

switch(format){
  case 1:

    sign= in_sign << 31;
    mant_bits_write= 24;
    expo_bits_write= 7;
    sign_mask= 0x80000000;
    expo_mask= 0x7F000000;
    mant_mask= 0x00FFFFFF;
    expo_bias= 64;
  break;


  case 2:

    sign= in_sign << 63;
    mant_bits_write= 56;
    expo_bits_write= 7;
    sign_mask= 0x8000000000000000LL;
    expo_mask= 0x7F00000000000000LL;
    mant_mask= 0x00FFFFFFFFFFFFFFLL;
    expo_bias= 64;
  break;

  case 3:

    sign= in_sign << 31;
    mant_bits_write= 23;
    expo_bits_write= 8;
    sign_mask= 0x80000000;
    expo_mask= 0x7F800000;
    mant_mask= 0x007FFFFF;
    expo_bias= 127;
  break;


  case 4:

    sign= in_sign << 63;
    mant_bits_write= 52;
    expo_bits_write= 11;
    sign_mask= 0x8000000000000000LL;
    expo_mask= 0x7FF0000000000000LL;
    mant_mask= 0x000FFFFFFFFFFFFFLL;
    expo_bias= 1023;
  break;


  case 6:

    sign= in_sign << 15;
    mant_bits_write= 15;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= 0x8000;
    mant_mask= 0x7FFF;
  break;

  case 7:

    sign= in_sign << 31;
    mant_bits_write= 31;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= 0x80000000;
    mant_mask= 0x7FFFFFFF;
  break;


  case 8:

    sign= in_sign << 63;
    mant_bits_write= 63;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= 0x8000000000000000LL;
    mant_mask= 0x7FFFFFFFFFFFFFFFLL;
  break;


}

if ( in_mant == 0 && in_expo == 0 ) {
  mant = 0 ;
  expo = 0 ;
}
else {


if(format == 1 || format == 2){




  if ( in_expo > -4 ) {
    expo = in_expo/4 + 1 ;

    if (in_expo < 0) {

      expo-- ;
    }

    expo_diff = in_expo - 4*expo ;
  }
  else {
    expo = in_expo/4 ;
    expo_diff = in_expo - 4*expo ;
    if (expo_diff == 0){
      expo++ ;
      expo_diff = in_expo - 4*expo ;
    }
  }
  mant = in_mant >> ( 0 - expo_diff ) ;



  mant_bits_read--;

  if ( expo >= expo_bias ) {

    expo = expo_bias - 1 ;
    mant = mant_mask ;
    expo_diff = 0 ;
    mant_bits_diff = 0 ;
  }
  else if ( expo < -expo_bias ) {

    expo = -expo_bias ;
    mant = 0 ;
    expo_diff = 0;
    mant_bits_diff = 0;
  }
  else {

    mant_bits_diff = mant_bits_write - mant_bits_read;
  }

  expo = expo + expo_bias;

  if (mant_bits_diff < 0) {
    mant = (mant >> (0 - mant_bits_diff - 1) ) ;


    add_one = 0 ;

    if ( (mant & 1) == 1 ) {
      add_one = 1 ;
      if ( (mant >> 1) == mant_mask ) {

        mant = (mant >> 1) + add_one ;
        mant = (mant >> 4) ;
        expo++ ;
      }
      else {
        mant = (mant >> 1) + add_one ;
      }
    }
    else {
      mant = (mant >> 1) ;
    }
  }
  else if(mant_bits_diff > 0) mant = ( mant << mant_bits_diff );
  else mant = mant;

}
else if (format == 3 || format == 4 ){





  one = 1 ;
  mant = in_mant ^ (one << (mant_bits_read - 1) );
  mant_bits_read--;

  if ( in_expo >= expo_bias ) {

    expo = 2*expo_bias + 1 ;
    mant = mant_mask + 1 ;
    expo_diff = 0 ;
    mant_bits_diff = 0 ;
  }
  else if ( in_expo < -expo_bias ) {

    expo = 0 ;
    mant = 0 ;
    expo_diff = 0;
    mant_bits_diff = 0;
  }
  else {

    expo = in_expo + expo_bias;
    mant_bits_diff = mant_bits_write - mant_bits_read;
  }
  if (mant_bits_diff < 0) {
    mant = (mant >> (0 - mant_bits_diff - 1) ) ;


    add_one = 0 ;
    if ( (mant & 1) == 1 ) {
      add_one = 1 ;
      if ( (mant >> 1) == mant_mask ) {

        mant = (mant >> 1) + add_one ;
        expo++ ;
      }
      else {
        mant = (mant >> 1) + add_one ;
      }
    }
    else {
      mant = (mant >> 1) ;
    }
  }
  else if(mant_bits_diff > 0) mant = (mant << mant_bits_diff );
  else mant = mant;

}

else if(format>5){







  first_bit=0;
  if (in_sign == 0){

    for (i=1;first_bit==0;i++){
      if ((( in_mant >> (mant_bits_read - i) ) & 1) == 1) first_bit = i;
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }


  else{

    temp = 0;
    for (i=1;first_bit==0;i++){
      if ( (in_mant >> (mant_bits_read - i)) < ( (temp<<1) + 1) )
          first_bit = i;
      temp = in_mant >> (mant_bits_read - i);
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }


  mant = in_mant & mant_mask;






  if(in_sign != 0 && mant_bits_write > mant_bits_read){
    one = 1 ;
    mant = mant | ( ( (one << (mant_bits_write - mant_bits_read)) - 1 )
                                << mant_bits_read ) ;


  }







  expo = 0;
}
}


sign= sign_mask & sign;
mant= mant_mask & mant;
expo= expo_mask & (expo << (mant_bits_write) );
conv_number= sign | expo | mant;
# 1383 "pio_data_conv.c"
  if (offset != 0 || size != (mant_bits_write + expo_bits_write + 1)
      || (word_length - offset - size) != 0) {

    one = 1;
    conv_number_mask = (one << size) - 1;
  }
  else {
    conv_number_mask = -1 ;
  }


  conv_number_mask = conv_number_mask << offset;
  conv_number = conv_number << offset;


  *out_number = conv_number | ( *out_number & (~conv_number_mask) );

}







void



movebytes_



(src,isb,num,dest,idb)
integer src[], *isb, *num, dest[], *idb;






{

  integer bytes, bytes_mask, dest_temp, one;
  integer nbytes_moved, nbytes_to_move, num_left ;
  integer src_start, dest_start ;
  integer word_length=sizeof(integer) ;
  int i ;

  nbytes_moved = 0 ;
  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;


  for (src_start = *isb ; src_start>word_length ;
                          src_start=src_start-word_length) {
    src++ ;
  }


  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {




    nbytes_to_move = word_length - src_start + 1;




    if (nbytes_to_move > num_left) nbytes_to_move = num_left ;






    if (nbytes_to_move > (word_length - dest_start + 1) ) {
       nbytes_to_move = word_length - dest_start + 1 ;
    }
# 1479 "pio_data_conv.c"
    if (nbytes_to_move == word_length) {

      bytes_mask = -1 ;
    }
    else{
      bytes_mask = 0;
      one = 1 ;
      bytes_mask = (one << (word_length * nbytes_to_move) ) - 1;
    }



    bytes =
         (*src >>
             word_length*(8 - (src_start-1) - nbytes_to_move ) )
          & bytes_mask;


    dest_temp = *dest &
                ( ~(bytes_mask << word_length*( 8 -
                   (dest_start-1) - nbytes_to_move ) ) );

    *dest = dest_temp |
        (bytes << word_length * ( 8 -
                    (dest_start-1) - nbytes_to_move ) );



    nbytes_moved = nbytes_moved + nbytes_to_move ;
    num_left = num_left - nbytes_to_move ;




    src_start = src_start + nbytes_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;

      src++ ;
    }




    dest_start = dest_start + nbytes_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;

      dest++ ;
    }
  }
}







void



movebits_



(src,isb,num,dest,idb)
integer src[], *isb, *num, dest[], *idb;






{


  integer bits, bits_mask, dest_temp, one;
  integer nbits_moved, nbits_to_move, num_left ;
  integer src_start, dest_start ;
  integer word_length=8*sizeof(integer) ;
  int i ;

  nbits_moved = 0 ;
  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;


  for (src_start = *isb; src_start>word_length ;
                         src_start=src_start-word_length) {
    src++ ;
  }


  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {




    nbits_to_move = word_length - src_start + 1;




    if (nbits_to_move > num_left) nbits_to_move = num_left ;






    if (nbits_to_move > (word_length - dest_start + 1) ) {
       nbits_to_move = word_length - dest_start + 1 ;
    }






    if (nbits_to_move == word_length) {

      bits_mask = -1 ;
    }
    else{
      bits_mask = 0;
      one = 1 ;
      bits_mask = (one << nbits_to_move) - 1;
    }


    bits = (*src >> ( word_length - src_start - nbits_to_move + 1) )
         & bits_mask;


    dest_temp = *dest
         & ( ~(bits_mask <<
              ( word_length - dest_start - nbits_to_move +1)) );
    *dest = dest_temp |
                ( bits <<
                ( word_length - dest_start - nbits_to_move +1) );



    nbits_moved = nbits_moved + nbits_to_move ;
    num_left = num_left - nbits_to_move ;




    src_start = src_start + nbits_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;

      src++ ;
    }




    dest_start = dest_start + nbits_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;

      dest++ ;
    }
  }
}
