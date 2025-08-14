##############################################################################
#
#                       H E R C U L E S . M 4
#
#               Hercules M4 macros for auto-configure
#
##############################################################################

# $Id$
#
# $Log$

#-----------------------------------------------------------------------------
#
# Macro:  HC_C99_FLEXIBLE_ARRAYS()
#
#  Checks whether or not the compiler supports C99 flexible arrays
#
# Parms:   None
# Input:   Nothing
# Output:  $hc_cv_c99_flexible_array,
#          AC_DEFINE() for 'C99_FLEXIBLE_ARRAYS' if yes.
# Note:    Since AC_DEFINE() might be issued, a corresponding
#          AH_TEMPLATE() for 'C99_FLEXIBLE_ARRAYS' is of course
#          also needed somewhere in configure.ac
#
# Credit:  Modeled after 'AC_C99_FLEXIBLE_ARRAY' macro copyright
#          by Erik de Castro Lopo <erikd AT mega-nerd DOT com>:
#
#             "Permission to use, copy, modify, distribute,
#              and sell this file for any purpose is hereby
#              granted without fee, provided that the above
#              copyright and this permission notice appear
#              in all copies.  No representations are made
#              about the suitability of this software for any
#              purpose.  It is provided "as is" without express
#              or implied warranty."
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_C99_FLEXIBLE_ARRAYS],
[
    AC_CACHE_CHECK( C99 struct flexible arrays support,

        [hc_cv_c99_flexible_array],
        [
            # Initialize to unknown
            hc_cv_c99_flexible_array=no

            AC_TRY_LINK(
                [
                    #include <stdlib.h>

                    typedef struct
                    {
                        int foo;
                        char bar[];
                    }
                    FOOBAR;
                ],
                [
                    int main(int argc, char *argv[])
                    {
                        FOOBAR* p = calloc( 1, sizeof(FOOBAR) + 16 );
                        return 0;
                    }
                ],
                [hc_cv_c99_flexible_array=yes],
                [hc_cv_c99_flexible_array=no]
            )
        ]
    )

    if test "$hc_cv_c99_flexible_array" = "yes"; then

        AC_DEFINE([C99_FLEXIBLE_ARRAYS])
    fi
])

#-----------------------------------------------------------------------------
#
# Macro:  HC_PROG_CC()       ((((( DEPRECATED )))))
#
#  Prevents undesired CFLAGS settings set by default AC_PROG_CC macro.
#
# Parms:   none
# Input:   $ac_env_CFLAGS_value
# Output:  CFLAGS initialized to our desired starting value
# Note:    use instead of the default AC_PROG_CC macro
#
#  AC_PROG_CC (actually _AC_PROC_CC_G) takes it upon itself to put "-g -O2"
#  in CFLAGS. While this may be good for most packages using autoconf, we
#  have our own "optimize" function that this interferes with. Thus this macro.
#
#  AC_BEFORE emits a warning if AC_PROG_CC was expanded prior to this macro
#  in case something gets put in configure.ac before us. AC_REQUIRE expands
#  AC_PROG_CC for us since we need it.
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_PROG_CC],
[
    AC_BEFORE(  [HC_PROG_CC], [AC_PROG_CC] )
    AC_REQUIRE( [AC_PROG_CC] )

    # (set CFLAGS to their initial value)

    CFLAGS=$ac_env_CFLAGS_value
])

#-----------------------------------------------------------------------------
#
# Macro:  HC_LD_DISALLOWDUPS()
#
#  Adds linker options to LDFLAGS variable to warn about duplicate symbols
#
# Parms:   none
# Input:   $lt_cv_prog_gnu_ld
# Output:  LDFLAGS variable modified as needed
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_LD_DISALLOWDUPS],
[
    AC_REQUIRE([AC_PROG_LIBTOOL])
    AC_REQUIRE([AC_PROG_LD_GNU])

    if test "$lt_cv_prog_gnu_ld" = "yes"; then

        LDFLAGS="$LDFLAGS -Wl,--warn-common"
    fi
])

#-----------------------------------------------------------------------------
#
# Macro:  HC_ARG_ENABLE_GETOPTWRAPPER()
#
#  Handles "--enable-getoptwrapper" configure option
#
# Parms:   none
# Input:   nothing
# Output:  $hc_cv_opt_getoptwrapper = "yes", "no" or "auto", as well as
#          other indirect side effects of _HC_CHECK_NEED_GETOPT_WRAPPER
#          helper macro being called
# Note:    calls _HC_CHECK_NEED_GETOPT_WRAPPER internal helper macro which
#          has its own requirements and side effects
#
#  Ultimately determines whether the 'getopt' wrapper kludge is necessary or
#  specifically requested. Issues AC_DEFINE for "NEED_GETOPT_WRAPPER" if true.
#  Refer to _HC_CHECK_NEED_GETOPT_WRAPPER helper macro for details.
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_ARG_ENABLE_GETOPTWRAPPER],
[
    AC_ARG_ENABLE( getoptwrapper,

        AC_HELP_STRING([--enable-getoptwrapper],

            [force use of the getopt wrapper kludge]
        ),
        [
            case "${enableval}" in
            yes)   hc_cv_opt_getoptwrapper=yes   ;;
            no)    hc_cv_opt_getoptwrapper=no    ;;
            auto)  hc_cv_opt_getoptwrapper=auto  ;;
            *)     hc_cv_opt_getoptwrapper=auto  ;;
            esac
        ],
        [hc_cv_opt_getoptwrapper=auto]
    )

    _HC_CHECK_NEED_GETOPT_WRAPPER($hc_cv_opt_getoptwrapper)
])

#------------------------------------------------------------------------------
#
# Macro:  HC_GET_C11_LOCK_FREE_VALUE( c11_lock_free_name )
#
#  Retrieve ATOMIC_*_LOCK_FREE macro compile-time / run-time value.
#
# Input:   nothing
#
# Parms:   c11_lock_free_name     the stdatomic.h pre-processor macro name
#                                 (e.g. ATOMIC_CHAR_LOCK_FREE) whose compile-
#                                 time or run-time value is to be returned.
#
# Output:  $c11_lock_free_value   the compile-time or run-time value of the
#                                 stdatomic.h ATOMIC_*_LOCK_FREE macro.
#
#------------------------------------------------------------------------------

AC_DEFUN( [HC_GET_C11_LOCK_FREE_VALUE],
[
    c11_lock_free_name="$1"

    if test "x$c11_lock_free_name" == "x"; then
        # PROGRAMMING NOTE: note special m4 escape syntax used in the below
        # message. Without it, our function name is (infinitely and recursively?)
        # expanded to our actual (as-yet not yet fully defined!) macro VALUE(!),
        # usually resulting in some type of weird m4 unmatched quoting error.
        AC_MSG_RESULT([LOGIC ERROR: no parm passed to ['HC_GET_C11_LOCK_FREE_VALUE'] function])
        hc_error=yes
        c11_lock_free_value="0"
    else
        AC_MSG_CHECKING([what '${c11_lock_free_name}'s value is])

        AC_TRY_RUN(
        [
            #if !defined( __STDC_NO_ATOMICS__ )
              #include <stdatomic.h>
            #endif
            int main( int argc, char* argv[] )
            {
            #if defined( __STDC_NO_ATOMICS__ ) || !defined( $c11_lock_free_name )
                return 0;               /* 0 = never atomic */
            #else
                return $c11_lock_free_name; /* 0=never, 1=sometimes or 2=always */
            #endif
            }
        ],
        [c11_lock_free_value="0"],
        [c11_lock_free_value="$?"],
        [c11_lock_free_value="0"] )

        AC_MSG_RESULT($c11_lock_free_value)
    fi
])

#-----------------------------------------------------------------------------
#
# Macro:  _HC_CHECK_NEED_GETOPT_WRAPPER( [opt = "auto"] )
#
#  Determines whether the 'getopt' wrapper kludge is necessary or requested
#
# Parms:   no       they specifically requested NOT to use the wrapper
#          yes      they specifically requested TO     use the wrapper
#          auto     we should determine ourselves if wrapper is needed
#
# Input:   nothing
# Output:  $hc_cv_need_getopt_wrapper = "yes" or "no",
#          AC_DEFINE for "NEED_GETOPT_WRAPPER" if yes.
#
# Note:    this macro is logically an internal subroutine of
#          the HC_ARG_ENABLE_GETOPTWRAPPER macro further above
#
#  If the passed value is "yes" or "no", we simply set the output variable
#  value to whatever value was passed. If the passed value is "auto" however,
#  then we determine for ourselves whether the wrapper is needed. We do this
#  by performing a test link of two small programs that each use getopt but
#  where one calls the other (thereby producing an interdependency between
#  the two) and then seeing if the link caused any duplicate symbol error.
#
#  If it is determined that the wrapper is needed (or if it was specifically
#  requested), an AC_DEFINE for "NEED_GETOPT_WRAPPER" is issued so Herc knows
#  to build/use its getopt wrapper kludge. Note that all AC_DEFINE's require
#  a corresponding AH_TEMPLATE statement somewhere in configure.ac.
#
#-----------------------------------------------------------------------------

AC_DEFUN([_HC_CHECK_NEED_GETOPT_WRAPPER],
[
    AC_REQUIRE([AC_PROG_LIBTOOL])
    AC_MSG_CHECKING([whether getopt wrapper kludge is necessary])

    if test "$1" != "auto"; then

        hc_cv_need_getopt_wrapper="$1"
        hc_cv_need_getopt_wrapper_result_msg="$1 (forced)"

    else

        if test $(./libtool --features | fgrep "enable shared libraries" | wc -l) -ne 1; then

            #  Libtool doesn't support shared libraries,
            #  and thus our wrapper kludge is not needed.

            hc_cv_need_getopt_wrapper=no
            hc_cv_need_getopt_wrapper_result_msg=no

        else

            rm -f libconftest*
            rm -f .libs/libconftest*

            cat > conftest1.c << DUPGETOPT1

                /*
                   Test program that needs getopt, called by
                   another program which itself needs getopt.
                   Will the linker complain about duplicate
                   symbols for getopt? We'll soon find out!
                */
		int getopt(int, char *const[], const char *);
                extern char *optarg;
                extern int optind;

                int test1()
                {
                    int i;
                    char *c;

                    i=optind;
                    c=optarg;

                    getopt(0,0,0);

                    return 0;
                }
DUPGETOPT1
            cat > conftest2.c << DUPGETOPT2

                /*
                   Test program that not only needs getopt,
                   but also calls another program which also
                   needs getopt. Will linker complain about
                   duplicate symbols for getopt? Let's see.
                */
		int getopt(int, char *const[], const char *);
                extern char *optarg;
                extern int optind;
                extern int test2();

                int test2()
                {
                    int i;
                    char *c;

                    i=optind;
                    c=optarg;

                    getopt(0,0,0);
                    test1();

                    return 0;
                }
DUPGETOPT2

            ./libtool --mode=compile ${CC-cc} conftest1.c -c -o conftest1.lo > /dev/null 2>&1
            ./libtool --mode=compile ${CC-cc} conftest2.c -c -o conftest2.lo > /dev/null 2>&1

            ./libtool --mode=link ${CC-cc} -shared -rpath /lib -no-undefined conftest1.lo                 -o libconftest1.la > /dev/null 2>&1
            ./libtool --mode=link ${CC-cc} -shared -rpath /lib -no-undefined conftest2.lo libconftest1.la -o libconftest2.la > /dev/null 2>&1

            if test $? = 0; then

                hc_cv_need_getopt_wrapper=no
                hc_cv_need_getopt_wrapper_result_msg=no
            else

                hc_cv_need_getopt_wrapper=yes
                hc_cv_need_getopt_wrapper_result_msg=yes
            fi

            rm -f *conftest*
            rm -f .libs/*conftest*
        fi
    fi

    AC_MSG_RESULT($hc_cv_need_getopt_wrapper_result_msg)

    if test "$hc_cv_need_getopt_wrapper" = "yes"; then

        AC_DEFINE([NEED_GETOPT_WRAPPER])
    fi
])

#-----------------------------------------------------------------------------
#
# Macro:  HC_CHECK_NEED_GETOPT_OPTRESET()
#
#  Checks whether or not 'optreset' needed for 'getopt' use
#
# Parms:   none
# Input:   nothing
# Output:  $hc_cv_need_getopt_optreset,
#          AC_DEFINE() for 'NEED_GETOPT_OPTRESET' if yes.
# Note:    since AC_DEFINE() might be issued, a corresponding AH_TEMPLATE()
#          for 'NEED_GETOPT_OPTRESET' is needed somewhere in configure.ac
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_CHECK_NEED_GETOPT_OPTRESET],
[
    AC_CACHE_CHECK( [whether 'optreset' needed for 'getopt' use],

        [hc_cv_need_getopt_optreset],
        [
            AC_TRY_LINK(
                [int getopt(int, char *const[], const char *);],
                [
                    extern int optreset;
                    optreset=1;
                    getopt(0,0,0);
                ],
                [hc_cv_need_getopt_optreset=yes],
                [hc_cv_need_getopt_optreset=no]
            )
        ]
    )

    if test "$hc_cv_need_getopt_optreset" = "yes"; then

        AC_DEFINE([NEED_GETOPT_OPTRESET])
    fi
])

#-----------------------------------------------------------------------------
#
# Macro:  HC_ADD_TO_CFLAGS_IF_SUPPORTED()
#
#  Checks whether the compiler supports 'option', adds to CFLAGS if yes
#
# Parms:   1. 'option' to be added
#          2. 'symbol' to be set
# Input:   nothing
# Output:  'option' is appended to CFLAGS if the compiler supports it;
#          'symbol' is set to "yes" or "no"
#
#-----------------------------------------------------------------------------

AC_DEFUN([HC_ADD_TO_CFLAGS_IF_SUPPORTED],
[
    AC_CACHE_CHECK( [whether compiler supports $1],
        [$2],
        [
            hc_temp="$CFLAGS"
            CFLAGS="$CFLAGS -Werror $1"
            AC_TRY_COMPILE([],
                [return(0);],
                $2="yes",
                $2="no")
            CFLAGS="$hc_temp"
        ])
    if test "${$2}" = "yes"; then
        CFLAGS="$CFLAGS $1"
    fi
])

###############################################################################
#   (end-of-file)
###############################################################################
