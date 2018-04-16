/* HCONSOLE.H   (c) Copyright "Fish" (David B. Trout), 2009          */
/*              Hercules console panel support functions header file */

//////////////////////////////////////////////////////////////////////////////////////////
// (c) Copyright "Fish" (David B. Trout), 2009. Released under the Q Public License
// (http://www.hercules-390.org/herclic.html) as modifications to Hercules.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef _HCONSOLE_H
#define _HCONSOLE_H

//-----------------------------------------------------------------------------
//
//                            VT100 User Guide
//
//                   Table 3-6 Cursor Control Key Codes
//
// Cursor Key   VT52       ANSI and Cursor Key Mode   ANSI and Cursor Key Mode
//  (Arrow)     Mode       Reset    (Normal mode)     Set  (Application mode)
//
//   Up         ESC A      ESC [ A                    ESC O A
//   Down       ESC B      ESC [ B                    ESC O B
//   Right      ESC C      ESC [ C                    ESC O C
//   Left       ESC D      ESC [ D                    ESC O D
//
//-----------------------------------------------------------------------------

#define ANSI_RESET_ALL_ATTRIBUTES  "\x1B[0m"

#define KBD_HOME                "\x1B[1~"
#define KBD_INSERT              "\x1B[2~"
#define KBD_DELETE              "\x1B[3~"
#define KBD_END                 "\x1B[4~"
#define KBD_PAGE_UP             "\x1B[5~"
#define KBD_PAGE_DOWN           "\x1B[6~"

#define KBD_UP_ARROW            "\x1B[A"
#define KBD_DOWN_ARROW          "\x1B[B"
#define KBD_RIGHT_ARROW         "\x1B[C"
#define KBD_LEFT_ARROW          "\x1B[D"

#define KBD_UP_ARROW2           "\x1BOA"
#define KBD_DOWN_ARROW2         "\x1BOB"
#define KBD_RIGHT_ARROW2        "\x1BOC"
#define KBD_LEFT_ARROW2         "\x1BOD"

// Does anyone know what the actual escape sequence that
// gets generated actually is on Linux for "Alt+UpArrow"
// and "Alt+DownArrow", etc?? Thanks! -- Fish

#define KBD_ALT_UP_ARROW        KBD_UP_ARROW2
#define KBD_ALT_DOWN_ARROW      KBD_DOWN_ARROW2
#define KBD_ALT_RIGHT_ARROW     KBD_RIGHT_ARROW2
#define KBD_ALT_LEFT_ARROW      KBD_LEFT_ARROW2

#define KBD_CTRL_HOME           "\x1B""w"           // (is this right??)
#define KBD_CTRL_END            "\x1B""u"           // (is this right??)

#define KBD_CTRL_UP_ARROW       "\x1B""D"           // (is this right??)
#define KBD_CTRL_DOWN_ARROW     "\x1B""M"           // (is this right??)
//efine KBD_CTRL_RIGHT_ARROW    "???????"           // (luckily we don't need it right now)
//efine KBD_CTRL_LEFT_ARROW     "???????"           // (luckily we don't need it right now)

#define KBD_ASK_CURSOR_POS      "\x1B[6n"           // Return value is the string "\x1B[n;mR"
                                                    // returned in the keyboard buffer where
                                                    // n = decimal row, m = decimal column.

// Hercules console color codes...

#define  COLOR_BLACK           0
#define  COLOR_RED             1
#define  COLOR_GREEN           2
#define  COLOR_BLUE            3
#define  COLOR_CYAN            4
#define  COLOR_MAGENTA         5
#define  COLOR_YELLOW          6
#define  COLOR_DARK_GREY       7
#define  COLOR_LIGHT_GREY      8
#define  COLOR_LIGHT_RED       9
#define  COLOR_LIGHT_GREEN     10
#define  COLOR_LIGHT_BLUE      11
#define  COLOR_LIGHT_CYAN      12
#define  COLOR_LIGHT_MAGENTA   13
#define  COLOR_LIGHT_YELLOW    14
#define  COLOR_WHITE           15
#define  COLOR_DEFAULT_FG      16
#define  COLOR_DEFAULT_BG      17
#define  COLOR_DEFAULT_LIGHT   18

extern  int   set_screen_color ( FILE* confp, short herc_fore,  short herc_back  );

// screen positions are 1-based; row 1 == top line; col 1 == leftmost column

extern  int   set_screen_pos   ( FILE* confp, short rowY1, short colX1 );

extern  int   clear_screen     ( FILE* confp );
extern  int   erase_to_eol     ( FILE* confp );

// 'save_and_set' = 1 --> just what it says; 0 --> restore from saved value.
extern  int   set_or_reset_console_mode ( int keybrd_fd, short save_and_set );

extern  void  translate_keystroke( char kbbuf[], int* pkblen );

extern  int   console_beep( FILE* confp );
extern  int   get_console_dim( FILE* confp, int* rows, int* cols );

/* Set by SIGWINCH handler if our window size changed */
extern  int  window_changed;

#ifdef OPTION_EXTCURS
extern  int   get_cursor_pos( int keybrd_fd, FILE* confp, short* row, short* col );
#endif // OPTION_EXTCURS
extern  int   set_console_cursor_shape( FILE* confp, int ins );

#if defined( _MSVC_ )
extern  int   w32_set_console_title( char* pszTitle );
#endif

#endif // _HCONSOLE_H
