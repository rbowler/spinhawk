/* $OpenBSD: rijndael.h,v 1.11 2005/05/25 05:47:53 markus Exp $ */
/* modified for use by dyncrypt */

// $Id$

/**
 * rijndael-alg-fst.h
 *
 * @version 3.0 (December 2000)
 *
 * Optimised ANSI C code for the Rijndael cipher (now AES)
 *
 * @author Vincent Rijmen <vincent.rijmen@esat.kuleuven.ac.be>
 * @author Antoon Bosselaers <antoon.bosselaers@esat.kuleuven.ac.be>
 * @author Paulo Barreto <paulo.barreto@terra.com.br>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// $Log$

#ifndef __RIJNDAEL_H
#define __RIJNDAEL_H

#define MAXKC   (256/32)
#define MAXKB   (256/8)
#define MAXNR   14

typedef u_int8_t        u8;
typedef u_int16_t       u16;
typedef u_int32_t       u32;

/*  The structure for key information */
typedef struct {
        int     enc_only;               /* context contains only encrypt schedule */
        int     Nr;                     /* key-length-dependent number of rounds */
        u32     ek[4*(MAXNR + 1)];      /* encrypt key schedule */
        u32     dk[4*(MAXNR + 1)];      /* decrypt key schedule */
} rijndael_ctx;

int      rijndael_set_key(rijndael_ctx *, u_char *, int);
int      rijndael_set_key_enc_only(rijndael_ctx *, u_char *, int);
void     rijndael_decrypt(rijndael_ctx *, u_char *, u_char *);
void     rijndael_encrypt(rijndael_ctx *, u_char *, u_char *);

int     rijndaelKeySetupEnc(unsigned int [], const unsigned char [], int);
int     rijndaelKeySetupDec(unsigned int [], const unsigned char [], int);
void    rijndaelEncrypt(const unsigned int [], int, const unsigned char [16],
            unsigned char [16]);

/* Additional definitions for dyncrypt */
typedef rijndael_ctx    aes_context;
#define aes_set_key     rijndael_set_key
#define aes_encrypt     rijndael_encrypt
#define aes_decrypt     rijndael_decrypt

#endif /* __RIJNDAEL_H */
