/***************************************************************************
                          gmx_rescue.c  -  description
                             -------------------
    begin                : Fri Apr  5 20:17:10 BST 2002
    origin               : from postings on Gromacs Mailing list by
                           Bert de Groot and Peter Tieleman
    modifications        : (C) 2003 by Marc Baaden
    email                : baaden@smplinux.de
 ***************************************************************************/
/*
 * $Id: gmx_rescue.c,v 1.4 2003/08/29 14:45:22 marc Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.0
 *
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>

// Debug libraries
//#include <libexplain/ferror.h>
#include <unistd.h>

// .xtc files are big-endian.
// Note: endianness is a CPU property so the original naming convention is not technically correct.
// Linux magic number. C on x86-style Linux defaults to little-endian.

#define MAGIC_LITTLE_ENDIAN -888733696
// SGI magic number. C on SGI defaults to Big-endian
#define MAGIC_BIG_ENDIAN 1995

#define BUF_SIZE 1024*32
#define VAR_SIZE 4

// Source: https://www.geeksforgeeks.org/write-an-efficient-c-program-to-reverse-bits-of-a-number
unsigned int reverse_bits(unsigned int num)
{
    unsigned int count = sizeof(num) * 8 - 1;
    unsigned int reverse_num = num;
      
    num >>= 1; 
    while(num)
    {
       reverse_num <<= 1;       
       reverse_num |= num & 1;
       num >>= 1;
       count--;
    }
    reverse_num <<= count;
    return reverse_num;
}

// source: https://codereview.stackexchange.com/questions/151049/endianness-conversion-in-cnsigned int reverse_endain
inline uint32_t reverse_endian(uint32_t value) 
{
    return (((value & 0x000000FF) << 24) |
            ((value & 0x0000FF00) <<  8) |
            ((value & 0x00FF0000) >>  8) |
            ((value & 0xFF000000) >> 24));
}


// Define simple struct for getting frame information.
struct frameHeader {
    unsigned int nAtoms ;
    int step ;
    float timeStamp;
    // float box[3][3] ;
} ;

void print_header(struct frameHeader *h) {
    printf("%-12s %-12s %-12s\n", "nAtoms", "step", "timeStamp");
    printf("%12d %12d %12e\n", h->nAtoms, h->step, h->timeStamp);
    return;
}

struct frameHeader read_header(FILE *fp, bool bReverse) {
    struct frameHeader h ;
    int temp;
    fread( &h.nAtoms, VAR_SIZE, 1 , fp );
    fread( &h.step,   VAR_SIZE, 1 , fp );
    fread( &temp,     VAR_SIZE, 1 , fp );
    if ( bReverse ) {
         h.nAtoms = reverse_endian( h.nAtoms ); 
         h.step   = reverse_endian( h.step   ); 
         temp     = reverse_endian( temp     ); 
    }
    h.timeStamp = *(float*)&temp ;
    return h ;
}

// source: https://stackoverflow.com/questions/111928/is-there-a-printf-converter-to-print-in-binary-format
// little endian
void print_bits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;   
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
        printf(" ");
    }
    //puts("");
}

void debug_stream(FILE * fp, int items, bool bReverse) {
    int i, t;
    for (i=0;i<items;i++){
        fread(&t, VAR_SIZE, 1, fp);
        if ( feof(fp) || ferror(fp) ) {
            return ;
        }
        if ( bReverse ) {
            t = reverse_endian( t );
        }
        print_bits(sizeof(t), &t);
        printf(": %12d ", t );
        printf(": %12e\n", *(float*)&t );
    }
}


void print_usage(char *argv[]) {
    fprintf(stderr, "usage: %s scan   <inputFile> [outputFile] : finds frames by magic numbers, optional output to text file.\n", argv[0]);
    fprintf(stderr, "       %s rescue <inputFile> <outputFile> <Offset>: print file starting from offset\n", argv[0]);
    fprintf(stderr, "       %s auto   <inputFile> <outputFile> : attempt automatic repair assuming first frame is OK\n", argv[0]);
    fprintf(stderr, "       %s query  <inputFile> <Offset> [context]: print raw-byte data around offset\n", argv[0]);
}

int main(int argc,  char *argv[]) {

    off64_t  f, out;
    FILE *fp, *fpOut;
    int value, magic, i, nAtoms = -1;
    unsigned int temp;
    struct frameHeader header;
    bool bReverse;
    off64_t res, res2, frames, offset, offsetPrev;
    char buf[BUF_SIZE];

    if (argc<3) {
        print_usage(argv);
        return 10;
    }

    // open xtc file readonly
    fp = fopen(argv[2], "rb");
    if ( fp == NULL ) {
        perror("= = Couldn't open file for reading");
        return 20;
    }
    // Check endianness of CPU versus the big-endian of the XTC.
    res = fread( &value, VAR_SIZE, 1 , fp ) ;
    if ( feof(fp) || ferror(fp) ) {
        perror("= = File is either empty or other error has occurred");
        return 20;
    }
    switch( value ) {
        case MAGIC_LITTLE_ENDIAN :
            magic = MAGIC_LITTLE_ENDIAN ;
            bReverse = true ;
            break;
        case MAGIC_BIG_ENDIAN :          
            magic = MAGIC_BIG_ENDIAN ;
            bReverse = false ;
            break;
        default:
            perror("= = First value read does not correspond to the magic number for file integrity check!");
            return 20;
    }
    // Reset.
    fseek(fp, 0, SEEK_SET);
    
    // Scan for magic number positions
    if ( strcmp(argv[1],"scan") == 0 ) {
        // Optional file output
        if ( argc==4 ) {
            fpOut = fopen( argv[3], "w" );
        } else {
            fpOut = stdout;
        }
        fprintf(fpOut, "%-12s %-12s %-12s %-12s %-12s %-s\n", "Frame", "numAtoms", "Step", "Time", "Delta", "Offset");
        offsetPrev=0 ; frames=0; res=1;
        while ( !( feof(fp) || ferror(fp) ) ) {
            // This steps seems to randomly catch a comparison that clearly shouldn't be equal.
            res = fread( &value, VAR_SIZE, 1 , fp ) ;
            if ( value == magic ) {
                //debug_stream(fp, 3);
                header = read_header( fp, bReverse ) ;
                offset = ftell(fp)/VAR_SIZE - 4 ;
                if ( nAtoms==-1 ) {
                    nAtoms = header.nAtoms ; 
                }
                if ( nAtoms == header.nAtoms ) {
                    fprintf(fpOut, "%12d %12d %12d %12e %12ld %ld\n",
                            frames, header.nAtoms, header.step, header.timeStamp, offset-offsetPrev , offset );
                    frames++ ;
                    offsetPrev = offset ;
                } 
            }
        }
        if ( ferror(fp) ) {
            fprintf(stderr, "...File read error encountered!\n");          
        }
        if ( feof(fp) ) {
            fprintf(stderr, "...End of file reached.\n");
        }
        fclose(fp);
        fclose(fpOut);
        return 0;
    }

    // Old rescue mechanic that prints everything after the offset.
    if ( strcmp(argv[1], "rescue") == 0 ) {
        char *ptr;
        off64_t offsetNew;
        offset = strtol(argv[4], &ptr, 10);
        fseek(fp, offset*VAR_SIZE, SEEK_SET);
        fprintf(stdout,"Debug: offset is %ld , ftell is %ld\n", offset, ftell(fp));
        // read up to 4 bytes from file into value until magic number is found
        while ( res=fread(&value, VAR_SIZE, 1, fp)> 0 && value != magic) {
            ;
        }
        if ( ferror(fp) ) {
            fprintf(stderr, "ERROR: ...File read error encountered!\n");
            fprintf(stderr, "     ftell position: %ld\n", ftell(fp) );
            return 40;
        }
        if ( feof(fp) ) {
            fprintf(stderr, "ERROR: ...End of file reached first!\n");
            fprintf(stderr, "     ftell position: %ld\n", ftell(fp) );
            return 40;
        }
        offsetNew = ftell(fp)/VAR_SIZE - 1 ;
        printf("Offset of next available frame: %ld\n", offsetNew );
        header = read_header( fp, bReverse ) ;
        print_header( &header ) ; 
        fseek(fp, offsetNew*VAR_SIZE, SEEK_SET);
        // Test if output file exists.
        if( access( argv[3], F_OK ) == 0 ) {
            printf("= = ERROR: Please delete the existing output file %s.\n", argv[3]);
            return 30;                                                               
        }
    
        fpOut = fopen( argv[3], "wb" );
        while ( !( feof(fp) || ferror(fp) ) ) {
            //Note: res retuen the number of chars read, so this naturally terminates when EOF is reached.
            res  = fread(  buf, 1, BUF_SIZE, fp ) ;
            res2 = fwrite( buf, 1, res, fpOut) ;
        }
        fclose(fp);
        fclose(fpOut);
        return 0;
    }

    // Attempt automated repair?
    if ( strcmp(argv[1], "auto") == 0 ) {
        if( access( argv[3], F_OK ) == 0 ) {
            printf("= = ERROR: Please delete the existing output file %s.\n", argv[3]);
            return 30;                                                               
        }
        // The plan is to parse the input file by frames and identify frames that would be sanely added.
        fpOut = fopen( argv[3], "wb" );

        frames=0; res=1;
        while ( !( feof(fp) || ferror(fp) ) ) {
            res = fread( &value, VAR_SIZE, 1 , fp ) ;
            if ( value == magic ) {
                //debug_stream(fp, 3);
                header = read_header( fp, bReverse ) ;
                //print_header( &header ) ; 
                frames++;
            }
        }

        
        fclose(fp);
        fclose(fpOut);
        return 0;
    }

    // Manual Query for advanced checking. In case there's a seomthing funny going on with the headers that hasn't been foreseen.
    if ( strcmp(argv[1], "query") == 0 ) {
        char *ptr;
        int range = 10;
        off64_t offPre ;
        if ( argc == 5 ) {
            range=atoi(argv[4]);
        }
        offset = strtol(argv[3], &ptr, 10);
        // reposition file offset by setting it to offset bytes
        offPre = (offset > range) ? offset - range : 0 ;
        fseek(fp, offPre * VAR_SIZE, SEEK_SET);
        debug_stream(fp, (offset-offPre), bReverse);
        puts("-------- -------- -------- --------");
        debug_stream(fp, 1, bReverse);
        puts("-------- -------- -------- --------");
        debug_stream(fp, range, bReverse);

        fclose(fp);
        return 0;
    }

}
