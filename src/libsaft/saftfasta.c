/* saftfasta.c
 * Copyright (C) 2008  Sylvain FORET
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *
 */

#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "safterror.h"
#include "saftfasta.h"

#define FASTA_CHUNK     256
#define READ_CHUNK     4096
#define NAME_INIT_SIZE  256
#define SEQ_INIT_SIZE  4096


typedef struct _SaftFastaParseData SaftFastaParseData;

struct _SaftFastaParseData
{
  SaftSequence **seqs;
  int            idx;
  int            alloc;
};

static int saft_fasta_append (SaftSequence       *seq,
                              SaftFastaParseData *data);


SaftSequence**
saft_fasta_read (const char   *filename,
                 unsigned int *n)
{
  SaftFastaParseData data;

  data.seqs  = malloc (FASTA_CHUNK * sizeof (*data.seqs));
  data.alloc =  FASTA_CHUNK;
  data.idx   = 0;

  saft_fasta_iter (filename,
                   (SaftFastaIterFunc)saft_fasta_append,
                   &data);

  data.seqs[data.idx] = NULL;
  if (n)
    *n = data.idx;

  return data.seqs;
}

static int
saft_fasta_append (SaftSequence       *seq,
                   SaftFastaParseData *data)
{
  SaftSequence *new_seq;

  new_seq = saft_sequence_copy (seq);

  /* Add +1 to make space for the terminal NULL */
  if (data->idx + 1 >= data->alloc)
    {
      data->alloc <<= 1;
      data->seqs    = realloc (data->seqs, data->alloc * sizeof (*data->seqs));
    }
  data->seqs[data->idx] = new_seq;
  data->idx++;

  return 1;
}

void
saft_fasta_iter (const char        *filename,
                 SaftFastaIterFunc  func,
                 void              *data)
{
  char           buffer[READ_CHUNK];
  SaftSequence  *seq              = NULL;
  int            in               = -1;
  int            status           = -1;
  char           in_header        =  0;
  char           started          =  0;

  if ((in = open (filename, O_RDONLY | O_NONBLOCK)) == -1)
    {
      saft_error ("Couldn't open `%s'", filename);
      return;
    }

  seq             = saft_sequence_new ();
  seq->name_alloc = NAME_INIT_SIZE;
  seq->seq_alloc  = SEQ_INIT_SIZE;
  seq->name       = malloc (seq->name_alloc);
  seq->seq        = malloc (seq->seq_alloc);

  while ((status = read (in, buffer, READ_CHUNK)) > 0)
    {
      char *max = buffer + status;
      char *start;

      start = buffer;
      if (!started)
        {
          start = memchr (buffer, '>', READ_CHUNK);
          if (!start)
            continue;

          started   = 1;
          in_header = 1;
          ++start;
        }

      while (start < max)
        {
          /* TODO handle '\r' */
          char   *end = memchr (start, '\n', max - start);
          size_t  size;
          int     has_eol = 0;

          if (end)
            has_eol = 1;
          else
            end = max;
          size = end - start;
          if (size == 0)
            {
              start = end + 1;
              continue;
            }

          if (in_header)
            {
              /* Check against size + 1 to make room for the terminating '\0' */
              if (seq->name_alloc < seq->name_length + size + 1)
                {
                  while (seq->name_alloc < seq->name_length + size)
                    seq->name_alloc <<= 1;
                  seq->name = realloc (seq->name, seq->name_alloc);
                }
              if (has_eol)
                  in_header = 0;
              memcpy (seq->name + seq->name_length, start, size);
              seq->name_length += size;
            }
          else
            {
              if (*start == '>')
                {
                  seq->name[seq->name_length] = '\0';
                  seq->seq[seq->seq_length]   = '\0';
                  if(!func (seq, data))
                    {
                      saft_sequence_free (seq);
                      close(in);
                      return;
                    }
                  seq->name_length = 0;
                  seq->seq_length  = 0;
                  ++start;
                  in_header = 1;
                  continue;
                }
              /* Check against size + 1 to make room for the terminating '\0' */
              if (seq->seq_alloc < seq->seq_length + size + 1)
                {
                  while (seq->seq_alloc < seq->seq_length + size)
                    seq->seq_alloc <<= 1;
                  seq->seq = realloc (seq->seq, seq->seq_alloc);
                }
              memcpy (seq->seq + seq->seq_length, start, size);
              seq->seq_length += size;
            }
          start = end + 1;
        }
    }
  if (status == -1)
    saft_error ("An IO error occured while reading `%s'", filename);
  if (started)
    {
      seq->name[seq->name_length] = '\0';
      seq->seq[seq->seq_length]   = '\0';
      func (seq, data);
      saft_sequence_free (seq);
    }
  close (in);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
