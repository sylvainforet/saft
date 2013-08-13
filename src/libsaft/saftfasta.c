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
  SaftFasta **seqs;
  int         idx;
  int         alloc;
};

static int saft_fasta_append (SaftFasta          *fasta,
                              SaftFastaParseData *data);


SaftFasta*
saft_fasta_new ()
{
  SaftFasta *fasta;

  fasta       = malloc (sizeof (*fasta));
  fasta->name = NULL;
  fasta->seq  = NULL;

  return fasta;
}

void
saft_fasta_free (SaftFasta *fasta)
{
  if (fasta)
    {
      if (fasta->name)
        free (fasta->name);
      if (fasta->seq)
        free (fasta->seq);
      free (fasta);
    }
}

SaftFasta*
saft_fasta_copy (SaftFasta *fasta)
{
  SaftFasta *new_fasta;

  new_fasta       = saft_fasta_new ();
  new_fasta->name = strcpy (new_fasta->name, fasta->name);
  new_fasta->seq  = strcpy (new_fasta->seq, fasta->seq);

  return new_fasta;
}

SaftFasta**
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
saft_fasta_append (SaftFasta          *fasta,
                   SaftFastaParseData *data)
{
  SaftFasta *new_fasta;

  new_fasta = saft_fasta_copy (fasta);

  /* Add +1 to make space for the terminal NULL */
  if (data->idx + 1 >= data->alloc)
    {
      data->alloc <<= 1;
      data->seqs    = realloc (data->seqs, data->alloc * sizeof (*data->seqs));
    }
  data->seqs[data->idx] = new_fasta;
  data->idx++;

  return 1;
}

void
saft_fasta_iter (const char        *filename,
                 SaftFastaIterFunc  func,
                 void              *data)
{
  char        buffer[READ_CHUNK];
  SaftFasta  *seq              = NULL;
  int         cur_name_alloc   =  0;
  int         cur_name_idx     = -1;
  int         cur_seq_alloc    =  0;
  int         cur_seq_idx      = -1;
  int         in               = -1;
  int         status           = -1;
  char        in_header        =  0;
  char        started          =  0;

  if ((in = open (filename, O_RDONLY | O_NONBLOCK)) == -1)
    {
      saft_error ("Couldn't open `%s'", filename);
      return;
    }

  seq = saft_fasta_new ();
  cur_name_alloc = NAME_INIT_SIZE;
  cur_seq_alloc  = SEQ_INIT_SIZE;
  seq->name      = malloc (cur_name_alloc);
  seq->seq       = malloc (cur_seq_alloc);

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
          char   *end = memchr (start, '\n', READ_CHUNK);
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
              if (__builtin_expect (cur_name_alloc < cur_name_idx + size + 1, 0))
                {
                  while (cur_name_alloc < cur_name_idx + size + 1)
                    cur_name_alloc <<= 1;
                  seq->name = realloc (seq->name, cur_name_alloc);
                }
              if (has_eol)
                  in_header = 0;
              memcpy (seq->name + cur_name_idx + 1, start, size);
              cur_name_idx += size;
            }
          else
            {
              if (*start == '>')
                {
                  seq->name[cur_name_idx + 1] = '\0';
                  seq->seq[cur_seq_idx + 1] = '\0';
                  if(!func (seq, data))
                    {
                      saft_fasta_free (seq);
                      close(in);
                      return;
                    }
                  ++start;
                  in_header = 1;
                  continue;
                }
              /* Check against size + 1 to make room for the terminating '\0' */
              if (__builtin_expect (cur_seq_alloc < cur_seq_idx + size + 1, 0))
                {
                  while (cur_seq_alloc < cur_seq_idx + size + 1)
                    cur_seq_alloc <<= 1;
                  seq->seq = realloc (seq->seq, cur_seq_alloc);
                }
              memcpy (seq->seq + cur_seq_idx + 1, start, size);
              cur_seq_idx += size;
            }
          start = end + 1;
        }
    }
  if (status == -1)
    saft_error ("An IO error occured while reading `%s'", filename);
  if (started)
    {
      seq->name[cur_name_idx + 1] = '\0';
      seq->seq[cur_seq_idx + 1] = '\0';
      func (seq, data);
      saft_fasta_free (seq);
    }
  close (in);
}

SaftSequence*
saft_fasta_to_seq (SaftFasta    *fasta,
                   SaftAlphabet *alphabet)
{
  SaftSequence  *seq;
  SaftSegment   *segment = NULL;
  unsigned char *tmp_f;
  SaftLetter    *tmp_s;
  unsigned int   in_segment;

  seq           = saft_sequence_new ();
  seq->alphabet = alphabet;
  seq->name     = strdup (fasta->name);
  seq->size     = strlen (fasta->seq);
  seq->seq      = malloc (seq->size * sizeof (*seq->seq));
  tmp_f         = (unsigned char *)fasta->seq;
  tmp_s         = seq->seq - 1;

  in_segment = 0;
  while (*tmp_f)
    {
      *++tmp_s = alphabet->codes[*tmp_f++];
      if (*tmp_s == 0)
        {
          /* FIXME only issue an error when the number of unknown letters is `large'
          saft_error ("Encountered symbol `%d' unknown in alphabet `%s'",
                      *tmp_s, alphabet->name);
          */
          if (in_segment)
            in_segment = 0;
        }
      else
        {
          if (!in_segment)
            {
              in_segment    = 1;
              segment       = saft_segment_new ();
              segment->seq  = tmp_s;
              segment->next = seq->segments;
              seq->segments = segment;
            }
          segment->size++;
        }
    }
  return seq;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
