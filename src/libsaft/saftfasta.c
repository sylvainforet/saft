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

#define READ_CHUNK       1024
#define SEQ_CHUNK        256
#define NAME_CHUNK       128
#define STRUCT_CHUNK     16
#define ENSURE(buf, buf_idx, chk_idx, nb_units) (buf_idx)++;                                              \
                                                if ((buf_idx) == (chk_idx))                               \
                                                  {                                                       \
                                                    (chk_idx) += (nb_units);                              \
                                                    (buf) = realloc ((buf), (chk_idx) * sizeof (*(buf))); \
                                                  }
#define CAP_STR(str, buf_idx, chk_idx) ENSURE (str, buf_idx, chk_idx, 1); \
                                       (str)[(buf_idx)] = '\0';

SaftFasta**
saft_fasta_read (const char   *filename,
                 unsigned int *n)
{
  SaftFastaParseData data;

  data.seqs  = NULL;
  data.idx   = -1;
  data.alloc =  0;

  saft_fasta_iter (filename,
                   (SaftFastaIterFunc)saft_fasta_append,
                   &data);

  ENSURE (data.seqs, data.idx, data.alloc, 1);
  data.seqs[data.idx] = NULL;
  if (n)
    *n = data.idx;

  return data.seqs;
}

static int
saft_fasta_append (SaftFasta          *fasta,
                   SaftFastaParseData *data)
{
  SaftFasta *seq_new;

  seq_new       = saft_fasta_new ();
  seq_new->name = fasta->name;
  seq_new->seq  = fasta->seq;
  fasta->name   = NULL;
  fasta->seq    = NULL;

  ENSURE (data->seqs, data->idx, data->alloc, STRUCT_CHUNK);
  data->seqs[data->idx] = seq_new;

  return 1;
}

void
saft_fasta_iter (const char        *filename,
                 SaftFastaIterFunc  func,
                 void              *data)
{
  char        buffer[READ_CHUNK];
  SaftFasta *seq              = NULL;
  int         cur_name_alloc   =  0;
  int         cur_name_idx     = -1;
  int         cur_seq_alloc    =  0;
  int         cur_seq_idx      = -1;
  int         in               = -1;
  int         status           = -1;
  char        in_header        =  0;

  if ((in = open (filename, O_RDONLY | O_NONBLOCK)) == -1)
    {
      saft_error ("Couldn't open `%s'", filename);
      return;
    }

  while ((status = read (in, buffer, READ_CHUNK)) > 0)
    {
      int i;

      for (i = 0 ; i < status; i++)
        {
          const char ch = buffer[i];

          if (ch == '>')
            {
              if (in_header)
                {
                  ENSURE (seq->name,
                          cur_name_idx, cur_name_alloc, NAME_CHUNK);
                  seq->name[cur_name_idx] = ch;
                  continue;
                }
              if (seq)
                {
                  CAP_STR (seq->name, cur_name_idx, cur_name_alloc);
                  CAP_STR (seq->seq, cur_seq_idx, cur_seq_alloc);
                  if (!func (seq, data))
                    {
                      saft_fasta_free (seq);
                      close(in);
                      return;
                    }
                  saft_fasta_free (seq);
                }
              cur_name_alloc =  0;
              cur_name_idx   = -1;
              cur_seq_alloc  =  0;
              cur_seq_idx    = -1;
              seq            = saft_fasta_new ();
              in_header      = 1;
            }
          else
            {
              if (!seq)
                continue;
              if (in_header)
                {
                  if ('\n' == ch ||
                      '\r' == ch)
                    {
                      in_header = 0;
                      continue;
                    }
                  ENSURE (seq->name,
                          cur_name_idx, cur_name_alloc, NAME_CHUNK);
                  seq->name[cur_name_idx] = ch;
                }
              else
                {
                  if ('\n' == ch ||
                      '\r' == ch ||
                      '\t' == ch ||
                      ' '  == ch)
                    continue;
                  ENSURE (seq->seq,
                          cur_seq_idx, cur_seq_alloc, SEQ_CHUNK);
                  seq->seq[cur_seq_idx] = ch;
                }
            }
        }
    }
  if (status == -1)
    saft_error ("An IO error occured while reading `%s'", filename);
  if (seq)
    {
      CAP_STR (seq->name, cur_name_idx, cur_name_alloc);
      CAP_STR (seq->seq, cur_seq_idx, cur_seq_alloc);
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
          saft_error ("Encountered symbol `%d' unknown in alphabet `%s'",
                      *tmp_s, alphabet->name);
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
