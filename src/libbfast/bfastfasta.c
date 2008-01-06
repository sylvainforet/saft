/**
 *
 */

#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bfasterror.h"
#include "bfastfasta.h"


typedef struct _BfastFastaParseData BfastFastaParseData;

struct _BfastFastaParseData
{
  BfastFasta **seqs;
  int          idx;
  int          alloc;
};

static int bfast_fasta_append (BfastFasta          *fasta,
                               BfastFastaParseData *data);


BfastFasta*
bfast_fasta_new ()
{
  BfastFasta *fasta;

  fasta       = malloc (sizeof (*fasta));
  fasta->name = NULL;
  fasta->seq  = NULL;

  return fasta;
}

void
bfast_fasta_free (BfastFasta *fasta)
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

BfastFasta**
bfast_fasta_read (const char   *filename,
                  unsigned int *n)
{
  BfastFastaParseData data;

  data.seqs  = NULL;
  data.idx   = -1;
  data.alloc =  0;

  bfast_fasta_iter (filename,
                    (BfastFastaIterFunc)bfast_fasta_append,
                    &data);

  ENSURE (data.seqs, data.idx, data.alloc, 1);
  data.seqs[data.idx] = NULL;
  if (n)
    *n = data.idx;

  return data.seqs;
}

static int
bfast_fasta_append (BfastFasta          *fasta,
                    BfastFastaParseData *data)
{
  BfastFasta *seq_new;

  seq_new       = bfast_fasta_new ();
  seq_new->name = fasta->name;
  seq_new->seq  = fasta->seq;
  fasta->name   = NULL;
  fasta->seq    = NULL;

  ENSURE (data->seqs, data->idx, data->alloc, STRUCT_CHUNK);
  data->seqs[data->idx] = seq_new;

  return 1;
}

void
bfast_fasta_iter (const char        *filename,
                  BfastFastaIterFunc func,
                  void              *data)
{
  char        buffer[READ_CHUNK];
  BfastFasta *seq              = NULL;
  int         cur_name_alloc   =  0;
  int         cur_name_idx     = -1;
  int         cur_seq_alloc    =  0;
  int         cur_seq_idx      = -1;
  int         in               = -1;
  int         status           = -1;
  char        in_header        =  0;

  if ((in = open (filename, O_RDONLY | O_NONBLOCK)) == -1)
    {
      bfast_error ("Couldn't open `%s'", filename);
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
                      bfast_fasta_free (seq);
                      close(in);
                      return;
                    }
                  bfast_fasta_free (seq);
                }
              cur_name_alloc =  0;
              cur_name_idx   = -1;
              cur_seq_alloc  =  0;
              cur_seq_idx    = -1;
              seq            = bfast_fasta_new ();
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
    bfast_error ("An IO error occured while reading `%s'", filename);
  if (seq)
    {
      CAP_STR (seq->name, cur_name_idx, cur_name_alloc);
      CAP_STR (seq->seq, cur_seq_idx, cur_seq_alloc);
      func (seq, data);
      bfast_fasta_free (seq);
    }
  close (in);
}

BfastSequence*
bfast_fasta_to_seq (BfastFasta    *fasta,
                    BfastAlphabet *alphabet)
{
  BfastSequence *seq;
  unsigned char *tmp_f;
  BfastLetter   *tmp_s;

  seq           = bfast_sequence_new ();
  seq->alphabet = alphabet;
  seq->name     = strdup (fasta->name);
  seq->size     = strlen (fasta->seq);
  seq->seq      = malloc (seq->size * sizeof (*seq->seq));
  tmp_f         = (unsigned char *)fasta->seq;
  tmp_s         = seq->seq - 1;

  while (*tmp_f)
    {
      *++tmp_s = alphabet->codes[*tmp_f++];
      if (*tmp_s == 0)
        bfast_error ("Encountered symbol `%d' unknown in alphabet `%s'",
                     *tmp_s, alphabet->name);
    }
  return seq;
}
