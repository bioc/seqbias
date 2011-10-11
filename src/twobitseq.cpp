/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include "twobitseq.hpp"
#include <algorithm>
#include <cstring>
#include <cstdlib>

kmer nuc_to_num(char c)
{
    switch( c ) {
        case 'a': case 'A':
        case 'u': case 'U': return 0;
        case 'c': case 'C': return 1;
        case 'g': case 'G': return 2;
        case 't': case 'T': return 3;
                  /* assume the N case is very rare, and just
                   * return a random nucleotide */
        default: return rand()%4;
    }
}



void num_to_nuc(char* dest, kmer K, int k)
{
    int i;
    /* read backwards, then reverse */
    for (i = 0; i < k; ++i) {
        switch (K & 0x3) {
            case 0: dest[i] = 'a'; break;
            case 1: dest[i] = 'c'; break;
            case 2: dest[i] = 'g'; break;
            case 3: dest[i] = 't'; break;
        }

        K >>= 2;
    }

    dest[i] = '\0';
    std::reverse(dest, dest + i);
}


const size_t twobitseq::max_kmer = 4 * sizeof(kmer);

twobitseq::twobitseq()
    : xs(NULL)
    , n(0)
{

}


twobitseq::twobitseq(const char* seq)
    : xs(NULL)
    , n(0)
{
    if (seq == NULL) return;
    n = strlen(seq);
    if (n == 0) return;

    xs = reinterpret_cast<kmer*>(
            malloc_or_die((n / max_kmer + 1) * sizeof(kmer)));
    memset(xs, 0, (n / max_kmer + 1) * sizeof(kmer));

    size_t i;
    size_t block, offset;
    for (i = 0; i < n; ++i) {
        block  = i / max_kmer;
        offset = i % max_kmer;

        xs[block] = xs[block] | (nuc_to_num(seq[i]) << (2*offset));
    }
}



twobitseq::twobitseq(const twobitseq& other)
    : xs(NULL)
    , n(0)
{
    n = other.n;
    if (n == 0) return;

    xs = reinterpret_cast<kmer*>(
            malloc_or_die((n / max_kmer + 1) * sizeof(kmer)));
    memcpy(xs, other.xs, (n / max_kmer + 1) * sizeof(kmer));
}



twobitseq::~twobitseq()
{
    free(xs);
}


void twobitseq::operator = (const twobitseq& other)
{
    n = other.n;
    xs = reinterpret_cast<kmer*>(
            realloc_or_die(xs, (n / max_kmer + 1) * sizeof(kmer)));
    memcpy(xs, other.xs, (n / max_kmer + 1) * sizeof(kmer));
}


void twobitseq::operator = (const char* seq)
{
    if (seq == NULL) {
        n = 0;
        free(xs);
        xs = NULL;
        return;
    }

    n = strlen(seq);
    xs = reinterpret_cast<kmer*>(
            realloc_or_die(xs, (n / max_kmer + 1) * sizeof(kmer)));
    memset(xs, 0, (n / max_kmer + 1) * sizeof(kmer));

    size_t i;
    size_t block, offset;
    for (i = 0; i < n; ++i) {
        block  = i / max_kmer;
        offset = i % max_kmer;

        xs[block] = xs[block] | (nuc_to_num(seq[i]) << (2*offset));
    }
}


kmer twobitseq::get_kmer(int k, pos_t pos)
{
    size_t block, offset;
    size_t i;
    kmer K = 0;

    for (i = 0; i < (size_t) k; ++i) {
        block  = (i + (size_t) (pos - (k - 1))) / max_kmer;
        offset = (i + (size_t) (pos - (k - 1))) % max_kmer;
        K = (K << 2) | ((xs[block] >> (2 * offset)) & 0x3);
    }

    return K;
}


int twobitseq::make_kmer(kmer& K, size_t pos, bool* mask, size_t mask_len) const
{
    int k = 0;
    size_t block, offset;
    K = 0;
    size_t i;
    for (i = 0; i < mask_len; ++i) {
        if (!mask[i]) continue;

        block  = (i + pos) / max_kmer;
        offset = (i + pos) % max_kmer;

        K = (K << 2) | ((xs[block] >> (2 * offset)) & 0x3);
        ++k;
    }

    return k;
}




