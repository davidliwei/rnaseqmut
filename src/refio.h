/* Rading and writing reference sequences */
#ifndef REFIO_H
#define REFIO_H

#include <cstdio>
#include <string>
#include <iostream>

using namespace std;


/* Initialize the reference sequence */
int refseq_init(string filename);

/* Finish loading the reference sequence */
int refseq_destroy();


/* Get the sequence content */
int refseq_getseq(string chrname, long a, int len, string& seq);

#endif
