/************************************************************/
/************************************************************/
/**             GABIOS-LIB 1.0 (1999)                      **/
/** A library for Greedy Alignment of BIOlogical Sequences **/
/**             Developed by Said Abdeddaim                **/
/**          Said.Abdeddaim@dir.univ-rouen.fr              **/
/************************************************************/
/************************************************************/

#ifndef _ALIG_GRAPH_CLOSURE_H
#define _ALIG_GRAPH_CLOSURE_H

typedef struct
{
    int* pos;
    int nbr;
} positionSet;

typedef struct
{
    int longueur;

    /// might map positions to eq class numbers?
    int *aligSetNbr;
    /// maybe this is next class information
    int *predAligSetPos;
    /// maybe this is next class information
    int *succAligSetPos;
} sequence;

typedef struct
{
    int seqNbr;
    sequence* seq;
    int maxLong;

    /// this stores information about which eq class is aligned with which 
    /// positions for each sequence
    /// is used to tell if a eq class is aligned with a seq (not aligned if val == 0)
    positionSet* aligSet;
    /// how many eq ckasses are there
    int nbrAligSets, oldNbrAligSets;

    /// these are prob used to save the actual frontiers,
    /// they are indexed with a eq class id and a sequence id
    int **predFrontier, **succFrontier;

    int* topolog;
    /// gauche and droite might be next class information?
    /// gauche is left
    /// and droite is right
    int *gauche1, *gauche2, *droite1, *droite2;
    /// indexed with 2 seq numbers
    int **pos_;

} CLOSURE;

CLOSURE* newAligGraphClosure(int nbreseq, int* longseq, int nbreancr,
    int** ancrages);

void freeAligGraphClosure(CLOSURE* clos);

int addAlignedPositions(CLOSURE* clos, int x, int i, int y, int j);

int alignablePositions(CLOSURE* clos, int x, int i, int y, int j);

int alignedPositions(CLOSURE* clos, int x, int i, int y, int j);

int addAlignedSegments(CLOSURE* clos, int x, int i, int y, int j, int l);

int alignableSegments(CLOSURE* clos, int x, int i, int y, int j, int l);

int alignedSegments(CLOSURE* clos, int x, int i, int y, int j, int l);

int predFrontier(CLOSURE* clos, int x, int i, int y);

int succFrontier(CLOSURE* clos, int x, int i, int y);

#endif /* _ALIG_GRAPH_CLOSURE_H */
