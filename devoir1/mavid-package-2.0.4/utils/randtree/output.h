
#ifndef MAVID_OUTPUT_H
#define MAVID_OUTPUT_H

enum { MAVID_BAD_ARGS, MAVID_BAD_TREE, MAVID_NON_BINARY, MAVID_NO_SEQ_FILE,
       MAVID_BAD_SEQ_FILE, MAVID_TREE_SEQ_DIFF, MAVID_NO_MASK_FILE,
       MAVID_BAD_MASK_FILE };

enum { SUCCESS, ERROR };

void splash();

void usage( );

void quit( int status );

#endif
