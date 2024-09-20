#ifndef SEQUENCE_AVL_TREE_H
#define SEQUENCE_AVL_TREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// Node structure for AVL tree
typedef struct Node
{
    char *sequence;  // The label or sequence name
    int start_coord; // Start position of the segment
    int end_coord;   // End position of the segment
    int length;      // Length of the segment
    char *type;      // Feature type (e.g., CDS, promoter)
    char *gene_name; // Gene name (optional)

    int height;   // Height of the AVL tree node
    int children; // Number of children nodes

    struct Node *left;   // Left child
    struct Node *right;  // Right child
    struct Node *parent; // Parent node

    int isNested; // Boolean flag to track if this segment is nested

} Node;

// BST structure for the AVL tree
typedef struct BST
{
    Node *root;            // Root node of the tree
    int sequenceLength;    // Total sequence length from the GenBank file

    int *labeled_starts;   // Dynamic array for CDS start positions
    int *unlabeled_starts; // Dynamic array for unlabeled start positions
    int labeled_count;     // Number of labeled sections
    int unlabeled_count;   // Number of unlabeled sections
} BST;

// BST General Functions
BST *newBST(); // Create a new binary search tree
void freeBST(BST *tree); // Free memory of tree
void freeNodes(Node *node); // Free memory of nodes in tree
Node *newNode(char *sequence, int start, int end, const char *type, char *gene_name); // Create a new BST node
Node *find(BST *tree, int start_coord); // Find a node in the tree by start coordinate

// AVL Tree Balance and Insertion
int height(Node *node);
int getBalance(Node *node);
Node *rightRotate(Node *node);
Node *leftRotate(Node *node);
Node *insert(Node *node, Node *new_node, BST *tree);

// Formatted Print Functions
void printStartPositions(BST *tree);
void printNodeDetailsInOrder(Node *root); // Prints tree segments in sorted order - used for range query output
void printAVLTreeSegments(Node *root);

// Misc. Helper Functions
int sumNodeLengths(Node *root);

// Node Structure for SeqHashTable - stores segments and their details
typedef struct SeqHashNode
{
    char *key;                // Feature type (e.g., CDS, promoter)
    int *startPositions;      // Array of start positions for each feature
    char **labels;            // Array of labels for each segment
    int startCount;           // Number of start positions
    int labelCount;           // Number of labels
    struct SeqHashNode *next; // For handling collisions in the hash table
} SeqHashNode;

// SeqHashTable Structure
typedef struct SeqHashTable
{
    SeqHashNode **table; // Array of pointers to SeqHashNodes (buckets)
    int size;            // Size of the hash table
} SeqHashTable;

// Parse Input and Insert Functions
void parseCSVAndInsertIntoAVL(BST *tree, SeqHashTable *hashTable, const char *csv_filename);

// SeqHashTable General Functions
unsigned int hash(const char *key, int table_size);
SeqHashTable *createSeqHashTable(int size);
void freeSeqHashTable(SeqHashTable *table);
void SeqHashInsert(SeqHashTable *table, const char *featureType, int startPosition, const char *label);
void insertUnlabeledSegmentToSeqHashTable(SeqHashTable *hashTable, int startPosition);
void printSeqHashTable(SeqHashTable *table, BST *tree);
SeqHashNode *searchSeqHashTable(SeqHashTable *table, const char *key);

// TreeHashNode for Managing AVL Trees
typedef struct TreeHashNode
{
    char *treeName;            // Name of the AVL tree
    BST *tree;                 // The AVL tree itself
    struct TreeHashNode *next; // For handling collisions
} TreeHashNode;

// TreeHashTable for Storing AVL Trees
typedef struct TreeHashTable
{
    TreeHashNode **table; // Array of pointers to TreeHashNode (buckets)
    int size;             // Size of the hash table
} TreeHashTable;

// TreeHashTable General Functinons
TreeHashTable *createTreeHashTable(int size);
void freeTreeHashTable(TreeHashTable *table);
void treeHashTableInsert(TreeHashTable *table, const char *treeName, BST *tree);
BST *searchTreeHashTable(TreeHashTable *table, const char *treeName);
void printTreeHashTable(TreeHashTable *table);

// Range Query Helper Functions
void collectAndPrintSegmentsInRange(Node *node, int x, int y);
void rangeQuery(BST *tree, SeqHashTable *table, BST *original);
int getLefttMostStartPosition(Node *root);
int getRightMostStartPosition(Node *root);

// Filter Tree Functions
BST *filterTreeByStartPositions(BST *originalTree, int *startPositionsToExclude, int excludeCount, SeqHashTable *seqHashTable, TreeHashTable *treeHashTable);
void filterInsertHelper(Node *currentNode, BST *newTree, int *startPositionsToExclude, int excludeCount, int **unlabeledStartPositions, int *unlabeledCount);
int isValidStartPosition(BST *tree, int startPosition);

// User Interface
void clearInputBuffer();
void promptForCommand(char *command);
void promptForRangeQueryBounds(SeqHashTable *hashTable, int *x, int *y, BST *tree);
int isInteger(const char *str);
int *promptForStartPositionsToExclude(int *excludeCount);
char *promptForUniqueTreeName(TreeHashTable *treeHashTable);
BST *promptUserForTreeSelection(TreeHashTable *treeHashTable);

#endif // SEQUENCE_AVL_TREE_H
