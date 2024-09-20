#include "sequenceAVLTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

extern int prev_start;
extern int prev_end;

////////////////////////////////////
// SeqHashTable General Functions //
////////////////////////////////////

unsigned int hash(const char *key, int table_size)
{
    unsigned long hash = 5381;
    int c;

    while ((c = *key++))
    {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }

    return hash % table_size;
}

SeqHashTable *createSeqHashTable(int size)
{
    SeqHashTable *newTable = (SeqHashTable *)malloc(sizeof(SeqHashTable));
    if (!newTable)
        return NULL;

    newTable->table = (SeqHashNode **)malloc(size * sizeof(SeqHashNode *));
    if (!newTable->table)
        return NULL;

    for (int i = 0; i < size; i++)
    {
        newTable->table[i] = NULL;
    }

    newTable->size = size;
    return newTable;
}

void freeSeqHashTable(SeqHashTable *table)
{
    if (table == NULL)
    {
        return;
    }

    // Loop through each slot in the hash table
    for (int i = 0; i < table->size; i++)
    {
        SeqHashNode *current = table->table[i]; // Start with the first node in the linked list

        // Traverse the linked list in each slot
        while (current != NULL)
        {
            SeqHashNode *temp = current; // Save current node

            // Free the feature type (key)
            if (temp->key != NULL)
            {
                free(temp->key);
            }

            // Free the array of start positions
            if (temp->startPositions != NULL)
            {
                free(temp->startPositions);
            }

            // Free the array of labels (if available)
            if (temp->labels != NULL)
            {
                for (int j = 0; j < temp->startCount; j++)
                {
                    if (temp->labels[j] != NULL)
                    {
                        free(temp->labels[j]); // Free each label string
                    }
                }
                free(temp->labels); // Free the labels array itself
            }

            current = current->next; // Move to the next node
            free(temp);              // Free the current node
        }
    }

    // Finally, free the table itself
    free(table->table);
    free(table);
}

void SeqHashInsert(SeqHashTable *table, const char *featureType, int startPosition, const char *label)
{
    int hashIndex = hash(featureType, table->size); // Use the featureType to calculate the hash index

    SeqHashNode *newNode = (SeqHashNode *)malloc(sizeof(SeqHashNode));
    newNode->key = strdup(featureType); // Store the feature type as the key

    // Allocate memory for start positions and labels
    newNode->startPositions = (int *)malloc(sizeof(int));
    newNode->startPositions[0] = startPosition; // Store the start position
    newNode->startCount = 1;                    // Only one start position for now

    if (label != NULL)
    {
        newNode->labels = (char **)malloc(sizeof(char *));
        newNode->labels[0] = strdup(label); // Store the label if provided
        newNode->labelCount = 1;
    }
    else
    {
        newNode->labels = NULL;
        newNode->labelCount = 0;
    }

    newNode->next = NULL;

    // Handle collisions (separate chaining)
    if (table->table[hashIndex] == NULL)
    {
        table->table[hashIndex] = newNode;
    }
    else
    {
        SeqHashNode *current = table->table[hashIndex];
        while (current->next != NULL)
        {
            current = current->next;
        }
        current->next = newNode;
    }
}

void insertUnlabeledSegmentToSeqHashTable(SeqHashTable *hashTable, int startPosition)
{
    SeqHashInsert(hashTable, "unlabeled", startPosition, NULL); // Using "unlabeled" as the feature type, with no label
}

void printSeqHashTable(SeqHashTable *table, BST *tree)
{
    for (int i = 0; i < table->size; i++)
    {
        SeqHashNode *node = table->table[i];

        // Use a set of keys (feature types) to avoid printing the same type multiple times
        while (node != NULL)
        {
            // Print feature type only once
            printf("Feature type: %s\n", node->key);
            printf("Segment labels and positions:\n");

            SeqHashNode *current = node;
            // Loop over each node in this chain (in case of collisions or multiple labels for the same feature type)
            while (current != NULL && strcmp(current->key, node->key) == 0)
            {
                // Print the corresponding labels and start positions
                for (int j = 0; j < current->startCount; j++)
                {
                    const char *label = (current->labels && current->labels[j]) ? current->labels[j] : "Unlabeled";
                    int end_coordinate = find(tree, current->startPositions[j])->end_coord;
                    printf("  %s - (%d, %d)\n", label, current->startPositions[j], end_coordinate);
                }
                current = current->next;
            }

            // Move to the next feature type
            node = current;
            printf("\n");
        }
    }
}

SeqHashNode *searchSeqHashTable(SeqHashTable *table, const char *key)
{
    int hashIndex = hash(key, table->size); // Compute the hash index for the key

    SeqHashNode *current = table->table[hashIndex]; // Get the first node at that index
    while (current != NULL)
    {
        if (strcmp(current->key, key) == 0) // Compare the keys to find a match
        {
            return current; // Return the node if found
        }
        current = current->next; // Continue searching in the chain if needed
    }

    return NULL; // Return NULL if no match is found
}

///////////////////////////
// Filter Tree Functions //
///////////////////////////

// Function to prompt the user for the start positions of segments to exclude
int *promptForStartPositionsToExclude(int *excludeCount)
{
    char input[1000]; // Buffer to store the user input
    int *startPositions = NULL;
    *excludeCount = 0;

    printf("Enter start positions of segments to exclude from input sequence - seperate with commas and without spaces(x,y,z): ");
    printf("\n");
    fgets(input, sizeof(input), stdin);
    input[strcspn(input, "\n")] = '\0'; // Remove newline character

    // Tokenize the string by commas
    char *token = strtok(input, ",");

    while (token != NULL)
    {
        if (!isInteger(token))
        {
            // printf("%s is not a valid start position and will be ignored.\n", token);
        }
        else
        {
            // Allocate or reallocate memory for the start positions array
            startPositions = realloc(startPositions, (*excludeCount + 1) * sizeof(int));

            // Convert the token to an integer and add it to the array
            startPositions[*excludeCount] = atoi(token);
            (*excludeCount)++;
        }

        // Get the next token
        token = strtok(NULL, ",");
    }

    return startPositions;
}

char *promptForUniqueTreeName(TreeHashTable *treeHashTable)
{
    // Allocate memory for the tree name buffer
    char *treeName = (char *)malloc(256 * sizeof(char));
    if (treeName == NULL)
    {
        printf("Memory allocation failed!\n");
        return NULL;
    }

    BST *existingTree = NULL;

    // Loop until a unique name is provided
    do
    {
        printf("\n\nEnter a name for the new filtered tree: ");
        fgets(treeName, 256, stdin); // Use fgets to read a line with spaces

        // Remove the trailing newline character that fgets captures
        treeName[strcspn(treeName, "\n")] = '\0';

        // Check if the tree name already exists in the tree hash table
        existingTree = searchTreeHashTable(treeHashTable, treeName);
        if (existingTree != NULL)
        {
            printf("Tree name already exists! Please enter a unique name.\n");
        }
    } while (existingTree != NULL); // Keep prompting until a unique name is given

    return treeName; // Return the valid tree name
}

// Helper function to check if a string is a valid integer
int isInteger(const char *str)
{
    // Check if the string is empty or null
    if (str == NULL || *str == '\0')
    {
        return 0;
    }

    // Check each character to ensure it's a digit
    for (int i = 0; str[i] != '\0'; i++)
    {
        if (!isdigit(str[i]))
        {
            return 0; // Return false if any character is not a digit
        }
    }

    return 1; // Return true if all characters are digits
}

// Function to check if the start position exists in the original tree
int isValidStartPosition(BST *tree, int startPosition)
{
    Node *node = find(tree, startPosition);
    if (node == NULL){
        // printf("%d is not a valid start position and will be ignored.\n", startPosition);
    }

    printf("\n");
    return node != NULL;
}

BST *filterTreeByStartPositions(BST *originalTree, int *startPositionsToExclude, int excludeCount, SeqHashTable *seqHashTable, TreeHashTable *treeHashTable)
{
    BST *filteredTree = newBST();
    int *unlabeledStartPositions = NULL;
    int unlabeledCount = 0;

    filterInsertHelper(originalTree->root, filteredTree, startPositionsToExclude, excludeCount, &unlabeledStartPositions, &unlabeledCount);

    // If no valid nodes were inserted into the new tree, return NULL
    if (filteredTree->root == NULL)
    {
        printf("Filtering failed: No valid nodes in the filtered tree.\n");
        free(filteredTree);
        free(unlabeledStartPositions);
        return NULL; // Indicate that filtering failed
    }

    free(unlabeledStartPositions);
    return filteredTree;
}

void filterInsertHelper(Node *currentNode, BST *newTree, int *startPositionsToExclude, int excludeCount, int **unlabeledStartPositions, int *unlabeledCount)
{
    if (currentNode == NULL)
    {
        return;
    }

    int exclude = 0;

    // Check if the current node's start position is in the exclude list
    for (int i = 0; i < excludeCount; i++)
    {
        if (currentNode->start_coord == startPositionsToExclude[i])
        {
            exclude = 1; // Mark as excluded
            printf("%d ", currentNode->start_coord);
            break;
        }
    }

    // If not excluded, insert into the new filtered tree
    if (!exclude)
    {
        Node *new_node = newNode(currentNode->sequence, currentNode->start_coord, currentNode->end_coord, currentNode->type, currentNode->gene_name);
        newTree->root = insert(newTree->root, new_node, newTree);
        // printf("inserted %d\n", new_node->start_coord);
    }

    // Traverse left and right subtrees
    filterInsertHelper(currentNode->left, newTree, startPositionsToExclude, excludeCount, unlabeledStartPositions, unlabeledCount);
    filterInsertHelper(currentNode->right, newTree, startPositionsToExclude, excludeCount, unlabeledStartPositions, unlabeledCount);
}

/////////////////////////////////////
// TreeHashTable General Functions //
/////////////////////////////////////

TreeHashTable *createTreeHashTable(int size)
{
    TreeHashTable *newTable = (TreeHashTable *)malloc(sizeof(TreeHashTable));
    if (!newTable)
        return NULL;

    newTable->table = (TreeHashNode **)malloc(size * sizeof(TreeHashNode *));
    if (!newTable->table)
    {
        free(newTable);
        return NULL;
    }

    for (int i = 0; i < size; i++)
    {
        newTable->table[i] = NULL;
    }

    newTable->size = size;
    return newTable;
}

void freeTreeHashTable(TreeHashTable *table)
{
    if (table == NULL)
    {
        return;
    }

    // Loop through each slot in the hash table
    for (int i = 0; i < table->size; i++)
    {
        TreeHashNode *current = table->table[i];

        // Traverse and free each linked node in the bucket (in case of collisions)
        while (current != NULL)
        {
            TreeHashNode *temp = current;

            // Free the tree stored in this node
            freeBST(temp->tree);

            // Free the tree name (key)
            if (temp->treeName != NULL)
            {
                free(temp->treeName);
            }

            current = current->next;
            free(temp); // Free the node itself
        }
    }

    // Finally, free the hash table array and the table structure itself
    free(table->table);
    free(table);
}

void treeHashTableInsert(TreeHashTable *table, const char *treeName, BST *tree)
{
    int hashIndex = hash(treeName, table->size); // Compute the hash index

    TreeHashNode *newNode = (TreeHashNode *)malloc(sizeof(TreeHashNode));
    if (!newNode)
    {
        printf("Memory allocation failed for tree hash node.\n");
        return;
    }

    // Initialize the new node
    newNode->treeName = strdup(treeName); // Duplicate the tree name
    newNode->tree = tree;                 // Assign the AVL tree
    newNode->next = NULL;                 // Next node is NULL (no collisions yet)

    // Handle collisions with separate chaining
    if (table->table[hashIndex] == NULL)
    {
        table->table[hashIndex] = newNode;
    }
    else
    {
        TreeHashNode *current = table->table[hashIndex];
        while (current->next != NULL)
        {
            current = current->next;
        }
        current->next = newNode;
    }
}

BST *searchTreeHashTable(TreeHashTable *table, const char *treeName)
{
    int hashIndex = hash(treeName, table->size); // Compute the hash index

    TreeHashNode *current = table->table[hashIndex];
    while (current != NULL)
    {
        if (strcmp(current->treeName, treeName) == 0)
        {
            return current->tree; // Return the AVL tree if the name matches
        }
        current = current->next; // Move to the next node in case of collisions
    }

    return NULL; // Tree with the given name not found
}

void printTreeHashTable(TreeHashTable *table)
{
    printf("Displaying all datasets:\n");

    for (int i = 0; i < table->size; i++)
    {
        TreeHashNode *node = table->table[i];

        while (node != NULL)
        {
            // Print the name of the tree
            printf("  %s\n", node->treeName);

            // Optionally, you could print the AVL tree structure or details here
            // printAVLTreeSegments(node->tree->root); // Example of printing tree segments

            node = node->next; // Move to the next node (collision handling)
        }
    }

    printf("\n");
}
