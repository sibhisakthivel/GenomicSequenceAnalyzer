#include "sequenceAVLTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// Used to check for duplicate segments 
int prev_start = -1;
int prev_end = -1;   

////////////////////////////////////////
// Input Sequence Retrieval Functions //
////////////////////////////////////////

void runPythonParser(const char *input_file, const char *output_file)
{
    char command[512]; // Buffer for the command string

    // Build the command to run the Python script with input and output files
    snprintf(command, sizeof(command), "python3 parse_genbank.py %s %s", input_file, output_file);

    // Execute the command
    int result = system(command);

    // Check if the command executed successfully
    if (result != 0)
    {
        printf("GenBank input parsing failed. Exiting program.\n");
        exit(1);
    }
}

void parseCSVAndInsertIntoAVL(BST *tree, SeqHashTable *hashTable, const char *csv_filename)
{
    FILE *csv_file = fopen(csv_filename, "r");
    if (csv_file == NULL)
    {
        printf("Error: Unable to open CSV file %s\n", csv_filename);
        return;
    }

    char line[256];
    char feature_name[50];
    char feature_type[50];
    int start, end;
    int prev_end = 0;            // Tracks the end of the last inserted segment
    int first_labeled_found = 0; // To track if we've found the first labeled segment
    Node *new_node;

    while (fgets(line, sizeof(line), csv_file))
    {
        // Parse the CSV line into feature name, type, start, and end
        sscanf(line, "%[^,],%[^,],%d,%d", feature_name, feature_type, &start, &end);

        // Skip "source" type segments, but store the sequence length
        if (strcmp(feature_type, "source") == 0)
        {
            tree->sequenceLength = end;
            continue;
        }

        // Insert an unlabeled segment if this is the first labeled segment and it doesn't start from base 1
        if (!first_labeled_found && start > 1)
        {
            new_node = newNode("Unlabeled", 1, start - 1, "Unlabeled", NULL);
            tree->root = insert(tree->root, new_node, tree);

            // Insert the unlabeled segment into the hash table
            insertUnlabeledSegmentToSeqHashTable(hashTable, 1);

            first_labeled_found = 1;
        }

        // Check for overlaps and handle them
        if (prev_end >= start)
        {
            start = prev_end + 1;
        }

        // Insert an unlabeled segment if there is a gap between the previous segment and this one
        if (prev_end > 0 && prev_end + 1 < start)
        {
            new_node = newNode("Unlabeled", prev_end + 1, start - 1, "Unlabeled", NULL);
            tree->root = insert(tree->root, new_node, tree);

            // Insert the unlabeled segment into the hash table
            insertUnlabeledSegmentToSeqHashTable(hashTable, prev_end + 1);
        }

        // Use the label as the name, fallback to feature_type if label is missing
        char *segment_name = (strcmp(feature_name, "") != 0) ? feature_name : feature_type;

        // Insert the labeled feature
        new_node = newNode(segment_name, start, end, feature_type, segment_name);
        new_node->isNested = 0; // For now, assume it is not nested
        tree->root = insert(tree->root, new_node, tree);

        // Insert the labeled segment into the hash table
        SeqHashInsert(hashTable, feature_type, start, segment_name);

        // Update the previous segment's end
        prev_end = end;
    }

    // Insert a final unlabeled segment if there is space at the end
    if (prev_end < tree->sequenceLength)
    {
        new_node = newNode("Unlabeled", prev_end + 1, tree->sequenceLength, "Unlabeled", NULL);
        tree->root = insert(tree->root, new_node, tree);

        // Insert the final unlabeled segment into the hash table
        insertUnlabeledSegmentToSeqHashTable(hashTable, prev_end + 1);
    }

    fclose(csv_file);
}

///////////////////////////
// BST General Functions //
///////////////////////////

BST *newBST()
{
    BST *new_bst = (BST *)malloc(sizeof(BST)); // Allocate memory for tree
    if (new_bst != NULL)
    {
        new_bst->root = NULL; // Tree is empty

        // Initialize dynamic arrays for CDS and NCDS starts
        new_bst->labeled_starts = (int *)malloc(0);   // Start with an empty array
        new_bst->unlabeled_starts = (int *)malloc(0); // Start with an empty array

        new_bst->labeled_count = 0; // Initially no CDS sections
        new_bst->unlabeled_count = 0; // Initially no NCDS sections
    }
    return new_bst;
}

void freeBST(BST *tree) 
{
    if (tree == NULL || tree->root == NULL)
    {
        return;
    }
    freeNodes(tree->root); // Assuming freeNodes is a helper function to free nodes recursively
    free(tree);            // Free the tree structure itself
}

Node *newNode(char *sequence, int start, int end, const char *type, char *gene_name)
{
    Node *new_node = (Node *)malloc(sizeof(Node)); // Allocate memory for new BST node

    if (new_node != NULL)
    {
        new_node->sequence = strdup(sequence);
        new_node->start_coord = start;
        new_node->end_coord = end;
        new_node->length = end - start + 1;
        new_node->type = strdup(type);

        if (gene_name != NULL)
        {
            new_node->gene_name = strdup(gene_name);
        }
        else
        {
            new_node->gene_name = NULL;
        }

        new_node->height = 1;
        new_node->children = 1;
        new_node->left = NULL;
        new_node->right = NULL;
        new_node->parent = NULL;

        new_node->isNested = 0; // Initialize isNested to 0
    }

    return new_node;
}

void freeNodes(Node *node)
{
    if (node == NULL)
    {
        return;
    }

    freeNodes(node->left);
    freeNodes(node->right);

    // Check and free memory safely
    if (node->sequence != NULL)
    {
        free(node->sequence);
    }
    if (node->type != NULL)
    {
        free(node->type);
    }
    if (node->gene_name != NULL)
    {
        free(node->gene_name);
    }
    free(node);
}

Node *find(BST *tree, int start_coord)
{
    if (tree->root == NULL)
    {
        return NULL; // Empty tree
    }

    Node *current = tree->root;
    while (current != NULL)
    {
        // printf("Checking node: Start: %d\n", current->start_coord); // Debug print
        if (start_coord == current->start_coord)
        {
            return current; // Found the node with the matching start coordinate
        }
        else if (start_coord < current->start_coord)
        {
            current = current->left; // Move left
        }
        else
        {
            current = current->right; // Move right
        }
    }

    return NULL; // Node with the specified start_coord not found
}

///////////////////////////
// Range Query Functions //
///////////////////////////

void eliminateLeftSubtrees(Node *current, int start_coord, int *total)
{
    if (current == NULL)
    {
        return;
    }

    if (current->start_coord < start_coord)
    { // If current node's start coordinate is less than the threshold
        if (current->left != NULL)
        {
            *total -= current->left->children; // Subtract the size of the left subtree
        }
        *total -= 1;                                               // Subtract for the current node itself
        eliminateLeftSubtrees(current->right, start_coord, total); // Move to the right subtree
    }
    else
    {
        eliminateLeftSubtrees(current->left, start_coord, total); // Continue searching in the left subtree
    }
}

void eliminateRightSubtrees(Node *current, int start_coord, int *total)
{
    if (current == NULL)
    {
        return;
    }

    if (current->start_coord > start_coord)
    { // If current node's start coordinate is greater than the threshold
        if (current->right != NULL)
        {
            *total -= current->right->children; // Subtract the size of the right subtree
        }
        *total -= 1;                                               // Subtract for the current node itself
        eliminateRightSubtrees(current->left, start_coord, total); // Move to the left subtree
    }
    else
    {
        eliminateRightSubtrees(current->right, start_coord, total); // Continue searching in the right subtree
    }
}

// Helper function to print segment details
void printSegment(Node *node)
{
    const char *label = node->gene_name ? node->gene_name : "Unlabeled";
    printf("%d, %d, %s\n", node->start_coord, node->end_coord, label);
}

// Recursive function to collect and print segments within the range
void collectAndPrintSegmentsInRange(Node *node, int x, int y)
{
    if (node == NULL)
    {
        return;
    }

    // Traverse the left subtree if the node's start position is greater than or equal to x
    if (node->start_coord >= x)
    {
        collectAndPrintSegmentsInRange(node->left, x, y);
    }

    // If the current node is within the range, print it
    if (node->start_coord >= x && node->start_coord <= y)
    {
        printSegment(node);
    }

    // Traverse the right subtree if the node's start position is less than or equal to y
    if (node->start_coord <= y)
    {
        collectAndPrintSegmentsInRange(node->right, x, y);
    }
}

BST *promptUserForTreeSelection(TreeHashTable *treeHashTable)
{
    char treeName[50];
    BST *selectedTree = NULL;

    // Print available trees
    printTreeHashTable(treeHashTable);

    // Prompt user to select a tree
    do
    {
        printf("Enter the name of the tree you want to perform a range query on: ");

        // Use fgets to allow tree names with spaces
        fgets(treeName, 256, stdin);

        // Remove the trailing newline character that fgets captures
        treeName[strcspn(treeName, "\n")] = '\0';

        void clearInputBuffer();

        // Search for the tree in the tree hash table
        selectedTree = searchTreeHashTable(treeHashTable, treeName);

        if (selectedTree == NULL)
        {
            printf("Tree not found! Please enter a valid tree name.\n");
        }
    } while (selectedTree == NULL);

    return selectedTree;
}

void promptForRangeQueryBounds(SeqHashTable *hashTable, int *x, int *y, BST *tree)
{
    // Print the sequence hash table so the user can see the label start positions
    printf("\nHere is the current SeqHashTable to help with selecting bounds:\n");
    printSeqHashTable(hashTable, tree);

    // Prompt the user for the lower bound
    printf("\nPlease enter the lower bound for the range query: ");
    scanf("%d", x);
    void clearInputBuffer();

    // Prompt the user for the upper bound
    printf("Please enter the upper bound for the range query: ");
    scanf("%d", y);
    void clearInputBuffer();

    // Confirm the selected bounds
    printf("\nYou selected the range: %d to %d\n\n", *x, *y);
}

int getLeftMostStartPosition(Node *root)
{
    if (root == NULL)
    {
        printf("The tree is empty.\n");
        return -1; // Return an error value if the tree is empty
    }

    Node *current = root;

    // Traverse to the left-most node
    while (current->left != NULL)
    {
        current = current->left;
    }

    // Return the start position of the left-most node
    return current->start_coord;
}

int getRightMostStartPosition(Node *root)
{
    if (root == NULL)
    {
        printf("The tree is empty.\n");
        return -1; // Return an error value if the tree is empty
    }

    Node *current = root;

    // Traverse to the right-most node
    while (current->right != NULL)
    {
        current = current->right;
    }

    // Return the start position of the right-most node
    return current->start_coord;
}

void rangeQuery(BST *tree, SeqHashTable *table, BST *original)
{
    int x, y;
    int firstStart = getLeftMostStartPosition(original->root);
    int lastStart = getRightMostStartPosition(original->root);

    if (firstStart == -1 || lastStart == -1){
        printf("This dataset is empty.\n");
        return;
    }

    printf("\nThis dataset's segments have start values ranging from position (%d) to (%d).\n\n", firstStart, lastStart);

    printf("Please review the start positions of all segments to select the positional bounds for the range query:\n\n");
    printSeqHashTable(table, original);
    printf("\n");

    // Prompt the user to enter the lower and upper bounds
    printf("Enter the lower bound position: ");
    scanf("%d", &x);
    void clearInputBuffer();
    printf("Enter the upper bound position: ");
    scanf("%d", &y);
    void clearInputBuffer();

    if (x > y){
        printf("Upper bound must be greater than lower bound.\n");
        return;
    }
    if (x < firstStart || y > lastStart){
        printf("Bound values must be within %d and %d.\n", firstStart, lastStart);
        return;
    }

    if (tree == NULL || tree->root == NULL)
    {
        printf("No segments found in the tree.\n");
        return;
    }

    // Print the inputted range
    printf("\nInputted Range: %d to %d\n", x, y);

    // int range = 0;  // Update with number of segments within bounds
    // printf("There are %d segments in this range.\n\n", range);

    // Collect and print segments in the range [x, y]
    collectAndPrintSegmentsInRange(tree->root, x, y);
    printf("\n");
}

//////////////////////////////////
// AVL Tree Balancing Functions //
//////////////////////////////////

int height(Node *node)
{
    if (node == NULL)
    {
        return 0;
    }
    else
    {
        return node->height; 
    }
}

int getBalance(Node *node)
{
    if (node == NULL)
    {
        return 0;
    }
    else
    {
        return height(node->left) - height(node->right);    // Balance denoted by difference in height
    }
}

Node *rightRotate(Node *node)
{
    Node *left = node->left;                                // Left subtree
    Node *right_sub = left->right;                          // Right sub-subtree

    // Perform rotation
    left->right = node;
    node->left = right_sub;

    // Update heights after rotation
    node->height = 1 + (height(node->left) > height(node->right) ? height(node->left) : height(node->right));
    left->height = 1 + (height(left->left) > height(left->right) ? height(left->left) : height(left->right));

    // Update children counts
    node->children = 1 + (node->left ? node->left->children : 0) + (node->right ? node->right->children : 0);
    left->children = 1 + (left->left ? left->left->children : 0) + (left->right ? left->right->children : 0);

    return left;                                            // Return the new head
}

Node *leftRotate(Node *node)
{
    Node *right = node->right;
    Node *left_sub = right->left;

    // Perform rotation
    right->left = node;
    node->right = left_sub;

    // Update heights after rotation
    node->height = 1 + (height(node->left) > height(node->right) ? height(node->left) : height(node->right));
    right->height = 1 + (height(right->left) > height(right->right) ? height(right->left) : height(right->right));

    // Update children counts
    node->children = 1 + (node->left ? node->left->children : 0) + (node->right ? node->right->children : 0);
    right->children = 1 + (right->left ? right->left->children : 0) + (right->right ? right->right->children : 0);

    return right;  // Return the new head
}

Node *insert(Node *node, Node *new_node, BST *tree)
{
    if (node == NULL)
    {
        // Check if the new node has the same start and end coordinates as the previous node
        if (new_node->start_coord == prev_start && new_node->end_coord == prev_end)
        {
            printf("Skipping duplicate segment with start %d and end %d.\n", new_node->start_coord, new_node->end_coord);
            return NULL; // Skip the insertion if it's a duplicate
        }

        // Update the prev_start and prev_end with the current node's coordinates
        prev_start = new_node->start_coord;
        prev_end = new_node->end_coord;

        if (strcmp(new_node->type, "Unlabeled") == 0)
        {
            tree->unlabeled_count++;
            tree->unlabeled_starts = realloc(tree->unlabeled_starts, tree->unlabeled_count * sizeof(int));
            tree->unlabeled_starts[tree->unlabeled_count - 1] = new_node->start_coord;
        }
        else
        {
            tree->labeled_count++;
            tree->labeled_starts = realloc(tree->labeled_starts, tree->labeled_count * sizeof(int));
            tree->labeled_starts[tree->labeled_count - 1] = new_node->start_coord;
        }
        return new_node;
    }

    // Check for nesting: if new_node's range is completely within the current node's range
    if (new_node->start_coord > node->start_coord && new_node->end_coord < node->end_coord)
    {
        new_node->isNested = 1; // Mark as nested
    }

    if (new_node->start_coord == node->start_coord)
    {
        return node; // Avoid inserting duplicates
    }
    else if (new_node->start_coord > node->start_coord)
    {
        node->right = insert(node->right, new_node, tree);
    }
    else
    {
        node->left = insert(node->left, new_node, tree);
    }

    // Update height
    node->height = 1 + (height(node->left) > height(node->right) ? height(node->left) : height(node->right));

    // Update children count
    node->children = 1;
    if (node->left != NULL)
    {
        node->children += node->left->children;
    }
    if (node->right != NULL)
    {
        node->children += node->right->children;
    }

    // Balance the tree (AVL operations)
    int balance = getBalance(node);

    if (balance > 1 && new_node->start_coord < node->left->start_coord)
    {
        return rightRotate(node);
    }
    if (balance < -1 && new_node->start_coord > node->right->start_coord)
    {
        return leftRotate(node);
    }
    if (balance > 1 && new_node->start_coord > node->left->start_coord)
    {
        node->left = leftRotate(node->left);
        return rightRotate(node);
    }
    if (balance < -1 && new_node->start_coord < node->right->start_coord)
    {
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }

    return node;
}

/////////////////////////////
// Various Print Functions //
/////////////////////////////

void printStartPositions(BST *tree)
{
    printf("Labeled Segment Start Positions:\n");
    for (int i = 0; i < tree->labeled_count; i++)
    {
        printf("%d\n", tree->labeled_starts[i]);
    }

    printf("Unlabeled Segment Start Positions:\n");
    for (int i = 0; i < tree->unlabeled_count; i++)
    {
        printf("%d\n", tree->unlabeled_starts[i]);
    }
}

void printNodeDetailsInOrder(Node *root)
{
    if (root == NULL)
    {
        return;
    }

    // Traverse the left subtree
    if (root->left != NULL)
    {
        printNodeDetailsInOrder(root->left);
    }

    // Print current node's details (name, start, end, length, type)
    printf("Name: %s, Start: %d, End: %d, Length: %d, Type: %s\n",
            root->sequence ? root->sequence : "Unnamed",
            root->start_coord, root->end_coord, root->length,
            root->type ? root->type : "Unknown");

    // Traverse the right subtree
    if (root->right != NULL)
    {
        printNodeDetailsInOrder(root->right);
    }
}

int sumNodeLengths(Node *root)
{
    if (root == NULL)
    {
        return 0;
    }

    int length = 0;
    if (!root->isNested)
    {
        length = root->length; // Only add length if it's not nested
    }

    return length + sumNodeLengths(root->left) + sumNodeLengths(root->right);
}

void printAVLTreeSegments(Node *root)
{
    if (root == NULL)
    {
        return;
    }

    // Traverse the left subtree
    printAVLTreeSegments(root->left);

    // Print current node's details (start, end, length)
    printf("Name: %s, Start: %d, End: %d, Length: %d, Type: %s\n",
           root->gene_name ? root->gene_name : "Unnamed",
           root->start_coord, root->end_coord, root->length, root->type);

    // Traverse the right subtree
    printAVLTreeSegments(root->right);
}

// UI Functions

void clearInputBuffer() // Cleans up input parsing
{
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF); // Only clear remaining input
}

void displayIntroduction()
{
    printf("\n");
    printf("========================================================\n");
    printf("        Welcome to the Genomic Sequence Analyzer        \n");
    printf("========================================================\n\n");

    printf("This tool organizes your sequence data and retrieves all segments within a chosen range.\n");
    printf("You can also choose to exclude specific segments, allowing for more selective analysis.\n");
    printf("This program is designed to help researchers highlight potential areas of interest within entire genomes.\n\n");

    printf("------------------------------------------------------------------------------------------------------------\n");
    printf("Please select a command from the following:\n\n");
    printf(" q - Perform a range query.\n");
    printf("     This shows all segments within a range. You will be prompted to enter the lower and upper bound values.\n\n");
    printf(" d - View your data.\n");
    printf("     This displays your sequence data split into labeled segments and grouped by feature type.\n\n");
    printf(" f - Filter your data.\n");
    printf("     You will be prompted to enter which feature type you want to exclude(use command 'd' for start positions)\n");
    printf("     Create a filtered dataset to perform a range query on.\n\n");
    printf(" c - Display command list.\n\n");
    printf(" e - Exit the program.\n\n");
    printf("------------------------------------------------------------------------------------------------------------\n");
}

void promptForCommand(char *command)
{
    printf("Enter your command: ");
    scanf(" %c", command); // Notice the space before %c to skip leading whitespaces
    clearInputBuffer();     // Optional: if needed to ensure buffer is cleared
}

// Main Program Command Loop

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s <input_genbank_file>\n", argv[0]);
        return 1;
    }

    const char *csv_filename = "parsed_features.csv";
    runPythonParser(argv[1], csv_filename);

    BST *tree = newBST();
    SeqHashTable *seqHashTable = createSeqHashTable(256);
    TreeHashTable *treeHashTable = createTreeHashTable(10);

    parseCSVAndInsertIntoAVL(tree, seqHashTable, csv_filename);
    treeHashTableInsert(treeHashTable, "original", tree);

    char command[100];
    displayIntroduction();

    while (1)
    {
        promptForCommand(command);
        printf("You selected: %s\n\n", command);

        if (strcmp(command, "q") == 0)
        {
            // Perform range query only after prompting for a tree
            BST *selectedTree = promptUserForTreeSelection(treeHashTable);
            rangeQuery(selectedTree, seqHashTable, tree);
        }
        else if (strcmp(command, "d") == 0)
        {
            int total = tree->labeled_count + tree->unlabeled_count;
            printf("The length of your input sequence is %d bases.\n", tree->sequenceLength);
            printf("Your sequence contains %d unique* segments: %d labeled and %d unlabeled.\n\n", total, tree->labeled_count, tree->unlabeled_count);
            printSeqHashTable(seqHashTable, tree);
        }
        else if (strcmp(command, "f") == 0)
        {
            while (true)
            {
                // Prompt user for start positions to exclude
                int excludeCount;
                int *startPositionsToExclude = promptForStartPositionsToExclude(&excludeCount);
                printf("Filtering out segments with start positions: ");

                // Filter the tree and create a new filtered tree
                BST *filteredTree = filterTreeByStartPositions(tree, startPositionsToExclude, excludeCount, seqHashTable, treeHashTable);

                if (!filteredTree || filteredTree->root->children == tree->root->children)
                {
                    // Filtering failed due to invalid start positions, re-prompt
                    printf("\nFiltering failed due to invalid start positions. Please try again.\n\n");
                    free(startPositionsToExclude);
                    continue; // Re-prompt for valid start positions
                }

                // Free the memory for start positions after use
                free(startPositionsToExclude);

                // Prompt the user for a name to identify the filtered tree
                char *treeName = promptForUniqueTreeName(treeHashTable);

                // Insert the filtered tree into the tree hash table
                treeHashTableInsert(treeHashTable, treeName, filteredTree);

                // Notify user that the filtered tree has been added
                printf("Filtered tree '%s' has been added to the TreeHashTable.\n\n", treeName);
                free(treeName);
                break; // Exit loop after successful filtering
            }
        }
        else if (strcmp(command, "c") == 0)
        {
            printf("------------------------------------------------------------------------------------------------------------\n");
            printf("Please select a command from the following:\n\n");
            printf(" q - Perform a range query.\n");
            printf("     This shows all segments within a range. You will be prompted to enter the lower and upper bound values.\n\n");
            printf(" d - View your data.\n");
            printf("     This displays your sequence data split into labeled segments and grouped by feature type.\n\n");
            printf(" f - Filter your data.\n");
            printf("     You will be prompted to enter which feature type you want to exclude(use command 'd' for start positions)\n");
            printf("     Create a filtered dataset to perform a range query on.\n\n");
            printf(" c - Display command list.\n\n");
            printf(" e - Exit the program.\n\n");
            printf("------------------------------------------------------------------------------------------------------------\n");
        }
        else if (strcmp(command, "e") == 0)
            {
                printf("Exiting the program.\n");

                // Clean up any other allocated resources
                freeBST(tree);
                freeSeqHashTable(seqHashTable);
                freeTreeHashTable(treeHashTable);
                return 0;
            }
        else
        {
            printf("Invalid command. Please try again.\n\n");
        }
    }

    return 0;
}
