#include <iostream>
#include <string>
#include <vector>

#include "CifFile.h"
#include "CifParserBase.h"
#include "ISTable.h"

int main(int argc, char **argv)
{
    // The name of the CIF file
    string cifFileName = argv[1];
    
    // A string to hold any parsing diagnostics
    string diagnostics;

    // Create CIF file and parser objects
    CifFile *cifFileP = new CifFile;
    CifParser *cifParserP = new CifParser(cifFileP);

    // Parse the CIF file
    cifParserP->Parse(cifFileName, diagnostics);

    // Delete the CIF parser, as it is no longer needed
    delete cifParserP;

    // Display any diagnostics
    if (!diagnostics.empty())
    {
	    std::cout << "Diagnostics: " << std::endl << diagnostics << std::endl;
    }

    // Get the first data block name in the CIF file 
    string firstBlockName = cifFileP->GetFirstBlockName();

    // Retrieve the first data block 
    Block &block = cifFileP->GetBlock(firstBlockName);

    // Retrieve the table corresponding to the atom_site category, which delineates atomic constituents
    ISTable& atom_site = block.GetTable("atom_site");

    // Will hold the atom_site row indices of any atoms fulfilling our search query
    vector<unsigned int> results;

    // Holds attribute names and their target values
    vector<string> colNames, targets;

    // We want alpha carbons in chain A
    colNames.push_back("label_atom_id"); 
    targets.push_back("CA");
    colNames.push_back("auth_asym_id");
    targets.push_back("A");

    // Perform the search, propagating the results vector with atom indices
    atom_site.Search(results, targets, colNames);

    // Retrieve and display the coordinates of every atom satisfying our query
    std::cout << results.size() << " atoms found: \n";
    vector<string> coords;
    for (unsigned int i = 0; i < results.size(); ++i)
    {
        atom_site.GetRow(coords, results[i], "Cartn_x", "Cartn_z");
        std::cout << "(" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
    }
    return 0;
}
