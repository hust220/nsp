#include "nsp.hpp"
#include "pdb.hpp"

BEGIN_JN

REGISTER_NSP_COMPONENT(fit) {
   Molecule m;
   mol_read(m, par.getv("global")[1]);
   for (auto && model : m) {
      for (auto && chain : model) {
         Residue r = chain[0];
         chain[0].clear();
         for (auto && atom : r) {
            if (atom.name != "P" && atom.name != "O1P" && atom.name != "O2P") {
               chain[0].push_back(std::move(atom));
            }
         }
      }
   }
   JN_OUT << m << std::endl;
}

END_JN

