#pragma once

#include "pdb.hpp"

BEGIN_JN

	class CompleteResidue {
	private:
		CompleteResidue();
		CompleteResidue(const CompleteResidue &);
		CompleteResidue &operator =(const CompleteResidue &);

	public:
		std::map<std::string, Residue> m_bb;
		std::map<std::string, Residue> m_base;

		static const CompleteResidue &instance();
		bool lack_atoms(const Residue &residue) const;
		void operator ()(Residue &residue) const;
	};
}