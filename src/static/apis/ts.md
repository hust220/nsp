## 3D structure

<h4 id='atom'>Object: Atom</h4>

*   Attibutes:

    *   name

            print(atom.name)

*   Methods:

    *   Atom.new

            atom = Atom.new("atom.pdb")

    *   operator: []

            print(atom[1])
            print(atom[2])
            print(atom[3])
            atom[1] = 1.25
            atom[2] = 5.74
            atom[3] = 4.31

    *   __tostring

            print(atom)

    *   Atom.dist

            print(Atom.dist(atom1, atom2))

    *   Atom.dist2

            print(Atom.dist2(atom1, atom2))

    *   Atom.ang

            print(Atom.ang(atom1, atom2, atom3))

    *   Atom.dih

            print(Atom.dih(atom1, atom2, atom3, atom4))

<h4 id='res'>Object: Residue</h4>

*   Attibutes:

    *   name

            print(res.name)

    *   len

            print(res.len)
            for i=1, res.len do
                print(res[i])
            end

*   Methods:

    *   Residue.new

            res = Residue.new("residue.pdb")

    *   operator: []

            print(res[1])
            res[1] = atom1

    *   __tostring

            print(res)

    *   Residue.iter

            for atom in Residue.iter(res):
                print(atom)

    *   Residue.push

            res.push(atom)

<h4 id='chain'>Object: Chain</h4>

*   Attibutes:

    *   name

            print(chain.name)

    *   len

            print(chain.len)

*   Methods:

    *   Chain.new

            chain = Chain.new("chain.pdb")

    *   operator: []

            print(chain[1])
            chain[1] = res

    *   __tostring

            print(chain)

<h4 id='model'>Object: Model</h4>

*   Attibutes:

    *   name

            print(model.name)

    *   len

            print(model.len)

*   Methods:

    *   Model.new

            model = Model.new("model.pdb")

    *   operator: []

            print(model[1])
            model[1] = chain

    *   __tostring

            print(model)

<h4 id='mol'>Object: Molecule</h4>

*   Attibutes:

    *   name

            print(mol.name)

    *   len

            print(mol.len)

*   Methods:

    *   Molecule.new

            mol = Molecule.new("mol.pdb")

    *   operator: []

            print(mol[1])
            mol[1] = model

    *   __tostring

            print(mol)

