#include "Constraints.hpp"
#include "jian/utils/file.hpp"
#include "jian/utils/ls.hpp"
#include "../dca.hpp"

BEGIN_JN

	template<typename T>
	static bool equal(T &&ls, int i, int j) {
		return ls.size() == 2 && (ls[0] == i && ls[1] == j || ls[0] == j && ls[1] == i);
	}

	template<typename T>
	static bool equal(T &&ls, int i, int j, int k) {
		return ls.size() == 3 && ls[1] == j && (ls[0] == i && ls[2] == k || ls[0] == k && ls[2] == i);
	}

	template<typename T>
	static bool equal(T &&ls, int i, int j, int k, int l) {
		return ls.size() == 4 && (ls[0] == i && ls[3] == l || ls[0] == l && ls[3] == i) &&
			(ls[1] == j && ls[2] == k || ls[1] == k && ls[2] == j);
	}

	Constraint make_contact(int i, int j, double weight) {
		Constraint c; append(c.key, i, j); c.weight = weight; return c;
	}

	Constraint make_distance(int i, int j, double value, double weight) {
		Constraint c; append(c.key, i, j); c.value = value; c.weight = weight; return c;
	}

	Constraint make_angle(int i, int j, int k, double value, double weight) {
		Constraint c; append(c.key, i, j, k); c.value = value; c.weight = weight; return c;
	}

	Constraint make_dihedral(int i, int j, int k, int l, double value, double weight) {
		Constraint c; append(c.key, i, j, k, l); c.value = value; c.weight = weight; return c;
	}

	void Constraints::add_contact(int i, int j, double weight) {
		contacts.push_back(make_contact(i, j, weight));
	}

	void Constraints::add_distance(int i, int j, double value, double weight) {
		distances.push_back(make_distance(i, j, value, weight));
	}

	void Constraints::add_angle(int i, int j, int k, double value, double weight) {
		angles.push_back(make_angle(i, j, k, value, weight));
	}

	void Constraints::add_dihedral(int i, int j, int k, int l, double value, double weight) {
		dihedrals.push_back(make_dihedral(i, j, k, l, value, weight));
	}

	bool Constraints::has_contact(int i, int j) {
		for (auto &&c : contacts) if (equal(c.key, i, j)) return true; return false;
	}

	bool Constraints::has_distance(int i, int j) {
		for (auto &&c : distances) if (equal(c.key, i, j)) return true; return false;
	}

	bool Constraints::has_angle(int i, int j, int k) {
		for (auto &&c : angles) if (equal(c.key, i, j, k)) return true; return false;
	}

	bool Constraints::has_dihedral(int i, int j, int k, int l) {
		for (auto &&c : dihedrals) if (equal(c.key, i, j, k, l)) return true; return false;
	}

	void Constraints::read_dca(const S &f, int size) {
		if (!f.empty()) {
			dca::tuples_t && tuples = dca::tuples_from_file(f, size);
			for (auto && tuple : tuples) {
				add_contact(tuple.a, tuple.b, tuple.c);
			}
		}
	}

	void Constraints::read_contacts(const S &f) {
		if (!f.empty()) {
			for (auto &&it : FileLines(f)) {
				if (size(it.arr) == 3) add_contact(JN_INT(it.arr[0]) - 1, JN_INT(it.arr[1]) - 1, JN_DBL(it.arr[2]));
				else if (size(it.arr) == 2) add_contact(JN_INT(it.arr[0]) - 1, JN_INT(it.arr[1]) - 1, 1.0);
			}
		}
	}

	void Constraints::read_distances(const S &f) {
		if (!f.empty()) {
			for (auto &&it : FileLines(f)) {
				if (size(it.arr) == 3) add_distance(JN_INT(it.arr[0]) - 1, JN_INT(it.arr[1]) - 1, JN_DBL(it.arr[2]));
				else if (size(it.arr) == 2) add_distance(JN_INT(it.arr[0]) - 1, JN_INT(it.arr[1]) - 1, 10.0);
			}
		}
	}

END_JN

