#pragma once

#include "../utils/ls.hpp"

BEGIN_JN

	class Constraint {
	public:
		std::vector<int> key;
		double value;
		double min, max;
		double weight;
	};

	Constraint make_contact(int i, int j, double weight = 1);
	Constraint make_distance(int i, int j, double value, double weight = 1);
	Constraint make_angle(int i, int j, int k, double value, double weight = 1);
	Constraint make_dihedral(int i, int j, int k, int l, double value, double weight = 1);

	class Constraints {
	public:
		std::deque<Constraint> contacts;
		std::deque<Constraint> distances;
		std::deque<Constraint> angles;
		std::deque<Constraint> dihedrals;

		void add_contact(int i, int j, double weight = 1);
		void add_distance(int i, int j, double value, double weight = 1);
		void add_angle(int i, int j, int k, double value, double weight = 1);
		void add_dihedral(int i, int j, int k, int l, double value, double weight = 1);
		bool has_contact(int i, int j);
		bool has_distance(int i, int j);
		bool has_angle(int i, int j, int k);
		bool has_dihedral(int i, int j, int k, int l);
		void read_dca(const S &f, int size);
		void read_contacts(const S &f);
		void read_distances(const S &f);
	};

END_JN

