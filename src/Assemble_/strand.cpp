#include "strand.h"

strand::strand() {
	head = tail = NULL;
}

int strand::getLen() {
	res *s;
	int len;

	for (len = 0, s = head; s != NULL; s = s->next) {
		len++;
	}
	return len;
}

double strand::getDist(RNA *mol) {
	int a, b;
	res *r;
	int i, j, k, flag;
	double x, y, z;

	a = head->num;
	for (r = head; r->next != NULL; r = r->next);
	b = r->num;

	x = 0; y = 0; z = 0;
	for (flag = 0, i = 0; i < (int) mol->chains.size(); i++) {
		for (j = 0; j < (int) mol->chains[i].residues.size(); j++) {
			flag++;
			if (flag == a) {
				for (k = 0; k < (int) mol->chains[i].residues[j].atoms.size(); k++) {
					if (mol->chains[i].residues[j].atoms[k].name == "O5*") {
						x = mol->chains[i].residues[j].atoms[k].x;
						y = mol->chains[i].residues[j].atoms[k].y;
						z = mol->chains[i].residues[j].atoms[k].z;
					}
				}
			}
			if (flag == b) {
				for (k = 0; k < (int) mol->chains[i].residues[j].atoms.size(); k++) {
					if (mol->chains[i].residues[j].atoms[k].name == "O3*") {
						x -= mol->chains[i].residues[j].atoms[k].x;
						y -= mol->chains[i].residues[j].atoms[k].y;
						z -= mol->chains[i].residues[j].atoms[k].z;
					}
				}
			}
		}
	}
	return sqrt(x * x + y * y + z * z);
}

strandInfoEx::strandInfoEx() {
	si = NULL;
	next = NULL;
	score = 0;
}

strandInfoEx::strandInfoEx(strandInfo *si, double score) {
	this->si = si;
	next = NULL;
	this->score = score;
}

strandInfoList::strandInfoList() {
	head = NULL;
}

int strandInfoList::getLen() {
	strandInfoEx *sie;
	int temp = 0;
	for (sie = head; sie != NULL; sie = sie->next) {
		temp++;
	}
	return temp;
}

void strandInfoList::add(strandInfo *si, double score) {
	strandInfoEx *sie = new strandInfoEx(si, score);
	if (getLen() == 0) {
		head = sie;
		return;
	}
	if (score > head->score) {
		sie->next = head;
		head = sie;
		return;
	}
	strandInfoEx *s;
	for (s = head; s->next != NULL; s = s->next) {
		if (score < s->next->score) {
			sie->next = s->next;
			s->next = sie;
			return;
		}
	}
	if (s->next == NULL) {
		s->next = sie;
	}
}

