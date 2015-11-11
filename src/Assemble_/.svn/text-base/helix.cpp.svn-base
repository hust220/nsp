#include "helix.h"

int helix::getLen() {
	bp *b;
	int temp;

	for (b = head, temp = 0; b != NULL; b = b->next) {
		temp++;
	}
	return temp;
}

string helix::seq() const {
	string str;
	int i = 0;
	for (auto b = head; b != NULL; b = b->next) {
		str.insert(i, string() + b->res1.name + b->res2.name);
		i++;
	}
	return str;
}

helixInfoEx::helixInfoEx() {
	hi = NULL;
	next = NULL;
	score = 0;
}

helixInfoEx::helixInfoEx(helixInfo *hi, double score) {
	this->hi = hi;
	next = NULL;
	this->score = score;
}

helixInfoList::helixInfoList() {
	head = NULL;
}

int helixInfoList::getLen() {
	helixInfoEx *hie;
	int temp = 0;
	for (hie = head; hie != NULL; hie = hie->next) {
		temp++;
	}
	return temp;
}

void helixInfoList::add(helixInfo *hi, double score) {
	helixInfoEx *hie = new helixInfoEx(hi, score);
	if (getLen() == 0) {
		head = hie;
		return;
	}
	if (score > head->score) {
		hie->next = head;
		head = hie;
		return;
	}
	helixInfoEx *h;
	for (h = head; h->next != NULL; h = h->next) {
		if (score > h->next->score) {
			hie->next = h->next;
			h->next = hie;
			return;
		}
	}
	if (h->next == NULL) {
		h->next = hie;
	}
}


