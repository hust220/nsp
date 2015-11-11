#include "loop.h"

int loop::getLen() {
	int temp;
	res *tempRes;

	for (temp = 0, tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
		if (tempRes->type == '(' || tempRes->type == ')' || tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
			temp++;
		}
	}
	return temp;
}

int loop::size() {
	int temp;
	res *tempRes;

	for (temp = 0, tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
		temp++;
	}
	return temp;
}

int loop::getType() {
	int temp = 0;
	for (res *tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
		if (tempRes->type == ')') {
			temp++;
		}
	}
	return temp;
}

int loop::isVirtual() {
	res *tempRes;
	int flag = 0;

	if (head != NULL) {
		for (tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
			if (tempRes->type == '&') {
				flag = 1;
				break;
			}
		}
		return flag;
	} else {
		return -1;
	}
}

int loop::getLoopCounts() {
	int counts = 0;
	if (son != NULL) {
		counts += son->getLoopCounts();
	}
	if (brother != NULL) {
		counts += brother->getLoopCounts();
	}
	if (head != NULL) {
		counts++;
	}
	return counts;
}

string loop::getFlag() {
	string flag;
	res *r;
	int n = 0;
	for (r = head; r->next != NULL; r = r->next);
	if (r->type == '&') {
		int temp = 0;
		stringstream ss;
		for (res *tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
			if (tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
				temp++;
			} else if (tempRes->type == ')') {
				if (n % 2 == 1) {
					ss << temp << '-';
					temp = 0;
				}
				n++;
			}
		}
		ss << temp;
		ss >> flag;
	} else {
		int temp = 0;
		stringstream ss;
		for (res *tempRes = head->next; tempRes->next != NULL; tempRes = tempRes->next) {
			if (tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
				temp++;
			} else if (tempRes->type == ')') {
				if (n % 2 == 1) {
					ss << temp << '-';
					temp = 0;
				}
				n++;
			}
		}
		ss << temp << '-';
		ss >> flag;
	}
	return flag;
}

string loop::getSS() {
	string ss;
	for (res *r = head; r != NULL; r = r->next) {
		ss += r->type;
	}

	return ss;
}

string loop::getSeq() {
	string seq;
	for (res *r = head; r != NULL; r = r->next) {
		if (r->type != '&') {
			seq += r->name;
		}
	}
	return seq;
}

void loop::getStrand(vector<strand> &v) {
	strand *s;
	res *r, *r1;
	int flag;
	
	v.clear();
	s = new strand;
	for (r = head; r != NULL; r = r->next) {
		if ((r->type == '(' || r->type == '&') && s->getLen() != 0) {
			for (flag = 0, r1 = s->head; r1 != NULL; r1 = r1->next) {
				if (r1->type == '[' || r1->type == ']') {
					flag = 1;
					break;
				}
			}
			if (flag == 0) {
				v.push_back(*s);
			}
			s = new strand;
		} else if (r->type != '(' && r->type != ')' && r->type != '&') {
			r1 = new res(r);
			if (s->getLen() == 0) {
				s->head = s->tail = r1;
			} else {
				s->tail->next = r1;
				s->tail = r1;
			}
		}
	}
	if (s->getLen() != 0) {
		for (flag = 0, r1 = s->head; r1 != NULL; r1 = r1->next) {
			if (r1->type == '[' || r1->type == ']') {
				flag = 1;
				break;
			}
		}
		if (flag == 0) {
			v.push_back(*s);
		}
	}
}

RNA *loop::create(string ss, string seq) {
	cout << "hello!" << endl;
	if (ss.size() != seq.size() ) {
		cout << "loop::create error!" << endl;
		exit(1);
	}
	int loopLen = ss.size();
	double *bondLens = new double[2 * loopLen];
	bondLens[0] = 17;
	bondLens[1] = 3.9;
	bondLens[2 * loopLen - 1] = 3.9;
	for (int i = 2; i < 2 * loopLen - 1; i++) {
		if (i % 2 == 1) {
			bondLens[i] = 3.9;
		} else {
			if (ss[i / 2] == ')') {
				bondLens[i] = 17;
			} else {
				bondLens[i] = 3.9;
			}
		}
	}
	for (int i = 0; i < 2 * loopLen; i++) {
		cout << bondLens[i] << ' ';
	}
	cout << endl;

	double *bondLenSums = new double[2 * loopLen];
	double sum = 0;
	for (int i = 2 * loopLen - 1; i >= 0; i--) {
		sum += bondLens[i];
		bondLenSums[i] = sum;
	}
	for (int i = 0; i < 2 * loopLen; i++) {
		cout << bondLenSums[i] << ' ';
	}
	cout << endl;

	double m1 = 0.001 * ((10 + 10 * sstmm()) % 1000);
	if (m1 < 0.001) m1 = 0.001;
	if (m1 > 0.999) m1 = 0.999;
	genrand(&m1);

	Point **pos = new Point *[2 * loopLen];
	pos[0] = new Point(0, 0, 0);
	pos[2 * loopLen - 1] = new Point(0, 17, 0);
	pos[1] = Point::grow(pos[2 * loopLen - 1], pos[0], 100, 360 * ran_uniform(), 3.9);
/*	while (pos[1]->dist(pos[2 * loopLen - 1]) > bondLenSums[2]) {
		delete pos[1];
		pos[1] = Point::grow(pos[2 * loopLen - 1], pos[0], 100, 360 * ran_uniform(), 3.9);
	}
*/

	for (int i = 2; i < 2 * loopLen - 1; i++) {
		pos[i] = Point::grow(pos[i - 2], pos[i - 1], 100, 360 * ran_uniform(), bondLens[i]);
/*		if (i == 2 * loopLen - 2) {
			break;
		}
		while (pos[i]->dist(pos[2 * loopLen - 1]) > bondLenSums[i + 1]) {
			delete pos[i];
			pos[i] = Point::grow(pos[i - 2], pos[i - 1], 100, 360 * ran_uniform(), bondLens[i]);
		}
*/
	}

	cout << pos[0]->x <<' ' << pos[0]->y << ' ' << pos[0]->z << ' ';
	cout << "distance:" << pos[0]->dist(pos[2 * loopLen - 1]) << endl;
	for (int i = 1; i < 2 * loopLen; i++) {
		cout << pos[i]->x << ' ' << pos[i]->y << ' ' << pos[i]->z << ' ';
		cout << "distance:" << pos[i - 1]->dist(pos[i]);
		cout << endl;
	}

	for(int n = 0; n < 10000; n++) {
		int i = int(ran_uniform() * (2 * loopLen - 2)) + 1;
		Point *temp = Point::rotate(pos[i - 1], pos[i], pos[i + 1], 360 * ran_uniform(), bondLens[i], bondLens[i + 1]);
		delete pos[i];
		pos[i] = temp;
	}

	cout << pos[0]->x <<' ' << pos[0]->y << ' ' << pos[0]->z << ' ';
	cout << "distance:" << pos[0]->dist(pos[2 * loopLen - 1]) << endl;
	for (int i = 1; i < 2 * loopLen; i++) {
		cout << pos[i]->x << ' ' << pos[i]->y << ' ' << pos[i]->z << ' ';
		cout << "distance:" << pos[i - 1]->dist(pos[i]);
		cout << endl;
	}

	return 0;
}

loopInfo::loopInfo() {
	n = 0;
	len = 0;
	score = 0;
	next = NULL;
}

loopInfo::loopInfo(string loopName, string srcName, string fam, loop *l) {
	name = loopName;
	src = srcName;
	score = 0;
	len = l->getLen();
	flag = l->getFlag();
	n = l->getType();
	ss = l->getSS();
	seq = l->getSeq();
	family = fam;
}

loopInfoList::loopInfoList() {
	head = NULL;
	len = 0;
}

loopInfoList::~loopInfoList() {
	for (loopInfo *tempLi = head; tempLi != NULL;) {
		loopInfo *next = tempLi->next;
		delete tempLi;
		tempLi = next;
	}
}

int loopInfoList::getLen() {
	return len;
}

void loopInfoList::add(loopInfo *li) {
	len++;
	if (head == NULL) {
		head = li;
		head->next = NULL;
		return;
	}

	if (li->score > head->score) {
		li->next = head;
		head = li;
		return;
	} 
	loopInfo *q1, *q2;
	for (q1 = head, q2 = head->next; q2 != NULL; q2 = q2->next) {
		if (li->score > q2->score) {
			q1->next = li;
			li->next = q2;
			return;
		}
		q1 = q2;
	}
	q1->next = li;
	li->next = q2;
}

