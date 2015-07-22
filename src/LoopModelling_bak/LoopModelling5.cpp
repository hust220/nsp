#include "LoopModelling5.h"

LoopModelling5::LoopModelling5(string seq, string ss_) {
	this->seq = seq;
	int temp = 0;
	for (int i = 0; i < ss_.size(); i++) {
		if (ss_[i] == '.' || ss_[i] == '(' || ss_[i] == ')' || ss_[i] == '[' || ss_[i] == ']') {
			ss += ss_[i];
			temp++;
		} else if (ss_[i] == '&' && i > 0 && temp < seq.size()) {
			brk.push_back(temp - 1);
		}
	}
	int len1 = seq.size();
	int len2 = ss.size();
	if (len1 != len2) {
		cerr << "Sequence's length is not equal to second structure's length!" << endl;
		exit(1);
	}
	resLen = ss.size();
	len = ss.size();
	type = new int[resLen];
	for (int i = 0; i < resLen; i++) {
		if (seq[i] == 'A') {
			type[i] = 0;
		} else if (seq[i] == 'U') {
			type[i] = 1;
		} else if (seq[i] == 'G') {
			type[i] = 2;
		} else {
			type[i] = 3;
		}
	}
	bound = NULL;

	init();
	dg = new DG(bound);
}

LoopModelling5::LoopModelling5(string filename) {
	string ss_;
	ifstream ifile(filename.c_str());
	while (ifile) {
		string line;
		getline(ifile, line);
		vector<string> splited_line;
		tokenize(line, splited_line, " ,:;");
		if (splited_line.size() < 2) continue;
		if (splited_line[0] == "sequence") {
			seq = splited_line[1];
		} else if (splited_line[0] == "secondary_structure") {
			ss_ = splited_line[1]; 
		} else if (splited_line[0] == "constraints") {
			int row = int((splited_line.size() - 1) / 3);
			int col = 3;
			constraints = new Matr_(row, col);
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < 3; j++) {
					constraints->data[i][j] = atof(splited_line[1 + i * 3 + j].c_str());
				}
			}
		}
	}
	ifile.close();

	int temp = 0;
	for (int i = 0; i < ss_.size(); i++) {
		if (ss_[i] == '.' || ss_[i] == '(' || ss_[i] == ')' || ss_[i] == '[' || ss_[i] == ']') {
			ss += ss_[i];
			temp++;
		} else if (ss_[i] == '&' && i > 0 && temp < seq.size()) {
			brk.push_back(temp - 1);
		}
	}
	int len1 = seq.size();
	int len2 = ss.size();
	if (len1 != len2) {
		cerr << "Sequence's length is not equal to second structure's length!" << endl;
		exit(1);
	}
	resLen = ss.size();
	len = ss.size();
	type = new int[resLen];
	for (int i = 0; i < resLen; i++) {
		if (seq[i] == 'A') {
			type[i] = 0;
		} else if (seq[i] == 'U') {
			type[i] = 1;
		} else if (seq[i] == 'G') {
			type[i] = 2;
		} else {
			type[i] = 3;
		}
	}
	bound = NULL;

	init();
	dg = new DG(bound);
}

void LoopModelling5::init() {
	delete bound;
	bound = new Matr_(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bound->data[i][j] = 0;
			} else {
				bound->data[j][i] = 7;
				bound->data[i][j] = 999;
			}
		}
	}
	for (int i = 0; i < resLen - 1; i++) {
		if (ss[i] == ss[i + 1] && ss[i] != '.') {
			bound->data[i][i + 1] = 6.2 + 0.001;
			bound->data[i + 1][i] = 6.2 - 0.001;
		} else {
			bound->data[i][i + 1] = 6.2 + 0.001;
			bound->data[i + 1][i] = 6.2 - 0.001;
		}
	}
	for (int i = 0; i < resLen - 2; i++) {
		if (ss[i] == ss[i + 1] && ss[i + 1] == ss[i + 2] && (ss[i] == '(' || ss[i] == ')')) {
			bound->data[i][i + 2] = 11.9 + 0.001;
			bound->data[i + 2][i] = 11.9 - 0.001;
		}
	}
	R2D r2d(seq, ss);
	setHelixPar(r2d._head);
	setHelixPar(r2d._pseudo_head);
/*
	for (int i = 2; i < resLen; i++) {
		if (ss[i] != '.') {
			if (ss[i - 2] == ss[i - 1] && ss[i - 1] == ss[i]) {
				int type1 = type[i - 2];
				int type2 = type[i];
				bound->data[i - 2][i] = nxtPrmt[type1 * 4 + type2]->data[n][1];
				bound->data[i][i - 2] = nxtPrmt[type1 * 4 + type2]->data[n][0];
			}
		}
	}
*/
	// break point
	for (int i = 0; i < (int) brk.size(); i++) {
		int a = brk[i];
		int b = a + 1;
		bound->data[a][b] = 999;
		bound->data[b][a] = 7;
	}

	// constraints
	if (constraints.get() != NULL) {
		for (int i = 0; i < constraints->row; i++) {
			int a = constraints->data[i][0] - 1;
			int b = constraints->data[i][1] - 1;
			int c = constraints->data[i][2];
			if (a > b) {
				int t = a;
				a = b; 
				b = t;
			}
			bound->data[a][b] = c + 0.001;
			bound->data[b][a] = c - 0.001;
		}
	}
}

void LoopModelling5::setHelixPar(Frag *head_module) {
	if (head_module == NULL) {
		return;
	}
	if (head_module->_type == 0) {
		int n = head_module->_len / 2;
		int num1[n];
		int num2[n];
		for (int i = 0; i < n; i++) {
			num1[i] = head_module->_num[i];
			num2[i] = head_module->_num[i + n];
		}

		string lib = getenv("RNA");
		string helix_par_file = lib + "helix.par";
		ifstream ifile(helix_par_file.c_str());
		int helix_par_dim;
		ifile >> helix_par_dim;
		Obj<Matr_> helix_par = new Matr_(helix_par_dim, helix_par_dim);
		for (int i = 0; i < helix_par_dim; i++) {
			for (int j = 0; j < helix_par_dim; j++) {
				ifile >> helix_par->data[i][j];
			}
		}
		ifile.close();
		// Obj<Matr_> temp_helix_par(2 * n, 2 * n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				/*
				cout << num1[i] << ' ' << num1[j] << ' ' << helix_par->data[i][j] << endl;
				cout << num1[i] << ' ' << num2[n - 1 - j] << endl;
				cout << num2[n - 1 - i] << ' ' << num1[j] << endl;
				cout << num2[n - 1 - i] << ' ' << num2[n - 1 - j] << endl;
				cout << "------" << endl;
				*/
				bound->data[num1[i]][num1[j]] /*= temp_helix_par->data[i][j]*/ = helix_par->data[i][j];
				bound->data[num1[i]][num2[j]] /*= temp_helix_par->data[i][j + n]*/ = helix_par->data[i][j + helix_par_dim - n];
				bound->data[num2[i]][num1[j]] /*= temp_helix_par->data[i + n][j]*/ = helix_par->data[i + helix_par_dim - n][j];
				bound->data[num2[i]][num2[j]] /*= temp_helix_par->data[i + n][j + n]*/ = helix_par->data[i + helix_par_dim - n][j + helix_par_dim - n];
			}
		}
		/*
		for (int i = 1; i < n - 1; i++) {
			bound->data[num1[i]][num2[i - 1]] = 17.5 + 0.001;
			bound->data[num2[i - 1]][num1[i]] = 17.5 - 0.001;
			bound->data[num1[i]][num2[i]] = 15.3 + 0.001;
			bound->data[num2[i]][num1[i]] = 15.3 - 0.001;
			bound->data[num1[i]][num2[i + 1]] = 13.2 + 0.001;
			bound->data[num2[i + 1]][num1[i]] = 13.2 - 0.001;
			bound->data[num1[i - 1]][num2[i]] = 17.5 + 0.001;
			bound->data[num2[i]][num1[i - 1]] = 17.5 - 0.001;
			bound->data[num1[i]][num2[i]] = 15.3 + 0.001;
			bound->data[num2[i]][num1[i]] = 15.3 - 0.001;
			bound->data[num1[i + 1]][num2[i]] = 13.2 + 0.001;
			bound->data[num2[i]][num1[i + 1]] = 13.2 - 0.001;
		}
		*/
	}
	Frag *son_module = head_module->_son;
	Frag *brother_module = head_module->_brother;
	setHelixPar(son_module);
	setHelixPar(brother_module);
}

Obj<RNA> LoopModelling5::run() {
	coord = dg->run();
	coord->print();
	mc();
	return c2a();
}

void LoopModelling5::to3pt() {
}

void LoopModelling5::mc() {
}

Obj<RNA> LoopModelling5::c2a() {
	return NULL;
}

void LoopModelling5::ss2ct() {
	int iTemp;
	vector<char> vcList;
	vector<int> viList;
	vector<char> vcList2;
	vector<int> viList2;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(' || ss[i] == '[') {
			iTemp++;
		}
	}
	if (iTemp == 0) return;
	delete ct;
	ct = new Matr_(iTemp, 2);
	int flag = 0;
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(') {
			vcList.push_back('(');
			viList.push_back(i);
		} else if (ss[i] == '[') {
			vcList2.push_back('[');
			viList2.push_back(i);
		} else if (ss[i] == ')') {
			if (vcList.size() == 0) {
				cerr << "wrong secondary structure:\n" << ss << endl;
				exit(1);
			}
			ct->data[flag][0] = viList.back();
			ct->data[flag][1] = i;
			flag++;
			vcList.pop_back();
			viList.pop_back();
		} else if (ss[i] == ']') {
			if (vcList2.size() == 0) {
				cerr << "wrong secondary structure:\n" << ss << endl;
				exit(1);
			}
			ct->data[flag][0] = viList2.back();
			ct->data[flag][1] = i;
			flag++;
			vcList2.pop_back();
			viList2.pop_back();
		}
	}
	if (vcList.size() != 0 || vcList2.size() != 0) {
		cerr << "wrong secondary structure:\n" << ss << endl;
		exit(1);
	}
}

double LoopModelling5::at(Matr_ *bound, int i, int j, int len) {
    if (j >= len) {
        j -= len;
        return at(bound, j, i, len);
    }
    if (i >= len) {
        i -= len;
        return at(bound, j, i, len);
    }
    return bound->data[i][j];
}

void LoopModelling5::assign(Matr_ *bound, int i, int j, double d, int len) {
    if (j >= len) {
        j -= len;
        assign(bound, j, i, d, len);
        return;
    }
    if (i >= len) {
        i -= len;
        assign(bound, j, i, d, len);
        return;
    }
    bound->data[i][j] = d;
}

