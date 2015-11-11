#include "R2D.h"

using namespace jian;

void R2D::operator ()(std::string seq, std::string second_struct, int view) {
	delTree(_head);
	delTree(_pseudo_head);

	// check sequence
	int len_seq = 0;
	for (int i = 0; i < seq.size(); i++) {
		if (std::set<char>{'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}.count(seq[i])) {
		/*
		if (seq[i] == 'A' || seq[i] == 'a' || 
		seq[i] == 'U' || seq[i] == 'u' || 
		seq[i] == 'G' || seq[i] == 'g' || 
		seq[i] == 'C' || seq[i] == 'c') {
		*/
			len_seq++;
		} else {
			cerr << "R2D::operator () error! RNA sequence has no residue named '" << seq[i] << "'" << endl;
			exit(1);
		}
	}
	_seq = seq;

	// check secondary structure
	int len_second_struct = 0;
	_ss = "";
	std::string second_struct1 = ""; // secondary structure which has no '[' and ']';
	std::string second_struct2 = ""; // secondary structure which has no '(' and ')';
	for (int i = 0; i < second_struct.size(); i++) {
		if (second_struct[i] == '.') {
			len_second_struct++;
			_ss += '.';
			second_struct1 += '.';
			second_struct2 += '.';
		} else if (second_struct[i] == '(' || second_struct[i] == ')') {
			len_second_struct++;
			_ss += second_struct[i];
			second_struct1 += second_struct[i];
			second_struct2 += '.';
		} else if (second_struct[i] == '[' || second_struct[i] == ']') {
			len_second_struct++;
			_ss += '.';
			second_struct1 += '.';
			if (second_struct[i] == '[') { 
				second_struct2 += '('; // transform '[' to '(' 
			} else {
				second_struct2 += ')'; // transform ']' to ')'
			}
		} else if (second_struct[i] == '&') {
		} else {
			cerr << "R2D::operator () error! Please input secondary structure of legal dot-bracket form!" << endl;
			exit(1);
		}
	}
	if (len_seq != len_second_struct) {
		cerr << "The length of the sequence and the length of the secondary structure must be equal!" << endl;
		exit(1);
	}
	_len = len_seq;

	_view = view;

	// set 2d structure tree
	_head = setTree(_seq, second_struct1);
	_pseudo_head = setTree(_seq, second_struct2);
}

Frag* R2D::setTree(std::string seq, std::string second_struct) {
	int len = seq.size();
	std::vector<Nuc> left_stack;
	std::vector<Nuc> right_stack;
	std::vector<Frag *> module_stack;
	Frag *head;

	for (int i = 0; i < second_struct.size(); i++) {
		Nuc nuc(seq[i], second_struct[i], i);
		left_stack.push_back(nuc);
		int left_num = 0;
		int right_num = 0;
		if (i == 0) {
			if (nuc._ss == ')') {
				cerr << "R2D::setTree error! Wrong secondary structure!" << endl;
				exit(1);
			} else if (nuc._ss == '(') {
				left_num++;
			}
		} else {
			if (nuc._ss == ')') {
				right_num++;
				if ((i + 1 < len && second_struct[i + 1] != ')') || i + 1 == len) {
					int left_bracket_number = 0;
					int right_bracket_number = 0;
					while (1) {
						Nuc temp_nuc = left_stack.back();
						right_stack.push_back(temp_nuc);
						left_stack.pop_back();
						if (temp_nuc._ss == '(') {
							left_bracket_number++;
						} else if (temp_nuc._ss == ')') {
							right_bracket_number++;
						}
						if (temp_nuc._ss == '(' && ((left_stack.back()._ss != '(') || left_stack.empty())) {
							std::vector<Nuc> temp_stack;
							int temp_left_number = 0;
							int temp_right_number = 0;
							int p1 = 0, p2, p3, p4;
							int loop_type = 0;
							int loop_length = 0;
							int j = 0;
							while (1) {
								Nuc temp_nuc2 = right_stack.back();
								temp_stack.push_back(temp_nuc2);
								right_stack.pop_back();
								if (temp_nuc2._ss == '(') {
									temp_left_number++;
									if (right_stack.back()._ss != '(') {
										p2 = j;
									}
								} else if (temp_nuc2._ss == ')') {
									temp_right_number++;
									if (temp_right_number == 1) {
										p3 = j;
									} else if (temp_right_number == temp_left_number) {
										p4 = j;
									}
								} else if (temp_nuc2._ss == '[') {
									loop_type++;
									loop_length++;
								} else if (temp_nuc2._ss == ']') {
									loop_length++;
								} else if (temp_nuc2._ss == '.') {
									loop_length++;
								}
								if (temp_left_number == temp_right_number) {
									int temp_length;
									std::string temp_seq;
									std::string temp_ss;
									int *temp_num;

									// construct helix
									temp_length = temp_left_number + temp_right_number;
									temp_num = new int[temp_length];
									temp_seq = "";
									temp_ss = "";
									for (int k = 0; k < temp_left_number; k++) {
										temp_seq += temp_stack[k]._seq;
										temp_ss += '(';
										temp_num[k] = temp_stack[k]._num;
									}
									for (int k = temp_left_number; k < temp_length; k++) {
										temp_seq += temp_stack[k + p3 - p2 - 1]._seq;
										temp_ss += ')';
										temp_num[k] = temp_stack[k + p3 - p2 - 1]._num;
									}
									Frag *temp_helix = new Frag(0, temp_length, temp_seq, temp_ss, temp_num);

									if (loop_length != 0) { // loop + helix
										// construct loop
										temp_seq = "";
										temp_ss = "";
										temp_length = p3 - p2 + 1;
										temp_num = new int[temp_length];
										for (int k = 0; k < temp_length; k++) {
											temp_seq += temp_stack[k + p2]._seq;
											temp_ss += temp_stack[k + p2]._ss;
											temp_num[k] = temp_stack[k + p2]._num;
										}
										Frag *temp_loop = new Frag(loop_type + 1, temp_length, temp_seq, temp_ss, temp_num);

										// push loop into the module stack
										if (loop_type > 0) { // non-hairpin loop
											Frag *temp_module = module_stack.back();
											temp_loop->_son = temp_module;
											module_stack.pop_back();
											for (int k = 0; k < loop_type - 1; k++) {
												Frag *new_temp_module = module_stack.back();
												temp_module->_brother = new_temp_module;
												module_stack.pop_back();
												temp_module = new_temp_module;
											}
											module_stack.push_back(temp_loop);
										} else { // hairpin loop
											module_stack.push_back(temp_loop);
										}

										// push helix into the module stack
										Frag *temp_module = module_stack.back();
										temp_helix->_son = temp_module;
										module_stack.pop_back();
										module_stack.push_back(temp_helix);
									
										// set head
										head = temp_helix;

									} else { // only helix

										// push helix into the module stack
										module_stack.push_back(temp_helix);

										// set head
										head = temp_helix;
									}

									/*
									if (!module_stack.empty()) {
										Frag *temp_module = module_stack.back();
										temp_helix->_son = temp_module;
										module_stack.pop_back();
									}
									module_stack.push_back(temp_helix);
									head = temp_helix;
									*/

									Nuc temp_nuc1(temp_nuc2._seq, ']', temp_nuc2._num);
									right_stack.push_back(temp_nuc1);
									Nuc temp_nuc2(temp_nuc._seq, '[', temp_nuc._num);
									right_stack.push_back(temp_nuc2);
									break;
								}
								j++;
							}
						}
						if (left_bracket_number == right_bracket_number) {
							Nuc temp_nuc1(temp_nuc._seq, '[', temp_nuc._num);
							left_stack.push_back(temp_nuc1);
							Nuc temp_nuc2(nuc._seq, ']', nuc._num);
							left_stack.push_back(temp_nuc2);
							break;
						}
						if (left_stack.empty() && left_bracket_number != right_bracket_number) {
							cerr << "Wrong secondary structure!" << endl;
							exit(1);
						}
					}
				}
			}
		}
	}
	if (!left_stack.empty()) {
		// construct loop
		int loop_type = 0;
		int temp_length = left_stack.size();
		std::string temp_seq;
		std::string temp_ss;
		int *temp_num = new int[temp_length];
		for (int i = 0; i < temp_length; i++) {
			if (left_stack[i]._ss == '(' || left_stack[i]._ss == ')') {
				cerr << "Wrong secondary structure!" << endl;
				exit(1);
			} else if (left_stack[i]._ss == '[') {
				loop_type++;
				temp_ss += '(';
			} else if (left_stack[i]._ss == ']') {
				temp_ss += ')';
			} else if (left_stack[i]._ss == '.') {
				temp_ss += '.';
			}
			temp_seq += left_stack[i]._seq;
			temp_num[i] = left_stack[i]._num;
		}
		Frag *temp_loop = new Frag(-loop_type - 1, temp_length, temp_seq, temp_ss, temp_num);
		if (temp_ss == "()") {
			return head;
		}
		if (!module_stack.empty()) {
			Frag *temp_module = module_stack.back();
			temp_loop->_son = temp_module;
			module_stack.pop_back();
			for (int i = 0; i < loop_type - 1; i++) {
				Frag *new_temp_module = module_stack.back();
				temp_module->_brother = new_temp_module;
				module_stack.pop_back();
				temp_module = new_temp_module;
			}
		}
		module_stack.push_back(temp_loop);
		head = temp_loop;
	}
	return head;
}

void R2D::delTree(Frag *head) {
	if (head == NULL) {
		return;
	}
	delTree(head->_son);
	delTree(head->_brother);
	delete head;
}

void R2D::print() {
	cout << "================= Secondary structure tree =====================" << endl;
	printTree(_head);
	cout << endl;
	cout << "================= Pseudoknot tree ======================" << endl;
	printTree(_pseudo_head);
	cout << endl;
}

void R2D::printTree(Frag *head) {
	if (head == NULL) {
		return;
	}

	cout << head->_seq << endl;
	cout << head->_ss << endl;
	for (int i = 0; i < head->_len; i++) {
		cout << head->_num[i] << ' ';
	}
	cout << endl;
	cout << endl;

	printTree(head->_son);
	printTree(head->_brother);
}

