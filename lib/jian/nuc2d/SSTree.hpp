#pragma once

#include <string>

namespace jian {

	class loop;
	class SSTreeImpl;

	std::pair<int, int> loop_head_tail(loop *l);
	loop *ss_tree(std::string seq, std::string ss, int hinge = 2);
	void ss_read_tree(std::string &ss, loop *l);
	void seq_read_tree(std::string &seq, loop *l);
	void free_ss_tree(loop *l);
	void print_ss_tree(loop *l);

	class SSTree {
	public:
		SSTree();
		~SSTree();
		loop *&head();
		bool empty() const;
		// make tree with no broken tag
		void make(const std::string &seq, const std::string &ss, int hinge = 2);
		// make tree with broken tag
		void make_b(const std::string &seq, const std::string &ss, int hinge = 2);
	private:
		SSTreeImpl *_impl;
	};

} // namespace jian

