#pragma once

#include <list>
#include "../utils/string.hpp"
#include "loop.hpp"

BEGIN_JN

	class loop;

	std::pair<int, int> loop_head_tail(loop *l);
	loop *ss_tree(Str seq, Str ss, int hinge = 2);
	void ss_read_tree(Str &ss, loop *l);
	void seq_read_tree(Str &seq, loop *l);
	void free_ss_tree(loop *l);
	void print_ss_tree(loop *l);

	class SSTree {
	public:
		using Path = std::list<loop *>;

        loop *m_head;

		SSTree();
		~SSTree();
		loop *&head();
		const loop *head() const;
		bool empty() const;
		// make tree with no broken tag
		void make(const Str &seq, const Str &ss, int hinge = 2);
		// make tree with broken tag
		void make_b(const Str &seq, const Str &ss, int hinge = 2);

        template<typename _Fn>
        loop *find(_Fn &&fn) const {
            loop *p = NULL;
            traverse([&p, &fn](auto &&q, auto &&path){
                if (fn(q, path)) {
                    p = q;
                    return true;
                }
                return false;
            });
            return p;
        }

		template<typename _Fn>
		void traverse(_Fn &&fn) const {
			Path ls;
			loop * L = m_head;
			while (true) {
				if (L == NULL) break;
				ls.push_back(L);
				if (fn(L, ls)) return;
				if (L->son != NULL) {
					L = L->son;
				}
				else if (L->brother != NULL) {
					L = L->brother;
					ls.pop_back();
				}
				else {
					while (true) {
						ls.pop_back();
						if (!ls.empty()) {
							L = ls.back()->brother;
							if (L == NULL) continue;
							else { ls.pop_back(); break; }
						}
						else { L = NULL; break; }
					}
				}
			}
		}

		template<typename _Fn>
		Path path(_Fn &&fn) const {
			Path ls;
			loop * L = m_head;
			while (true) {
				if (L == NULL) break;
				ls.push_back(L);
				if (fn(L)) return ls;
				if (L->son != NULL) {
					L = L->son;
				}
				else if (L->brother != NULL) {
					L = L->brother;
					ls.pop_back();
				}
				else {
					while (true) {
						ls.pop_back();
						if (!ls.empty()) {
							L = ls.back()->brother;
							if (L == NULL) continue;
							else { ls.pop_back(); break; }
						}
						else { L = NULL; break; }
					}
				}
			}
			return ls;
		}

	};

END_JN

