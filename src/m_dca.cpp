#include <jian/dca/Dca.hpp>
#include <jian/dca/ss_pairs.hpp>
#include <jian/nuc2d/SSTree.hpp>
#include <jian/nuc2d/SSE.hpp>
#include "nsp.hpp"

BEGIN_JN
	namespace dca {

		TreePath<SSE> direct_path(const TreePath<SSE> &p1, const TreePath<SSE> &p2) {
			auto pos = std::adjacent_find(p1.begin(), p1.end(), [&p2](const SSE *a1, const SSE *a2) {
				return std::adjacent_find(p2.begin(), p2.end(), [&a1, &a2](const SSE *b1, const SSE *b2) {
					return a1 == b1 && a2 != b2;
				}) != p2.end() || size(p2) == 1;
			});
			if (pos == p1.end() && size(p1) == 1) pos = p1.begin();
			TreePath<SSE> path;
			for (auto it = std::next(p1.rbegin()); it != p1.rend(); it++) {
				if (*it == *pos) {
					path.push_back(*it);
					break;
				}
				else {
					path.push_back(*it);
				}
			}
			bool flag = false;
			for (auto it = p2.begin(); std::next(it) != p2.end(); it++) {
				if (*it == *pos) {
					flag = true;
					continue;
				}
				if (flag) path.push_back(*it);
			}
			return path;
		}

		void print_path(const TreePath<SSE> &path) {
			for (auto && l : path) {
				LOG << l << ' ';
			}
			LOG << std::endl;
		}

		bool connect_by_ml(int a, int b, const SSTree &sst) {
			TreePath<SSE> path_a, path_b;
			for (auto &&sse : sst) {
				if (sse.has(a + 1)) path_a.push_back(&sse);
				if (sse.has(b + 1)) path_b.push_back(&sse);
			}
			auto path = direct_path(path_a, path_b);
			LOG << a << ' ' << b << std::endl;
			print_path(path_a);
			print_path(path_b);
			print_path(path);
			return std::any_of(path.begin(), path.end(), [](const SSE *l) {return l->is_ml(); });
		}

		void trim_di(const Par &par) {
			Str ss = par.get("ss");
            int l = size(ss);
			Str seq;
			seq.resize(l, 'A');
			Str difile = par.get("di");
			tuples_t tuples = tuples_from_file(difile, ss.size());
			tuples_t ls;
			SSTree sst(seq, ss, 1);
			//sst.head()->print_tree();
			//print_tuples(tuples);
            V<SSE *> v;
            v.resize(l);
			for (Int i = 0; i < l; i++) {
				v[i] = STD_ find_if(sst.begin(), sst.end(), [&i](auto &&sse) {return sse.has(i+1); }).val;
			}
            for (auto && l : v) LOG << l << ' '; LOG << std::endl;
			for (auto && tuple : tuples) {
				//if (connect_by_ml(tuple.a, tuple.b, sst)) {
				//if (v[tuple.a] != v[tuple.b] && ss[tuple.a] == '.' && ss[tuple.b] == '.') {
				if (v[tuple.a] != v[tuple.b]) {
					ls.push_back(tuple);
				}
			}
			print_tuples(ls);
		}

		REGISTER_NSP_COMPONENT(dca) {
			int n = 1;
			float step = 1;
			float pw = 0.5;
			S mol_type = "RNA";
			S method = "mf";
			Par::pars_t global = par.getv("global");
			if (size(global) == 2 && global[1] == "trim_di") {
				trim_di(par);
				return;
			}

			S out_file = par.get("o", "out");
			S fa_file = par.get("i", "in");

			par.set(method, "m", "method");
			par.set(n, "n");
			par.set(step, "step");
			par.set(mol_type, "t", "type");
			par.set(pw, "w");

			Dca *dca = FacDca::create(method, mol_type, pw);
			if (method == "mp") dca->set_step(step);
			dca->run(fa_file, out_file, n - 1);
			delete dca;
		}

	}
END_JN

