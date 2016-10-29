#pragma once

namespace jian {

	template<typename Derived>
	class CG {
	public:
		static Chain chain(const Chain &chain) {
			Chain c;
			c.name = chain.name;
			c.model_name = chain.model_name;
			for (auto && r : chain) {
				Residue res = r;
				c.push_back(res.cg<Derived>());
			}
			return c;
		}
	};

}
