#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <pdb/util.h>

namespace jian {

namespace nuc3d {

class Transform {
public:
	Transform() {}
	Transform(const Model &model) {
		_model = model;
	}
    Transform(Model &&model) {
        std::swap(_model, model);
    }

	Model operator() (string, string);
	Model to_rna(string);
	Model to_dna(string);
	Model &model() {
		return _model;
	}
	const Model &model() const {
		return _model;
	}
	void model(const Model &model) {
		_model = model;
	}
    Convert cvt;

private:
	Model _model;
};


}


}





#endif 

