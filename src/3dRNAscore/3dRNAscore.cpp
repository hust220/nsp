#include "../Score"

int main(int argc, char **argv) {
	if (!strcmp(argv[1], "-score:list")) {
		string str;
		Obj<Score> score;
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(argv[3]);
		}

		ifstream ifile(argv[2]);
		while (ifile >> str) {
			Obj<RNA> rna = new RNA(str);
			cout << str << ' ' << score->run(rna) << endl;
		}
		ifile.close();
	} else if (!strcmp(argv[1], "-score")) {
		Obj<Score> score;
		Obj<RNA> rna = new RNA(argv[2]);
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(argv[3]);
		}
		cout << score->run(rna) << endl;
	} else if (!strcmp(argv[1], "-train:distance")) {
		Obj<Train> train;
		if (argc == 3) {
			train = new Train;
		} else if (argc == 4) {
			train = new Train(argv[3]);
		}
		train->dist(argv[2]);
	} else if (!strcmp(argv[1], "-train:dihedral")) {
		Train train;
		train.dih(argv[2]);
	}
	return 0;
}
