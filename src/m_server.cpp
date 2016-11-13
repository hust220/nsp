#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/socket.hpp>

namespace jian {

	REGISTER_NSP_COMPONENT(server) {
		int port = 80;
		Socket socket;

		par.set(port, "port");

		auto s = socket.open(SOCK_STREAM);
		socket.bind_local(s, port);
		socket.listen(s);

		while (true) {
			auto rt = socket.accept(s);
			LOG << socket.recv_all(rt.socket) << std::endl;
		}

		socket.close(s);
	}

} // namespace jian

