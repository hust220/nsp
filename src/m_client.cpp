#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/socket.hpp>

namespace jian {

	REGISTER_NSP_COMPONENT(client) {
		std::string ip = "127.0.0.1";
		int port = 80;
		std::string content = "hi";
		Socket socket;

		par.set(ip, "ip");
		par.set(port, "port");
		par.set(content, "content");

		auto s = socket.open(SOCK_STREAM);
		socket.connect(s, ip, port);
		socket.send(s, content);
		socket.close(s);
	}

} // namespace jian

