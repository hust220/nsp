#include "nsp.hpp"
#include "pdb.hpp"
#include "socket.hpp"

BEGIN_JN

	REGISTER_NSP_COMPONENT(client) {
		S ip = "127.0.0.1";
		int port = 80;
		S content = "hi";
		Socket socket;

		par.set(ip, "ip");
		par.set(port, "port");
		par.set(content, "content");

		auto s = socket.open(SOCK_STREAM);
		socket.connect(s, ip, port);
		socket.send(s, content);
		socket.close(s);
	}

END_JN

