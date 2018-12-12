#include "nsp.hpp"
#include "pdb.hpp"
#include "socket.hpp"

namespace jian {
	namespace {

		class Server {
		public:
			int m_port;

			Server(int port = 80) {
				m_port = port;
			}

		};

		class Request {
		public:
			using head_t = std::map<std::string, std::string>;

			S m_method;
			S m_url;
			S m_version;
			head_t m_head;
			S m_content;

			Request(const S &s) {
				std::istringstream stream;
				stream.str(s);
				tokenize_v v, w;
				int n = 0;
				for (S line; std::getline(stream, line); ) {
					tokenize(line, v, " ");
					//LOG << v.size() << std::endl;
					//for (auto && s : v) LOG << s << ','; LOG << std::endl;
					if (v.empty()) break;
					if (n == 0) {
						m_method = v[0];
						m_url = v[1];
						m_version = v[2];
					}
					else {
						tokenize(v[0], w, ":");
						m_head[w[0]] = w[1];
					}
					n++;
				}
				for (S line; std::getline(stream, line); ) {
					m_content += line;
				}
			}

			friend std::ostream &operator <<(std::ostream &stream, const Request &request) {
				stream << "Method: " << request.m_method << std::endl;
				stream << "URL: " << request.m_url << std::endl;
				stream << "Version: " << request.m_version << std::endl;
				for (auto && pair : request.m_head) {
					stream << pair.first << ": " << pair.second << std::endl;
				}
				stream << request.m_content << std::endl;
				return stream;
			}
		};

		class Response {
		public:
			using head_t = typename Request::head_t;

			S m_version;
			int m_status_code;
			S m_status_identifier;
			head_t m_head;
			S m_content;

			Response(const Request &request, const S &c = "") {
				m_version = request.m_version;
				m_status_code = 200;
				m_status_identifier = "OK";
				m_content = c;
				m_head["Server"] = "Jian";
				content("");
				m_head["Content-Type"] = "text/html";
			}

			Response &content(const S &c) {
				m_content = c;
				m_head["Content-Length"] = JN_STR(c.size());
				return *this;
			}

			friend std::ostream &operator <<(std::ostream &stream, const Response &response) {
				stream
					<< response.m_version << ' ' << response.m_status_code << ' '
					<< response.m_status_identifier << '\r' << '\n';
				for (auto && pair : response.m_head) {
					stream << pair.first << ':' << pair.second << '\r' << '\n';
				}
				stream << '\r' << '\n';
				stream << response.m_content;
				return stream;
			}
		};

		class HttpServer : public Server {
		public:
			std::shared_ptr<Socket> m_socket;

			HttpServer(int port = 80) : Server(port) {}

			template<typename Fn>
			void run(Fn &&fn) {
				m_socket.reset(new Socket);

				auto s = m_socket->open(SOCK_STREAM);
				m_socket->bind_local(s, m_port);
				m_socket->listen(s);

				while (true) {
					auto rt = m_socket->accept(s);
					S r;
					while (true) {
						try {
							S s = m_socket->recv(rt.socket);
							for (int i = 0; i + 3 < s.size(); i++) {
								if (s[i] == '\r' && s[i + 1] == '\n' && s[i + 2] == '\r' && s[i + 3] == '\n') {
									r += s.substr(0, i);
									Request request(r);
									Response response(request);
									fn(request, response);
									send_response(rt.socket, response);
									r = "";
									s = s.substr(i+4);
									break;
								}
							}
							r += s;
						}
						catch (...) {
							break;
						}
					}
				}

				m_socket->close(s);
			}

			void send_response(Socket::socket_t s, const Response &r) {
				std::ostringstream stream;
				stream << r;
				m_socket->send(s, stream.str());
			}

		};

		REGISTER_NSP_COMPONENT(server) {
			int port = 80;

			par.set(port, "port");
			HttpServer server(port);
			server.run([](const Request &request, Response &response) {
				LOG << request << std::endl;
				response.content("Hi, this is the Jian Server!");
			});
		}
	}
}

