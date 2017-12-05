#include "rtsp_sample_loop.hpp"
#include "rna_mc.hpp"

BEGIN_JN

void SampleLoop::init(const Chain &chain, S ss) {
	m_mc = STD_ make_shared<DHMC>();
	m_mc->init(Par("seq", seq(chain))("ss", ss)("name", "")("queue", "samc:30000:200"));
	m_mc->_pred_chain = chain;
}

Chain SampleLoop::operator ()() {
	if (!m_mc) throw "jian::SampleLoop error! Please initialize first!";
	m_mc->run();
	return m_mc->_pred_chain;
}

END_JN
