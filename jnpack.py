"""========================================
jnpack     A C++ package manager
Author:    Jian Wang
Email:     wj_hust08@hust.edu.cn
Usage:     python jnpack.py <install|update>
========================================"""

import json
import os
import sys
import glob
import fnmatch

SOLUTION_PATH = os.getcwd()
JNP_PATH = SOLUTION_PATH + '/.jnp'
PROJS_PATH = JNP_PATH + '/projs'
JNPACKS_PATH = SOLUTION_PATH + '/jnpacks'

def pars_to_cmake(pars, src):
    print pars
    str = ""
    str += "cmake_minimum_required(VERSION %s)\n" % pars['cmake_version']
    str += "project(%s)\n" % pars['project']
    str += "set(CMAKE_INSTALL_PREFIX %s)\n" % pars['install_prefix']
    str += "set(CMAKE_BUILD_TYPE %s)\n" % pars['build_type']
    for i in pars['defs']:
        str += "add_definitions(-D%s)\n" % i
    str += "set(CMAKE_CXX_COMPILER %s)\n" % pars['cxx_compiler']
    str += "set(CMAKE_VERBOSE_MAKEFILE on)\n"
    str += "set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} %s\")\n" % ' '.join(pars['cxx_flags'])
    for i in glob.glob('%s/*' % JNPACKS_PATH):
        basename = os.path.basename(i)
        if basename != 'CMakeLists.txt':
            str += "include_directories(%s)\n" % i
    str += "add_executable(%s %s)\n" % (pars['target'], ' '.join(src))
    str += "install(TARGETS %s DESTINATION bin)" % pars['target']
    return str

def read_config(dir):
    path = dir + '/jnpack.json'
    if not os.path.exists(path): return {}

    f = open(path)
    content = f.read()
    f.close()

    return json.loads(content)

def add_packs(ls, n):
    if len(ls) == 0: return

    packs = []
    for pack in ls:
        pack_path = JNPACKS_PATH + '/' + pack
        if not os.path.exists(pack_path) or n == 0:
            packs.append(pack)
        if not os.path.exists(pack_path):
            os.system('git clone https://github.com/hust220/%s.git' % pack)

    ls = []
    for pack in packs:
        pack_path = JNPACKS_PATH + '/' + pack
        ls.extend(read_config(pack_path)['depends'])

    add_packs(ls, n + 1)

def find_cpps(dir):
    matches = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in fnmatch.filter(filenames, '*.cpp'):
            matches.append(os.path.join(root, filename))
    return matches

class Jnpack(object):
    def __init__(self):
        self.config = read_config(SOLUTION_PATH)
        self.depends = self.config['depends']
        init_pars = {\
            "cmake_version": "2.8.7",\
            "install_prefix": JNP_PATH,\
            "build_type": "Release"\
            }
        self.projs = self.config['projs']
        self.projs_elim_extends()
        self.projs_merge(init_pars)

    def projs_merge(self, init_pars):
        for proj_name, proj in self.projs.iteritems():
            proj['project'] = proj_name
            for k in init_pars:
                if k not in proj:
                    proj[k] = init_pars[k]

    def projs_elim_extends(self):
        n = -1
        while n != 0:
            n = 0
            for proj_name, proj in self.projs.iteritems():
                for k, v in proj.iteritems():
                    if k == 'extends':
                        n += 1
                        del proj[k]
                        for k2, v2 in self.projs[v].iteritems(): 
                            if k2 not in proj:
                                proj[k2] = v2
                        break

    def get_jnpacks_cpps(self):
        ls = []
        for it in glob.glob('%s/*' % JNPACKS_PATH):
            basename = os.path.basename(it)
            if basename != "CMakeLists.txt":
                ls.extend(find_cpps(it + '/' + basename))
        return ls
        #return "set(jnpacks_sources %s)\n" % (' '.join(ls),)

    def get_dirs_cpps(self, dir):
        return find_cpps(SOLUTION_PATH + '/' + dir)

    def get_proj_src(self, proj_name):
        ls = []
        ls.extend(self.get_jnpacks_cpps())
        #print self.projs[proj_name]
        for dir in self.projs[proj_name]['dirs']:
            ls.extend(self.get_dirs_cpps(dir))
        return ls

    def gen_cmake(self, proj_name):
        if not os.path.exists(JNP_PATH): os.mkdir(JNP_PATH)
        if not os.path.exists(PROJS_PATH): os.mkdir(PROJS_PATH)
        proj_path = PROJS_PATH + '/' + proj_name
        if not os.path.exists(proj_path): os.mkdir(proj_path)
        os.chdir(proj_path)
        f = open('CMakeLists.txt', 'w')
        f.write(pars_to_cmake(self.projs[proj_name], self.get_proj_src(proj_name)))
        f.close()

    def build(self, proj_name):
        self.gen_cmake(proj_name)

        proj_path = PROJS_PATH + '/' + proj_name
        os.chdir(proj_path)
        os.system('cmake .')
        os.system('make -j %d' % self.projs[proj_name]['threads'])

    def install(self, proj_name):
        proj_path = PROJS_PATH + '/' + proj_name
        self.build(proj_name)
        os.chdir(proj_path)
        os.system('make install')

    def update(self):
        if not os.path.exists(JNPACKS_PATH): os.mkdir(JNPACKS_PATH)
        os.chdir(JNPACKS_PATH)
        add_packs(self.depends, 0)
        os.chdir(SOLUTION_PATH)

    def upgrade(self):
        self.update()

        os.chdir(JNPACKS_PATH)
        for it in os.listdir('.'):
            if it != 'CMakeLists.txt':
                print it
                os.chdir(it)
                os.system('git pull')
                os.chdir('..')
        os.chdir(SOLUTION_PATH)
        os.system('git pull')

    def help(self):
        print __doc__

    def init(self):
        print __doc__

if __name__ == '__main__':
    jnpack = Jnpack()
    if len(sys.argv) == 1:
        jnpack.help()
    else:
        func = sys.argv[1]
        if func == 'build':
            if len(sys.argv) == 2: jnpack.build('default')
            else: jnpack.build(sys.argv[2])
        elif func == 'gen-cmake':
            jnpack.gen_cmake()
        elif func == 'install':
            if len(sys.argv) == 2: jnpack.install('default')
            else: jnpack.install(sys.argv[2])
        elif func == 'init':
            jnpack.init()
        elif func == 'upgrade':
            jnpack.upgrade()
        elif func == 'update':
            jnpack.update()
        else:
            jnpack.help()

